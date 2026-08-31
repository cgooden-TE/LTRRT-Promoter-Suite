#!/usr/bin/env python3
"""
Classifier for transposable-element and gene isoforms with read/TSS summaries across replicates.

This classifier is intended for use with high quality long reads (PacBio HiFi or ONT Q20+) that are
mapped to a reference genome and can be accurately assigned a locus and an isoform category. It is
not intended for use with short reads or low quality long reads.

Feature model
-------------
Every annotated feature that can originate a transcript is loaded as an "element" with a class:

  LTR_structural : LTR retrotransposon with a Parent attribute and both LTRs resolved
                   (EDTA Method=structural). TSS is assigned from the LTR the read starts in,
                   so these get the full LTR-aware isoform vocabulary.
  LTR_fragment   : LTR retrotransposon called by homology (EDTA Method=homology). No usable LTR
                   pair, so TSS is assigned from the read coordinate within the element body.
  TIR, Helitron,
  LINE, SINE     : DNA transposons and non-LTR retrotransposons. Body-coordinate TSS.
  Gene           : protein-coding gene. Body-coordinate TSS, with full-length judged against the
                   union of annotated exons rather than the gene span.

Structural bookkeeping rows (repeat_region, target_site_duplication) and non-transposon repeats
(knob, centromeric_repeat, subtelomere, rDNA_intergenic_spacer_element, low_complexity) are dropped.

Assignment
----------
A read is assigned to the *innermost* feature containing its TSS: of every feature whose span
contains the 5' end, the one with the smallest span wins. Every other containing feature is
reported as a nesting parent, so a read starting inside a TE that sits inside a gene is credited
to the TE and flagged as nested in the gene. Genes and TEs compete in the same index, so no
gene-versus-TE precedence rule is needed.

Orientation
-----------
Categories are computed in the *read's* frame: the LTR a read initiates in is its 5' LTR whether
or not the read agrees with the element's annotated strand. A record is labelled sense when the
read strand matches the element strand and antisense when it does not, and the two are written to
separate isoform files. Features annotated with an unknown strand ('.' or '?') are resolved after
classification by majority read strand.

TSS agreement
-------------
For each feature and orientation the modal 5' end is taken as the peak, and the fraction of that
feature's reads falling in the window centred on it is reported. A peak holding at least
--tss-min-frac of the reads is flagged TSS_Called=1. This is reported, never used to discard reads.

Inputs:
    - GFF with TE annotations (EDTA-style) and genes
    - Optional gene GFF3 with exon records, for exon-aware gene full-length calls
    - One or more BAM files with mapped reads
    - Minimum MAPQ for filtering reads

Outputs:
    - <prefix>_sense.isoforms.tsv / <prefix>_antisense.isoforms.tsv
    - TSS, cleavage, and density summaries carrying an Orientation column
    - Gene read count and TSS summary carrying an Orientation column
    - Per-read soft-clip and exon-statistics tables
    - Optional U3/LTR and promoter FASTAs

Example Usage:
    python3 IsoClassifier.py \
        --gff TEs_LTRs_Genes_UpdStr.gff \
        --gene-gff Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.gff3 \
        --bam sample1.bam sample2.bam \
        --min_mapq 30 \
        --threads 4 \
        --output ltr_isoforms \
        --tss_out ltr_tss_summary.tsv \
        --gene_out gene_summary.tsv
"""
import argparse
import logging
import re
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
import pandas as pd
import pysam
import numpy as np
from intervaltree import IntervalTree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Feature classes and isoform vocabularies
# ---------------------------------------------------------------------------

# Feature-type patterns mapped to element classes. Order matters: first match wins.
# LTR retrotransposons are split into structural/fragment later, by Method and LTR availability.
CLASS_PATTERNS = [
    (re.compile(r'LTR_retrotransposon$'),        'LTR'),
    (re.compile(r'_TIR_transposon$'),            'TIR'),
    (re.compile(r'^helitron$', re.I),            'Helitron'),
    (re.compile(r'LINE'),                        'LINE'),
    (re.compile(r'SINE'),                        'SINE'),
    (re.compile(r'^gene$'),                      'Gene'),
]

# Rows that describe annotation bookkeeping or non-transposon repeats. Never become elements.
DROP_FEATURES = {
    'repeat_region', 'target_site_duplication', 'knob', 'centromeric_repeat',
    'subtelomere', 'rDNA_intergenic_spacer_element', 'low_complexity',
    'chromosome', 'scaffold', 'contig',
}

# Isoform vocabulary for structural LTR-RTs, named in the read's frame.
LTR_CATS = ['ltr5_contained', 'ltr3_contained', 'spanning',
            'readout_5ltr', 'readout_3ltr', 'readout_internal',
            'spliced_ltr5', 'spliced_ltr3', 'spliced_spanning', 'partial']

# Isoform vocabulary for every other class, and for structural elements whose TSS falls in the
# internal domain rather than an LTR.
SIMPLE_CATS = ['full_length', 'spliced', 'readout', 'partial']

# Union used for the wide isoform table; a feature only ever populates its own vocabulary.
ALL_CATS = LTR_CATS + [c for c in SIMPLE_CATS if c not in LTR_CATS]

# Categories whose membership test already fixes the splice counter, so emitting it would just
# restate another column: these require a real intron, so n_spliced always equals <cat>_reads.
SPLICED_BY_DEF = {'spliced_ltr5', 'spliced_ltr3', 'spliced_spanning', 'spliced'}

# ...and these require the absence of one, so n_spliced and unique_juncts are always 0.
UNSPLICED_BY_DEF = {'ltr5_contained', 'ltr3_contained'}

# Classes eligible for the LTR-aware vocabulary.
LTR_CLASSES = {'LTR_structural'}

# Tie-break rank when two containing features have identical spans. Lower wins.
CLASS_RANK = {'LTR_structural': 0, 'LTR_fragment': 1, 'TIR': 2, 'Helitron': 3,
              'LINE': 4, 'SINE': 5, 'Gene': 6}

DEFAULT_TE_CLASSES = 'LTR_structural,LTR_fragment,TIR,Helitron,LINE,SINE,Gene'


def parse_args():
    parser = argparse.ArgumentParser(
        description='Classify TE and gene isoforms with sense/antisense read and TSS summaries.'
    )
    parser.add_argument('--gff', required=True,
                        help='GFF with TE annotations (EDTA-style), LTRs, and genes')
    parser.add_argument('--gene-gff', dest='gene_gff', default=None,
                        help='Optional GFF3 with exon records (e.g. the Ensembl gene annotation). '
                             'Exons are keyed to genes by gene ID, so seqid naming need not match '
                             'the main GFF, but coordinates must be on the same assembly. Without '
                             'it, gene full-length falls back to gene-span coverage, which almost '
                             'never fires for spliced transcripts.')
    parser.add_argument('--canonical-exons-only', dest='canonical_exons_only',
                        action='store_true',
                        help='Use only the Ensembl_canonical transcript when building each gene\'s '
                             'exon union (default: union across all transcripts)')
    parser.add_argument('--bam', nargs='+', required=True,
                        help='One or more BAM files')
    parser.add_argument('--te-classes', dest='te_classes', default=DEFAULT_TE_CLASSES,
                        help=f'Comma-separated element classes to classify (default: {DEFAULT_TE_CLASSES})')
    parser.add_argument('--trust-st-tag', dest='trust_st_tag', action='store_true',
                        help='Trust the PyChopper/AccuMap ST tag as a genomic strand call. '
                             'OFF by default: PyChopper reorients reads, so ST records a '
                             'transformation already applied and is ~50%% accurate (no better '
                             'than chance). Only enable for BAMs whose ST was written by '
                             'something that did NOT reorient reads.')
    parser.add_argument('--min_mapq', type=int, default=30,
                        help='Minimum MAPQ for filtering reads')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of parallel processes for BAM classification (default: 1)')
    parser.add_argument('--min-intron-len', dest='min_intron_len', type=int, default=69,
                        help='Minimum CIGAR N length counted as a real intron (default: 69)')
    parser.add_argument('--full-length-frac', dest='full_length_frac', type=float, default=0.8,
                        help='Covered fraction of an element (or of a gene exon union) at or above '
                             'which a contained read is called full_length (default: 0.8)')
    parser.add_argument('--spliced-cov-frac', dest='spliced_cov_frac', type=float, default=0.5,
                        help='Covered fraction below which a spliced read is called spliced rather '
                             'than full_length/spanning (default: 0.5)')
    parser.add_argument('--tss-window', dest='tss_window', type=int, default=3,
                        help='Width in bp of the TSS agreement window, centred on the modal 5\' end '
                             '(default: 3, i.e. peak-1..peak+1)')
    parser.add_argument('--tss-min-frac', dest='tss_min_frac', type=float, default=0.5,
                        help='Fraction of a feature\'s reads that must fall in the TSS window for '
                             'TSS_Called=1 (default: 0.5). Reported only; never discards reads.')
    parser.add_argument('--progress-every', dest='progress_every', type=int, default=1_000_000,
                        help='Print a progress line every N reads assessed, with an ETA read from '
                             'the BAM indexes (default: 1000000)')
    parser.add_argument('--include-partial', dest='include_partial', action='store_true',
                        help='Include the partial category in the isoform tables. Partial reads are '
                             'always written to the per-read tables with Nested_In and '
                             'Terminates_In populated, regardless of this flag.')
    parser.add_argument('--output', required=True,
                        help='Output prefix for isoform TSVs and summaries')
    parser.add_argument('--tss_out', required=True,
                        help='Output file for the isoform TSS summary')
    parser.add_argument('--gene_out', required=True,
                        help='Output file for gene read counts and TSS summary')
    parser.add_argument('--genome-fasta', dest='genome_fasta', default=None,
                        help='Genome FASTA for U3/promoter sequence extraction '
                             '(enables u3_seq_extraction after classification)')
    return parser.parse_args()


# Helper function to parse GFF attributes into a dictionary
def parse_attributes(attrs):
    attr_dict = {}
    for pair in str(attrs).split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attr_dict[key.strip()] = value.strip()
    return attr_dict


def _base_class(feature):
    """Map a GFF feature type to a coarse element class, or None if it is not an element."""
    if feature in DROP_FEATURES:
        return None
    for pat, klass in CLASS_PATTERNS:
        if pat.search(feature):
            return klass
    return None


def _merge_intervals(ivs):
    """Merge a list of (start, end) half-open intervals; returns sorted, non-overlapping tuples."""
    if not ivs:
        return ()
    ivs = sorted(ivs)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return tuple((s, e) for s, e in merged)


def load_gene_exons(gene_gff, canonical_only=False):
    """
    Build gene_id -> (merged_exon_intervals, exonic_length) from a GFF3 carrying exon records.

    Exons are keyed to genes through their transcript Parent, so only IDs need to agree with the
    main GFF; the seqid column is ignored entirely. That matters because the Ensembl annotation
    names contigs '1'..'10' while the combined TE/gene reference names them 'B73_chr1'..'B73_chr10'
    on identical coordinates.
    """
    print(f"Reading exon annotation from {gene_gff} ...")
    cols = ['chrom', 'src', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attrs']
    df = pd.read_csv(gene_gff, sep='\t', comment='#', header=None, names=cols,
                     dtype={'attrs': str}, low_memory=False)

    tx2gene = {}
    n_canon = 0
    tx_rows = df[df['feature'].isin(['mRNA', 'transcript'])]
    for feature, attrs in zip(tx_rows['feature'], tx_rows['attrs']):
        a = parse_attributes(attrs)
        tid = a.get('ID', '')
        gid = a.get('Parent', '')
        if not tid or not gid:
            continue
        if canonical_only and 'Ensembl_canonical' not in str(attrs):
            continue
        # Strip Ensembl 'transcript:' / 'gene:' prefixes so IDs match the main GFF's bare form
        tx2gene[tid.split(':', 1)[-1] if ':' in tid else tid] = gid.split(':', 1)[-1] if ':' in gid else gid
        n_canon += 1

    per_gene = defaultdict(list)
    ex_rows = df[df['feature'] == 'exon']
    for s, e, attrs in zip(ex_rows['start'], ex_rows['end'], ex_rows['attrs']):
        a = parse_attributes(attrs)
        tid = a.get('Parent', '')
        if not tid:
            continue
        tid = tid.split(':', 1)[-1] if ':' in tid else tid
        gid = tx2gene.get(tid)
        if gid is None:
            continue
        per_gene[gid].append((int(s) - 1, int(e)))

    gene_exons = {}
    for gid, ivs in per_gene.items():
        merged = _merge_intervals(ivs)
        gene_exons[gid] = (merged, sum(e - s for s, e in merged))

    print(f"  {len(gene_exons)} genes with exons from {n_canon} transcripts "
          f"({'canonical only' if canonical_only else 'all transcripts'}).")
    return gene_exons


def load_elements_and_ranges(gff_path, te_classes, gene_gff=None, canonical_exons_only=False):
    """
    Load every transcribable feature as an element and index bodies in one interval tree.

    Returns
    -------
    feat_info : dict fid -> element record
    body_tree : dict chrom -> IntervalTree, data = (fid, class_rank)
    """
    print("Reading GFF and parsing elements...")
    cols = ['chrom', 'src', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attrs']
    df = pd.read_csv(gff_path, sep='\t', comment='#', header=None, names=cols,
                     dtype={'attrs': str}, low_memory=False)
    df = df.dropna(subset=['feature', 'start', 'end'])

    keep = set(c.strip() for c in te_classes.split(',') if c.strip())

    feat_info = {}
    ltr_rows = []          # (parent_id, ltr_id, start0, end)
    class_counts = Counter()
    skipped = Counter()

    for row in df.itertuples(index=False):
        feature = str(row.feature)

        # LTR sub-features are consumed into their structural parent, never elements themselves
        if feature == 'long_terminal_repeat':
            a = parse_attributes(row.attrs)
            parent = a.get('Parent')
            if parent:
                ltr_rows.append((parent, a.get('ID', ''), int(row.start) - 1, int(row.end)))
            continue

        base = _base_class(feature)
        if base is None:
            skipped[feature] += 1
            continue

        a = parse_attributes(row.attrs)

        if base == 'LTR':
            method = a.get('Method', '')
            parent = a.get('Parent')
            if parent and method == 'structural':
                # Structural elements are keyed on Parent so their long_terminal_repeat
                # children (which carry the same Parent) resolve to the same record.
                fid, klass = parent, 'LTR_structural'
            else:
                fid, klass = a.get('ID'), 'LTR_fragment'
        else:
            fid, klass = a.get('ID'), base

        if not fid:
            skipped[feature] += 1
            continue

        feat_info[fid] = {
            'chrom': row.chrom,
            'start': int(row.start) - 1,
            'end': int(row.end),
            'strand': row.strand if row.strand in ('+', '-') else '.',
            'strand_known': row.strand in ('+', '-'),
            'name': a.get('Name', ''),
            'attrs': row.attrs,
            'class': klass,
            'ltr_left': None,
            'ltr_right': None,
            'exons': None,
            'exon_len': 0,
        }
        class_counts[klass] += 1

    # Attach LTR coordinates to their structural parents
    n_ltr = 0
    for parent, lid, s, e in ltr_rows:
        rec = feat_info.get(parent)
        if rec is None or rec['class'] != 'LTR_structural':
            continue
        if lid.startswith('l'):
            rec['ltr_left'] = (s, e)
            n_ltr += 1
        elif lid.startswith('r'):
            rec['ltr_right'] = (s, e)
            n_ltr += 1

    # A structural element missing either LTR cannot support LTR-based TSS assignment,
    # so it is demoted to the body-coordinate vocabulary alongside homology fragments.
    n_demoted = 0
    for fid, rec in feat_info.items():
        if rec['class'] == 'LTR_structural' and not (rec['ltr_left'] and rec['ltr_right']):
            rec['class'] = 'LTR_fragment'
            rec['ltr_left'] = rec['ltr_right'] = None
            n_demoted += 1
    if n_demoted:
        class_counts['LTR_structural'] -= n_demoted
        class_counts['LTR_fragment'] += n_demoted

    # Drop classes the user excluded
    if keep:
        dropped = [fid for fid, rec in feat_info.items() if rec['class'] not in keep]
        for fid in dropped:
            del feat_info[fid]

    # Exon annotation for gene full-length calls
    if gene_gff:
        gene_exons = load_gene_exons(gene_gff, canonical_exons_only)
        n_hit = 0
        n_bad = 0
        for fid, rec in feat_info.items():
            if rec['class'] != 'Gene':
                continue
            ex = gene_exons.get(fid)
            if not ex:
                continue
            merged, exlen = ex
            # Sanity check: exons must sit inside the gene span from the main GFF, otherwise the
            # two annotations are on different assemblies and the coverage call would be garbage.
            if merged[0][0] < rec['start'] or merged[-1][1] > rec['end']:
                n_bad += 1
                continue
            rec['exons'] = merged
            rec['exon_len'] = exlen
            n_hit += 1
        n_genes = class_counts.get('Gene', 0)
        print(f"  Attached exons to {n_hit}/{n_genes} genes"
              + (f"; {n_bad} rejected for falling outside the gene span" if n_bad else ""))
        if n_bad > n_hit:
            logging.warning("Most exon records fall outside their gene span -- the two "
                            "annotations are probably on different assemblies. Gene full_length "
                            "calls will fall back to gene-span coverage.")

    print("Element classes loaded: " + ", ".join(f"{k}={v}" for k, v in sorted(class_counts.items())))
    if skipped:
        top = ", ".join(f"{k}={v}" for k, v in skipped.most_common(6))
        print(f"  Skipped non-element rows: {top}")
    print(f"  {n_ltr} long_terminal_repeat records attached; "
          f"{n_demoted} structural elements demoted to LTR_fragment for missing an LTR.")

    print("Indexing element bodies with intervaltree...")
    body_tree = defaultdict(IntervalTree)
    for fid, rec in feat_info.items():
        if rec['end'] <= rec['start']:
            continue
        body_tree[rec['chrom']][rec['start']:rec['end']] = (fid, CLASS_RANK.get(rec['class'], 9))
    print(f"Indexed {len(feat_info)} elements across {len(body_tree)} contigs.")

    return feat_info, dict(body_tree)


# Get the strand for the read according to transcription direction determined by PyChopper if available.
# Fallbacks are in descending priority.
def infer_read_strand(read, trust_st_tag=False):
    """
    Return '+', '-', or '.' using fallbacks:
      1) alignment orientation (is_reverse)  -- PRIMARY, see note below
      2) ST tag (PyChopper)  -- opt-in only, see below
      3) ts tag (minimap2 transcript strand)
      4) XS tag
      5) jM tag (minimap2 splice motif)

    Why alignment orientation is primary
    -------------------------------------
    PyChopper REVERSE-COMPLEMENTS every read it calls '-' before it ever reaches the
    mapper, so every read classified in this pipeline is already in mRNA-sense
    (5'->3') orientation -- the read's sequence *is* the transcript's sequence. That
    makes the genomic strand it aligns to the strand it was transcribed from: forward
    alignment of a sense-oriented read is a '+' transcript, reverse alignment is '-'.
    This is what `origin`/`endpos` (and therefore every TSS/isoform call below) are
    computed from, so it needs to be right, not just plausible.

    Measured on the Dec2025 ONT libraries (CT1/CT2/CT3), scoring against intron motifs
    and, independently, against annotated gene strand on unspliced reads:

        strand = ST                    50.10%   (n ~ 600k, genome-wide)
        strand = ST XOR is_reverse     49.41%
        strand = alignment orientation 99.98%   (99.83% spliced / 99.44% unspliced)

    ST records which way PyChopper flipped a read *before* mapping (provenance, not
    post-mapping genomic strand) and is statistically independent of true strand.
    `ts`/`jM` are splice-motif calls that only add independent information when
    minimap2 is allowed to search both strands (-ub); under the old -uf mapping
    setting they were a forced constant ('ts:A:+' on every read) that happened to
    equal is_reverse's complement for every '-'-strand read, silently overriding it.
    Because ST/ts sat above alignment orientation in the old cascade, and 100% of ONT
    primary alignments carried one or the other, read strand -- and therefore every
    ONT TSS call -- was effectively a coin flip between the read's two ends.

    ST/ts/XS/jM are now fallbacks only, for the case alignment orientation is somehow
    unavailable (it should always be available for a mapped, non-secondary,
    non-supplementary read, which is the only kind this is called on). Set
    trust_st_tag=True only for BAMs whose ST was written by something that did NOT
    reorient the reads.

    Ref: splice_ESE/st_tag_diagnosis/FINDINGS_ST_TAG_20260804.md
    """
    # 1) alignment orientation -- primary, see docstring above.
    if not read.is_unmapped:
        try:
            return "-" if read.is_reverse else "+"
        except Exception:
            pass

    # 2) PyChopper ST tag -- opt-in only (see docstring)
    if trust_st_tag and read.has_tag("ST"):
        st = read.get_tag("ST")
        if st in {"+", "-"}:
            return st

    # 3) minimap2 transcript strand
    if read.has_tag("ts"):
        # For some PacBio alignments, minimap2 assumes the query is the coding strand
        # and emits ts:A:+ statically, ignoring whether it mapped forward or reverse.
        # We need to XOR it with the alignment orientation just like the ST tag.
        ts = read.get_tag("ts")
        if ts in {"+", "-"}:
            if read.is_reverse:
                return "-" if ts == "+" else "+"
            else:
                return ts

    # 4) XS tag (common in some spliced alignments)
    if read.has_tag("XS"):
        xs = read.get_tag("XS")
        if xs in {"+", "-"}:
            return xs

    # 5) minimap2 intron motif-based inference
    # jM codes: 1=GT-AG(+), 2=CT-AC(-), 3=GC-AG(+), 4=CT-GC(-)
    if read.has_tag("jM"):
        try:
            jm = read.get_tag("jM")
            # jM can be a list/tuple of codes; take the first informative one
            if isinstance(jm, (list, tuple)):
                for code in jm:
                    try:
                        c = int(code)
                    except Exception:
                        continue
                    if c in (1, 3):
                        return "+"
                    if c in (2, 4):
                        return "-"
            else:
                # Sometimes it may be encoded differently; try a single value
                try:
                    c = int(jm)
                    if c in (1, 3):
                        return "+"
                    if c in (2, 4):
                        return "-"
                except Exception:
                    pass
        except Exception:
            pass

    # 6) unknown
    return "."

        
# Collects soft clipping data from 3' ends
def softclip_3prime(read):
    """
    Return soft-clipped bases at the biological 3' end (CIGAR op 4).
    Uses read.is_reverse to decide which end is 3'.
    """
    if not read.cigartuples:
        return 0
    op, length = (read.cigartuples[0] if read.is_reverse else read.cigartuples[-1])
    return length if op == 4 else 0  # 4 = S (soft clip)

# Collects intron and exon statistics from CIGAR tuples, with flexible length modes and intron definitions
def exon_intron_stats_from_cigar(cigartuples, min_intron_len=69, length_mode="query"):
    """
    Parse CIGAR tuples into exon + intron statistics.
    length_mode:
      - "ref": exon length in reference coords (M,=,X,D)
      - "query": exon length in query/read bases (M,=,X,I)

    Introns are derived from N ops (op==3) with length >= min_intron_len.
    Returns:
      exon_count, intron_count,
      exon_lens, intron_lens,
      exon_total, intron_total,
      saw_real_intron
    """
    if not cigartuples:
        return 0, 0, [], [], 0, 0, False

    if length_mode == "ref":
        exon_ops = {0, 2, 7, 8}   # M, D, =, X
    elif length_mode == "query":
        exon_ops = {0, 1, 7, 8}   # M, I, =, X
    else:
        raise ValueError("length_mode must be 'ref' or 'query'")

    exon_lens = []
    intron_lens = []
    cur_exon = 0
    saw_real_intron = False

    for op, length in cigartuples:
        # real intron boundary
        if op == 3 and length >= min_intron_len:  # N
            saw_real_intron = True
            exon_lens.append(cur_exon)
            intron_lens.append(length)
            cur_exon = 0
            continue

        # accumulate exon length
        if op in exon_ops:
            cur_exon += length

        # ignore S/H and also ignore N below threshold

    # last exon
    exon_lens.append(cur_exon)

    # If no real intron, treat as single-exon (and no introns)
    if not saw_real_intron:
        exon_lens = [exon_lens[0]]
        return 1, 0, exon_lens, [], exon_lens[0], 0, False

    # Drop 0-length exons (can occur if cigar begins/ends with N)
    exon_lens = [e for e in exon_lens if e > 0]

    # If dropped exons, intron count may no longer equal exon_count-1.
    # Usually this only happens with pathological CIGARs (starting/ending with N).
    # Optionally you can also trim introns to match exon_count-1:
    if len(intron_lens) > max(0, len(exon_lens) - 1):
        intron_lens = intron_lens[:len(exon_lens) - 1]

    exon_total = sum(exon_lens)
    intron_total = sum(intron_lens)
    return len(exon_lens), len(intron_lens), exon_lens, intron_lens, exon_total, intron_total, True

# Builds the alternating exon/intron fields for output, with padding up to max_exons and intron definition based on min_intron_len
def exon_intron_row_fields(read, min_intron_len=69, max_exons=5, length_mode="query"):
    """
    Returns:
      exon_count, intron_count, exon_total, intron_total, padded_fields

    padded_fields are alternating: exon1, intron1, exon2, intron2, ... up to max_exons
    (so intron slots are max_exons-1).
    """
    (ex_n, in_n,
     ex_lens, in_lens,
     ex_total, in_total,
     saw_real) = exon_intron_stats_from_cigar(read.cigartuples, min_intron_len, length_mode)

    # Build alternating exon/intron columns
    # For max_exons=5 => fields: exon1 intron1 exon2 intron2 exon3 intron3 exon4 intron4 exon5
    fields = []
    ex_keep = ex_lens[:max_exons]
    in_keep = in_lens[:max_exons - 1]

    for i in range(max_exons):
        fields.append(ex_keep[i] if i < len(ex_keep) else "NA")
        if i < max_exons - 1:
            fields.append(in_keep[i] if i < len(in_keep) else "NA")

    return ex_n, in_n, ex_total, in_total, fields, saw_real


# ---------------------------------------------------------------------------
# Assignment helpers
# ---------------------------------------------------------------------------

def junctions_from_read(read, min_intron_len):
    """Reference-coordinate (donor, acceptor) pairs for every real intron in a read."""
    juncs = set()
    if not read.cigartuples:
        return juncs
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 2, 7, 8):        # M, D, =, X consume reference
            pos += length
        elif op == 3:                 # N
            if length >= min_intron_len:
                juncs.add((pos, pos + length))
            pos += length
    return juncs


def covered_bases(blocks, intervals):
    """Aligned read bases falling inside a sorted set of half-open reference intervals."""
    total = 0
    for b1, b2 in blocks:
        for s, e in intervals:
            if e <= b1:
                continue
            if s >= b2:
                break
            total += min(b2, e) - max(b1, s)
    return total


def innermost_at(tree_chrom, pos):
    """
    Resolve the feature a position belongs to: of every feature whose span contains pos, the one
    with the smallest span wins, breaking ties by class rank then feature ID.

    Returns (winner_fid, [other containing fids]). This single rule implements both the nested-TE
    assignment and gene-versus-TE precedence: a read starting inside a TE that sits inside a gene
    goes to the TE, with the gene reported as a nesting parent.
    """
    if tree_chrom is None:
        return None, []
    hits = tree_chrom.at(pos)
    if not hits:
        return None, []
    best = None
    best_key = None
    for iv in hits:
        fid, rank = iv.data
        key = (iv.end - iv.begin, rank, fid)
        if best_key is None or key < best_key:
            best_key, best = key, fid
    others = [iv.data[0] for iv in hits if iv.data[0] != best]
    return best, others


def classify_ltr(rec, origin, endpos, strand, blocks, has_intron, spliced_cov_frac):
    """
    Isoform category for a structural LTR-RT, computed in the read's frame.

    The LTR the read initiates in is its 5' LTR regardless of the element's annotated strand, so
    an antisense read on a '+' element is evaluated against the element's right LTR. That is the
    whole of the sense/antisense inversion; nothing else needs a second code path.

    A consequence worth stating, because it is the only way one of these categories arises: a read
    that starts in the element's annotated 3' LTR and ends in its annotated 5' LTR is not a
    backwards transcript, it is an antisense one. In its own frame it started in its 5' LTR and
    ended in its 3' LTR, so it is scored as spanning and lands in the antisense table. The
    genomically-backwards case cannot occur for a sense read and is not represented.

    The three spliced categories are kept distinct because they are different transcripts:
      spliced_ltr5 / spliced_ltr3 : contained within the LTR the read initiated in, carrying a
        real intron. This is the dominant observed phenotype, so containment alone is not enough
        to call a read ltr5_contained -- that category means contained *and* unspliced.
      spliced_spanning : spans both LTRs but skipped most of the internal domain (coverage below
        spliced_cov_frac). Formerly folded into spliced_ltr5, which made the pair asymmetric --
        a spanning read always initiates in its 5' LTR, so the demotion could never yield
        spliced_ltr3 and the contained phenotype could not be counted on its own.
    """
    ltr_l, ltr_r = rec['ltr_left'], rec['ltr_right']
    ltr5, ltr3 = (ltr_l, ltr_r) if strand == '+' else (ltr_r, ltr_l)

    in5 = ltr5[0] <= origin < ltr5[1]
    in3 = ltr3[0] <= origin < ltr3[1]
    past3 = endpos >= ltr3[1] if strand == '+' else endpos < ltr3[0]

    if in5:
        if ltr5[0] <= endpos < ltr5[1]:
            # Started and finished inside the LTR it initiated in. A real intron makes this a
            # spliced isoform rather than a plain LTR-contained read.
            return 'spliced_ltr5' if has_intron else 'ltr5_contained'
        if ltr3[0] <= endpos < ltr3[1]:
            cat = 'spanning'
        elif past3:
            return 'readout_5ltr'
        else:
            return 'partial'
    elif in3:
        if ltr3[0] <= endpos < ltr3[1]:
            return 'spliced_ltr3' if has_intron else 'ltr3_contained'
        if past3:
            return 'readout_3ltr'
        return 'partial'
    else:
        # TSS sits in the internal domain, so this is not an LTR-driven transcript. It gets its
        # own read-out category rather than the body vocabulary's plain 'readout', which would
        # have put three different meanings of "ran past the 3' end" in one class's columns.
        return 'readout_internal' if past3 else 'partial'

    # A spanning read that skipped most of the internal domain is a spliced isoform, not coding.
    # It keeps its own category rather than being folded in with the LTR-contained spliced reads.
    coding_start, coding_end = ltr_l[1], ltr_r[0]
    coding_len = coding_end - coding_start
    if coding_len > 0 and has_intron:
        if covered_bases(blocks, ((coding_start, coding_end),)) / coding_len < spliced_cov_frac:
            return 'spliced_spanning'
    return cat


def classify_simple(rec, origin, endpos, strand, blocks, has_intron,
                    full_length_frac, spliced_cov_frac):
    """
    Isoform category from body coordinates, for fragments, DNA TEs, LINEs/SINEs, and genes.

    For genes with exon annotation the covered fraction is measured against the union of annotated
    exons rather than the gene span. Measuring against the span would put essentially every spliced
    transcript below the full-length threshold, since introns contribute nothing to the alignment.
    """
    b0, b1 = rec['start'], rec['end']
    past_end = endpos >= b1 if strand == '+' else endpos < b0
    if past_end:
        return 'readout'

    if rec['exons']:
        intervals, denom = rec['exons'], rec['exon_len']
    else:
        intervals, denom = ((b0, b1),), b1 - b0
    if denom <= 0:
        return 'partial'

    frac = covered_bases(blocks, intervals) / denom
    if frac >= full_length_frac:
        return 'full_length'
    if has_intron and frac < spliced_cov_frac:
        return 'spliced'
    return 'partial'


# ---------------------------------------------------------------------------
# Per-chunk worker for multiprocessing
# ---------------------------------------------------------------------------

# Globals for worker inheritance (zero-copy on Linux via fork)
_GLOBAL_FEAT_INFO = {}
_GLOBAL_BODY_TREE = {}


def init_worker(feat_info, body_tree):
    """Initializer to populate globals for spawned workers (Windows/Mac)."""
    global _GLOBAL_FEAT_INFO, _GLOBAL_BODY_TREE
    _GLOBAL_FEAT_INFO = feat_info
    _GLOBAL_BODY_TREE = body_tree


def _empty_result(chrom):
    return {'chrom': chrom, 'n_reads': 0, 'n_assigned': 0,
            'strand_counts': Counter(), 'stats': {},
            'tss_positions': {}, 'end_positions': {},
            'clip_rows_splice': [], 'clip_rows_nonsplice': [],
            'elem_exon_rows': [], 'gene_exon_rows': []}


def _classify_chunk(chrom, claim_start, claim_end, fetch_start, fetch_end, bam_paths, min_mapq,
                    trust_st_tag, min_intron_len, full_length_frac, spliced_cov_frac):
    """Process primary reads whose alignment starts in [claim_start, claim_end) on *chrom*.

    Reads are fetched from [fetch_start, fetch_end) (which may extend past the claim bounds by the
    overlap margin so reads straddling a chunk boundary are still seen), but only reads with
    reference_start inside [claim_start, claim_end) are counted -- that is the deduplication rule
    across adjacent chunks. Annotation lookups use the whole-chromosome tree, so a read is never
    truncated by chunk boundaries.
    """
    feat_info = _GLOBAL_FEAT_INFO
    body_tree_chrom = _GLOBAL_BODY_TREE.get(chrom, None)

    if not body_tree_chrom:
        return _empty_result(chrom)

    n_reads = 0
    n_assigned = 0
    strand_counts = Counter()

    stats = {}            # (fid, orient) -> per-feature stats
    tss_positions = {}    # (fid, orient) -> {cat: [int, ...]}
    end_positions = {}    # (fid, orient) -> {cat: [int, ...]}

    clip_rows_splice = []
    clip_rows_nonsplice = []
    elem_exon_rows = []
    gene_exon_rows = []

    def _ensure(key):
        if key not in stats:
            stats[key] = {'total': 0,
                          'counts': {c: 0 for c in ALL_CATS},
                          'lengths': {c: 0 for c in ALL_CATS},
                          'spliced': {c: 0 for c in ALL_CATS},
                          'junctions': {c: set() for c in ALL_CATS},
                          'strands': {c: [] for c in ALL_CATS}}
            tss_positions[key] = {c: [] for c in ALL_CATS}
            end_positions[key] = {c: [] for c in ALL_CATS}

    for path in bam_paths:
        bf = pysam.AlignmentFile(path, 'rb')
        try:
            read_iter = bf.fetch(contig=chrom, start=fetch_start, end=fetch_end)
        except ValueError:
            bf.close()
            continue

        for read in read_iter:
            rs = read.reference_start
            if rs is None or rs < claim_start or rs >= claim_end:
                continue

            n_reads += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            strand = infer_read_strand(read, trust_st_tag=trust_st_tag)
            if strand not in ("+", "-"):
                continue
            strand_counts[strand] += 1

            origin = read.reference_start if strand == '+' else read.reference_end - 1
            endpos = read.reference_end - 1 if strand == '+' else read.reference_start

            # Innermost containing feature wins; the rest become nesting parents.
            winner, others = innermost_at(body_tree_chrom, origin)
            if winner is None:
                continue
            rec = feat_info.get(winner)
            if rec is None:
                continue

            klass = rec['class']

            # Orientation. Features with a known strand are labelled directly. Features annotated
            # '.' or '?' are bucketed by read strand and resolved to sense/antisense after the
            # run, by majority read support -- the old code voted on strand but applied the result
            # only after classification had already finished, so the vote never took effect.
            if rec['strand_known']:
                orient = 'sense' if strand == rec['strand'] else 'antisense'
            else:
                orient = 'rs+' if strand == '+' else 'rs-'

            blocks = read.get_blocks()
            juncs = junctions_from_read(read, min_intron_len)
            has_intron = bool(juncs)

            if klass == 'LTR_structural':
                cat = classify_ltr(rec, origin, endpos, strand, blocks, has_intron,
                                   spliced_cov_frac)
            else:
                cat = classify_simple(rec, origin, endpos, strand, blocks, has_intron,
                                      full_length_frac, spliced_cov_frac)
            if not cat:
                continue

            # Feature the 3' end lands in, when it differs from the assigned feature. Lets a
            # transcript that starts in a gene and terminates inside a nested TE (or starts in a
            # TE and runs out through a gene) be picked out for chimerism follow-up.
            term_fid, _ = innermost_at(body_tree_chrom, endpos)
            terminates_in = term_fid if (term_fid and term_fid != winner) else ''

            key = (winner, orient)
            _ensure(key)
            srec = stats[key]
            srec['total'] += 1
            srec['counts'][cat] += 1
            srec['lengths'][cat] += (read.query_length or read.infer_query_length() or 0)
            srec['strands'][cat].append(strand)
            srec['junctions'][cat] |= juncs
            tss_positions[key][cat].append(origin)
            end_positions[key][cat].append(endpos)
            if has_intron:
                srec['spliced'][cat] += 1
            n_assigned += 1

            nested_in = ";".join(others) if others else ""
            clip3S = softclip_3prime(read)

            n_ex, n_in, ex_total, in_total, fields, saw_real = exon_intron_row_fields(
                read, min_intron_len=min_intron_len, max_exons=5, length_mode="query")

            row = (winner, klass, orient, chrom, origin, cat,
                   read.query_name, read.mapping_quality, int(read.is_reverse),
                   clip3S, read.reference_start, read.reference_end,
                   nested_in, terminates_in)
            if saw_real:
                clip_rows_splice.append(row)
                exon_row = (winner, klass, orient, chrom, origin, cat,
                            read.query_name, read.mapping_quality, int(read.is_reverse),
                            n_ex, n_in, ex_total, in_total, *fields,
                            nested_in, terminates_in)
                if klass == 'Gene':
                    gene_exon_rows.append(exon_row)
                else:
                    elem_exon_rows.append(exon_row)
            else:
                clip_rows_nonsplice.append(row)

        bf.close()

    return {
        'chrom': chrom,
        'n_reads': n_reads,
        'n_assigned': n_assigned,
        'strand_counts': strand_counts,
        'stats': stats,
        'tss_positions': tss_positions,
        'end_positions': end_positions,
        'clip_rows_splice': clip_rows_splice,
        'clip_rows_nonsplice': clip_rows_nonsplice,
        'elem_exon_rows': elem_exon_rows,
        'gene_exon_rows': gene_exon_rows,
    }


# ---------------------------------------------------------------------------
# Dispatcher: shards by chunk, merges results, writes per-read files
# ---------------------------------------------------------------------------

def _classify_chunk_star(a):
    """imap_unordered passes a single argument, so unpack the chunk tuple here."""
    return _classify_chunk(*a)


def _fmt_hms(seconds):
    seconds = int(max(0, seconds))
    h, rem = divmod(seconds, 3600)
    m, sec = divmod(rem, 60)
    return f"{h}:{m:02d}:{sec:02d}" if h else f"{m}:{sec:02d}"


def _expected_read_total(bam_paths, contigs):
    """
    Mapped reads on the contigs that will actually be scanned, read from the BAM indexes.

    Used only to turn elapsed time into an ETA. Returns 0 if any index is unusable, in which case
    progress is reported without a projection rather than with a wrong one.
    """
    total = 0
    for path in bam_paths:
        try:
            bf = pysam.AlignmentFile(path, 'rb')
            try:
                for st in bf.get_index_statistics():
                    if st.contig in contigs:
                        total += st.mapped
            finally:
                bf.close()
        except Exception:
            return 0
    return total


def _new_stat_rec():
    return {'total': 0,
            'counts': {c: 0 for c in ALL_CATS},
            'lengths': {c: 0 for c in ALL_CATS},
            'spliced': {c: 0 for c in ALL_CATS},
            'junctions': {c: set() for c in ALL_CATS},
            'strands': {c: [] for c in ALL_CATS}}


def compute_tss_peak(origins, window, min_frac):
    """
    Modal 5' end for a feature, and how much of its read support sits on that peak.

    The window is centred on the peak, so the default width of 3 covers peak-1..peak+1. Reported
    only -- TSS_Called never gates which reads are counted.
    """
    if not origins:
        return ('NA', 0, 0.0, 0)
    ctr = Counter(origins)
    total = len(origins)
    peak, _ = ctr.most_common(1)[0]
    half = max(0, (window - 1) // 2)
    cnt = sum(v for pos, v in ctr.items() if abs(pos - peak) <= half)
    frac = cnt / total
    return (peak, cnt, frac, int(frac >= min_frac))


def nesting_parents(fid, feat_info, body_tree):
    """Features whose span strictly contains this one, using the same rank tiebreak as assignment."""
    rec = feat_info[fid]
    tree = body_tree.get(rec['chrom'])
    if tree is None:
        return []
    span = rec['end'] - rec['start']
    rank = CLASS_RANK.get(rec['class'], 9)
    out = []
    for iv in tree.overlap(rec['start'], rec['end']):
        other, orank = iv.data
        if other == fid:
            continue
        if iv.begin <= rec['start'] and iv.end >= rec['end']:
            ospan = iv.end - iv.begin
            if ospan > span or (ospan == span and orank < rank):
                out.append(other)
    return sorted(out)


def classify_multiple_bams(feat_info, body_tree, bam_paths, min_mapq, args,
                           clip_out_splice=None, clip_out_nonsplice=None,
                           elem_exon_out=None, gene_exon_out=None,
                           n_threads=1, chunk_size=1_000_000, overlap=10_000):
    # Discover chromosome sizes from BAM headers (union across all BAMs, max length wins).
    chrom_sizes = {}
    for p in bam_paths:
        bf = pysam.AlignmentFile(p, 'rb')
        for ref, length in zip(bf.references, bf.lengths):
            chrom_sizes[ref] = max(chrom_sizes.get(ref, 0), length)
        bf.close()

    chunks = []
    for chrom in sorted(chrom_sizes):
        length = chrom_sizes[chrom]
        if length <= 0 or chrom not in body_tree:
            continue
        pos = 0
        while pos < length:
            claim_start = pos
            claim_end = min(pos + chunk_size, length)
            chunks.append((chrom, claim_start, claim_end,
                           max(0, claim_start - overlap), min(length, claim_end + overlap),
                           bam_paths, min_mapq, args.trust_st_tag, args.min_intron_len,
                           args.full_length_frac, args.spliced_cov_frac))
            pos = claim_end

    print(f"Dispatching {len(chunks)} chunks ({chunk_size // 1000} kb each, "
          f"{overlap // 1000} kb overlap) across {n_threads} process(es)...")

    start_time = time.time()

    global _GLOBAL_FEAT_INFO, _GLOBAL_BODY_TREE
    _GLOBAL_FEAT_INFO = feat_info
    _GLOBAL_BODY_TREE = body_tree

    # Progress is reported on reads assessed, aggregated across workers, rather than one line per
    # chunk -- 1 Mb chunks meant thousands of lines per run.
    step = max(1, args.progress_every)
    expected = _expected_read_total(bam_paths, {c[0] for c in chunks})
    results = []
    seen = 0
    assigned = 0
    next_mark = step

    def _tick(res):
        nonlocal seen, assigned, next_mark
        results.append(res)
        seen += res['n_reads']
        assigned += res['n_assigned']
        if seen < next_mark:
            return
        elapsed = time.time() - start_time
        msg = f"  {seen:,} reads assessed, {assigned:,} assigned | {_fmt_hms(elapsed)} elapsed"
        if expected and seen < expected and elapsed > 0:
            msg += (f" | {100.0 * seen / expected:.1f}% of {expected:,}"
                    f" | ETA {_fmt_hms(elapsed * (expected - seen) / seen)}")
        print(msg, flush=True)
        next_mark = (seen // step + 1) * step

    if n_threads > 1:
        with Pool(processes=n_threads, initializer=init_worker,
                  initargs=(feat_info, body_tree)) as pool:
            for res in pool.imap_unordered(_classify_chunk_star, chunks, chunksize=1):
                _tick(res)
    else:
        for a in chunks:
            _tick(_classify_chunk(*a))

    # Merge. Accumulators are sparse -- only features that actually received reads get a record.
    # Pre-allocating every feature was affordable at 52k structural elements but is not at the
    # ~1.1M features the full TE annotation brings in.
    strand_counts = Counter()
    stats = {}
    tss_positions = {}
    end_positions = {}
    clip_rows_splice = []
    clip_rows_nonsplice = []
    elem_exon_rows = []
    gene_exon_rows = []

    for res in results:
        strand_counts += res['strand_counts']

        for key, r in res['stats'].items():
            s = stats.get(key)
            if s is None:
                s = stats[key] = _new_stat_rec()
                tss_positions[key] = {c: [] for c in ALL_CATS}
                end_positions[key] = {c: [] for c in ALL_CATS}
            s['total'] += r['total']
            for c in ALL_CATS:
                if r['counts'][c]:
                    s['counts'][c] += r['counts'][c]
                    s['lengths'][c] += r['lengths'][c]
                    s['spliced'][c] += r['spliced'][c]
                    s['strands'][c].extend(r['strands'][c])
                if r['junctions'][c]:
                    s['junctions'][c] |= r['junctions'][c]

        for key, pos_dict in res['tss_positions'].items():
            dst = tss_positions[key]
            for c, lst in pos_dict.items():
                if lst:
                    dst[c].extend(lst)
        for key, pos_dict in res['end_positions'].items():
            dst = end_positions[key]
            for c, lst in pos_dict.items():
                if lst:
                    dst[c].extend(lst)

        clip_rows_splice.extend(res['clip_rows_splice'])
        clip_rows_nonsplice.extend(res['clip_rows_nonsplice'])
        elem_exon_rows.extend(res['elem_exon_rows'])
        gene_exon_rows.extend(res['gene_exon_rows'])

    # Resolve orientation for features annotated with an unknown strand. Reads were bucketed by
    # their own strand; the better-supported bucket defines the element strand, and therefore
    # which bucket is sense.
    unknown = defaultdict(dict)
    for (fid, orient) in list(stats):
        if orient in ('rs+', 'rs-'):
            unknown[fid][orient] = stats[(fid, orient)]['total']

    n_inferred = 0
    for fid, buckets in unknown.items():
        plus = buckets.get('rs+', 0)
        minus = buckets.get('rs-', 0)
        inferred = '+' if plus >= minus else '-'
        feat_info[fid]['strand'] = inferred
        feat_info[fid]['strand_source'] = 'inferred'
        sense_key = 'rs+' if inferred == '+' else 'rs-'
        for orient in ('rs+', 'rs-'):
            old = (fid, orient)
            if old not in stats:
                continue
            new = (fid, 'sense' if orient == sense_key else 'antisense')
            stats[new] = stats.pop(old)
            tss_positions[new] = tss_positions.pop(old)
            end_positions[new] = end_positions.pop(old)
        n_inferred += 1

    print(f"Classification complete in {(time.time() - start_time) / 60:.1f} minutes.")
    print(f"Strand counts across all reads: {dict(strand_counts)}")
    print(f"Features with reads: {len({f for f, _ in stats})} "
          f"({n_inferred} had their strand inferred from read support)")
    obs = Counter(feat_info[f]['class'] for f in {f for f, _ in stats})
    print("  by class: " + ", ".join(f"{k}={v}" for k, v in sorted(obs.items())))
    orient_totals = Counter()
    for (fid, orient), rec in stats.items():
        orient_totals[orient] += rec['total']
    print(f"  reads by orientation: {dict(orient_totals)}")

    # ---- per-read outputs ----
    read_hdr = ["Feature", "Class", "Orientation", "Chrom", "TSS", "Category",
                "Read", "MAPQ", "Aln_Reverse", "softclip_3p", "aln_start", "aln_end",
                "Nested_In", "Terminates_In"]
    exon_hdr = ["Feature", "Class", "Orientation", "Chrom", "TSS", "Category",
                "Read", "MAPQ", "Aln_Reverse",
                "exon_count", "intron_count", "exon_len_combined", "intron_len_combined",
                "exon1", "intron1", "exon2", "intron2", "exon3", "intron3",
                "exon4", "intron4", "exon5", "Nested_In", "Terminates_In"]

    for path, rows, hdr, label in (
            (clip_out_splice, clip_rows_splice, read_hdr, "spliced reads"),
            (clip_out_nonsplice, clip_rows_nonsplice, read_hdr, "non-spliced reads"),
            (elem_exon_out, elem_exon_rows, exon_hdr, "TE exon statistics"),
            (gene_exon_out, gene_exon_rows, exon_hdr, "gene exon statistics")):
        if not path:
            continue
        with open(path, "w") as out:
            out.write("\t".join(hdr) + "\n")
            for row in rows:
                out.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote {len(rows)} {label} records to {path}")

    return stats, tss_positions, end_positions


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _report_cats(include_partial):
    return [c for c in ALL_CATS if include_partial or c != 'partial']


def _pooled_origins(tss_positions, key):
    out = []
    for lst in tss_positions[key].values():
        out.extend(lst)
    return out


def _read_frame_strand(rec, orient):
    """Transcription strand of the reads in a record: flipped from the element for antisense."""
    s = rec.get('strand', '.')
    if orient == 'antisense':
        return '-' if s == '+' else ('+' if s == '-' else '.')
    return s


def write_isoform_tables(feat_info, body_tree, stats, tss_positions, prefix, args):
    """
    Write per-feature isoform tables, split into sense and antisense files.

    Every feature reports against the union vocabulary; categories outside a feature's own
    vocabulary stay at zero (a gene never has an ltr5_contained read, a fragment never has a
    spanning one). Nesting parents are computed only for features that received reads.
    """
    cats = _report_cats(args.include_partial)
    parent_cache = {}
    rows_by_orient = defaultdict(list)

    for (fid, orient), rec in stats.items():
        info = feat_info[fid]
        total = rec['total']
        if total == 0:
            continue

        if fid not in parent_cache:
            parents = nesting_parents(fid, feat_info, body_tree)
            parent_cache[fid] = (parents,
                                 [feat_info[p]['class'] for p in parents])
        parents, parent_classes = parent_cache[fid]

        peak, _peak_cnt, peak_frac, called = compute_tss_peak(
            _pooled_origins(tss_positions, (fid, orient)),
            args.tss_window, args.tss_min_frac)

        row = {
            'Feature': fid,
            'Chrom': info['chrom'],
            'Start': info['start'],
            'End': info['end'],
            'Class': info['class'],
            'Name': info['name'],
            'Strand': info['strand'],
            'Strand_Source': info.get('strand_source', 'annotated'),
            'Orientation': orient,
            'Total_Reads': total,
            'Nested_In': ";".join(parents) if parents else "NA",
            'Nested_In_Class': ";".join(parent_classes) if parents else "NA",
            'TSS_Called': called,
            'TSS_Peak': peak,
            'TSS_Peak_Frac': f"{peak_frac:.4f}",
        }
        for cat in cats:
            count = rec['counts'][cat]
            row[f'{cat}_reads'] = count
            row[f'mean_len_{cat}'] = f"{rec['lengths'][cat] / count:.1f}" if count else "0.0"
            # 'n_spliced_' rather than 'spliced_': several categories are themselves named
            # spliced_*, and 'spliced_spanning' as a counter prefix collided with the
            # spliced_spanning category's own columns.
            if cat not in SPLICED_BY_DEF and cat not in UNSPLICED_BY_DEF:
                row[f'n_spliced_{cat}'] = rec['spliced'][cat]
            if cat not in UNSPLICED_BY_DEF:
                row[f'unique_juncts_{cat}'] = len(rec['junctions'][cat])
        row['Attrs'] = info['attrs']
        rows_by_orient[orient].append(row)

    for orient in ('sense', 'antisense'):
        out_path = f"{prefix}_{orient}.isoforms.tsv"
        rows = rows_by_orient.get(orient, [])
        if rows:
            df = pd.DataFrame(rows).sort_values('Total_Reads', ascending=False)
        else:
            df = pd.DataFrame(columns=['Feature', 'Chrom', 'Start', 'End', 'Class', 'Orientation',
                                       'Total_Reads'])
        df.to_csv(out_path, sep='\t', index=False)
        print(f"Wrote {len(rows)} {orient} isoform records to {out_path}")


def write_isoform_tss_summary(feat_info, stats, tss_positions, out_path, args, n=2,
                              orient_filter=None):
    """
    Top isoform and top-n TSS positions per feature and orientation.

    orient_filter restricts output to one orientation. Downstream tools that key a lookup on
    Feature (WindowScrubber builds tss_map[Feature] from TSS1) need one row per feature, so the
    orientation-split copies rather than the combined table are what should be fed to them.
    """
    cats = _report_cats(args.include_partial)
    with open(out_path, 'w') as out:
        hdr = ["Feature", "Class", "Orientation", "Top_Isoform", "Strand", "Total_Reads",
               "TSS_Called", "TSS_Peak_Frac"]
        for i in range(1, n + 1):
            hdr += [f"TSS{i}", f"Count{i}"]
        out.write("\t".join(hdr) + "\n")

        for (fid, orient), rec in stats.items():
            total = rec['total']
            if total == 0 or (orient_filter and orient != orient_filter):
                continue
            info = feat_info[fid]
            top_cat = max(cats, key=lambda c: rec['counts'][c])
            if rec['counts'][top_cat] == 0:
                continue

            origins = _pooled_origins(tss_positions, (fid, orient))
            peak, peak_cnt, peak_frac, called = compute_tss_peak(
                origins, args.tss_window, args.tss_min_frac)

            # Pooled across categories, matching TSS_Peak in the isoform table. Ranking only
            # the top category made TSS1 here a different quantity from TSS_Peak there.
            common = Counter(origins).most_common(n)
            common += [("NA", 0)] * (n - len(common))

            row = [fid, info['class'], orient, top_cat, _read_frame_strand(info, orient),
                   str(total), str(called), f"{peak_frac:.4f}"]
            for tss, cnt in common:
                row += [str(tss), str(cnt)]
            out.write("\t".join(row) + "\n")
    print(f"Wrote isoform TSS summary to {out_path}")


def write_cleavage_summary(feat_info, stats, end_positions, out_path, args, n=10, window=5):
    """3'-end (cleavage) distribution per feature and orientation, pooled across categories."""
    with open(out_path, "w") as out:
        hdr = ["Feature", "Class", "Orientation", "Strand", "Total_Reads",
               "CleavageSitePeak", "CleavageSiteCnt", "CleavageSiteFrac",
               f"PeakWindow{window}_Cnt", f"PeakWindow{window}_Frac"]
        for i in range(1, n + 1):
            hdr += [f"End{i}", f"Count{i}"]
        out.write("\t".join(hdr) + "\n")

        for (fid, orient), rec in stats.items():
            all_ends = []
            for cat, cnt in rec["counts"].items():
                if cnt:
                    all_ends.extend(end_positions[(fid, orient)][cat])
            if not all_ends:
                continue

            ctr = Counter(all_ends)
            total = sum(ctr.values())
            common = ctr.most_common()
            top_n = common[:n] + [("NA", 0)] * (n - len(common[:n]))

            peak_site, peak_cnt = common[0]
            peak_window_cnt = sum(v for pos, v in ctr.items() if abs(pos - peak_site) <= window)

            info = feat_info[fid]
            row = [fid, info['class'], orient, _read_frame_strand(info, orient), str(total),
                   str(peak_site), str(peak_cnt), f"{peak_cnt / total:.4f}",
                   str(peak_window_cnt), f"{peak_window_cnt / total:.4f}"]
            for end, cnt in top_n:
                row += [str(end), str(cnt)]
            out.write("\t".join(row) + "\n")
    print(f"Wrote top-{n} cleavage summary to {out_path}")


def write_gene_summary(feat_info, stats, tss_positions, out_path, args, n=2):
    """Gene read counts and top-n TSS positions, one row per gene and orientation."""
    with open(out_path, 'w') as out:
        hdr = ["Gene", "Orientation", "Total_Reads", "TSS_Called", "TSS_Peak_Frac"]
        for i in range(1, n + 1):
            hdr += [f"TSS{i}", f"Count{i}"]
        out.write("\t".join(hdr) + "\n")

        n_rows = 0
        for (fid, orient), rec in stats.items():
            if feat_info[fid]['class'] != 'Gene' or rec['total'] == 0:
                continue
            origins = _pooled_origins(tss_positions, (fid, orient))
            _, _, peak_frac, called = compute_tss_peak(origins, args.tss_window, args.tss_min_frac)
            common = Counter(origins).most_common(n)
            common += [("NA", 0)] * (n - len(common))
            row = [fid, orient, str(rec['total']), str(called), f"{peak_frac:.4f}"]
            for tss, cnt in common:
                row += [str(tss), str(cnt)]
            out.write("\t".join(row) + "\n")
            n_rows += 1
    print(f"Wrote {n_rows} gene summary records to {out_path}")


def compute_tss_density(stats, tss_positions, feat_info, window=10, min_reads=1,
                        include_tss_reads=True, klass=None):
    """
    Per feature and orientation, histograms of read 5' ends around the primary and secondary TSS.

    Distances are measured in the read's frame, so positive is always downstream of the TSS in the
    direction of transcription -- for antisense records that is the opposite genomic direction from
    the element's annotated strand.
    """
    density = {}
    size = 2 * window + 1

    for (fid, orient), rec in stats.items():
        if rec.get('total', 0) < min_reads:
            continue
        info = feat_info[fid]
        if klass is not None and info['class'] != klass:
            continue
        if klass is None and info['class'] == 'Gene':
            continue

        origins = _pooled_origins(tss_positions, (fid, orient))
        if not origins:
            continue
        common = Counter(origins).most_common(2)
        while len(common) < 2:
            common.append((None, 0))
        (tss1, _), (tss2, _) = common
        if tss1 is None:
            continue

        strand = _read_frame_strand(info, orient)
        hist1 = np.zeros(size, dtype=int)
        hist2 = np.zeros(size, dtype=int)

        for target, hist in ((tss1, hist1), (tss2, hist2)):
            if target is None:
                continue
            for o in origins:
                if not include_tss_reads and o == target:
                    continue
                d = o - target
                if strand == '-':
                    d = -d
                if -window <= d <= window:
                    hist[d + window] += 1

        density[(fid, orient)] = {'primary': hist1, 'secondary': hist2}

    return density


def write_density(density, out_path, which, window=10, label="Feature"):
    with open(out_path, 'w') as out:
        out.write(f"{label}\tOrientation\tDistance\tCount\n")
        for (fid, orient), hists in density.items():
            for i, cnt in enumerate(hists[which]):
                out.write(f"{fid}\t{orient}\t{i - window}\t{cnt}\n")
    print(f"Wrote {which} TSS densities to {out_path}")


# ---------------------------------------------------------------------------
# U3 / Promoter Sequence Extraction
# ---------------------------------------------------------------------------

def _extract_region(fa, chrom, start1, end1, strand):
    """Extract 1-based inclusive [start1, end1] and reverse-complement if strand == '-'."""
    seq = Seq(fa.fetch(chrom, start1 - 1, end1))
    return seq.reverse_complement() if strand == '-' else seq


def u3_seq_extraction(genome_fasta, feat_info, stats, tss_positions, output_prefix, args):
    """
    Extract U3/LTR sequences for structural LTR-RTs and promoter regions for genes.

    U3 extraction is restricted to LTR_structural: fragments have no resolved LTR pair and DNA
    transposons have no U3 at all. Both orientations are emitted so antisense initiation can be
    compared against the canonical U3-driven TSS downstream. The upstream LTR boundary and the
    reverse-complement are both taken in the read's frame, so an antisense record measures back
    to the edge of the LTR the antisense transcript actually initiated in.
    """
    print("Opening genome FASTA for U3/promoter extraction...")
    fa = pysam.FastaFile(genome_fasta)
    contigs = set(fa.references)

    records = {('u3', 'sense'): [], ('u3', 'antisense'): [],
               ('ltr', 'sense'): [], ('ltr', 'antisense'): []}
    n_skipped = Counter()

    for (fid, orient), rec in stats.items():
        if rec.get('total', 0) == 0:
            continue
        info = feat_info[fid]
        if info['class'] != 'LTR_structural':
            continue

        chrom = info['chrom']
        if chrom not in contigs:
            n_skipped['contig missing from genome'] += 1
            continue

        tstrand = _read_frame_strand(info, orient)
        if tstrand not in ('+', '-'):
            n_skipped['unresolved strand'] += 1
            continue

        origins = _pooled_origins(tss_positions, (fid, orient))
        peak, _, _, _ = compute_tss_peak(origins, args.tss_window, args.tss_min_frac)
        if peak == 'NA':
            continue
        tss0 = peak
        tss1 = tss0 + 1

        ltr_left, ltr_right = info['ltr_left'], info['ltr_right']
        if not ltr_left or not ltr_right:
            continue

        if ltr_left[0] <= tss0 < ltr_left[1]:
            chosen = ltr_left
        elif ltr_right[0] <= tss0 < ltr_right[1]:
            chosen = ltr_right
        else:
            n_skipped['TSS outside both LTRs'] += 1
            continue

        # Upstream LTR edge in the read's frame
        boundary1 = (chosen[0] + 1) if tstrand == '+' else chosen[1]
        u3_start1, u3_end1 = sorted((tss1, boundary1))

        u3_seq = _extract_region(fa, chrom, u3_start1, u3_end1, tstrand)
        records[('u3', orient)].append(SeqRecord(
            u3_seq, id=f"{fid}|{chrom}:{u3_start1}-{u3_end1}({tstrand})|{orient}", description=""))

        full_start1, full_end1 = chosen[0] + 1, chosen[1]
        records[('ltr', orient)].append(SeqRecord(
            _extract_region(fa, chrom, full_start1, full_end1, tstrand),
            id=f"{fid}|{chrom}:{full_start1}-{full_end1}({tstrand})|{orient}", description=""))

    for (kind, orient), recs in records.items():
        path = f"{output_prefix}_{kind}_seqs_{orient}.fa"
        SeqIO.write(recs, path, 'fasta')
        print(f"[U3-LTR] Wrote {len(recs)} {kind} sequences to {path}")
    if n_skipped:
        print("  skipped: " + ", ".join(f"{k}={v}" for k, v in n_skipped.items()))

    # ---- Gene promoters (sense only; the antisense 5' end is not the annotated promoter) ----
    promoter_records = []
    gene_u3_records = []
    bed_path = output_prefix + '_gene_2kb_proms.bed'
    prom_fa_path = output_prefix + '_gene_2kb_proms.fa'
    gene_u3_path = output_prefix + '_gene_dummy_u3.fa'

    with open(bed_path, 'w') as bed_out:
        for (fid, orient), rec in stats.items():
            if orient != 'sense' or rec.get('total', 0) == 0:
                continue
            info = feat_info[fid]
            if info['class'] != 'Gene':
                continue
            chrom = info['chrom']
            strand = info['strand']
            if chrom not in contigs or strand not in ('+', '-'):
                continue

            origins = _pooled_origins(tss_positions, (fid, orient))
            peak, _, _, _ = compute_tss_peak(origins, args.tss_window, args.tss_min_frac)
            if peak == 'NA':
                continue
            tss1 = peak + 1
            chrom_len = fa.get_reference_length(chrom)

            prom_start0 = max(0, tss1 - 1000)
            prom_end1 = min(tss1 + 1000, chrom_len)
            bed_out.write(f"{chrom}\t{prom_start0}\t{prom_end1}\t{fid}\t0\t{strand}\n")
            promoter_records.append(SeqRecord(
                _extract_region(fa, chrom, prom_start0 + 1, prom_end1, strand),
                id=f"{fid}|{chrom}:{prom_start0}-{prom_end1}({strand})", description=""))

            if strand == '+':
                u3_start0, u3_end1_g = max(0, tss1 - 1000), tss1
            else:
                u3_start0, u3_end1_g = tss1, min(tss1 + 1000, chrom_len)
            gene_u3_records.append(SeqRecord(
                Seq("NNNNN"), id=f"{fid}|{chrom}:{u3_start0}-{u3_end1_g}({strand})", description=""))

    SeqIO.write(promoter_records, prom_fa_path, 'fasta')
    print(f"[Gene] Wrote {len(promoter_records)} promoter entries to {bed_path} and {prom_fa_path}")
    SeqIO.write(gene_u3_records, gene_u3_path, 'fasta')
    print(f"[Gene] Wrote {len(gene_u3_records)} dummy U3 sequences to {gene_u3_path}")
    fa.close()


def main():
    args = parse_args()

    clip_out_splice = args.output + "_3p_softclip_per_read_spliced.tsv"
    clip_out_nonsplice = args.output + "_3p_softclip_per_read_nonspliced.tsv"
    elem_exon_out = args.output + "_te_exon_stats_per_read.tsv"
    gene_exon_out = args.output + "_gene_exon_stats_per_read.tsv"

    feat_info, body_tree = load_elements_and_ranges(
        args.gff, args.te_classes, args.gene_gff, args.canonical_exons_only)

    stats, tss_positions, end_positions = classify_multiple_bams(
        feat_info, body_tree, args.bam, args.min_mapq, args,
        clip_out_splice=clip_out_splice,
        clip_out_nonsplice=clip_out_nonsplice,
        elem_exon_out=elem_exon_out,
        gene_exon_out=gene_exon_out,
        n_threads=args.threads)

    # Isoform tables, split sense/antisense
    write_isoform_tables(feat_info, body_tree, stats, tss_positions, args.output, args)

    # TSS and cleavage summaries (Orientation column)
    write_isoform_tss_summary(feat_info, stats, tss_positions, args.tss_out, args, n=2)
    write_isoform_tss_summary(feat_info, stats, tss_positions,
                              args.output + "_10site.tss_summary.tsv", args, n=10)
    # Single-orientation copies: WindowScrubber keys tss_map on Feature, so it needs exactly one
    # row per feature. Pair these with the matching _u3_seqs_<orient>.fa / _ltr_seqs_<orient>.fa.
    for _o in ('sense', 'antisense'):
        write_isoform_tss_summary(feat_info, stats, tss_positions,
                                  f"{args.output}_{_o}.tss_summary.tsv", args, n=10,
                                  orient_filter=_o)
    write_cleavage_summary(feat_info, stats, end_positions,
                           args.output + "_10site.cleavage_summary.tsv", args, n=10)

    # Gene summaries (Orientation column)
    write_gene_summary(feat_info, stats, tss_positions, args.gene_out, args, n=2)
    write_gene_summary(feat_info, stats, tss_positions,
                       args.output + "_10site.gene_summary.tsv", args, n=10)

    # TSS densities
    te_dens = compute_tss_density(stats, tss_positions, feat_info, window=10, min_reads=6)
    write_density(te_dens, args.output + '_primary_tss_density.tsv', 'primary')
    write_density(te_dens, args.output + '_secondary_tss_density.tsv', 'secondary')

    gene_dens = compute_tss_density(stats, tss_positions, feat_info, window=10, min_reads=7,
                                    klass='Gene')
    write_density(gene_dens, args.output + '_gene_primary_density.tsv', 'primary', label="Gene")
    write_density(gene_dens, args.output + '_gene_secondary_density.tsv', 'secondary', label="Gene")

    if args.genome_fasta:
        u3_seq_extraction(args.genome_fasta, feat_info, stats, tss_positions, args.output, args)


if __name__ == '__main__':
    main()
