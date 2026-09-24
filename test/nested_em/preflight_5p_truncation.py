#!/usr/bin/env python3
"""
Pre-flight check: how often do reads in this library miss their 5' end?

A read is "5'-truncated" when its first base lies downstream of where the RNA actually started.
Reverse transcription starts at the poly(A) tail and copies toward the 5' end; if it falls off, or
the RNA was already broken, the read begins part-way along the molecule. Strand-switching kits
(e.g. ONT PCB114) add the switch oligo wherever RT stopped, so these reads still pass PyChopper as
"full-length". This is different from premature termination, which shortens the 3' end and never
moves the 5' end.

5'-truncated reads are what let a gene's (or a host TE's) reads land their first base inside a
nested TE and be credited to it as a TSS. This script measures the rate on genes that contain or
overlap no other annotated feature, so no nesting can confuse the answer:

  - for each such gene, the reference TSS is the gene's dominant 5' end within its first exon (or
    the first --search-bp), because annotated TSSs are often tens of bp off; if there is no clear
    peak the annotated TSS is used;
  - sense reads whose 5' end lies in the gene are classed as starting at the TSS (within --offset
    bp along the transcript), downstream in an exon, or in an intron;
  - results are reported per BAM (libraries differ) and by transcript length.

How to read the result:
  - median per-gene 5'-truncated fraction < ~5%: the gene-into-TE problem is minor for you; EM
    mainly matters for TE-in-TE nesting.
  - 20-50% (typical for ONT cDNA): contested reads are common; use --assign em.
  - fraction rising with transcript length: truncation is RT-processivity-like, which is exactly
    what the EM's read-length survival model assumes.
  - flat across lengths: truncation is mostly degradation/fragmentation, or annotated TSSs are off.
    Genes whose reference TSS came from the annotation (Ref_TSS_Source = annotated) are the
    least reliable; check a few by eye in a browser.

Uses IsoClassifier's own GFF loader and strand logic, so run it with the same --gff / --gene-gff.
"""
import argparse
import bisect
import os
import statistics
import sys
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from IsoClassifier import (DEFAULT_TE_CLASSES, infer_read_strand,  # noqa: E402
                           load_elements_and_ranges)

LEN_BINS = [(0, 1000, '<1 kb'), (1000, 2000, '1-2 kb'), (2000, 4000, '2-4 kb'),
            (4000, 8000, '4-8 kb'), (8000, 10 ** 12, '>8 kb')]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--gff', required=True, help='Same --gff as IsoClassifier')
    ap.add_argument('--gene-gff', dest='gene_gff', default=None,
                    help='Same --gene-gff as IsoClassifier (exons make distances transcript-based)')
    ap.add_argument('--bam', nargs='+', required=True)
    ap.add_argument('--min_mapq', type=int, default=30)
    ap.add_argument('--offset', type=int, default=50,
                    help='A read starting more than this many transcript bases downstream of the '
                         'annotated TSS counts as 5\'-truncated (default 50)')
    ap.add_argument('--search-bp', dest='search_bp', type=int, default=300,
                    help='Window at the gene 5\' end (plus the first exon) searched for the '
                         'dominant 5\' end used as reference TSS (default 300)')
    ap.add_argument('--tss-window', dest='tss_window', type=int, default=10,
                    help='+/- bp pooled when finding the dominant 5\' end (default 10)')
    ap.add_argument('--min-reads', dest='min_reads', type=int, default=10,
                    help='Genes need at least this many sense reads per BAM to be scored (default 10)')
    ap.add_argument('--te-classes', dest='te_classes', default=DEFAULT_TE_CLASSES)
    ap.add_argument('--trust-st-tag', dest='trust_st_tag', action='store_true')
    ap.add_argument('--canonical-exons-only', dest='canonical_exons_only', action='store_true')
    ap.add_argument('--out', default='preflight_5p_truncation.tsv', help='Per-gene output table')
    return ap.parse_args()


def isolated_genes(feat_info, body_tree):
    """Genes with a known strand whose span overlaps no other annotated feature."""
    out = []
    for fid, rec in feat_info.items():
        if rec['class'] != 'Gene' or rec['strand'] not in ('+', '-'):
            continue
        tree = body_tree.get(rec['chrom'])
        if tree is not None and len(tree.overlap(rec['start'], rec['end'])) == 1:
            out.append(fid)
    return out


def reference_tss(rec, origins, a):
    """Dominant 5' end in the gene's promoter window (first exon + search_bp), else annotated."""
    ann = rec['start'] if rec['strand'] == '+' else rec['end'] - 1
    lo, hi = (rec['start'], rec['start'] + a.search_bp) if rec['strand'] == '+' else \
             (rec['end'] - a.search_bp, rec['end'])
    exons = rec.get('exons')
    if exons:
        first = exons[0] if rec['strand'] == '+' else exons[-1]
        lo, hi = min(lo, first[0]), max(hi, first[1])
    pos = sorted(o for o in origins if lo <= o < hi)
    if not pos:
        return ann, 'annotated'
    w = a.tss_window
    best, best_n = None, 0
    for x in sorted(set(pos)):
        n = bisect.bisect_right(pos, x + w) - bisect.bisect_left(pos, x - w)
        if n > best_n:
            best, best_n = x, n
    # a peak needs real support: at least 3 reads and 10% of the gene's reads
    if best_n >= 3 and best_n >= 0.10 * len(origins):
        return best, 'observed'
    return ann, 'annotated'


def classify_start(rec, origin, tss, offset):
    """'tss', 'exon' (downstream, exonic) or 'intron' for a sense read's 5' end."""
    downstream = origin > tss if rec['strand'] == '+' else origin < tss
    exons = rec.get('exons')
    if exons:
        starts = [s for s, _ in exons]
        i = bisect.bisect_right(starts, origin) - 1
        if i < 0 or origin >= exons[i][1]:
            return 'intron'
    if not downstream:
        return 'tss'           # at or upstream of the reference TSS: not truncated
    lo, hi = min(tss, origin), max(tss, origin) + 1
    if exons:
        dist = sum(max(0, min(e, hi) - max(s, lo)) for s, e in exons) - 1
    else:
        dist = hi - lo - 1
    return 'tss' if dist <= offset else 'exon'


def main():
    a = parse_args()
    feat_info, body_tree = load_elements_and_ranges(a.gff, a.te_classes, a.gene_gff,
                                                    a.canonical_exons_only)
    genes = isolated_genes(feat_info, body_tree)
    n_genes = sum(1 for r in feat_info.values() if r['class'] == 'Gene')
    print(f"Pre-flight: {len(genes):,} of {n_genes:,} genes overlap no other feature and are used.")
    if not any(feat_info[g].get('exons') for g in genes):
        print("  (no exon models attached: distances are genomic, and intronic starts cannot be "
              "separated; pass --gene-gff for a sharper estimate)")

    base = [os.path.basename(b) for b in a.bam]
    labels = [b if base.count(os.path.basename(b)) == 1 else p for b, p in zip(base, a.bam)]
    rows = []
    for bam_path, sample in zip(a.bam, labels):
        bf = pysam.AlignmentFile(bam_path, 'rb')
        contigs = set(bf.references)
        for fid in genes:
            rec = feat_info[fid]
            if rec['chrom'] not in contigs:
                continue
            origins = []
            for read in bf.fetch(rec['chrom'], rec['start'], rec['end']):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < a.min_mapq:
                    continue
                strand = infer_read_strand(read, trust_st_tag=a.trust_st_tag)
                if strand != rec['strand']:
                    continue
                origin = read.reference_start if strand == '+' else read.reference_end - 1
                if rec['start'] <= origin < rec['end']:
                    origins.append(origin)
            if len(origins) < a.min_reads:
                continue
            tss, src = reference_tss(rec, origins, a)
            counts = defaultdict(int)
            for o in origins:
                counts[classify_start(rec, o, tss, a.offset)] += 1
            tx_len = rec.get('exon_len') or (rec['end'] - rec['start'])
            rows.append((sample, fid, tx_len, len(origins), counts['tss'], counts['exon'],
                         counts['intron'], tss, src))
        bf.close()

    with open(a.out, 'w') as out:
        out.write("Sample\tGene\tTx_Length\tReads\tStart_At_TSS\tStart_Downstream_Exon\t"
                  "Start_Intron\tFrac_5p_Truncated\tRef_TSS\tRef_TSS_Source\n")
        for s, g, L, n, t, e, i, tss, src in rows:
            out.write(f"{s}\t{g}\t{L}\t{n}\t{t}\t{e}\t{i}\t{(e + i) / n:.4f}\t{tss}\t{src}\n")
    print(f"Wrote {len(rows):,} gene x sample rows to {a.out}")

    for sample in dict.fromkeys(r[0] for r in rows):
        sr = [r for r in rows if r[0] == sample]
        reads = sum(r[3] for r in sr)
        trunc = sum(r[5] + r[6] for r in sr)
        intr = sum(r[6] for r in sr)
        med = statistics.median((r[5] + r[6]) / r[3] for r in sr)
        print(f"\n== {sample}: {len(sr):,} genes, {reads:,} sense reads ==")
        n_ann = sum(1 for r in sr if r[8] == 'annotated')
        print(f"  reference TSS observed for {len(sr) - n_ann:,} genes, annotated for {n_ann:,}")
        print(f"  5'-truncated (start > {a.offset} bp downstream of the reference TSS): "
              f"{trunc / reads:.1%} of reads; median per gene {med:.1%}")
        print(f"  of which start in an intron (pre-mRNA / retained intron): {intr / reads:.1%} of reads")
        by_bin = []
        for lo, hi, label in LEN_BINS:
            b = [r for r in sr if lo <= r[2] < hi]
            if not b:
                continue
            fr = sum(r[5] + r[6] for r in b) / sum(r[3] for r in b)
            by_bin.append(fr)
            print(f"    transcripts {label:>7s}: {fr:6.1%} truncated  ({len(b):,} genes)")
        if med < 0.05:
            verdict = "minor: nesting errors from 5'-truncated reads will be rare"
        elif med < 0.20:
            verdict = "moderate: --assign em worthwhile around nested features"
        else:
            verdict = "substantial: use --assign em; innermost-rule counts in nested loci are unreliable"
        print(f"  verdict: {verdict}")
        if len(by_bin) >= 2 and by_bin[0] > 0:
            trend = by_bin[-1] / by_bin[0]
            kind = ("rises with length (RT-processivity-like; matches the EM survival model)"
                    if trend >= 1.5 else
                    "roughly flat with length (degradation/fragmentation, or annotated TSSs are off)")
            print(f"  length trend: longest/shortest bin = {trend:.2f}x, {kind}")


if __name__ == '__main__':
    main()
