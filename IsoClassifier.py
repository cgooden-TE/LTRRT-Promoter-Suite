#!/usr/bin/env python3
"""
Classifier for LTR retrotransposon isoforms and gene read/TSS summary across multiple replicates
with nested-element reporting using intervaltree

Inputs: 
    - GFF with LTR elements, genes, and nested features
    - One or more BAM files with mapped reads
    - Minimum MAPQ for filtering reads

Outputs:
    - TSV of LTR isoforms with read counts, lengths, splicing, and junctions
    - TSS summary for LTR isoforms
    - TSV of gene read counts and TSS summary (total reads + top TSS positions)

Example Usage:
    python3 IsoClassifier.py \
        --gff annotations.gff \
        --bam sample1.bam sample2.bam \
        --min_mapq 30 \
        --output ltr_isoforms \
        --tss_out ltr_tss_summary.tsv \
        --gene_out gene_summary.tsv
"""
import argparse
import time
from collections import Counter, defaultdict
import pandas as pd
import pysam
import numpy as np
from intervaltree import IntervalTree

def parse_args():
    parser = argparse.ArgumentParser(
        description='Classify LTR isoforms and gene read/TSS summary. All inputs required.'
    )
    parser.add_argument('--gff', required=True,
                        help='GFF with elements, LTRs, genes, and nested features')
    parser.add_argument('--bam', nargs='+', required=True,
                        help='One or more BAM files')
    parser.add_argument('--min_mapq', type=int, default=30,
                        help='Minimum MAPQ for filtering reads')
    parser.add_argument('--output', required=True,
                        help='Output prefix for LTR TSV and summaries')
    parser.add_argument('--tss_out', required=True,
                        help='Output file for LTR isoform TSS summary')
    parser.add_argument('--gene_out', required=True,
                        help='Output file for gene read counts and TSS summary')
    return parser.parse_args()

# Helper function to parse GFF attributes into a dictionary
def parse_attributes(attrs):
    attr_dict = {}
    for pair in attrs.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attr_dict[key] = value
    return attr_dict

# Load GFF and build data structures for elements, genes, and nested features for removal
def load_elements_and_ranges(gff_path):
    print("Reading GFF and parsing elements...")
    cols = ['chrom','src','feature','start','end','score','strand','frame','attrs']
    df = pd.read_csv(gff_path, sep='\t', comment='#', header=None, names=cols)

    full_df   = df[df['feature'].str.contains('LTR_retrotransposon', na=False)]
    ltr_df    = df[df['feature']=='long_terminal_repeat']
    gene_df   = df[df['feature']=='gene']
    nested_df = df.drop(full_df.index).drop(ltr_df.index).drop(gene_df.index)

    element_info = {}
    gene_info = {}

    # Parse LTR elements
    for _, row in full_df.iterrows():
        attrs = parse_attributes(row['attrs'])
        eid = attrs.get('Parent')
        if not eid:
            continue
        element_info[eid] = {
            'chrom': row['chrom'],
            'start': int(row['start'])-1,
            'end': int(row['end']),
            'name': attrs.get('Name',''),
            'attrs': row['attrs'],
            'ltr_left': None,
            'ltr_right': None,
            'strand': row['strand']
        }

    # Parse genes
    for _, row in gene_df.iterrows():
        attrs = parse_attributes(row['attrs'])
        gid = attrs.get('ID')
        if not gid:
            continue
        gene_info[gid] = {
            'chrom': row['chrom'],
            'start': int(row['start'])-1,
            'end': int(row['end']),
            'name': attrs.get('Name',''),
            'attrs': row['attrs'],
            'strand': row['strand']
        }

    # Parse LTR coordinates
    for _, row in ltr_df.iterrows():
        attrs = parse_attributes(row['attrs'])
        eid = attrs.get('Parent')
        if eid in element_info:
            s, e = int(row['start'])-1, int(row['end'])
            lid = attrs.get('ID','')
            if lid.startswith('l'):
                element_info[eid]['ltr_left'] = (s, e)
            elif lid.startswith('r'):
                element_info[eid]['ltr_right'] = (s, e)

    print(f"Found {len(element_info)} LTR elements; "
          f"{sum(1 for v in element_info.values() if v['ltr_left'] and v['ltr_right'])} have both LTRs.")
    print("Indexing intervals with intervaltree...")

    full_tree      = defaultdict(IntervalTree)
    nested_tree    = defaultdict(IntervalTree)
    gene_tree      = defaultdict(IntervalTree)
    nested_info    = {}

    # Index LTR elements
    for eid, info in element_info.items():
        chrom = info['chrom']
        if info['ltr_left']:
            full_tree[chrom][info['ltr_left'][0]:info['ltr_left'][1]] = eid
        if info['ltr_right']:
            full_tree[chrom][info['ltr_right'][0]:info['ltr_right'][1]] = eid

    # Index genes
    for gid, info in gene_info.items():
        chrom = info['chrom']
        gene_tree[chrom][info['start']:info['end']] = gid

    # Index nested features
    for _, row in nested_df.iterrows():
        chrom = row['chrom']
        s, e = int(row['start'])-1, int(row['end'])
        nid = parse_attributes(row['attrs']).get('ID','nested')
        parents = {iv.data for iv in full_tree[chrom].overlap(s, e)}
        if parents:
            nested_info[nid] = {'parent': list(parents), 'chrom': chrom, 'start': s, 'end': e}
            nested_tree[chrom][s:e] = nid

    print(f"Indexed {len(nested_info)} nested features and {len(gene_info)} genes.")

    return element_info, gene_info, full_tree, nested_tree, nested_info, gene_tree

# Get the strand for the read according to transcription direction determined by PyChopper if available.
# Fallbacks are in descending priority.
def infer_read_strand(read):
    """
    Return '+', '-', or '.' using fallbacks:
      1) ST tag (PyChopper)
      2) ts tag (minimap2 transcript strand)
      3) XS tag
      4) jM tag (minimap2 splice motif)
      5) alignment orientation (is_reverse)
    """
    # 1) PyChopper ST tag
    if read.has_tag("ST"):
        st = read.get_tag("ST")
        if st in {"+", "-"}:
            return st

    # 2) minimap2 transcript strand
    if read.has_tag("ts"):
        ts = read.get_tag("ts")
        # pysam usually returns '+' or '-' for ts:A:+ / ts:A:-
        if ts in {"+", "-"}:
            return ts

    # 3) XS tag (common in some spliced alignments)
    if read.has_tag("XS"):
        xs = read.get_tag("XS")
        if xs in {"+", "-"}:
            return xs

    # 4) minimap2 intron motif-based inference
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

    # 5) fallback: alignment orientation
    try:
        return "-" if read.is_reverse else "+"
    except Exception:
        return "."

# Get the strand for element by aggregating read evidence across all BAMs, only for those without a valid strand in the GFF.
def infer_element_strands(elem_info, full_tree, bam_paths, min_mapq):
    # only track those without a valid '+' or '-'
    to_infer = {
        eid: {'+': 0, '-': 0}
        for eid, e in elem_info.items()
        if e.get('strand') not in ('+','-')
    }

    for path in bam_paths:
        bf = pysam.AlignmentFile(path, 'rb')
        for read in bf.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            strand = infer_read_strand(read)
            if strand not in ("+", "-"):
                continue  # optional: skip unknowns instead of forcing them

            origin = read.reference_start if strand == '+' else read.reference_end - 1
            chrom  = bf.get_reference_name(read.reference_id)

            # for each element overlapping this read-origin
            for iv in full_tree[chrom].at(origin):
                eid = iv.data
                if eid in to_infer:
                    to_infer[eid][strand] += 1

    # assign the majority
    for eid, counts in to_infer.items():
        elem_info[eid]['strand'] = '+' if counts['+'] >= counts['-'] else '-'
        
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

# Reads in multiple bam files if desired, aggregates all reads overlapping each element, 
# and classifies into categories based on start/end positions relative to LTRs and splicing patterns.
def classify_multiple_bams(elem_info, gene_info, full_tree, nested_tree, nested_info,
                           gene_tree, bam_paths, min_mapq, clip_out_splice=None, clip_out_nonsplice=None,
                           ltr_exon_out=None, gene_exon_out=None):
    strand_counts = Counter()
    cats = ['ltr_left','ltr_right','spanning','ro5','ro3', 'spliced_ltr_left','spliced_ltr_right']

    stats = {eid: {'total':0,
                   'counts':{c:0 for c in cats},
                   'lengths':{c:0 for c in cats},
                   'spliced':{c:0 for c in cats},
                   'junctions':{c:set() for c in cats},
                   'strands':{c:[] for c in cats}}
             for eid in elem_info}

    gene_stats = {gid:{'total':0} for gid in gene_info}
    gene_tss = {gid:[] for gid in gene_info}
    tss_positions = {eid: {c: [] for c in cats} for eid in elem_info}
    end_positions = {eid: {c: [] for c in cats} for eid in elem_info}

    # Store per-read clipping records
    clip_rows_splice = []
    clip_rows_nonsplice = []
    ltr_exon_rows = []
    gene_exon_rows = []

    total_reads = sum(
        sum(s.mapped for s in pysam.AlignmentFile(p,'rb').get_index_statistics())
        for p in bam_paths
    )
    print(f"Total reads across {len(bam_paths)} BAMs: {total_reads}")
    print("Filtering reads originating in nested features...")

    processed, start_time = 0, time.time()
    interval = max(total_reads//10, 10000)

    for path in bam_paths:
        print(f"Processing {path}...")
        bf = pysam.AlignmentFile(path,'rb')
        for read in bf.fetch(until_eof=True):
            processed += 1
            if processed % interval == 0:
                elapsed = time.time() - start_time
                est = elapsed / processed * total_reads
                print(f"{processed}/{total_reads} reads; elapsed {elapsed:.1f}s; ~{(est-elapsed)/60:.1f}m remaining")

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            strand = infer_read_strand(read)
            if strand not in ("+", "-"):
                continue
            strand_counts[strand] += 1

            origin = read.reference_start if strand=='+' else read.reference_end-1
            endpos = read.reference_end-1 if strand=='+' else read.reference_start
            chrom = bf.get_reference_name(read.reference_id)

            if nested_tree[chrom].at(origin):
                continue

            L = read.reference_start
            R = read.reference_end  # half-open
            gene_hits = gene_tree[chrom].overlap(L, R)
            if gene_hits:       
                MIN_INTRON_LEN = 69
                ex_n, in_n, ex_total, in_total, exin_fields, saw_real = exon_intron_row_fields(
                read,
                min_intron_len=MIN_INTRON_LEN,
                max_exons=5,
                length_mode="query"
                )
                #n_real_N = sum(1 for op,ln in (read.cigartuples or []) if op==3 and ln>=MIN_INTRON_LEN)
                #if n_real_N >= 4:
                    #print("5+ exon candidate:", read.query_name, "real_N=", n_real_N, "cigar=", read.cigarstring)
                # Only write rows when spliced
                if saw_real:
                    for iv in gene_hits:
                        gid = iv.data
                        gstrand = gene_info[gid]['strand']
                        if strand != gstrand:
                            continue
                        gene_stats[gid]['total'] += 1
                        gene_tss[gid].append(origin)
                        gene_exon_rows.append((
                            gid,
                            chrom,
                            origin,
                            read.query_name,
                            read.mapping_quality,
                            int(read.is_reverse),
                            ex_n,
                            in_n,
                            ex_total,
                            in_total,
                            *exin_fields
                        ))
                else:
                    for iv in gene_hits:
                        gid = iv.data
                        gstrand = gene_info[gid]['strand']
                        if strand != gstrand:
                            continue
                        gene_stats[gid]['total'] += 1
                        gene_tss[gid].append(origin)

                continue

            # compute once per read (cheap)
            clip3S = softclip_3prime(read)

            for eid in {iv.data for iv in full_tree[chrom].at(origin)}:
                element = elem_info[eid]
                elem_strand = element['strand']
                ltr5 = element['ltr_left'] if elem_strand=='+' else element['ltr_right']
                ltr3 = element['ltr_right'] if elem_strand=='+' else element['ltr_left']
                if not ltr5 or not ltr3:
                    continue

                cat = None
                same_strand = (strand == elem_strand)

                # Starts in the 5′ LTR
                if ltr5[0] <= origin < ltr5[1] and same_strand:
                    if ltr5[0] <= endpos < ltr5[1]:
                        cat = 'ltr_left' if elem_strand == '+' else 'ltr_right'
                    elif ltr3[0] <= endpos < ltr3[1]:
                        cat = 'spanning'
                    elif ((elem_strand == '+' and endpos >= ltr3[1]) or
                          (elem_strand == '-' and endpos <  ltr3[0])):
                        cat = 'ro5'

                # Starts in the 3′ LTR
                elif ltr3[0] <= origin < ltr3[1] and same_strand:
                    if ltr3[0] <= endpos < ltr3[1]:
                        cat = 'ltr_right' if elem_strand == '+' else 'ltr_left'
                    elif ((elem_strand == '+' and endpos >= ltr3[1]) or
                          (elem_strand == '-' and endpos <  ltr3[0])):
                        cat = 'ro3'
                    elif ltr5[0] <= endpos < ltr5[1]:
                        cat = 'spanning'

                MIN_INTRON_LEN = 69 
                has_real_intron = False
                if read.cigartuples:
                    for op, length in read.cigartuples:
                        if op == 3 and length >= MIN_INTRON_LEN:
                            has_real_intron = True
                            break
                        
                if cat == 'spanning':
                    coding_start = element['ltr_left'][1]
                    coding_end   = element['ltr_right'][0]
                    coding_len   = coding_end - coding_start
                    if coding_len > 0:
                        covered = sum(
                            max(0, min(b2, coding_end) - max(b1, coding_start))
                            for b1, b2 in read.get_blocks()
                        )
                        if covered/coding_len < 0.5 and has_real_intron:
                            if ltr5[0] <= origin < ltr5[1]:
                                cat = 'spliced_ltr_left' if elem_strand=='+' else 'spliced_ltr_right'
                            elif ltr3[0] <= origin < ltr3[1]:
                                cat = 'spliced_ltr_right' if elem_strand=='+' else 'spliced_ltr_left'
                            else:
                                cat = None

                if not cat:
                    continue

                # ---- record isoform stats (existing) ----
                rec = stats[eid]
                rec['total'] += 1
                rec['counts'][cat] += 1
                rec['lengths'][cat] += read.query_length
                rec['strands'][cat].append(strand)
                tss_positions[eid][cat].append(origin)
                end_positions[eid][cat].append(endpos)
                if has_real_intron:
                    rec['spliced'][cat] += 1

                # record of 3' soft-clipping for spliced reads indicating if ltr-contained ----
                stayed_in_ltr5 = (ltr5[0] <= origin < ltr5[1]) and (ltr5[0] <= endpos < ltr5[1])
                stayed_in_ltr3 = (ltr3[0] <= origin < ltr3[1]) and (ltr3[0] <= endpos < ltr3[1])

                n_exons, n_introns, ex_total, in_total, fields, saw_real = exon_intron_row_fields(
                    read,
                    min_intron_len=MIN_INTRON_LEN,
                    max_exons=5,
                    length_mode="query"
                )

                if saw_real:
                    # Checking exons only of potential lncRNA loci
                    if (stayed_in_ltr5 or stayed_in_ltr3):
                        ltr_exon_rows.append((
                            eid, chrom, origin, cat, read.query_name, read.mapping_quality,
                            int(read.is_reverse),
                            n_exons, n_introns, ex_total, in_total,
                            *fields
                        ))

                    clip_rows_splice.append((
                        eid,
                        chrom,
                        origin,
                        "ltr5c_spliced" if stayed_in_ltr5 else "ltr3c_spliced" if stayed_in_ltr3 else cat,
                        cat,
                        read.query_name,
                        read.mapping_quality,
                        int(read.is_reverse),
                        clip3S,
                        read.reference_start,
                        read.reference_end
                    ))
                else:
                    clip_rows_nonsplice.append((
                        eid,
                        chrom,
                        origin,
                        "ltr5c_nonspliced" if stayed_in_ltr5 else "ltr3c_nonspliced" if stayed_in_ltr3 else cat,
                        cat,
                        read.query_name,
                        read.mapping_quality,
                        int(read.is_reverse),
                        clip3S,
                        read.reference_start,
                        read.reference_end
                    ))

        bf.close()

    print(f"Classification complete in {(time.time()-start_time)/60:.1f} minutes.")
    print(f"Strand counts across all reads:", strand_counts)

    # Write clipping table if requested
    if clip_out_splice:
        with open(clip_out_splice, "w") as out:
            out.write("\t".join([
                "EID","Chrom","Origin","LTR_side","Category","Read","MAPQ","is_reverse",
                "softclip_3p","aln_start","aln_end"
            ]) + "\n")
            for row in clip_rows_splice:
                out.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote 3' soft-clipping records for spliced reads to {clip_out_splice}")

    if clip_out_nonsplice:
        with open(clip_out_nonsplice, "w") as out:
            out.write("\t".join([
                "EID","Chrom","Origin","LTR_side","Category","Read","MAPQ","is_reverse",
                "softclip_3p","aln_start","aln_end"
            ]) + "\n")
            for row in clip_rows_nonsplice:
                out.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote 3' soft-clipping records for non-spliced reads to {clip_out_nonsplice}")
        
    if ltr_exon_out:
        with open(ltr_exon_out, "w") as out:
            out.write("\t".join([
                "EID","Chrom","Origin","Category","Read","MAPQ","is_reverse",
                "exon_count","intron_count","exon_len_combined","intron_len_combined",
                "exon1","intron1","exon2","intron2","exon3","intron3","exon4","intron4","exon5"
            ]) + "\n")
            for row in ltr_exon_rows:
                out.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote LTR-RT exon statistics to {ltr_exon_out}")

    if gene_exon_out:
        with open(gene_exon_out, "w") as out:
            out.write("\t".join([
                "GeneID","Chrom","Origin","Read","MAPQ","is_reverse",
                "exon_count","intron_count","exon_len_combined","intron_len_combined",
                "exon1","intron1","exon2","intron2","exon3","intron3","exon4","intron4","exon5"
            ]) + "\n")
            for row in gene_exon_rows:
                out.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote Gene exon statistics to {gene_exon_out}")

    return stats, tss_positions, end_positions, gene_stats, gene_tss
    
def write_tsv(elem_info, stats, prefix):
    cats = ['ltr_left','ltr_right','spanning','ro5','ro3']
    rows = []
    for eid, info in elem_info.items():
        rec = stats[eid]
        total = rec['total']
        row = {
            'chrom': info['chrom'], 'name': info['name'],
            'start': info['start'], 'end': info['end'],
            'attrs': info['attrs'], 'total_reads': total
        }
        for cat in cats:
            count = rec['counts'][cat]
            pct = count/total*100 if total else 0
            mean_len = rec['lengths'][cat]/count if count else 0
            spc = rec['spliced'][cat]
            junc_count = len(rec['junctions'][cat])
            row.update({
                f'{cat}_reads': count,
                f'{cat}_pct': f"{pct:.2f}",
                f'mean_len_{cat}': f"{mean_len:.1f}",
                f'spliced_{cat}': spc,
                f'unique_juncts_{cat}': junc_count
            })
        rows.append(row)
    df = pd.DataFrame(rows).sort_values('total_reads', ascending=False)
    out_path = prefix + '.tsv'
    df.to_csv(out_path, sep='\t', index=False)
    print(f"Wrote LTR summary to {out_path}")


def write_isoform_tss_summary(stats, tss_positions, out_path):
    with open(out_path, 'w') as out:
        out.write("Feature\tTop_Isoform\tTop_Isoform_Strand\tTotal_Reads\tTSS1\tCount1\tTSS2\tCount2\n")
        for eid, rec in stats.items():
            total = rec['total']
            if total == 0:
                continue
            top_cat = max(rec['counts'], key=rec['counts'].get)
            common = Counter(tss_positions[eid][top_cat]).most_common(2)
            (tss1,c1),(tss2,c2) = (common + [( 'NA',0),( 'NA',0) ])[:2]

            strand_counts = Counter([s for s in rec.get('strands', {}).get(top_cat, []) if s in ('+','-')])
            strand = strand_counts.most_common(1)[0][0] if strand_counts else 'NA'

            out.write(f"{eid}\t{top_cat}\t{strand}\t{total}\t{tss1}\t{c1}\t{tss2}\t{c2}\n")
    print(f"Wrote LTR TSS summary to {out_path}")
    

def write_isoform_tss_summary_top_n(stats, tss_positions, end_positions, tss_out_path, cleave_out_path, n=10):
    """
    Writes Feature, Strand, Total_Reads, then TSS1,Count1,...,TSSn,Countn
    up to n sites (pads with NA if fewer). Counts are drawn from the full
    distribution of TSS positions (not just the top n).
    """
     # ---------- TSS ----------
    with open(tss_out_path, "w") as out:
        hdr = ["Feature", "Strand", "Total_Reads"]
        for i in range(1, n + 1):
            hdr += [f"TSS{i}", f"Count{i}"]
        out.write("\t".join(hdr) + "\n")

        for eid, rec in stats.items():
            total = rec.get("total", 0)
            if total == 0:
                continue

            top_cat = max(rec["counts"], key=rec["counts"].get)

            ctr = Counter(tss_positions[eid][top_cat])
            top_n = ctr.most_common(n)
            top_n += [("NA", 0)] * (n - len(top_n))

            strands = rec.get("strands", {}).get(top_cat, [])
            strand = Counter(s for s in strands if s in ("+", "-")).most_common(1)
            strand = strand[0][0] if strand else "NA"

            row = [eid, strand, str(total)]
            for tss, cnt in top_n:
                row += [str(tss), str(cnt)]
            out.write("\t".join(row) + "\n")

    # ---------- Cleavage ----------
    window = 5
    with open(cleave_out_path, "w") as out:
        hdr = ["Feature", "Strand", "Total_Reads",
            "CleavageSitePeak", "CleavageSiteCnt", "CleavageSiteFrac",
            f"PeakWindow{window}_Cnt", f"PeakWindow{window}_Frac"]
        for i in range(1, n + 1):
            hdr += [f"End{i}", f"Count{i}"]
        out.write("\t".join(hdr) + "\n")

        for eid, rec in stats.items():
            # collect ALL end positions across ALL categories
            all_ends = []
            all_strands = []
            for cat, cnt in rec["counts"].items():
                if cnt == 0:
                    continue
                all_ends.extend(end_positions[eid][cat])
                all_strands.extend(rec.get("strands", {}).get(cat, []))

            ctr = Counter(all_ends)
            total = sum(ctr.values())  # safest denominator

            if total == 0:
                continue

            common = ctr.most_common()
            top_n = common[:n] + [("NA", 0)] * (n - len(common[:n]))

            peak_site, peak_cnt = common[0]
            cleave_frac = peak_cnt / total

            peak_window_cnt = sum(v for pos, v in ctr.items()
                                if abs(pos - peak_site) <= window)
            peak_window_frac = peak_window_cnt / total

            strand = Counter(s for s in all_strands if s in ("+", "-")).most_common(1)
            strand = strand[0][0] if strand else "NA"

            row = [eid, strand, str(total),
                str(peak_site), str(peak_cnt), f"{cleave_frac:.4f}",
                str(peak_window_cnt), f"{peak_window_frac:.4f}"]
            for end, cnt in top_n:
                row += [str(end), str(cnt)]

            out.write("\t".join(row) + "\n")

        print(f"Wrote top-{n} LTR TSS summary to {tss_out_path}")
        print(f"Wrote top-{n} LTR cleavage summary to {cleave_out_path}")

def write_gene_summary(gene_stats, gene_tss, out_path):
    with open(out_path, 'w') as out:
        out.write("Gene\tTotal_Reads\tTSS1\tCount1\tTSS2\tCount2\n")
        for gid, rec in gene_stats.items():
            total = rec['total']
            if total == 0:
                continue
            common = Counter(gene_tss[gid]).most_common(2)
            (tss1,c1),(tss2,c2) = (common + [( 'NA',0),( 'NA',0) ])[:2]
            out.write(f"{gid}\t{total}\t{tss1}\t{c1}\t{tss2}\t{c2}\n")
    print(f"Wrote gene summary to {out_path}")

def write_gene_summary_top_n(gene_stats, gene_tss, out_path, n=10):
    """
    Writes a tab-delimited file with columns:
      Gene, Total_Reads,
      TSS1, Count1, ..., TSSn, Countn
    Up to n sites. Counts are drawn from the full
    distribution of gene_tss positions (not just the top n).
    """
    with open(out_path, 'w') as out:
        # header
        headers = ["Gene", "Total_Reads"]
        for i in range(1, n+1):
            headers += [f"TSS{i}", f"Count{i}"]
        out.write("\t".join(headers) + "\n")

        for gid, rec in gene_stats.items():
            total = rec.get('total', 0)
            if total == 0:
                continue

            # build full Counter over gene-TSS positions
            ctr = Counter(gene_tss.get(gid, []))
            common_full = ctr.most_common()
            # slice to top-n and pad
            top_n = common_full[:n]
            top_n += [("NA", 0)] * (n - len(top_n))

            # assemble row
            row = [gid, str(total)]
            for tss, cnt in top_n:
                row += [str(tss), str(cnt)]

            out.write("\t".join(row) + "\n")

    print(f"Wrote top-{n} gene TSS summary to {out_path}")
    
# Computes TSS density for histograms around the primary and secondary TSS for each feature, 
# with strand-aware distances and optional inclusion of reads at the TSS.
def compute_tss_density_separate(stats, tss_positions, elem_info,
                                 window=10, min_reads=1, include_tss_reads=True):
    """
    For each feature in stats:
      - find the two most common TSS (primary, secondary)
      - collect *all* read‐end origins
      - build TWO histograms (±window):
         • hist1 around primary TSS, dropping only reads == primary
         • hist2 around secondary TSS, dropping only reads == secondary
      - strand‐aware: on '-' features flip the sign of (o - TSS)
      - If include_tss_reads=True, reads at the TSS (distance==0) are INCLUDED.
    Returns:
      density: dict mapping feature_id ->
        { 'primary': np.array(length=2*window+1),
          'secondary': np.array(length=2*window+1) }
    """
    density = {}
    size = 2 * window + 1

    for fid, rec in stats.items():
        total = rec.get('total', 0)
        if total < min_reads:
            continue

        # pick primary & secondary TSS
        top_cat = max(rec['counts'], key=rec['counts'].get)
        common = Counter(tss_positions[fid][top_cat]).most_common(2)
        while len(common) < 2:
            common.append((None, 0))
        (tss1, _), (tss2, _) = common
        if tss1 is None:
            continue

        # collect all read origins for this feature
        origins = []
        for lst in tss_positions[fid].values():
            origins.extend(lst)

        # strand-aware distance
        strand = elem_info[fid].get('strand', '+')
        def stranded_dist(o, t):
            d = o - t
            return -d if strand == '-' else d

        hist1 = np.zeros(size, dtype=int)
        hist2 = np.zeros(size, dtype=int)

        # primary window (now INCLUDING reads at TSS1 if include_tss_reads=True)
        for o in origins:
            if not include_tss_reads and o == tss1:
                continue
            d = stranded_dist(o, tss1)
            if -window <= d <= window:
                hist1[d + window] += 1

        # secondary window (same inclusion behavior)
        if tss2 is not None:
            for o in origins:
                if not include_tss_reads and o == tss2:
                    continue
                d = stranded_dist(o, tss2)
                if -window <= d <= window:
                    hist2[d + window] += 1

        density[fid] = {'primary': hist1, 'secondary': hist2}

    return density

def compute_gene_tss_density_separate(gene_stats, gene_tss, gene_info,
                                      window=10, min_reads=5, include_tss_reads=True):
    """
    For each gene in gene_stats:
      - Identify the two most common TSS positions from gene_tss[gid]
      - Build two histograms (length 2*window+1):
         • hist1: counts of ALL read-ends within ±window of primary TSS
         • hist2: counts of ALL read-ends within ±window of secondary TSS
        (If include_tss_reads=True, reads exactly at the TSS are INCLUDED.)
      - Strand-aware: distances on '-' genes are flipped so upstream is positive.

    Returns:
      dict: gid -> {'primary': hist1, 'secondary': hist2}
    """
    density = {}
    size = 2 * window + 1

    for gid, rec in gene_stats.items():
        total = rec.get('total', 0)
        if total < min_reads:
            continue

        origins = gene_tss.get(gid, [])
        if not origins:
            continue

        # primary & secondary TSS
        common = Counter(origins).most_common(2)
        while len(common) < 2:
            common.append((None, 0))
        (tss1, _), (tss2, _) = common
        if tss1 is None:
            continue

        # strand-aware distance
        strand = gene_info.get(gid, {}).get('strand', '+')
        def stranded_dist(o, t):
            d = o - t
            return -d if strand == '-' else d

        hist1 = np.zeros(size, dtype=int)
        hist2 = np.zeros(size, dtype=int)

        # primary window
        for o in origins:
            if not include_tss_reads and o == tss1:
                continue
            d = stranded_dist(o, tss1)
            if -window <= d <= window:
                hist1[d + window] += 1

        # secondary window
        if tss2 is not None:
            for o in origins:
                if not include_tss_reads and o == tss2:
                    continue
                d = stranded_dist(o, tss2)
                if -window <= d <= window:
                    hist2[d + window] += 1

        density[gid] = {'primary': hist1, 'secondary': hist2}

    return density

def main():
    args = parse_args()
    clip_out_splice = args.output + "_ltr_3p_softclip_per_read_spliced.tsv"
    clip_out_nonsplice = args.output + "_ltr_3p_softclip_per_read_nonspliced.tsv"
    ltr_exon_out = args.output + "_ltr_exon_stats_per_read.tsv"
    gene_exon_out = args.output + "_gene_exon_stats_per_read.tsv"
    elem_info, gene_info, full_tree, nested_tree, nested_info, gene_tree = load_elements_and_ranges(args.gff)
    stats, tss_positions, end_positions, gene_stats, gene_tss = classify_multiple_bams(
        elem_info, gene_info, full_tree, nested_tree, nested_info, gene_tree,
        args.bam, args.min_mapq,
        clip_out_splice,
        clip_out_nonsplice,
        ltr_exon_out,
        gene_exon_out
    )
    infer_element_strands(elem_info, full_tree, args.bam, args.min_mapq)
    densities = compute_tss_density_separate(stats,
                                         tss_positions,
                                         elem_info,
                                         window=10,
                                         min_reads=6)
    with open(args.output + '_primary_tss_density.tsv','w') as out:
        out.write("Feature\tDistance\tCount\n")
        for fid, twohists in densities.items():
            hist = twohists['primary']
            for i, cnt in enumerate(hist):
                dist = i - 10
                out.write(f"{fid}\t{dist}\t{cnt}\n")
    print("Wrote primary TSS densities to", args.output + '_primary_tss_density.tsv')

    # write SECONDARY‐only densities
    with open(args.output + '_secondary_tss_density.tsv','w') as out:
        out.write("Feature\tDistance\tCount\n")
        for fid, twohists in densities.items():
            hist = twohists['secondary']
            for i, cnt in enumerate(hist):
                dist = i - 10
                out.write(f"{fid}\t{dist}\t{cnt}\n")
    print("Wrote secondary TSS densities to", args.output + '_secondary_tss_density.tsv')
    
    gene_dens = compute_gene_tss_density_separate(gene_stats,
                                                  gene_tss,
                                                  gene_info,
                                                  window=10,
                                                  min_reads=7)

    # write PRIMARY gene densities
    with open(args.output + '_gene_primary_density.tsv','w') as out:
        out.write("Gene\tDistance\tCount\n")
        for gid, h in gene_dens.items():
            hist = h['primary']
            for i, cnt in enumerate(hist):
                dist = i - 10
                out.write(f"{gid}\t{dist}\t{cnt}\n")
    print("Wrote gene primary densities to", args.output + '_gene_primary_density.tsv')

    # write SECONDARY gene densities
    with open(args.output + '_gene_secondary_density.tsv','w') as out:
        out.write("Gene\tDistance\tCount\n")
        for gid, h in gene_dens.items():
            hist = h['secondary']
            for i, cnt in enumerate(hist):
                dist = i - 10
                out.write(f"{gid}\t{dist}\t{cnt}\n")
    print("Wrote gene secondary densities to", args.output + '_gene_secondary_density.tsv')

    write_tsv(elem_info, stats, args.output)
    write_isoform_tss_summary(stats, tss_positions, args.tss_out)
    write_isoform_tss_summary_top_n(stats, tss_positions,end_positions,
                                    cleave_out_path="10site.ltr_cleavage_summary.tsv", tss_out_path="10site.ltr_tss_summary.tsv",n=10)
    write_gene_summary(gene_stats, gene_tss, args.gene_out)
    write_gene_summary_top_n(gene_stats,gene_tss,out_path="10site.gene_summary.tsv",n=10)

if __name__ == '__main__':
    main()