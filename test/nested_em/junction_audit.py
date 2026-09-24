#!/usr/bin/env python3
"""
Audit --assign em against evidence the EM never sees: splice junctions, alignment direction
and mapping quality.

The EM decides among candidate features from 5'-end position, 3'-end position, read length and
a spliced/unspliced summary. It never compares a read's junctions to the annotated intron chains
of the candidates. That makes junctions an independent ground truth wherever a read carries them:

  - a read whose junctions match the annotated introns of exactly one candidate gene, and no
    others, was transcribed from that gene, whatever its 5' end suggests;
  - the gene's strand then fixes the orientation the read must have been assigned;
  - reads carrying junctions that match no candidate are dropped (mismapping, chimeras, unannotated
    isoforms), as are low-MAPQ reads.

Accuracy is reported for the EM and for the innermost rule on the same labelled subset, so the two
are compared on identical reads. Junction-bearing reads are not a random sample of contested reads
(they are longer and better anchored), so the numbers are a diagnostic of relative behaviour on
resolvable reads, not an unbiased error rate over all reads.

Run it on the outputs of a --assign em run over a whole chromosome.
"""
import argparse
import sys
from collections import defaultdict

import pysam


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--em-reads', dest='em_reads', required=True,
                    help='<output>_em_reads.tsv from a --assign em run')
    ap.add_argument('--bam', required=True, help='The same BAM that run used')
    ap.add_argument('--gene-gff', dest='gene_gff', required=True,
                    help='The same --gene-gff, for annotated intron chains')
    ap.add_argument('--min-mapq', dest='min_mapq', type=int, default=30)
    ap.add_argument('--min-intron-len', dest='min_intron_len', type=int, default=69,
                    help='Matches IsoClassifier --min-intron-len (default 69)')
    ap.add_argument('--tol', type=int, default=10,
                    help='Junction is called matching within this many bp at both ends (default 10)')
    ap.add_argument('--out', default='junction_audit.tsv', help='Per-read audit table')
    return ap.parse_args()


def load_gene_introns(path):
    """Union of annotated introns per gene, over every transcript, as 0-based half-open pairs."""
    tx_exons = defaultdict(list)
    tx_gene, gene_strand = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            attrs = dict(kv.split('=', 1) for kv in f[8].rstrip(';').split(';') if '=' in kv)
            if f[2] == 'exon':
                parent = attrs.get('Parent', '').split(',')[0]
                if parent:
                    tx_exons[parent].append((int(f[3]) - 1, int(f[4])))
            elif f[2] in ('mRNA', 'transcript'):
                tid, gid = attrs.get('ID'), attrs.get('Parent', '').split(',')[0]
                if tid and gid:
                    tx_gene[tid] = gid
            elif f[2] == 'gene':
                gid = attrs.get('ID')
                if gid:
                    gene_strand[gid] = f[6]
    introns = defaultdict(set)
    for tid, exons in tx_exons.items():
        gid = tx_gene.get(tid)
        if not gid:
            continue
        exons.sort()
        for i in range(len(exons) - 1):
            if exons[i + 1][0] > exons[i][1]:
                introns[gid].add((exons[i][1], exons[i + 1][0]))
    return introns, gene_strand


def load_contested(path):
    """Read name -> (chrom, innermost, assigned, assigned_orient, posterior, candidates)."""
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            cands = []
            for item in f[idx['Candidates']].split(';'):
                parts = item.split('|')
                if len(parts) >= 2:
                    cands.append((parts[0], parts[1]))
            out[f[idx['Read']]] = dict(
                chrom=f[idx['Chrom']],
                innermost=f[idx['Innermost_Feature']],
                innermost_or=f[idx['Innermost_Orientation']],
                assigned=f[idx['Assigned_Feature']],
                assigned_or=f[idx['Assigned_Orientation']],
                posterior=float(f[idx['Posterior']]),
                tss_prob=float(f[idx['TSS_Init_Prob']]),
                candidates=cands)
    return out


def junctions_of(aln, min_intron_len):
    """0-based half-open intron coordinates implied by N ops of at least min_intron_len."""
    pos, out = aln.reference_start, []
    for op, ln in aln.cigartuples or []:
        if op in (0, 2, 7, 8):        # M D = X consume reference
            pos += ln
        elif op == 3:                 # N
            if ln >= min_intron_len:
                out.append((pos, pos + ln))
            pos += ln
    return out


def main():
    a = parse_args()
    introns, gene_strand = load_gene_introns(a.gene_gff)
    print(f'Annotated intron chains for {len(introns)} genes', file=sys.stderr)
    contested = load_contested(a.em_reads)
    chroms = {r['chrom'] for r in contested.values()}
    print(f'{len(contested)} contested reads on {len(chroms)} contig(s)', file=sys.stderr)

    def matches(j, gid):
        ref = introns.get(gid)
        if not ref:
            return False
        return any(abs(j[0] - s) <= a.tol and abs(j[1] - e) <= a.tol for s, e in ref)

    bam = pysam.AlignmentFile(a.bam)
    seen = set()
    counts = defaultdict(int)
    rows = []
    for chrom in sorted(chroms):
        for aln in bam.fetch(chrom):
            name = aln.query_name
            if name in seen or name not in contested:
                continue
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            seen.add(name)
            rec = contested[name]
            counts['fetched'] += 1
            if aln.mapping_quality < a.min_mapq:
                counts['drop_mapq'] += 1
                continue
            jn = junctions_of(aln, a.min_intron_len)
            if not jn:
                counts['drop_unspliced'] += 1
                continue
            hits = {gid: sum(matches(j, gid) for j in jn)
                    for gid, _ in rec['candidates'] if gid in introns}
            supported = [g for g, n in hits.items() if n > 0]
            total_matched = sum(hits.values())
            if len(supported) != 1:
                counts['drop_ambiguous' if supported else 'drop_no_candidate_match'] += 1
                continue
            truth = supported[0]
            if total_matched != hits[truth]:
                counts['drop_ambiguous'] += 1
                continue
            counts['labelled'] += 1
            strand = gene_strand.get(truth)
            read_strand = '-' if aln.is_reverse else '+'
            st = aln.get_tag('ST') if aln.has_tag('ST') else 'NA'
            rows.append(dict(
                read=name, truth=truth, truth_strand=strand or 'NA',
                em=rec['assigned'], em_or=rec['assigned_or'],
                innermost=rec['innermost'], innermost_or=rec['innermost_or'],
                posterior=rec['posterior'], tss_prob=rec['tss_prob'],
                mapq=aln.mapping_quality, njunc=len(jn), nmatch=hits[truth],
                novel_junc=len(jn) - total_matched, read_strand=read_strand, st=st,
                ncand=len(rec['candidates']),
                em_ok=int(rec['assigned'] == truth),
                inn_ok=int(rec['innermost'] == truth),
                flipped=int(rec['assigned_or'] != rec['innermost_or'])))

    with open(a.out, 'w') as fh:
        cols = list(rows[0].keys()) if rows else []
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    n = len(rows)
    print(f'\n== junction audit: {a.em_reads} ==')
    for k in ('fetched', 'drop_mapq', 'drop_unspliced', 'drop_no_candidate_match',
              'drop_ambiguous', 'labelled'):
        print(f'  {k:26s} {counts[k]:8d}')
    if not n:
        print('  no labelled reads; nothing to score')
        return
    em = sum(r['em_ok'] for r in rows)
    inn = sum(r['inn_ok'] for r in rows)
    print(f'\n  reads with an independent junction label: {n}')
    print(f'  EM correct       : {em:6d}  ({100 * em / n:5.1f}%)')
    print(f'  innermost correct: {inn:6d}  ({100 * inn / n:5.1f}%)')
    agree = sum(1 for r in rows if r['em'] == r['innermost'])
    print(f'  the two rules agree on {agree} reads ({100 * agree / n:.1f}%); '
          f'on the {n - agree} they differ, EM is right '
          f'{sum(r["em_ok"] for r in rows if r["em"] != r["innermost"])} times, innermost '
          f'{sum(r["inn_ok"] for r in rows if r["em"] != r["innermost"])} times')

    def breakdown(title, keyfn):
        agg = defaultdict(lambda: [0, 0, 0])
        for r in rows:
            c = agg[keyfn(r)]
            c[0] += 1
            c[1] += r['em_ok']
            c[2] += r['inn_ok']
        print(f'\n  -- {title} --')
        print(f'    {"group":28s} {"reads":>7s} {"EM":>8s} {"innermost":>10s}')
        for k in sorted(agg, key=lambda x: -agg[x][0]):
            t, e, i = agg[k]
            print(f'    {str(k):28s} {t:7d} {100 * e / t:7.1f}% {100 * i / t:9.1f}%')

    breakdown('EM changed the orientation', lambda r: 'flipped' if r['flipped'] else 'same orientation')
    breakdown('EM posterior', lambda r: ('>=0.99' if r['posterior'] >= 0.99 else
                                         '0.9-0.99' if r['posterior'] >= 0.9 else
                                         '0.6-0.9' if r['posterior'] >= 0.6 else '<0.6'))
    breakdown('junctions matched', lambda r: f'{min(r["nmatch"], 4)}{"+" if r["nmatch"] > 4 else ""}')
    breakdown('MAPQ', lambda r: '60' if r['mapq'] >= 60 else '50-59' if r['mapq'] >= 50 else '<50')
    breakdown('novel junctions on the read', lambda r: 'none' if r['novel_junc'] == 0 else 'some')

    # Orientation check: the gene's strand fixes what the orientation should have been.
    ok = sum(1 for r in rows if r['truth_strand'] != 'NA' and r['em'] == r['truth'])
    st_ok = sum(1 for r in rows if r['st'] == r['truth_strand'])
    st_known = sum(1 for r in rows if r['st'] in '+-')
    print(f'\n  reads whose EM feature is the junction-supported gene: {ok}')
    if st_known:
        print(f'  ST tag agrees with the junction-supported gene strand on '
              f'{st_ok}/{st_known} ({100 * st_ok / st_known:.1f}%) reads with a tag '
              f'-- the ceiling on any direction-based rule here')
    fwd = sum(1 for r in rows if r['read_strand'] == '+')
    print(f'  alignment direction: {fwd}/{n} reads map forward '
          f'({100 * fwd / n:.1f}%); reads are oriented before mapping, so this tracks '
          f'transcript direction only as well as that orientation step does')


if __name__ == '__main__':
    main()
