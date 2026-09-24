#!/usr/bin/env python3
"""
Simulate long reads from nested TE / TE-in-gene loci with known origins.

Writes a GFF (EDTA-style TEs + genes), a gene GFF3 with exons, and a sorted, indexed BAM whose
read names carry the true origin: <scenario>|<true_feature>|<true_orient>|<kind>|<i>
kind is 'full' (5' end at the TSS) or 'trunc' (the read lost its 5' end to RT drop-off or RNA
degradation, so its first base lies downstream of the TSS; its 3' end is unaffected).

Scenarios (0-based half-open coordinates):
  A  same-strand LTR-RT in LTR-RT          outer 200 molecules, inner 60
  B  opposite-strand LTR-RT in LTR-RT      outer 150, inner 40 (innermost rule invents antisense)
  C  TE in a gene intron, same strand      gene 300 (25% pre-mRNA), TE 25 + 15 TE-initiated
                                           chimeric reads that splice into the gene's exons
  E  silent LTR-RT in an expressed LTR-RT  outer 300, inner 0
"""
import argparse
import os
import random

import pysam

CHROM, CHROM_LEN = 'chrS', 145_000


def ltr_rt(n, start, end, strand, lltr, rltr):
    rr = f'repeat_region_{n}'
    rows = [
        ('Copia_LTR_retrotransposon', start, end, strand,
         f'ID=LTRRT_{n};Parent={rr};Name=sim_{n};Classification=LTR/Copia;Method=structural'),
        ('long_terminal_repeat', lltr[0], lltr[1], strand,
         f'ID=lLTR_{n};Parent={rr};Name=sim_{n};Method=structural'),
        ('long_terminal_repeat', rltr[0], rltr[1], strand,
         f'ID=rLTR_{n};Parent={rr};Name=sim_{n};Method=structural'),
    ]
    return rr, rows


# ---- annotation ------------------------------------------------------------------------------
FEATS = []
O1, r = ltr_rt(1, 10000, 20000, '+', (10000, 11000), (19000, 20000)); FEATS += r
I1, r = ltr_rt(2, 14000, 16000, '+', (14000, 14400), (15600, 16000)); FEATS += r
O2, r = ltr_rt(3, 30000, 40000, '+', (30000, 31000), (39000, 40000)); FEATS += r
I2, r = ltr_rt(4, 34000, 36000, '-', (34000, 34400), (35600, 36000)); FEATS += r
O3, r = ltr_rt(6, 70000, 80000, '+', (70000, 71000), (79000, 80000)); FEATS += r
O4, r = ltr_rt(8, 90000, 100000, '+', (90000, 91000), (99000, 100000)); FEATS += r
I4, r = ltr_rt(9, 94000, 96000, '+', (94000, 94400), (95600, 96000)); FEATS += r
I3, r = ltr_rt(7, 74000, 76000, '+', (74000, 74400), (75600, 76000)); FEATS += r
G1 = 'Gene1'
FEATS.append(('gene', 50000, 62000, '-', f'ID={G1};biotype=protein_coding'))
T1 = 'TE_homo_5'
FEATS.append(('LTR_retrotransposon', 58000, 59000, '-',
              f'ID={T1};Name=frag;Classification=LTR/unknown;Method=homology'))
GENE_EXONS = [(61000, 62000), (56000, 56500), (50000, 51000)]

# Scenarios G and H: a short TE annotated across a gene's own TSS. This is the geometry that lets
# the child-span subtraction cut a gene's promoter out of its own zone. G's transcripts are
# spliced, so junction identity can break the tie; H's are single-exon, so nothing can.
G2 = 'Gene2'
FEATS.append(('gene', 110000, 118000, '+', f'ID={G2};biotype=protein_coding'))
P2 = 'TE_homo_7'
FEATS.append(('LTR_retrotransposon', 110050, 110250, '+',
              f'ID={P2};Name=frag;Classification=LTR/unknown;Method=homology'))
GENE2_EXONS = [(110000, 111000), (114000, 114500), (117000, 118000)]

G3 = 'Gene3'
FEATS.append(('gene', 130000, 138000, '+', f'ID={G3};biotype=protein_coding'))
P3 = 'TE_homo_8'
FEATS.append(('LTR_retrotransposon', 130050, 130250, '+',
              f'ID={P3};Name=frag;Classification=LTR/unknown;Method=homology'))
GENE3_EXONS = [(130000, 138000)]

# ---- transcript models: (scenario, true_feature, strand, TSS, TES, exon_chain, n, trunc_frac) --
# exon chains are listed 5'->3' in transcription order, genomic half-open intervals.
MODELS = [
    ('A', O1, '+', 10300, 19700, [(10300, 12000), (13000, 19700)], 60, 0.55),   # spliced isoform
    ('A', O1, '+', 10300, 19700, [(10300, 19700)], 140, 0.55),
    ('A', I1, '+', 14150, 15850, [(14150, 15850)], 60, 0.50),
    ('B', O2, '+', 30300, 39700, [(30300, 39700)], 150, 0.55),
    ('B', I2, '-', 35850, 34150, [(34150, 35850)], 40, 0.50),
    ('C', G1, '-', 61950, 50050, [(61000, 61950), (56000, 56500), (50050, 51000)], 225, 0.50),
    ('C', G1, '-', 61950, 50050, [(50050, 61950)], 75, 0.50),                   # pre-mRNA
    ('C', T1, '-', 58950, 58050, [(58050, 58950)], 25, 0.40),
    ('C', T1, '-', 58950, 50050, [(58100, 58950), (56000, 56500), (50050, 51000)], 15, 0.30),
    ('E', O3, '+', 70300, 79700, [(70300, 79700)], 300, 0.55),
    ('G', G2, '+', 110100, 117900, [(110100, 111000), (114000, 114500), (117000, 117900)], 200, 0.55),
    ('H', G3, '+', 130100, 137900, [(130100, 137900)], 200, 0.55),
]

# Harder variant: inner read-out, outer transcripts that terminate early at the inner poly(A)
# signal (a 3'-end event), and a barely expressed inner element. Used with --hard, which also
# makes reads of long transcripts lose their 5' end more often (as in ONT cDNA) and adds 5' jitter.
MODELS_HARD = [
    ('A', O1, '+', 10300, 19700, [(10300, 12000), (13000, 19700)], 50, 0.55),
    ('A', O1, '+', 10300, 19700, [(10300, 19700)], 110, 0.55),
    ('A', O1, '+', 10300, 15850, [(10300, 15850)], 40, 0.55),   # stops at inner 3' LTR
    ('A', I1, '+', 14150, 15850, [(14150, 15850)], 42, 0.50),
    ('A', I1, '+', 14150, 19700, [(14150, 19700)], 18, 0.50),   # inner read-out to outer pA
] + [m for m in MODELS if m[0] in ('B', 'C', 'E')] + [
    ('F', O4, '+', 90300, 99700, [(90300, 99700)], 200, 0.55),
    ('F', I4, '+', 94150, 95850, [(94150, 95850)], 6, 0.50),     # barely expressed inner
    ('G', G2, '+', 110100, 117900, [(110100, 111000), (114000, 114500), (117000, 117900)], 200, 0.55),
    ('H', G3, '+', 130100, 137900, [(130100, 137900)], 200, 0.55),
]


def tx_positions(chain, strand):
    """Genomic coordinates of the mature transcript, 5'->3' (chain is listed 5'->3')."""
    out = []
    for s, e in chain:
        out.extend(range(s, e) if strand == '+' else range(e - 1, s - 1, -1))
    return out


def blocks_for(chain, five, three, strand):
    """Genomic aligned blocks of a read covering transcript from 5' end `five` to 3' end `three`."""
    lo, hi = (five, three + 1) if strand == '+' else (three, five + 1)
    blocks = []
    for s, e in sorted(chain):
        s2, e2 = max(s, lo), min(e, hi)
        if s2 < e2:
            blocks.append((s2, e2))
    return blocks


def cigar_from_blocks(blocks):
    cig, pos = [], blocks[0][0]
    for s, e in blocks:
        if s > pos:
            cig.append((3, s - pos))       # N
        cig.append((0, e - s))             # M
        pos = e
    return cig


def simulate(outdir, seed=7, hard=False):
    random.seed(seed)
    models = MODELS_HARD if hard else MODELS
    os.makedirs(outdir, exist_ok=True)

    with open(os.path.join(outdir, 'sim.gff'), 'w') as fh:
        for feature, s, e, strand, attrs in FEATS:
            fh.write(f'{CHROM}\tSIM\t{feature}\t{s + 1}\t{e}\t.\t{strand}\t.\t{attrs}\n')
    with open(os.path.join(outdir, 'sim_genes.gff3'), 'w') as fh:
        fh.write('##gff-version 3\n')
        for gid, gs, ge, gstrand, gexons in ((G1, 50000, 62000, '-', GENE_EXONS),
                                             (G2, 110000, 118000, '+', GENE2_EXONS),
                                             (G3, 130000, 138000, '+', GENE3_EXONS)):
            fh.write(f'{CHROM}\tSIM\tgene\t{gs + 1}\t{ge}\t.\t{gstrand}\t.\tID=gene:{gid}\n')
            fh.write(f'{CHROM}\tSIM\tmRNA\t{gs + 1}\t{ge}\t.\t{gstrand}\t.\t'
                     f'ID=transcript:{gid}_T1;Parent=gene:{gid}\n')
            for s, e in gexons:
                fh.write(f'{CHROM}\tSIM\texon\t{s + 1}\t{e}\t.\t{gstrand}\t.\t'
                         f'Parent=transcript:{gid}_T1\n')

    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': [{'SN': CHROM, 'LN': CHROM_LEN}]}
    tmp = os.path.join(outdir, 'sim.unsorted.bam')
    n_written = 0
    with pysam.AlignmentFile(tmp, 'wb', header=header) as bam:
        for scen, feat, strand, tss, tes, chain, n, trunc in models:
            pos = tx_positions(chain, strand)
            for i in range(n):
                kind = 'trunc' if random.random() < trunc else 'full'
                # 3' end: poly(A) site with a little jitter, clamped to the model
                tlen = len(pos)
                j3 = min(tlen - 1, max(tlen - 12, tlen - 1 + round(random.gauss(0, 4))))
                if kind == 'full':
                    sd, cap = (4, 15) if hard else (1.5, 5)
                    j5 = max(0, min(round(abs(random.gauss(0, sd))), cap))
                elif hard:
                    # 5' truncation: the read keeps an exponential length back from its
                    # 3' end, so long transcripts lose their 5' end more often (ONT cDNA)
                    length = max(150, min(j3 - 10, int(random.expovariate(1 / 1500))))
                    j5 = j3 - length
                else:
                    j5 = random.randrange(10, max(11, j3 - 150))
                five, three = pos[j5], pos[j3]
                blocks = blocks_for(chain, five, three, strand)
                qlen = sum(e - s for s, e in blocks)
                a = pysam.AlignedSegment()
                a.query_name = f'{scen}|{feat}|sense|{kind}|{i}'
                a.query_sequence = 'A' * qlen
                a.flag = 16 if strand == '-' else 0
                a.reference_id = 0
                a.reference_start = blocks[0][0]
                a.mapping_quality = 60
                a.cigartuples = cigar_from_blocks(blocks)
                bam.write(a)
                n_written += 1
    out_bam = os.path.join(outdir, 'sim.bam')
    pysam.sort('-o', out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    print(f'Wrote {n_written} reads to {out_bam}')
    return out_bam


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='sim_out')
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--hard', action='store_true')
    a = ap.parse_args()
    simulate(a.outdir, a.seed, a.hard)
