#!/usr/bin/env python3
"""Unit tests for tss_consensus.py. Plain asserts; run with python3, exits non-zero on failure."""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import numpy as np  # noqa: E402
import tss_consensus as tc  # noqa: E402

HW = {'ont': 8, 'ccs': 3, 'ont2': 8, 'tiny': 8, 'a': 8, 'b': 8}
failures = []


def check(name, cond, detail=''):
    print(f"{'ok  ' if cond else 'FAIL'} {name}{'' if cond else ': ' + str(detail)}")
    if not cond:
        failures.append(name)


def call(libs, strand='+', window=3, min_frac=0.5, **kw):
    return tc.call_tss(libs, {k: HW[k] for k in libs}, strand, window, min_frac, **kw)


# Depth invariance: scaling one library's counts leaves the call and score unchanged.
a = Counter({1000: 30, 1003: 10, 1200: 25})
b = Counter({1200: 20, 1201: 8, 1000: 5})
base = call({'ont': a, 'ont2': b})
scaled = call({'ont': Counter({k: v * 10 for k, v in a.items()}), 'ont2': b})
check('depth invariance: TSS1', base.sites[0][0] == scaled.sites[0][0],
      (base.sites[0], scaled.sites[0]))

# A 3-read library at a wrong site cannot move TSS1 against a 500-read library.
deep = Counter({5000: 300, 5002: 100, 5100: 100})
tiny = Counter({5100: 3})
r = call({'ont': deep, 'tiny': tiny})
check('tiny library cannot move TSS1', r.sites[0][0] in (5000, 5002), r.sites[0])

# Raw pooling would let a deep library win; the balanced call lets two agreeing libraries win.
r = call({'ont': Counter({100: 1000, 400: 900}), 'a': Counter({400: 50}), 'b': Counter({400: 40})})
check('agreeing libraries outvote one deep library', r.sites[0][0] == 400, r.sites[0])
check('Samples_Supporting counts agreeing libraries', (r.n_qualifying, r.n_supporting) == (3, 2),
      (r.n_qualifying, r.n_supporting))

# Ties resolve deterministically, upstream in the read's frame.
tie = Counter({200: 10, 300: 10})
check('tie on + strand -> lower coordinate', call({'ont': tie}).sites[0][0] == 200)
check('tie on - strand -> higher coordinate', call({'ont': tie}, strand='-').sites[0][0] == 300)

# TSS2 is a separate site, outside +/- h of TSS1.
r = call({'ont': Counter({1000: 20, 1001: 15, 1002: 12, 1500: 5})})
check('TSS2 outside +/- h of TSS1', abs(r.sites[1][0] - r.sites[0][0]) > 8, r.sites[:2])

# The reported site is a summit, not the centre of the densest window: two piles 10 bp apart put
# the best window between them, but TSS1 must land on the taller pile.
r = call({'ont': Counter({1000: 40, 1010: 30})})
check('TSS1 sits on a pile, not between piles', r.sites[0][0] == 1000 and r.sites[0][1] == 40,
      r.sites[0])

# One library: peak fraction equals the plain fraction of reads within the window.
c = Counter({50: 6, 51: 2, 60: 2})
r = call({'ont': c}, window=3)
check('single-library peak frac', abs(r.peak_frac - 8 / 10) < 1e-12, r.peak_frac)

# Mode method reproduces Counter.most_common, ties in insertion order.
pooled = Counter()
for x in [7, 3, 3, 7, 9]:
    pooled[x] += 1
m = tc.call_tss_mode(pooled, {'ont': pooled}, {'ont': 8}, '+', 3, 0.5)
check('mode: top site and tie order', [s[0] for s in m.sites] == [7, 3, 9], m.sites)

# Narrow platforms smooth less: ccs keeps two sites 5 bp apart distinct.
r = tc.call_tss({'ccs': Counter({700: 10, 705: 9})}, {'ccs': 3}, '+', 3, 0.5)
check('ccs half-width keeps sites 5 bp apart', [s[0] for s in r.sites[:2]] == [700, 705], r.sites)

# Low agreement flag.
r = call({'a': Counter({100: 20}), 'b': Counter({900: 20}), 'ont': Counter({900: 20})})
check('low agreement flagged', r.low_agreement == 0 and r.n_supporting == 2, r)
r = call({'a': Counter({100: 20}), 'b': Counter({500: 20}), 'ont': Counter({900: 21})})
check('low agreement when libraries disagree', r.low_agreement == 1, r)

check('halfwidth parsing', tc.parse_halfwidths('ont=5,nanopore_drna=12') ==
      {'ont': 5, 'pacbio_ccs': 3, 'pacbio_clr': 8, 'nanopore_drna': 12})

# --- external evidence ---
def track(name, pts):
    return tc.EvidenceTrack(name, [('c1', st, pos, w) for st, pos, w in pts])


# Reads pile hardest at 5000 (an RT stall, say); a smaller pile at 4000 is where Smar2C2 tags are.
reads = {'ont': Counter({5000: 100, 4000: 6, 4500: 1})}
smar = track('smar2c2', [('+', 4002, 20)])
cage = track('cage', [('+', 4500, 50), ('+', 5000, 50)])
r = call(reads)
check('no evidence: reads decide', r.sites[0][0] == 5000 and r.source == 'reads', r.sites[0])
r = call(reads, evidence=[smar], chrom='c1')
check('evidence moves TSS1 to the supported read pile', r.sites[0][0] == 4000
      and r.source == 'smar2c2' and r.reads_tss1 == 5000, (r.sites[0], r.source))
check('reads-only TSS1 kept as a later site', 5000 in [s[0] for s in r.sites], r.sites)
r = call(reads, evidence=[smar, cage], chrom='c1')
check('first track in priority order wins', r.source == 'smar2c2' and r.sites[0][0] == 4000, r)
r = call(reads, evidence=[track('smar2c2', [('+', 9000, 20)]), cage], chrom='c1')
check('falls through to the next track when the first has no support near reads',
      r.source == 'cage' and r.sites[0][0] == 5000, (r.source, r.sites[0]))
r = call(reads, evidence=[track('smar2c2', [('+', 4500, 20)])], chrom='c1')
check('a single-read position cannot be chosen', r.source == 'reads' and r.sites[0][0] == 5000,
      (r.source, r.sites[0]))
r = call(reads, evidence=[track('smar2c2', [('-', 4000, 20)])], chrom='c1')
check('evidence on the other strand is ignored', r.source == 'reads', r.source)
r = call(reads, evidence=[track('smar2c2', [('+', 4000, 2)])], chrom='c1')
check('evidence below ev_min is ignored', r.source == 'reads', r.source)
m = tc.call_tss_mode(reads['ont'], reads, {'ont': 8}, '+', 3, 0.5, evidence=[smar], chrom='c1')
check('mode method arbitrates too', m.sites[0][0] == 4000 and m.source == 'smar2c2', m.sites[:2])

import tempfile  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    open(f'{d}/t.txt', 'w').write('seq\tTSS\tstrand\tnTAGs\tisreal\nc1\t4001\t+\t7\tTRUE\n')
    open(f'{d}/c.gff3', 'w').write('c1\tlab\tTSS_cluster\t3990\t4010\t9.5\t-\t4003\tID=x\n')
    open(f'{d}/b.bed', 'w').write('c1\t4000\t4001\t+\n')
    t = tc.load_evidence(f'smar={d}/t.txt')
    check('TSS table parsed 1-based -> 0-based', t.support('c1', '+', np.array([4000]), 0)[0] == 7)
    g = tc.load_evidence(f'cage={d}/c.gff3')
    check('GFF dominant TSS and score', g.support('c1', '-', np.array([4002]), 0)[0] == 9.5)
    b = tc.load_evidence(f'bed={d}/t.txt,{d}/b.bed')
    check('paths pool, BED4 weight 1', b.support('c1', '+', np.array([4000]), 0)[0] == 8)

print(f"\n{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
