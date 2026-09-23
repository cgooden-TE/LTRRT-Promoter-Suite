#!/usr/bin/env python3
"""
Library-balanced consensus TSS calling for IsoClassifier.

A feature's TSS used to be the modal 5' end of every read pooled across every BAM, which lets the
deepest library decide it and breaks ties by read order. Here each library contributes the share of
the feature's reads starting near a position, weighted by how much evidence the library has for the
feature, so library depth cancels out and neither a tiny nor a very deep library can dominate.

For one feature/orientation with read-start counts c_s(x) per library s:

    c~_s(x) = starts of s within +/- h_s of x          (h_s set by platform: 5' end precision)
    p_s(x)  = c~_s(x) / n_s                            (n_s = the library's reads for the feature)
    w_s     = n_s / (n_s + k)                          (evidence weight, saturating at 1)
    Score(x) = sum_s w_s p_s(x) / sum_s w_s

The window with the highest Score picks the region; the site reported is that window's summit, the
position with the largest library-weighted share of exact 5' ends, so a TSS always sits on a pile of
reads rather than between two piles or on the flank of a broad one. Further sites are taken greedily
the same way, each summit outside +/- max(h_s) of every site already chosen, so TSS2 is a second
promoter rather than TSS1 +/- 1. Ties go to the higher raw pooled count, then to the most upstream
position in the read's frame (5' truncation only moves starts downstream).

Agreement is reported alongside: each library's own peak (argmax of c~_s), how many libraries with
enough reads have that peak within h_s of TSS1, and a flag when fewer than half do.

External TSS evidence (Smar2C2, CAGE, ...) can then arbitrate. Long-read 5' ends usually include a
pile at the true TSS but cannot say which pile it is: 5' truncation and RT stalls build piles too.
Evidence tracks are tried in priority order; the first one with enough support near any long-read
start position holding at least ev_min_reads reads picks, among those positions, the one it supports
most, and that becomes TSS1. With no track qualifying, the long-read call stands. TSS1 therefore
always sits where reads in these libraries start; evidence only chooses among those sites.

call_tss_mode() reproduces the original modal call exactly, for --tss-method mode.
"""
import gzip
import os
from collections import Counter, namedtuple

import numpy as np

PLATFORM_HALFWIDTH = {'ont': 8, 'pacbio_ccs': 3, 'pacbio_clr': 8}

# sites: [(position, raw pooled starts at that position, score or None), ...] best first
# per_lib: [(library, n_s, weight, own_peak, within_h), ...] in library order
# source: 'reads', or the name of the evidence track that chose TSS1; support: that track's weight
# within the evidence window of TSS1; reads_tss1: TSS1 as the reads alone would have called it
TSSCall = namedtuple('TSSCall', [
    'sites', 'peak_frac', 'called', 'n_qualifying', 'n_supporting', 'support_frac',
    'low_agreement', 'per_lib', 'source', 'support', 'reads_tss1'])


# ---------------------------------------------------------------------------------------------
# External TSS evidence
# ---------------------------------------------------------------------------------------------

class EvidenceTrack:
    """Weighted TSS positions (0-based) per (chrom, strand), with prefix sums for window support."""

    def __init__(self, name, points):
        self.name = name
        self.tracks = {}
        by = {}
        for chrom, strand, pos, w in points:
            by.setdefault((chrom, strand), []).append((pos, w))
        for key, lst in by.items():
            lst.sort()
            pos = np.fromiter((p for p, _ in lst), dtype=np.int64, count=len(lst))
            w = np.fromiter((x for _, x in lst), dtype=float, count=len(lst))
            self.tracks[key] = (pos, np.concatenate(([0.0], np.cumsum(w))))

    def support(self, chrom, strand, xs, half):
        t = self.tracks.get((chrom, strand))
        if t is None:
            return np.zeros(len(xs))
        pos, cum = t
        lo = np.searchsorted(pos, xs - half, side='left')
        hi = np.searchsorted(pos, xs + half, side='right')
        return cum[hi] - cum[lo]

    def __len__(self):
        return sum(len(p) for p, _c in self.tracks.values())


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _points_from_file(path):
    """
    TSS points from one file, format chosen by extension:
      .gff/.gff3 : CAGE-style clusters; column 8 is the dominant TSS (1-based) when numeric, else
                   the cluster midpoint; weight is column 6 (1 when '.')
      .bed       : 0-based start; strand is the first '+'/'-' among columns 4-6; weight is column 5
                   of a BED6 when numeric, else 1
      otherwise  : a TSS table with a header naming seq/chrom, TSS/pos, strand and optionally
                   nTAGs/count/score (TSRexplorer TSS sets), positions 1-based
    """
    base = path[:-3] if path.endswith('.gz') else path
    ext = os.path.splitext(base)[1].lower()
    with _open(path) as fh:
        if ext in ('.gff', '.gff3'):
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 8 or c[6] not in '+-':
                    continue
                pos = int(c[7]) - 1 if c[7].isdigit() else (int(c[3]) + int(c[4])) // 2 - 1
                w = float(c[5]) if c[5] not in ('.', '') else 1.0
                yield c[0], c[6], pos, w
        elif ext == '.bed':
            for line in fh:
                if line.startswith(('#', 'track', 'browser')) or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                strand = next((x for x in c[3:6] if x in ('+', '-')), None)
                if strand is None:
                    continue
                w = 1.0
                if len(c) >= 6:
                    try:
                        w = float(c[4])
                    except ValueError:
                        pass
                yield c[0], strand, int(c[1]), w
        else:
            header = fh.readline().split()
            low = [h.lower() for h in header]

            def col(*names):
                for n in names:
                    if n in low:
                        return low.index(n)
                return None
            ci, pi, si = col('seq', 'chrom', 'chr'), col('tss', 'pos', 'position'), col('strand')
            wi = col('ntags', 'count', 'score', 'tags')
            if None in (ci, pi, si):
                raise ValueError(f"{path}: TSS table needs seq/chrom, TSS/pos and strand columns")
            for line in fh:
                c = line.split()
                if len(c) <= max(ci, pi, si) or c[si] not in ('+', '-'):
                    continue
                yield c[ci], c[si], int(c[pi]) - 1, float(c[wi]) if wi is not None else 1.0


def load_evidence(spec):
    """'name=path[,path...]' -> EvidenceTrack; several paths (e.g. CAGE shoot and root) pool."""
    name, _, paths = spec.partition('=')
    if not paths:
        raise ValueError(f"--tss-evidence {spec!r}: expected name=path[,path...]")
    pts = []
    for path in paths.split(','):
        pts.extend(_points_from_file(path))
    return EvidenceTrack(name.strip(), pts)


def _arbitrate(sites, xs, raw, score, strand, excl, n, evidence, chrom, ev_window, ev_min,
               ev_min_reads):
    """
    Let evidence tracks, in priority order, pick TSS1 among start positions with at least
    ev_min_reads reads. Returns (sites, source, support). score is None for the mode call.
    """
    if not evidence or strand not in ('+', '-'):
        return sites, 'reads', None
    cand = np.nonzero(raw >= ev_min_reads)[0]
    if len(cand) == 0:
        return sites, 'reads', None
    cx = xs[cand]
    for track in evidence:
        sup = track.support(chrom, strand, cx, ev_window)
        if sup.max() < ev_min:
            continue
        j = min(range(len(cand)), key=lambda t: (
            -round(float(sup[t]), 9),
            -round(float(score[cand[t]]), 12) if score is not None else 0.0,
            -int(raw[cand[t]]), int(cx[t]) if strand != '-' else -int(cx[t])))
        i = cand[j]
        x = int(xs[i])
        first = (x, int(raw[i]), float(score[i]) if score is not None else None)
        rest = [s for s in sites if abs(s[0] - x) > excl]
        return [first] + rest[:n - 1], track.name, float(sup[j])
    return sites, 'reads', None


def parse_halfwidths(spec):
    """'ont=8,pacbio_ccs=3' -> dict, layered over the defaults."""
    out = dict(PLATFORM_HALFWIDTH)
    if spec:
        for item in spec.split(','):
            item = item.strip()
            if not item:
                continue
            name, _, val = item.partition('=')
            out[name.strip()] = int(val)
    return out


class _Lib:
    """Sorted start positions of one library for one feature, with prefix sums for window counts."""
    __slots__ = ('name', 'h', 'pos', 'cnt', 'cum', 'n')

    def __init__(self, name, h, counter):
        self.name = name
        self.h = h
        items = sorted(counter.items())
        self.pos = np.fromiter((p for p, _ in items), dtype=np.int64, count=len(items))
        self.cnt = np.fromiter((c for _, c in items), dtype=np.int64, count=len(items))
        self.cum = np.concatenate(([0], np.cumsum(self.cnt)))
        self.n = int(self.cum[-1])

    def window(self, xs, half):
        """Starts within +/- half of each x in xs."""
        lo = np.searchsorted(self.pos, xs - half, side='left')
        hi = np.searchsorted(self.pos, xs + half, side='right')
        return self.cum[hi] - self.cum[lo]


def _order_key(score, raw, x, strand):
    # Rounded so mathematically equal scores built from different terms still tie exactly.
    return (-round(score, 12), -raw, x if strand != '-' else -x)


def _own_peak(lib, strand):
    sm = lib.window(lib.pos, lib.h)
    best = min(range(len(lib.pos)),
               key=lambda i: (-int(sm[i]), -int(lib.cnt[i]),
                              int(lib.pos[i]) if strand != '-' else -int(lib.pos[i])))
    return int(lib.pos[best])


def _agreement(libs, tss1, strand, k, min_lib_reads):
    per_lib = []
    n_qual = n_sup = 0
    for lib in libs:
        own = _own_peak(lib, strand)
        within = abs(own - tss1) <= lib.h
        per_lib.append((lib.name, lib.n, lib.n / (lib.n + k), own, int(within)))
        if lib.n >= min_lib_reads:
            n_qual += 1
            n_sup += int(within)
    frac = n_sup / n_qual if n_qual else None
    low = int(n_qual >= 2 and frac < 0.5)
    return n_qual, n_sup, frac, low, per_lib


def _build_libs(lib_counts, lib_hw):
    libs = []
    for name in sorted(lib_counts):
        c = lib_counts[name]
        if c and sum(c.values()) > 0:
            libs.append(_Lib(name, lib_hw[name], c))
    return libs


def call_tss(lib_counts, lib_hw, strand, window, min_frac, n=10, k=10.0, min_lib_reads=5,
             evidence=None, chrom=None, ev_window=10, ev_min=3.0, ev_min_reads=2):
    """
    lib_counts : {library: Counter(position -> read starts)} for one feature/orientation
    lib_hw     : {library: half-width in bp}
    strand     : transcription strand of the reads ('+', '-' or '.'), for the upstream tie-break
    window, min_frac : TSS agreement window (bp, centred) and fraction for TSS_Called
    evidence   : EvidenceTracks in priority order (optional), with chrom and the ev_* thresholds
    Returns a TSSCall, or None when there are no reads.
    """
    libs = _build_libs(lib_counts, lib_hw)
    if not libs:
        return None
    xs = np.unique(np.concatenate([lib.pos for lib in libs]))
    raw = np.zeros(len(xs), dtype=np.int64)
    score = np.zeros(len(xs), dtype=float)     # smoothed: share within +/- h_s
    exact = np.zeros(len(xs), dtype=float)     # unsmoothed: share at exactly x
    wsum = 0.0
    for lib in libs:
        w = lib.n / (lib.n + k)
        wsum += w
        score += w * lib.window(xs, lib.h) / lib.n
        at_x = lib.window(xs, 0)
        exact += w * at_x / lib.n
        raw += at_x
    score /= wsum
    exact /= wsum

    order = sorted(range(len(xs)),
                   key=lambda i: _order_key(float(score[i]), int(raw[i]), int(xs[i]), strand))
    excl = max(lib.h for lib in libs)
    sites = []
    taken = set()
    for i in order:
        x = int(xs[i])
        # Summit of this window: the tallest exact position within +/- excl of its centre.
        lo = int(np.searchsorted(xs, x - excl, side='left'))
        hi = int(np.searchsorted(xs, x + excl, side='right'))
        j = min(range(lo, hi),
                key=lambda t: _order_key(float(exact[t]), int(raw[t]), int(xs[t]), strand))
        y = int(xs[j])
        if y in taken or any(abs(y - s) <= excl for s, _r, _sc in sites):
            continue
        taken.add(y)
        sites.append((y, int(raw[j]), float(score[j])))
        if len(sites) >= n:
            break

    reads_tss1 = sites[0][0]
    sites, source, support = _arbitrate(sites, xs, raw, score, strand, excl, n, evidence, chrom,
                                        ev_window, ev_min, ev_min_reads)
    tss1 = sites[0][0]
    half = max(0, (window - 1) // 2)
    x1 = np.array([tss1], dtype=np.int64)
    frac = sum((lib.n / (lib.n + k)) * int(lib.window(x1, half)[0]) / lib.n
               for lib in libs) / wsum
    n_qual, n_sup, sfrac, low, per_lib = _agreement(libs, tss1, strand, k, min_lib_reads)
    return TSSCall(sites, frac, int(frac >= min_frac), n_qual, n_sup, sfrac, low, per_lib,
                   source, support, reads_tss1)


def call_tss_mode(pooled, lib_counts, lib_hw, strand, window, min_frac, n=10, k=10.0,
                  min_lib_reads=5, evidence=None, chrom=None, ev_window=10, ev_min=3.0,
                  ev_min_reads=2):
    """
    The original call: modal 5' end of the pooled reads, top-n by Counter.most_common (ties in
    insertion order), and the unweighted fraction of reads within the window around the mode.
    pooled must be built in the same order the original read lists were, for identical ties.
    Agreement fields are still computed from the per-library profiles, as a diagnostic. Evidence,
    when given, arbitrates exactly as in call_tss.
    """
    if not pooled:
        return None
    total = sum(pooled.values())
    common = pooled.most_common(n)
    sites = [(p, c, None) for p, c in common]
    reads_tss1 = sites[0][0]
    if evidence:
        xs = np.fromiter(pooled.keys(), dtype=np.int64, count=len(pooled))
        raw = np.fromiter(pooled.values(), dtype=np.int64, count=len(pooled))
        order = np.argsort(xs, kind='stable')
        sites, source, support = _arbitrate(sites, xs[order], raw[order], None, strand, 0, n,
                                            evidence, chrom, ev_window, ev_min, ev_min_reads)
    else:
        source, support = 'reads', None
    peak = sites[0][0]
    half = max(0, (window - 1) // 2)
    cnt = sum(v for pos, v in pooled.items() if abs(pos - peak) <= half)
    frac = cnt / total
    libs = _build_libs(lib_counts, lib_hw)
    n_qual, n_sup, sfrac, low, per_lib = _agreement(libs, peak, strand, k, min_lib_reads)
    return TSSCall(sites, frac, int(frac >= min_frac), n_qual, n_sup, sfrac, low, per_lib,
                   source, support, reads_tss1)
