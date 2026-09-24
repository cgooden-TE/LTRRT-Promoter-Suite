#!/usr/bin/env python3
"""
Probabilistic origin assignment for reads whose 5' end lies inside more than one annotated feature
(a TE nested in a TE, a TE inside a gene, overlapping annotations). Used by IsoClassifier.py when
run with --assign em; the innermost-feature rule remains the default.

Model
-----
A transcription unit (TU) is a feature in one orientation, keyed exactly like IsoClassifier's
stats: (feature_id, 'sense' | 'antisense' | 'rs+' | 'rs-'). The candidate TUs for a read are the
features whose span contains the read's 5' end, i.e. the hits innermost_at() already sees.

Each TU u has an abundance N_u and generates a read with 5' end o and 3' end e with likelihood

    L(r | u) = P3_u(e) * P5_u(o | e) * P_splice(r | u) * P_junc(r | u)

P3_u(e):   where u's molecules end. A kernel-smoothed histogram of 3' ends (poly(A) sites) learned
           per TU, shrunk toward a base density that puts most mass on u's termination zone (3'
           LTR in the read's frame, the last exon for genes, the last term_bp otherwise), some on
           the rest of the body, and some on read-out past the end. Ends inside u's own sequence
           are learned from all of u's reads; ends beyond it (read-out) only from reads u
           demonstrably initiated, so a nested element cannot adopt its host's poly(A) site.
P5_u(o|e): where the read starts, given where it ends. A molecule initiated at TSS t and ending at
           e gives a full-length read (5' end at t) with probability S(l_te), the survival of the
           library's truncation process (RT drop-off, degradation) over the transcript length;
           otherwise the read is truncated and its 5' end sits at aligned length L upstream of e
           with density f(L). So

               P5_u(o | e) = sum_t w_ut * [ S(l_te) * K_t(o) + f(L) * 1(o downstream of t) ]

           over the TSS peaks t in u's promoter zone, with learned usage w_ut. S and f are shared
           by every unit and fitted by Kaplan-Meier (full-length reads are censored at their
           length), so initiation versus truncation is decided by read length and geometry, not
           by a free per-unit parameter that could drift. A useful consequence: a truncated read
           has the same 5'-end density f(L) under an inner and an outer element, so it is shared
           out by how many of each unit's molecules end at its 3' end, which is the true ratio.
TSS peaks: called from all 5' ends at the locus before assignment (Poisson test against local
           background), so "is there a TSS here?" does not depend on assignment. A unit may only
           initiate at peaks inside its promoter zone: both LTRs for LTR_structural, the first exon
           plus a window at the annotated 5' end for genes, the whole span for other classes. Spans
           of features nested inside u are removed from its zone, so an outer element can never
           claim a nested element's promoter. That exclusion is what makes nesting identifiable.
           Units with no peak in their zone fall back to truncation from an unknown start.
Genes:     truncated 5' ends in introns (pre-mRNA) are scaled by a learned intronic share eta_u.

Splicing and junctions: a learned spliced fraction per TU, and a learned junction-usage
distribution per TU whose prior puts mass on annotated introns (genes) and a little on everything
else, so a TU can learn its own recurrent junctions.

Fitting: EM over the TUs of each connected cluster. Reads whose 5' end lies in exactly one feature
of the cluster anchor that TU's abundance and shape; reads in several features are shared out.
TSS usage, 3'-end histograms, splicing and junction terms are evaluated leave-one-out (a read never
supports itself). S and f are re-estimated by Kaplan-Meier between passes, weighting each read as
censored (full-length) or observed (truncated) by its EM probability of being initiated. Clusters
still moving after a few plain iterations are accelerated with SQUAREM (SqS3), which extrapolates
the per-read state along the last two steps and keeps the jump only if the map's residual shrinks.

Outputs, per read: posterior over candidate TUs, and the probability that its 5' end is a TSS of
the assigned TU (rather than a truncation). Per TU: expected reads, expected TSS-initiated reads,
EM TSS peak, and the learned 3'/splicing parameters.
"""
import bisect
import math
from collections import Counter, defaultdict
from dataclasses import dataclass

TERM, BODY, READOUT = 0, 1, 2


@dataclass
class EMParams:
    tss_halfwidth: int = 8          # TSS peak half-width (bp); ONT 5' ends jitter several bp
    tss_bg_bp: int = 250            # flank used to estimate local 5'-end background for peak calls
    tss_min_reads: int = 3          # a TSS peak needs at least this many 5' ends in its window...
    tss_max_p: float = 1e-3         # ...and a Poisson tail probability below this vs background
    tss_bg_floor: float = 0.02      # floor on expected background reads per window
    tss_trunc_ratio: float = 20.0   # a peak must hold this many times the 5' ends that read-through
                                    # transcription plus RT drop-off would leave in its window
    tes_halfwidth: int = 10         # 3' kernel half-width (bp); poly(A) sites are fuzzier
    gene_promoter_bp: int = 300     # gene promoter window at the annotated 5' end
    term_bp: int = 500              # 3' window used as the termination zone for non-LTR classes
    term_ext_bp: int = 200          # termination zone extension past the annotated 3' end
    readout_bp: int = 10_000        # read-out region length for the 3' base density
    far_penalty: float = 1e-3       # extra factor for 3' ends beyond the read-out region
    base_3p: tuple = (0.6, 0.2, 0.2)   # base 3'-end mass: termination zone, body, read-out
    junc_tol: int = 5               # bp tolerance for matching annotated introns / pooling junctions
    max_iter: int = 500             # cap on E-steps (map evaluations), SQUAREM steps included
    tol: float = 1e-3               # stop when no TU's expected read count moves more than this
    accel: str = 'squarem'          # 'squarem' extrapolates the fixed-point map; 'none' is plain EM
    accel_warmup: int = 5           # plain E-steps before extrapolating; most clusters finish first
    accel_step_max0: float = 1.0    # SQUAREM step-length bounds (Varadhan & Roland 2008, SqS3)
    accel_mstep: float = 4.0
    trunc_min_sample: int = 2000    # a BAM needs this many reads for its own survival fit
    trunc_passes: int = 2           # EM passes; S and f are re-fitted (Kaplan-Meier) between passes
    eps_unknown_tss: float = 0.05   # share of a unit's reads allowed to start from an uncalled TSS
    eps_evidence_k: float = 5.0     # ...shrunk toward eps_min_frac of that by how little the unit
    eps_min_frac: float = 1e-3      # has shown it transcribes (certain reads + called peaks)
    eps_peak_weight: float = 5.0    # a called peak survived the truncation test, so it counts for
                                    # more than one certain read in that evidence
    alpha_prior: float = 0.5        # Dirichlet prior on a unit's abundance, so absorbed reads do
                                    # not make a unit more attractive to the next read
    prior_peak: float = 0.5         # pseudo-reads per called peak in a unit's TSS-usage weights
    prior_3p: float = 1.0           # pseudo-reads spread over the base 3' density
    prior_splice: tuple = (1.0, 1.0)
    prior_eta: tuple = (1.0, 9.0)      # gene 5'-truncation falling in introns (pre-mRNA)
    junc_prior_annot: float = 1.0
    junc_prior_novel: float = 0.5
    junc_prior_foreign: float = 0.05  # a junction another candidate annotates; cannot be learned
    junc_init_frac: float = 0.5       # a read may use a unit's learned read-out junctions only if
                                      # the model gives it at least this share of initiation there
    antisense_prior: float = 0.25     # prior odds that a read arose from antisense transcription
                                      # rather than sense at the same locus
    share_promoter_peaks: bool = True  # a child overlapping the host's own promoter (gene 5' end,
                                       # LTRs) does not get exclusive claim to peaks there
    junc_novel_space: float = 1e4
    alpha_floor: float = 1e-9


# ---------------------------------------------------------------------------------------------
# Interval helpers (half-open, sorted, non-overlapping tuples)
# ---------------------------------------------------------------------------------------------

def _merge(ivs):
    ivs = sorted((s, e) for s, e in ivs if e > s)
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return tuple((s, e) for s, e in out)


def _subtract(ivs, cuts):
    out = []
    cuts = _merge(cuts)
    for s, e in ivs:
        cur = s
        for cs, ce in cuts:
            if ce <= cur or cs >= e:
                continue
            if cs > cur:
                out.append((cur, cs))
            cur = max(cur, ce)
        if cur < e:
            out.append((cur, e))
    return _merge(out)


def _total(ivs):
    return sum(e - s for s, e in ivs)


def _contains(ivs, starts, x):
    i = bisect.bisect_right(starts, x) - 1
    return i >= 0 and ivs[i][0] <= x < ivs[i][1]


def _covered(ivs, lo, hi):
    """Bases of ivs inside [lo, hi)."""
    return sum(max(0, min(e, hi) - max(s, lo)) for s, e in ivs)


# ---------------------------------------------------------------------------------------------
# Global truncation-length density
# ---------------------------------------------------------------------------------------------

class TruncModel:
    """
    Library-wide read-length survival: S(L) = P(a read extends at least L aligned bases upstream
    of its 3' end) and its density f(L). A molecule whose TSS lies l bases upstream of its 3' end
    gives a full-length read with probability S(l); otherwise the read is truncated at length L
    with density f(L). Fitted by Kaplan-Meier in log-spaced bins: truncated reads are observed
    events at their length, full-length reads are censored at their length (all we know is that
    the library would have let them run further).
    """

    def __init__(self, edges, hazard):
        self.edges = edges
        self.h = hazard
        self.S_edges = [1.0]
        for h in hazard:
            self.S_edges.append(self.S_edges[-1] * (1.0 - h))

    @classmethod
    def fit_km(cls, lengths, w_event, w_cens, lo=50, hi=500_000, nbins=64, a0=1.0):
        step = (math.log(hi) - math.log(lo)) / nbins
        edges = [0.0] + [lo * math.exp(i * step) for i in range(nbins + 1)]
        nb = len(edges) - 1
        d = [0.0] * nb
        tot = [0.0] * nb
        for L, we, wc in zip(lengths, w_event, w_cens):
            i = min(nb - 1, max(0, bisect.bisect_right(edges, L) - 1))
            d[i] += we
            tot[i] += we + wc
        at_risk = [0.0] * nb
        acc = 0.0
        for i in range(nb - 1, -1, -1):
            acc += tot[i]
            at_risk[i] = acc
        hbar = sum(d) / max(sum(at_risk), 1e-12)
        hazard = [min(0.999, max(1e-6, (d[i] + a0 * hbar) / (at_risk[i] + a0))) for i in range(nb)]
        return cls(edges, hazard)

    def _bin(self, L):
        return min(len(self.h) - 1, max(0, bisect.bisect_right(self.edges, L) - 1))

    def pdf(self, L):
        i = self._bin(L)
        return self.S_edges[i] * self.h[i] / (self.edges[i + 1] - self.edges[i])

    def sf(self, L):
        i = self._bin(L)
        return max(self.S_edges[i + 1], self.S_edges[i] - self.pdf(L) * (L - self.edges[i]))

    def cdf(self, L):
        return 1.0 - self.sf(L)

    def median(self):
        for i in range(len(self.h)):
            if self.S_edges[i + 1] <= 0.5:
                f = self.S_edges[i] * self.h[i] / (self.edges[i + 1] - self.edges[i])
                return self.edges[i] + (self.S_edges[i] - 0.5) / f
        return self.edges[-1]

    def table(self):
        return [(self.edges[i], self.edges[i + 1], self.S_edges[i], self.h[i])
                for i in range(len(self.h))]

    def hazard_per_bp(self):
        """(lo, hi, per-base hazard) for each fitted bin. Used to ask how often RT would stop
        inside a given window, for reads that were transcribed through it."""
        return [(self.edges[i], self.edges[i + 1],
                 self.h[i] / max(1e-9, self.edges[i + 1] - self.edges[i]))
                for i in range(len(self.h))]


class TruncSet:
    """The library-wide survival model, plus one per BAM when several are pooled and each has
    enough reads to fit. Peak calling uses the pooled model -- a 5'-end pile-up is a property of
    the locus, not of one library -- while the per-read 5' terms use the read's own sample, so
    libraries that differ in RT processivity are not forced to share one curve."""

    def __init__(self, pooled, by_sample=None):
        self.pooled = pooled
        self.by_sample = by_sample or {}

    def get(self, sample):
        return self.by_sample.get(sample, self.pooled)

    def table(self):
        return self.pooled.table()

    def median(self):
        return self.pooled.median()


# ---------------------------------------------------------------------------------------------
# Transcription-unit specification
# ---------------------------------------------------------------------------------------------

def unit_strand(rec, orient):
    if orient == 'sense':
        return rec['strand']
    if orient == 'antisense':
        return '-' if rec['strand'] == '+' else '+'
    return '+' if orient == 'rs+' else '-'


def build_unit(key, rec, child_spans, gene_introns, p):
    """Static description of one TU. child_spans: spans of features strictly nested in rec."""
    fid, orient = key
    strand = unit_strand(rec, orient)
    start, end = rec['start'], rec['end']
    span = ((start, end),)
    klass = rec['class']
    exons = rec.get('exons') if klass == 'Gene' else None
    sense_like = orient in ('sense', 'rs+', 'rs-') and strand in ('+', '-')

    # promoter (initiation) zone
    if klass == 'LTR_structural':
        zone = _merge([rec['ltr_left'], rec['ltr_right']])
    elif klass == 'Gene' and orient == 'sense':
        w = p.gene_promoter_bp
        win = (start, min(end, start + w)) if strand == '+' else (max(start, end - w), end)
        first = [exons[0] if strand == '+' else exons[-1]] if exons else []
        zone = _merge([win] + first)
    else:
        zone = span
    # A child inside the host's BODY gets exclusive claim to its own promoter; that exclusion is
    # what makes nesting identifiable. But where the host has a promoter of its own -- a gene's 5'
    # end, an element's LTRs -- a child overlapping THAT is as likely to be an annotation artefact
    # sitting on the host's TSS as a real nested promoter, and cutting it out leaves the host with
    # no promoter at all. There both keep the peak and the rest of the evidence decides.
    if not (p.share_promoter_peaks and zone != span):
        zone = _subtract(zone, child_spans)

    # termination zone (for the base 3' density)
    ext = p.term_ext_bp
    if klass == 'LTR_structural':
        term = (rec['ltr_right'] if strand == '+' else rec['ltr_left'],)
    elif exons and orient == 'sense':
        term = ((exons[-1][0], end + ext),) if strand == '+' else ((start - ext, exons[0][1]),)
    else:
        w = p.term_bp
        term = ((max(start, end - w), end + ext),) if strand == '+' else \
               ((start - ext, min(end, start + w)),)
    term = _merge(term)
    in_span_term = _covered(term, start, end)
    exon_len = _total(exons) if exons else 0
    return {
        'key': key, 'fid': fid, 'orient': orient, 'klass': klass, 'strand': strand,
        'start': start, 'end': end, 'span_len': max(1, end - start),
        'zone': zone, 'zone_starts': [s for s, _ in zone], 'zone_len': _total(zone),
        'term': term, 'term_starts': [s for s, _ in term], 'term_len': max(1, _total(term)),
        'body_len': max(1, (end - start) - in_span_term),
        'readout_hi': max(end, term[-1][1]), 'readout_lo': min(start, term[0][0]),
        'exons': exons, 'exon_starts': [s for s, _ in exons] if exons else None,
        'exon_len': exon_len,
        # A gene's annotated introns belong to the unit transcribed in the gene's own direction.
        # sense_like used to include rs+/rs-, so the unit reading a + gene backwards was handed
        # that gene's intron chain and could score its own reads as annotated.
        'annot': (frozenset(gene_introns or ())
                  if (klass == 'Gene' and sense_like and strand == rec['strand'])
                  else frozenset()),
    }


# ---------------------------------------------------------------------------------------------
# Locus-level TSS peak calling
# ---------------------------------------------------------------------------------------------

def _poisson_sf(n, mu):
    """P(X >= n) for X ~ Poisson(mu)."""
    if n <= 0:
        return 1.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, n):
        term *= mu / i
        cdf += term
    return max(0.0, 1.0 - cdf)


def _truncation_expectation(cands, reads, strand, p, trunc):
    """
    Keep only candidate peaks that read-through transcription cannot account for.

    A read whose 5' end lies upstream of a window and whose 3' end lies downstream of it was
    transcribed through that position, so the library's own hazard says how often reverse
    transcription would have stopped inside the window instead of running on. Summed over those
    reads, that is the number of 5' ends the window is expected to collect for free. A genuine
    TSS clears it by orders of magnitude; an RT stall site inside a transcribed feature does not,
    however sharp it looks against a flat background.

    This is what keeps a truncation hotspot from being sold to whichever unit happens to have the
    position in its promoter zone -- typically an antisense unit of an overlapping gene, or a TE
    sitting in an intron, neither of which transcribed anything.
    """
    kw = 2 * p.tss_halfwidth + 1
    bins = trunc.hazard_per_bp()
    fwd = strand == '+'
    near = p.tss_halfwidth + 1
    # Walk candidates along the direction of transcription, accumulating the 3' ends of reads
    # already known to start upstream of the current window.
    order = sorted(cands, key=lambda t: t[1], reverse=not fwd)
    by_origin = sorted(reads, key=lambda r: r[0], reverse=not fwd)
    ends, i, keep = [], 0, []
    for w, x in order:
        edge = x - p.tss_halfwidth if fwd else x + p.tss_halfwidth
        while i < len(by_origin) and ((by_origin[i][0] < edge) if fwd else (by_origin[i][0] > edge)):
            bisect.insort(ends, by_origin[i][1])
            i += 1
        mu = 0.0
        for lo, hi, h in bins:
            if fwd:
                n = (bisect.bisect_right(ends, x + hi)
                     - bisect.bisect_left(ends, x + max(lo, near)))
            else:
                n = (bisect.bisect_right(ends, x - max(lo, near))
                     - bisect.bisect_left(ends, x - hi))
            if n:
                mu += n * h
        if w >= p.tss_trunc_ratio * mu * kw:
            keep.append((w, x))
    return keep


def call_tss_peaks(reads, p, strand='+', trunc=None):
    """
    Call TSS peaks from the (5' end, 3' end) pairs of the reads on one strand, regardless of which
    unit they belong to. A peak is a window of +/- tss_halfwidth holding >= tss_min_reads 5' ends,
    significantly more than the local background (Poisson), and -- once a TruncModel is available
    -- more than read-through transcription plus RT drop-off would put there anyway. Returns a
    sorted list of (lo, hi, center) windows.

    Peaks are called before assignment so that "is there a TSS here?" is answered by the reads
    alone, and "whose TSS is it?" by which unit's promoter zone contains it. Without this, two
    stray 5' ends that happen to coincide would support each other as a TSS in EM.
    """
    if not reads:
        return []
    k, bg = p.tss_halfwidth, p.tss_bg_bp
    kw = 2 * k + 1
    pos = sorted(o for o, _e in reads)
    uniq = sorted(set(pos))

    def count(lo, hi):
        return bisect.bisect_right(pos, hi) - bisect.bisect_left(pos, lo)

    cands = []
    for x in uniq:
        w = count(x - k, x + k)
        if w < p.tss_min_reads:
            continue
        b = count(x - bg, x + bg) - w
        mu = max(p.tss_bg_floor, b / max(1, (2 * bg + 1 - kw)) * kw)
        if _poisson_sf(w, mu) < p.tss_max_p:
            cands.append((w, x))
    if trunc is not None and cands:
        cands = _truncation_expectation(cands, reads, strand, p, trunc)
    peaks = []
    for w, x in sorted(cands, key=lambda t: (-t[0], t[1])):
        if all(abs(x - c) >= kw for _lo, _hi, c in peaks):
            peaks.append((x - k, x + k, x))
    return sorted(peaks)


# ---------------------------------------------------------------------------------------------
# Per (read, unit) static features
# ---------------------------------------------------------------------------------------------

def _upstream(strand, a, b):
    """True if position a is at or upstream of b in the direction of transcription."""
    return a <= b if strand == '+' else a >= b


def _tx_len(u, a, b):
    """Transcript length between genomic positions a and b (inclusive) within unit u."""
    lo, hi = min(a, b), max(a, b) + 1
    if u['exons']:
        beyond = max(0, hi - u['end']) + max(0, u['start'] - lo)
        return max(1, _covered(u['exons'], lo, hi) + beyond)
    return hi - lo


def _pair_features(u, origin, endpos, alen, juncs, peaks, p, cand_annot=frozenset()):
    """Static (assignment-independent) description of one read against one candidate unit.
    cand_annot: annotated introns of every candidate for this read, used to tell a junction this
    unit simply has not seen before from one that demonstrably belongs to a rival."""
    k = p.tss_halfwidth
    st = u['strand']
    # TSS peaks this unit could have initiated the read from: called peak, centre in the unit's
    # promoter zone (nested features already removed), and upstream of the read's 3' end.
    usable = []
    for t, (lo, hi, c) in enumerate(peaks):
        if u['zone_len'] == 0 or not _contains(u['zone'], u['zone_starts'], c):
            continue
        if not _upstream(st, c, endpos):
            continue
        in_window = lo <= origin <= hi
        downstream = _upstream(st, c - k if st == '+' else c + k, origin)
        if not (in_window or downstream):
            continue           # read starts upstream of this TSS: it cannot come from it
        l_te = _tx_len(u, c, endpos)
        if l_te < alen:        # read carries sequence the spliced transcript lacks (pre-mRNA)
            l_te = abs(endpos - c) + 1
        usable.append((t, in_window, downstream, l_te))
    intronic = bool(u['exons']) and not _contains(u['exons'], u['exon_starts'], origin)

    # base 3' density at endpos, and whether endpos lies within the unit's own sequence
    b_term, b_body, b_out = p.base_3p
    in_own = u['start'] <= endpos < u['end'] or _contains(u['term'], u['term_starts'], endpos)
    if _contains(u['term'], u['term_starts'], endpos):
        cat, base = TERM, b_term / u['term_len']
    elif (st == '+' and endpos >= u['readout_hi']) or (st == '-' and endpos < u['readout_lo']):
        dist = endpos - u['readout_hi'] if st == '+' else u['readout_lo'] - endpos
        cat, base = READOUT, b_out / p.readout_bp * (p.far_penalty if dist > p.readout_bp else 1.0)
    else:
        cat, base = BODY, b_body / u['body_len']

    # longest read that fits between the unit's 5' boundary and endpos (unknown-TSS fallback)
    bound = u['start'] if st == '+' else u['end'] - 1
    lmax = max(abs(endpos - bound) + 1, alen)

    jkeys = []
    tol = p.junc_tol
    for d, a in juncs:
        hit = None
        for ad, aa in u['annot']:
            if abs(ad - d) <= tol and abs(aa - a) <= tol:
                hit = (ad, aa)
                break
        if hit:
            jkeys.append((hit, 'annot'))
            continue
        # A junction another candidate annotates is evidence for that candidate. Tagging it keeps
        # this unit from learning it (M-step) and from pricing it as merely unfamiliar (E-step).
        foreign = any(abs(ad - d) <= tol and abs(aa - a) <= tol for ad, aa in cand_annot)
        jkeys.append(((d // tol, a // tol), 'foreign' if foreign else 'novel'))
    return tuple(usable), intronic, cat, base, in_own, lmax, tuple(jkeys)


def _cluster_peaks(U, idx, records, p, trunc=None):
    by_strand = defaultdict(list)
    for rec in records:
        by_strand[U[idx[rec[4][0]]]['strand']].append((rec[0], rec[1]))
    return {st: call_tss_peaks(r, p, st, trunc) for st, r in by_strand.items()}


def initial_censoring(units, records, p):
    """Pass-0 guess for Kaplan-Meier: a read is 'full-length' if it starts in a called TSS peak
    that lies in one of its candidates' promoter zones."""
    keys = sorted(units)
    idx = {k: i for i, k in enumerate(keys)}
    U = [units[k] for k in keys]
    peaks = _cluster_peaks(U, idx, records, p)
    out = []
    for origin, endpos, _sp, _j, cands, alen, _sa in records:
        u0 = U[idx[cands[0]]]
        pk = peaks.get(u0['strand'], [])
        full = False
        for lo, hi, c in pk:
            if lo <= origin <= hi and any(
                    U[idx[cc]]['zone_len'] and _contains(U[idx[cc]]['zone'],
                                                         U[idx[cc]]['zone_starts'], c)
                    for cc in cands):
                full = True
                break
        out.append(1.0 if full else 0.0)
    return out


# ---------------------------------------------------------------------------------------------
# EM on one cluster
# ---------------------------------------------------------------------------------------------

def em_cluster(units, records, p, trunc):
    """
    units   : dict key -> unit spec (build_unit)
    records : list of (origin, endpos, spliced, juncs, cand_keys, aligned_len)
    trunc   : TruncSet (pooled S and f, plus per-sample models when several BAMs are pooled)
    Returns (per_record, summary):
      per_record[i] = [(key, posterior, posterior_full_length), ...] over the read's candidates
      summary[key]  = dict of fitted unit parameters
    """
    keys = sorted(units)
    idx = {k: i for i, k in enumerate(keys)}
    U = [units[k] for k in keys]
    n_u = len(U)
    k3 = p.tes_halfwidth
    kw5, kw3 = 2 * p.tss_halfwidth + 1, 2 * k3 + 1
    a_s, b_s = p.prior_splice
    a_e, b_e = p.prior_eta
    peaks = _cluster_peaks(U, idx, records, p, trunc.pooled)
    n_peaks = [0] * n_u
    for j, u in enumerate(U):
        n_peaks[j] = sum(1 for _lo, _hi, c in peaks.get(u['strand'], [])
                         if u['zone_len'] > 0 and _contains(u['zone'], u['zone_starts'], c))

    # Reads with a single candidate anchor a unit independently of anything the EM decides.
    n_single = [0] * n_u
    for _o, _e, _sp, _ju, cands, _al, _sa in records:
        if len(cands) == 1:
            n_single[idx[cands[0]]] += 1

    # Every unit may claim reads that started at a TSS it has no called peak for, but only in
    # proportion to its evidence of transcribing at all. Without this a unit with no peaks and no
    # certain reads gets the widest unknown-start channel of any unit in the cluster, which is
    # backwards: it is exactly the unit that has shown nothing.
    eps_u = []
    for j in range(n_u):
        ev = n_single[j] + p.eps_peak_weight * n_peaks[j]
        eps_u.append(p.eps_unknown_tss * max(ev / (ev + p.eps_evidence_k), p.eps_min_frac))

    # per pair: static features plus the S/f numbers they imply under this pass's TruncModel
    R = []
    annot_cache = {}
    for origin, endpos, spliced, juncs, cands, alen, sample in records:
        ci = [idx[c] for c in cands]
        st = U[ci[0]]['strand']
        tr = trunc.get(sample)
        f_len = tr.pdf(alen)
        ck = tuple(ci)
        cand_annot = annot_cache.get(ck)
        if cand_annot is None:
            cand_annot = frozenset().union(*(U[j]['annot'] for j in ci))
            annot_cache[ck] = cand_annot
        pairs = []
        for j in ci:
            usable, intronic, cat, base3, in_own, lmax, jkeys = _pair_features(
                U[j], origin, endpos, alen, juncs, peaks.get(st, []), p, cand_annot)
            ab = tuple((t, (tr.sf(l_te) / kw5) if inw else 0.0, f_len if dn else 0.0)
                       for t, inw, dn, l_te in usable)
            fb = f_len / max(tr.cdf(lmax), 1e-12)
            pairs.append((ab, intronic, cat, base3, in_own, fb, jkeys))
        R.append((origin, endpos, bool(spliced), ci, pairs))

    # state per (read, candidate): posterior g, full-length part gi, per-peak responsibilities
    G = [[1.0 / len(r[3])] * len(r[3]) for r in R]
    # Warm start: the posterior initiation gate (junc_init_frac) has nothing to read on the first
    # iteration, and starting everything at "not initiated" makes genuine read-out transcripts fail
    # the junction test before the EM has seen anything. Seed it from geometry -- a 5' end inside one
    # of the unit's peak windows -- and let the first E-step correct it.
    GI = [[(g if any(aa > 0 for _t, aa, _bb in pr[0]) else 0.0)
           for pr, g in zip(r[4], gs)] for r, gs in zip(R, G)]
    RESP = [[tuple(1.0 / len(pr[0]) for _ in pr[0]) if pr[0] else () for pr in r[4]] for r in R]

    def m_step(G, GI, RESP):
        S = {'N': [0.0] * n_u, 'I': [0.0] * n_u, 'Intr': [0.0] * n_u, 'Sp': [0.0] * n_u,
             'JC': [0.0] * n_u, 'M3': [0.0] * n_u, 'W': [0.0] * n_u,
             'C3': [[0.0, 0.0, 0.0] for _ in range(n_u)],
             'Jc': [defaultdict(float) for _ in range(n_u)],
             'Wpk': [defaultdict(float) for _ in range(n_u)],
             'Ipk': [defaultdict(float) for _ in range(n_u)],
             'H3': [defaultdict(float) for _ in range(n_u)]}
        for (origin, endpos, spliced, ci, pairs), gs, gis, rs in zip(R, G, GI, RESP):
            for j, pr, g, gi, resp in zip(ci, pairs, gs, gis, rs):
                ab, intronic, cat, _b, in_own, _fb, jkeys = pr
                S['N'][j] += g
                S['I'][j] += gi
                for (t, a, bb), rr in zip(ab, resp):
                    S['Wpk'][j][t] += g * rr
                    S['W'][j] += g * rr
                if gi > 0 and ab:
                    tot_a = sum(a for _t, a, _bb in ab) or 1.0
                    for t, a, _bb in ab:
                        if a:
                            S['Ipk'][j][t] += gi * a / tot_a
                # 3' ends inside the unit's own sequence are learned from all its reads; ends
                # beyond it (read-out) only from reads it demonstrably initiated
                w3 = g if in_own else gi
                if w3 > 0:
                    S['H3'][j][endpos] += w3
                    S['M3'][j] += w3
                S['C3'][j][cat] += g
                if intronic:
                    S['Intr'][j] += g
                if spliced:
                    S['Sp'][j] += g
                    # Same rule as the E-step: a unit may learn a rival's intron only from reads it
                    # demonstrably initiated (read-out), never from reads it explains as truncations.
                    init_here = (any(a > 0 for _t, a, _bb in ab)
                                 and g > 0.0 and gi >= p.junc_init_frac * g)
                    n_learn = 0
                    for jk, tag in jkeys:
                        if tag == 'foreign' and not init_here:
                            continue
                        S['Jc'][j][jk] += g
                        n_learn += 1
                    S['JC'][j] += g * n_learn
        return S

    def window(H, x, k):
        return sum(H.get(x + d, 0.0) for d in range(-k, k + 1))

    def pair_terms(S, j, origin, endpos, spliced, pr, g, gi, resp, full):
        """Leave-one-out log-likelihood of the read under unit j, P(full-length | read, j), and
        the read's responsibilities over j's usable TSS peaks."""
        u = U[j]
        eps = eps_u[j]
        ab, intronic, _cat, base3, in_own, fb, jkeys = pr
        N_ = max(S['N'][j] - g, 0.0)
        # pre-mRNA share for genes scales truncated 5' ends that fall in introns
        if u['exons']:
            eta = (max(S['Intr'][j] - (g if intronic else 0.0), 0.0) + a_e) / (N_ + a_e + b_e)
            gf = eta if intronic else 1.0
        else:
            gf = 1.0
        if n_peaks[j]:
            own = dict(zip((t for t, _a, _b in ab), resp))
            W_ = max(S['W'][j] - g * sum(resp), 0.0)
            denom = W_ + p.prior_peak * n_peaks[j]
            terms = []
            for t, a, bb in ab:
                w = (max(S['Wpk'][j].get(t, 0.0) - g * own.get(t, 0.0), 0.0) + p.prior_peak) / denom
                terms.append((w * a, w * bb * gf))
            s_full = sum(x for x, _y in terms)
            s_all = sum(x + y for x, y in terms)
            p5 = (1.0 - eps) * s_all + eps * fb * gf
            frac = (1.0 - eps) * s_full / p5 if p5 > 0 else 0.0
            new_resp = tuple((x + y) / s_all for x, y in terms) if s_all > 0 else \
                tuple(0.0 for _ in terms)
            # responsibilities are over peaks; the unknown-TSS share is left out of them
            share = (1.0 - eps) * s_all / p5 if p5 > 0 else 0.0
            new_resp = tuple(r * share for r in new_resp)
        else:
            # No called peak in this unit's zone: every read it claims has to come through the
            # unknown-start channel, which is only as wide as the unit's evidence. Giving it the
            # undiscounted density here made having no evidence worth 1/eps over having some.
            p5 = eps * fb * gf
            frac, new_resp = 0.0, ()
        if not full:
            return 0.0, frac, new_resp
        ll = math.log(max(p5, 1e-300))
        own3 = g if in_own else gi
        M_ = max(S['M3'][j] - own3, 0.0)
        win3 = max(window(S['H3'][j], endpos, k3) - own3, 0.0)
        p3 = (win3 / kw3 + p.prior_3p * base3) / (M_ + p.prior_3p)
        ll += math.log(max(p3, 1e-300))
        sp = (max(S['Sp'][j] - (g if spliced else 0.0), 0.0) + a_s) / (N_ + a_s + b_s)
        ll += math.log(sp if spliced else 1.0 - sp)
        if spliced and jkeys:
            n_ann = len(u['annot'])
            # A read that initiated here is reading out, so junctions it picks up downstream say
            # nothing about its origin. "Initiated" has to be the model's own verdict (gi/g), not
            # merely a 5' end landing in a peak window: a truncation hotspot collects 5' ends too,
            # and geometry alone let those reads teach the unit its neighbour's introns.
            has_init = (any(a > 0 for _t, a, _bb in ab)
                        and g > 0.0 and gi >= p.junc_init_frac * g)
            mult = Counter(jk for jk, tag in jkeys if tag != 'foreign' or has_init)
            denom_j = max(S['JC'][j] - g * sum(mult.values()), 0.0) + p.junc_prior_novel + \
                (p.junc_prior_annot if n_ann else 0.0)
            Jc = S['Jc'][j]
            for jk, tag in jkeys:
                # ...and a unit that has already shown read-out through this junction, on reads it
                # did initiate, keeps it for its truncated reads too. A unit with no called peak
                # never learns one, so a phantom can never buy its way out of the floor.
                if tag == 'foreign' and not has_init and Jc.get(jk, 0.0) <= 0.0:
                    # A rival candidate annotates this intron. Fixed floor, never learned away.
                    ll += math.log(p.junc_prior_foreign / p.junc_novel_space / denom_j)
                    continue
                base = (p.junc_prior_annot / n_ann if tag == 'annot'
                        else p.junc_prior_novel / p.junc_novel_space)
                ll += math.log((max(Jc.get(jk, 0.0) - g * mult[jk], 0.0) + base) / denom_j)
        # Antisense transcription is real but rarer than sense at the same locus, and nothing else
        # in the likelihood says so -- an antisense unit competes on equal terms for every truncated
        # read whose 5' end it happens to cover.
        if u['orient'] == 'antisense':
            ll += math.log(p.antisense_prior)
        # Leave-one-out and Dirichlet-smoothed, like every other learned term: a unit's own
        # absorbed mass must not be what makes it attractive to the next read.
        ll += math.log(max(S['N'][j] - g, 0.0) + p.alpha_prior)
        return ll, frac, new_resp

    def e_step(G, GI, RESP, S):
        newG, newGI, newR = [], [], []
        for (origin, endpos, spliced, ci, pairs), gs, gis, rs in zip(R, G, GI, RESP):
            if len(ci) == 1:
                _, fi, rr = pair_terms(S, ci[0], origin, endpos, spliced, pairs[0],
                                       gs[0], gis[0], rs[0], False)
                newG.append([1.0])
                newGI.append([fi])
                newR.append([rr])
                continue
            lls, fis, rrs = [], [], []
            for j, pr, g, gi, resp in zip(ci, pairs, gs, gis, rs):
                ll, fi, rr = pair_terms(S, j, origin, endpos, spliced, pr, g, gi, resp, True)
                lls.append(ll)
                fis.append(fi)
                rrs.append(rr)
            m = max(lls)
            w = [math.exp(x - m) for x in lls]
            z = sum(w)
            post = [x / z for x in w]
            newG.append(post)
            newGI.append([pp * fi for pp, fi in zip(post, fis)])
            newR.append(rrs)
        return newG, newGI, newR

    def n_delta(Sa, Sb):
        return max(abs(a - b) for a, b in zip(Sa['N'], Sb['N'])) if n_u else 0.0

    def extrapolate(s0, s1, s2, step_max):
        """SqS3 step on the per-read state: x0 + 2a*r + a^2*v, r = x1 - x0, v = x2 - 2*x1 + x0,
        a = |r|/|v| bounded to [1, step_max] (a = 1 gives back x2). The result is projected back
        onto valid posteriors: each read's G on the simplex, GI within [0, G], peak
        responsibilities non-negative and summing to at most 1."""
        sr = sv = 0.0
        for c in range(3):
            for r0, r1, r2 in zip(s0[c], s1[c], s2[c]):
                for a0, a1, a2 in zip(r0, r1, r2):
                    if c == 2:
                        for b0, b1, b2 in zip(a0, a1, a2):
                            sr += (b1 - b0) ** 2
                            sv += (b2 - 2 * b1 + b0) ** 2
                    else:
                        sr += (a1 - a0) ** 2
                        sv += (a2 - 2 * a1 + a0) ** 2
        if sv <= 0.0 or sr <= 0.0:
            return None, 1.0
        a = min(max(math.sqrt(sr / sv), 1.0), step_max)
        c1, c2 = 2.0 * a, a * a

        def ext(x0, x1, x2):
            return x0 + c1 * (x1 - x0) + c2 * (x2 - 2 * x1 + x0)

        G, GI, RESP = [], [], []
        for g0, g1, g2, i0, i1, i2, q0, q1, q2 in zip(s0[0], s1[0], s2[0], s0[1], s1[1], s2[1],
                                                      s0[2], s1[2], s2[2]):
            if len(g0) == 1:
                g = [1.0]
            else:
                g = [min(max(ext(x0, x1, x2), 0.0), 1.0) for x0, x1, x2 in zip(g0, g1, g2)]
                z = sum(g)
                g = [x / z for x in g] if z > 0 else [1.0 / len(g)] * len(g)
            gi = [min(max(ext(x0, x1, x2), 0.0), gg) for x0, x1, x2, gg in zip(i0, i1, i2, g)]
            rs = []
            for p0, p1, p2 in zip(q0, q1, q2):
                rr = [min(max(ext(x0, x1, x2), 0.0), 1.0) for x0, x1, x2 in zip(p0, p1, p2)]
                z = sum(rr)
                rs.append(tuple(x / z for x in rr) if z > 1.0 else tuple(rr))
            G.append(g)
            GI.append(gi)
            RESP.append(rs)
        return (G, GI, RESP), a

    state = (G, GI, RESP)
    S = m_step(*state)
    n_iter = 0
    step_max = p.accel_step_max0
    while n_iter < p.max_iter:
        s1 = e_step(*state, S)
        S1 = m_step(*s1)
        n_iter += 1
        d1 = n_delta(S1, S)
        if d1 < p.tol or p.accel != 'squarem' or n_iter < p.accel_warmup \
                or n_iter >= p.max_iter:
            state, S = s1, S1
            if d1 < p.tol:
                break
            continue
        # SQUAREM cycle: two map evaluations, an extrapolated jump, and a stabilising map step.
        # The map (leave-one-out E-step + M-step) has no objective to check a jump against, so the
        # merit is its residual: a jump is kept only if the step taken from it moves N less than
        # the first plain step of the cycle did; otherwise the cycle falls back to x2.
        s2 = e_step(*s1, S1)
        S2 = m_step(*s2)
        n_iter += 1
        d2 = n_delta(S2, S1)
        if d2 < p.tol or n_iter >= p.max_iter:
            state, S = s2, S2
            if d2 < p.tol:
                break
            continue
        sx, a = extrapolate(state, s1, s2, step_max)
        if sx is None or a <= 1.0:
            # a = 1 is the plain step x2; at the bound, allow a longer jump next cycle
            state, S = s2, S2
            if sx is not None and a >= step_max:
                step_max *= p.accel_mstep
            continue
        Sx = m_step(*sx)
        s4 = e_step(*sx, Sx)
        S4 = m_step(*s4)
        n_iter += 1
        d4 = n_delta(S4, Sx)
        if math.isfinite(d4) and d4 <= d1:
            state, S = s4, S4
            if a >= step_max:
                step_max *= p.accel_mstep
            if d4 < p.tol:
                break
        else:
            state, S = s2, S2
            step_max = max(p.accel_step_max0, step_max / p.accel_mstep)
    G, GI, RESP = state

    per_record = [[(keys[j], g, gi) for j, g, gi in zip(r[3], gs, gis)]
                  for r, gs, gis in zip(R, G, GI)]

    summary = {}
    for j, u in enumerate(U):
        N, I = S['N'][j], S['I'][j]
        Ipk, H3 = S['Ipk'][j], S['H3'][j]
        pk = peaks.get(u['strand'], [])
        peak, peak_mass, tes, tes_mass = 'NA', 0.0, 'NA', 0.0
        if Ipk and I >= 0.5:
            best = max(Ipk, key=lambda t: Ipk[t])
            peak, peak_mass = pk[best][2], Ipk[best] / I
        if H3 and S['M3'][j] >= 0.5:
            tes = max(H3, key=lambda x: (window(H3, x, k3), -x))
            tes_mass = window(H3, tes, k3) / S['M3'][j]
        c3 = S['C3'][j]
        summary[u['key']] = {
            'EM_Reads': N, 'EM_TSS_Reads': I, 'Unique_Reads': n_single[j],
            'TSS_Peaks_In_Zone': n_peaks[j],
            'Full_Length_Frac': I / N if N > 0 else 0.0,
            'EM_TSS_Peak': peak, 'EM_TSS_Peak_Frac': peak_mass,
            # Where this unit's promoter zone came from, and how far the reads disagree with it.
            # A gene whose annotated 5' end is wrong puts its zone in the wrong place, and every
            # 5' end it should have explained looks like a truncation instead.
            'Promoter_Anchor': 'observed' if n_peaks[j] else 'annotated',
            'TSS_Offset': ((peak - (u['start'] if u['strand'] == '+' else u['end'] - 1))
                           * (1 if u['strand'] == '+' else -1)) if peak != 'NA' else 'NA',
            'EM_End_Peak': tes, 'EM_End_Peak_Frac': tes_mass,
            'Term_Frac': c3[TERM] / N if N > 0 else 0.0,
            'Body_End_Frac': c3[BODY] / N if N > 0 else 0.0,
            'Readout_Frac': c3[READOUT] / N if N > 0 else 0.0,
            'Spliced_Frac': (S['Sp'][j] + a_s) / (N + a_s + b_s),
            'Intronic_5p_Frac': ((S['Intr'][j] + a_e) / (N + a_e + b_e)) if u['exons'] else 'NA',
            'Iterations': n_iter,
        }
    return per_record, summary


# ---------------------------------------------------------------------------------------------
# Clustering and driver
# ---------------------------------------------------------------------------------------------

def clusters_from_records(records):
    """Union TUs that share a multi-candidate read; returns list of (unit_keys, record_indices)."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for rec in records:
        cands = rec[4]
        if len(cands) < 2:
            continue
        for c in cands:
            parent.setdefault(c, c)
        root = find(cands[0])
        for c in cands[1:]:
            rc = find(c)
            if rc != root:
                parent[rc] = root
    groups = defaultdict(lambda: [set(), []])
    for i, rec in enumerate(records):
        cands = rec[4]
        if cands[0] not in parent:
            continue          # single-candidate read of a TU with no shared reads: nothing to do
        root = find(cands[0])
        groups[root][0].update(cands)
        groups[root][1].append(i)
    return [(sorted(ks), ri) for ks, ri in groups.values()]


def _run_one(args):
    units, recs, p, trunc = args
    return em_cluster(units, recs, p, trunc)


def _run_init(args):
    units, recs, p = args
    return initial_censoring(units, recs, p)


def _fit_trunc(L, WE, WC, SA, p, log, label, n_full=None):
    """Pooled Kaplan-Meier fit, plus one per sample when more than one BAM contributed enough
    reads. A sample below trunc_min_sample falls back to the pooled curve."""
    pooled = TruncModel.fit_km(L, WE, WC)
    by_sample = {}
    counts = Counter(SA)
    if len(counts) > 1:
        for sid, n in sorted(counts.items()):
            if n < p.trunc_min_sample:
                log(f"EM: sample {sid}: {n:,} reads, below trunc_min_sample; using the pooled curve")
                continue
            sl = [(l, we, wc) for l, we, wc, sa in zip(L, WE, WC, SA) if sa == sid]
            m = TruncModel.fit_km([x[0] for x in sl], [x[1] for x in sl], [x[2] for x in sl])
            by_sample[sid] = m
            log(f"EM: sample {sid}: read-length survival median {m.median():,.0f} bp ({n:,} reads)")
    tail = f", {n_full:,.0f} full-length reads" if n_full is not None else ""
    log(f"EM: {label} read-length survival: median {pooled.median():,.0f} bp{tail}")
    return TruncSet(pooled, by_sample)


def resolve(records, unit_builder, p, pool=None, log=print):
    """
    records      : list of (origin, endpos, spliced, juncs, cand_keys, aligned_len, sample);
                   cand_keys[0] is the innermost feature, sample indexes the BAM the read came from
    unit_builder : callable key -> unit spec
    Returns (per_record, summary, cluster_of_record, trunc_model); per_record[i] is None for
    records outside any cluster (single-candidate reads of units that never share a read).
    """
    clusters = clusters_from_records(records)
    unit_cache = {}
    cl_units = []
    for keys, _ridx in clusters:
        units = {}
        for kk in keys:
            if kk not in unit_cache:
                unit_cache[kk] = unit_builder(kk)
            units[kk] = unit_cache[kk]
        cl_units.append(units)

    mapper = (lambda f, jobs: pool.map(f, jobs, chunksize=8)) if pool else \
        (lambda f, jobs: [f(j) for j in jobs])

    # pass 0: reads starting in a called TSS peak are treated as full-length (censored)
    init = mapper(_run_init, [(units, [records[i] for i in ridx], p)
                              for units, (_k, ridx) in zip(cl_units, clusters)])
    L, WE, WC, SA = [], [], [], []
    for (_k, ridx), flags in zip(clusters, init):
        for i, fl in zip(ridx, flags):
            L.append(records[i][5])
            WE.append(1.0 - fl)
            WC.append(fl)
            SA.append(records[i][6])
    trunc = _fit_trunc(L, WE, WC, SA, p, log, 'initial')

    results = None
    for pas in range(max(1, p.trunc_passes)):
        jobs = [(units, [records[i] for i in ridx], p, trunc)
                for units, (_k, ridx) in zip(cl_units, clusters)]
        results = mapper(_run_one, jobs)
        if pas + 1 < p.trunc_passes:
            L, WE, WC, SA = [], [], [], []
            for (_k, ridx), (pr, _s) in zip(clusters, results):
                for i, res in zip(ridx, pr):
                    full = sum(gi for _kk, _g, gi in res)
                    L.append(records[i][5])
                    WE.append(max(0.0, 1.0 - full))
                    WC.append(full)
                    SA.append(records[i][6])
            trunc = _fit_trunc(L, WE, WC, SA, p, log, f'pass {pas + 1}',
                               n_full=sum(WC))

    per_record = [None] * len(records)
    cluster_of = [None] * len(records)
    summary = {}
    for cid, ((keys, ridx), (pr, summ)) in enumerate(zip(clusters, results)):
        for i, res in zip(ridx, pr):
            per_record[i] = res
            cluster_of[i] = cid
        for kk, v in summ.items():
            v['Cluster'] = cid
            v['Cluster_Units'] = len(keys)
            summary[kk] = v
    return per_record, summary, cluster_of, trunc


def _weighted_median(xs, ws):
    pairs = sorted(zip(xs, ws))
    tot = sum(ws)
    if tot <= 0:
        return 0.0
    acc = 0.0
    for x, w in pairs:
        acc += w
        if acc >= tot / 2:
            return x
    return pairs[-1][0]
