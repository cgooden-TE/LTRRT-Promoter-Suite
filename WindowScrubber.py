#!/usr/bin/env python3
"""
Database-enabled motif scanner: stores ALL motif hits for flexible post-analysis filtering.

Database Schema:
  elements: Feature metadata and TSS
  ta_regions: TA-rich regions per element
  ca_runs: All significant CA runs with p-values
  motif_hits: All motif matches above threshold
  
Inputs:
    - LTR FASTA with headers: >Feature|chr:start-end(strand)
    - U3 FASTA with same headers, sequences = U3 regions
    - TSS summary TSV: Feature, TSS1 absolute coordinates
    - DB Name for output SQLite database

Usage Example:
  python3 WindowScrubber.py \
    -l full_ltr_sequences.fa \
    -u3 u3_regions.fa \
    -t tss_summary.tsv \
    -db motif_hits.db \
    --tata_mismatch 1 --ccaat_mismatch 0
    
Query examples (see query_db.py companion script):
  - Best TATA per element by score
  - All TATA hits within TA-rich regions
  - Motifs within distance ranges
  - Statistical summaries
"""
import argparse
import csv
import json
import logging
import math
import re
import sqlite3
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from statistics import mean, pstdev, median
from scipy import stats
import random
import bisect


from Bio import SeqIO
from Bio.Seq import Seq

# ------------------------------
# Utilities: intervals & windows
# ------------------------------

def sliding_windows(seq: str, w: int, step: int, start: int = 0, end: Optional[int] = None):
    """Yield (s, e, subseq) for half-open windows [s,e) of length w within [start,end)."""
    if end is None:
        end = len(seq)
    stop = max(start, min(end, len(seq))) - w + 1
    for i in range(max(0, start), max(0, stop), step):
        if i + w <= end:
            yield i, i + w, seq[i:i + w]


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge half-open intervals given as (start, end). Returns half-open intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


# ------------------------------
# PWM scoring / helpers
# ------------------------------

Base = str
PWM = Dict[Base, List[float]]


def log2(x: float) -> float:
    return math.log(x, 2) if x > 0 else float('-inf')


def pwm_len(pwm: PWM) -> int:
    return len(next(iter(pwm.values())))


def pwm_logodds_score(window: str, pwm: PWM, bg: Dict[Base, float]) -> float:
    """Compute log-odds score sum_i log2(pwm[b][i] / bg[b])."""
    L = pwm_len(pwm)
    if len(window) != L:
        return float('-inf')
    score = 0.0
    for i, b in enumerate(window):
        if b not in 'ACGT':
            return float('-inf')
        pb = pwm.get(b, [0.0] * L)[i]
        qb = bg.get(b, 0.0)
        pb = max(pb, 1e-6)
        qb = max(qb, 1e-6)
        score += log2(pb / qb)
    return score


def pwm_scores_matrix(pwm: PWM, bg: Dict[Base, float]) -> Dict[Base, List[float]]:
    """Precompute per-base, per-position log2 odds."""
    L = pwm_len(pwm)
    mat: Dict[Base, List[float]] = {b: [0.0] * L for b in 'ACGT'}
    for b in 'ACGT':
        for i in range(L):
            pb = max(pwm.get(b, [0.0]*L)[i], 1e-6)
            qb = max(bg.get(b, 0.0), 1e-6)
            mat[b][i] = log2(pb / qb)
    return mat


def pwm_perfect_score(pwm: PWM, bg: Dict[Base, float]) -> float:
    mat = pwm_scores_matrix(pwm, bg)
    L = pwm_len(pwm)
    s = 0.0
    for i in range(L):
        s += max(mat[b][i] for b in 'ACGT')
    return s


def pwm_position_drops(pwm: PWM, bg: Dict[Base, float]) -> List[float]:
    """For each position, drop = (best base score - second-best base score)."""
    mat = pwm_scores_matrix(pwm, bg)
    L = pwm_len(pwm)
    drops = []
    for i in range(L):
        vals = sorted([mat[b][i] for b in 'ACGT'], reverse=True)
        drops.append(max(0.0, vals[0] - vals[1]))
    return drops


def threshold_from_mismatches(pwm: PWM, bg: Dict[Base, float], mismatches: int) -> float:
    """Convert mismatch allowance to bits threshold."""
    perf = pwm_perfect_score(pwm, bg)
    drops = sorted(pwm_position_drops(pwm, bg))
    m = max(0, min(mismatches, len(drops)))
    return perf - sum(drops[:m])


def scan_all_pwm_hits(seq: str, pwm: PWM, bg: Dict[Base, float], 
                      start: int, end: int, threshold: float) -> List[Tuple[int, int, float, str]]:
    """Return ALL hits above threshold as (rel_start, rel_end, score, sequence).
    No filtering - just raw hits.
    """
    L = pwm_len(pwm)
    hits = []
    for s, e, sub in sliding_windows(seq, L, 1, start, end):
        sc = pwm_logodds_score(sub, pwm, bg)
        if sc >= threshold:
            hits.append((s, e, sc, sub))
    return hits

# ------------------------------
# Threshold helpers
# ------------------------------

def pwm_score_bounds(pwm: PWM, bg: Dict[Base, float]) -> Tuple[float, float]:
    """
    Returns (min_possible_score, max_possible_score) for this PWM under bg.
    Implemented as sum over positions of the min and max LLR.
    """
    mat = pwm_scores_matrix(pwm, bg)
    L = pwm_len(pwm)
    s_min = 0.0
    s_max = 0.0
    for i in range(L):
        vals = [mat[b][i] for b in 'ACGT']
        s_min += min(vals)
        s_max += max(vals)
    return float(s_min), float(s_max)


def relative_threshold(pwm: PWM, bg: Dict[Base, float], rel: float = 0.85) -> float:
    """
    Relative cutoff: min + rel * (max - min). rel in [0,1].
    """
    s_min, s_max = pwm_score_bounds(pwm, bg)
    return s_min + rel * (s_max - s_min)


def threshold_from_mismatches_anywhere(pwm: PWM, bg: Dict[Base, float], m: int, epsilon: float = 1e-9) -> float:
    """
    Worst-case conversion from 'allow up to m mismatches anywhere' to a score cutoff:
    subtract the sum of the m LARGEST per-position deltas from the perfect score.
    This guarantees acceptance for ANY configuration of m mismatches.
    """
    perf = pwm_perfect_score(pwm, bg)
    drops = sorted(pwm_position_drops(pwm, bg), reverse=True)
    m = max(0, min(m, len(drops)))
    penalty = sum(drops[:m])
    return perf - penalty - epsilon  # epsilon avoids equality traps at S_max

# ------------------------------
# LLR CDF and p-value estimation
# ------------------------------
def sample_bg_scores_iid(pwm: PWM, bg: Dict[Base, float], n: int, seed: int) -> List[float]:
    """Monte-Carlo null: sample n iid A/C/G/T kmers of motif length, score with PWM LLR."""
    rng = random.Random(seed)
    L = pwm_len(pwm)
    weights = [bg['A'], bg['C'], bg['G'], bg['T']]
    cum = [weights[0], weights[0]+weights[1], weights[0]+weights[1]+weights[2], 1.0]
    def draw_base():
        x = rng.random()
        if x < cum[0]: return 'A'
        if x < cum[1]: return 'C'
        if x < cum[2]: return 'G'
        return 'T'
    scores = []
    for _ in range(n):
        kmer = ''.join(draw_base() for _ in range(L))
        scores.append(pwm_logodds_score(kmer, pwm, bg))
    scores.sort()
    return scores

def score_to_pvalue(sorted_scores: List[float], s: float) -> float:
    """Right-tail p = P_bg(LLR >= s) using empirical CDF with a +1/(N+1) continuity correction."""
    i = bisect.bisect_left(sorted_scores, s)  # first index >= s
    tail = len(sorted_scores) - i
    return (tail + 1.0) / (len(sorted_scores) + 1.0)

def cutoff_from_bg(sorted_scores: List[float], alpha: float) -> float:
    """Return LLR cutoff so that P_bg(LLR >= cutoff) <= alpha (1 - alpha quantile)."""
    N = len(sorted_scores)
    idx = max(0, min(N-1, int(math.ceil((1.0 - alpha) * N) - 1)))
    return float(sorted_scores[idx])

# ------------------------------
# Default PWMs
# ------------------------------

def y_bias():
    return {"A": 0.05, "C": 0.475, "G": 0.05, "T": 0.425}


def n_uniform():
    return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}


def cons(base: str):
    p = {"A": 0.033, "C": 0.033, "G": 0.033, "T": 0.033}
    p[base] = 0.901
    return p


def dual(b1: str, b2: str):
    p = {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.10}
    p[b1] = 0.40
    p[b2] = 0.40
    return p

def mix(col, bg, alpha=0.6):
    # convex blend: alpha*bg + (1-alpha)*col
    out = {}
    for b in "ACGT":
        out[b] = alpha * bg[b] + (1 - alpha) * col[b]
    return out

bg_uniform = {"A":0.25,"C":0.25,"G":0.25,"T":0.25}

# Core: strong T, A, T, A
core_cols = [cons("T"), cons("A"), cons("T"), cons("A")]

# Tail: start from original intents, then blend toward background to weaken them
tail_cols = [
    mix(dual("A","T"), bg_uniform, alpha=0.6),  # position 5
    mix(cons("A"),      bg_uniform, alpha=0.6), # position 6
    mix(dual("A","T"), bg_uniform, alpha=0.6),  # position 7
    mix(dual("A","G"), bg_uniform, alpha=0.6),  # position 8
]

_cols = core_cols + tail_cols
TATA_PWM = {b: [col[b] for col in _cols] for b in "ACGT"}

DEFAULT_PWMS: Dict[str, PWM] = {
    "TATA": TATA_PWM,
    "CCAAT": {
        "A": [cons("C")["A"], cons("C")["A"], cons("A")["A"], cons("A")["A"], cons("T")["A"]],
        "C": [cons("C")["C"], cons("C")["C"], cons("A")["C"], cons("A")["C"], cons("T")["C"]],
        "G": [cons("C")["G"], cons("C")["G"], cons("A")["G"], cons("A")["G"], cons("T")["G"]],
        "T": [cons("C")["T"], cons("C")["T"], cons("A")["T"], cons("A")["T"], cons("T")["T"]],
    },
    "YPATCH": {b: [y_bias()[b]] * 8 for b in 'ACGT'},
    "INR7": {
        "A": [dual("C", "T")["A"], cons("T")["A"], cons("C")["A"], cons("A")["A"], n_uniform()["A"], dual("C", "T")["A"], dual("C", "T")["A"]],
        "C": [dual("C", "T")["C"], cons("T")["C"], cons("C")["C"], cons("A")["C"], n_uniform()["C"], dual("C", "T")["C"], dual("C", "T")["C"]],
        "G": [dual("C", "T")["G"], cons("T")["G"], cons("C")["G"], cons("A")["G"], n_uniform()["G"], dual("C", "T")["G"], dual("C", "T")["G"]],
        "T": [dual("C", "T")["T"], cons("T")["T"], cons("C")["T"], cons("A")["T"], n_uniform()["T"], dual("C", "T")["T"], dual("C", "T")["T"]],
    },
    "DPE7": {
        "A": [dual("A", "G")["A"], cons("G")["A"], dual("A", "T")["A"], cons("C")["A"], cons("G")["A"], cons("T")["A"], cons("G")["A"]],
        "C": [dual("A", "G")["C"], cons("G")["C"], dual("A", "T")["C"], cons("C")["C"], cons("G")["C"], cons("T")["C"], cons("G")["C"]],
        "G": [dual("A", "G")["G"], cons("G")["G"], dual("A", "T")["G"], cons("C")["G"], cons("G")["G"], cons("T")["G"], cons("G")["G"]],
        "T": [dual("A", "G")["T"], cons("G")["T"], dual("A", "T")["T"], cons("C")["T"], cons("G")["T"], cons("T")["T"], cons("G")["T"]],
    },
}


# ------------------------------
# CA run detection
# ------------------------------

def find_alternating_ca_runs(seq: str, min_length: int = 10) -> List[Tuple[int, int]]:
    """Find all maximal alternating CA runs of length >= min_length."""
    runs = []
    i = 0
    while i < len(seq) - 1:
        if (seq[i] == 'C' and seq[i+1] == 'A') or (seq[i] == 'A' and seq[i+1] == 'C'):
            start = i
            if seq[i] == 'C':
                j = i
                while j < len(seq) - 1 and seq[j] == 'C' and seq[j+1] == 'A':
                    j += 2
                end = j
            else:
                j = i
                while j < len(seq) - 1 and seq[j] == 'A' and seq[j+1] == 'C':
                    j += 2
                end = j
            
            run_length = end - start
            if run_length >= min_length:
                runs.append((start, end))
            i = end
        else:
            i += 1
    return runs


def binomial_test_ca_run(run_length: int, p_c: float, p_a: float) -> float:
    """Test if an alternating CA run is significant."""
    if run_length < 2:
        return 1.0
    n_dinucs = run_length // 2
    p_ca = p_c * p_a
    p_ac = p_a * p_c
    p_either = p_ca + p_ac
    p_alternating = (p_either) ** n_dinucs
    return p_alternating


# ------------------------------
# Database schema and setup
# ------------------------------

def init_database(db_path: str) -> sqlite3.Connection:
    """Initialize SQLite database with schema."""
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    c.execute('''
        CREATE TABLE IF NOT EXISTS thresholds (
            motif_type TEXT,
            method TEXT,
            param TEXT,
            cutoff REAL,
            PRIMARY KEY (motif_type, method, param)
        )
    ''')
    
    # Elements table
    c.execute('''
        CREATE TABLE IF NOT EXISTS elements (
            feature TEXT PRIMARY KEY,
            chrom TEXT,
            ltr_start INTEGER,
            ltr_end INTEGER,
            strand TEXT,
            u3_start INTEGER,
            u3_end INTEGER,
            tss_abs INTEGER,
            tss_rel INTEGER,
            sequence TEXT,
            sequence_length INTEGER
        )
    ''')
    
    # TA-rich regions
    c.execute('''
        CREATE TABLE IF NOT EXISTS ta_regions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            feature TEXT,
            start_rel INTEGER,
            end_rel INTEGER,
            start_abs INTEGER,
            end_abs INTEGER,
            FOREIGN KEY(feature) REFERENCES elements(feature)
        )
    ''')
    
    # CA runs
    c.execute('''
        CREATE TABLE IF NOT EXISTS ca_runs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            feature TEXT,
            start_rel INTEGER,
            end_rel INTEGER,
            start_abs INTEGER,
            end_abs INTEGER,
            length INTEGER,
            p_value REAL,
            is_significant INTEGER,
            FOREIGN KEY(feature) REFERENCES elements(feature)
        )
    ''')
    
    # Motif hits storage
    c.execute('''
        CREATE TABLE IF NOT EXISTS motif_hits (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            feature TEXT,
            motif_type TEXT,
            motif_length INTEGER,
            start_rel INTEGER,
            end_rel INTEGER,
            start_abs INTEGER,
            end_abs INTEGER,
            score REAL,
            p_value REAL,
            sequence TEXT,
            dist_to_tss INTEGER,
            in_ta_region INTEGER,
            y_count INTEGER,
            search_window_start INTEGER,
            search_window_end INTEGER,
            FOREIGN KEY(feature) REFERENCES elements(feature)
        )
    ''')
    
    # Create indices for common queries
    c.executescript('''
        CREATE INDEX IF NOT EXISTS idx_mh_type ON motif_hits(motif_type);
        CREATE INDEX IF NOT EXISTS idx_mh_feature ON motif_hits(feature);
        CREATE INDEX IF NOT EXISTS idx_mh_type_feature ON motif_hits(motif_type, feature);
        CREATE INDEX IF NOT EXISTS idx_mh_score ON motif_hits(score DESC);
        CREATE INDEX IF NOT EXISTS idx_mh_dist ON motif_hits(dist_to_tss);
        CREATE INDEX IF NOT EXISTS idx_mh_pval ON motif_hits(p_value);

        -- Filter for "by ranges per motif"
        CREATE INDEX IF NOT EXISTS idx_mh_type_dist_score
        ON motif_hits(motif_type, dist_to_tss, score DESC);

        -- On ta_regions
        CREATE INDEX IF NOT EXISTS idx_ta_feature ON ta_regions(feature);

        -- On thresholds
        CREATE INDEX IF NOT EXISTS idx_thr_lookup ON thresholds(method, param, motif_type);
    ''')
    
    conn.commit()
    return conn

# ------------------------------
# Parsing
# ------------------------------

def parse_header(header: str) -> Tuple[str, str, int, int, str]:
    """Parse header and return (feature, chrom, start, end, strand)."""
    feat, rest = header.split('|', 1)
    m = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)", rest)
    if not m:
        raise ValueError(f"Bad header: {header}")
    chrom = m.group(1)
    return feat, chrom, int(m.group(2)), int(m.group(3)), m.group(4)


class OrientedCoords:
    """Map between oriented relative indices (0-based) and absolute genomic (1-based inclusive)."""
    def __init__(self, ltr_s: int, ltr_e: int, strand: str):
        self.ltr_s = ltr_s
        self.ltr_e = ltr_e
        self.strand = strand

    def abs_to_rel(self, abs_pos: int) -> int:
        """Map absolute 1-based position to oriented 0-based index."""
        if self.strand == '+':
            return abs_pos - self.ltr_s
        else:
            return self.ltr_e - abs_pos

    def rel_to_abs_interval(self, rel_start: int, k: int) -> Tuple[int, int]:
        """Map oriented [rel_start, rel_start+k) to absolute 1-based inclusive (start<=end)."""
        if self.strand == '+':
            a = self.ltr_s + rel_start
            b = a + k - 1
            return a, b
        else:
            end_abs = self.ltr_e - rel_start
            start_abs = end_abs - (k - 1)
            return min(start_abs, end_abs), max(start_abs, end_abs)


def frac_bases(seq: str, allowed: set) -> float:
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in allowed) / len(seq)


# ------------------------------
# Troubleshooting Functions
# ------------------------------
def at_fraction(s: str) -> float:
    if not s:
        return 0.0
    s = s.upper()
    return (s.count("A") + s.count("T")) / len(s)

def ta_fraction(s: str) -> float:
    if len(s) < 2:
        return 0.0
    s = s.upper()
    ta = sum(1 for a, b in zip(s[:-1], s[1:]) if a == "T" and b == "A")
    return ta / (len(s) - 1)

def mono_shuffle(win: str, rng: random.Random) -> str:
    s = list(win)
    rng.shuffle(s)
    return "".join(s)

def dinuc_shuffle(win: str, rng: random.Random) -> str:
    """
    Dinucleotide-preserving shuffle for A/C/G/T-only strings.
    Uses randomized Hierholzer traversal on the implied multigraph.
    """
    if len(win) < 2:
        return win
    win = win.upper()
    # Guard: you said no ambiguous bases
    for ch in win:
        if ch not in "ACGT":
            raise ValueError(f"dinuc_shuffle: non-ACGT base found: {ch}")

    adj = defaultdict(list)
    for a, b in zip(win[:-1], win[1:]):
        adj[a].append(b)
    for a in adj:
        rng.shuffle(adj[a])

    start = win[0]
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if adj[v]:
            stack.append(adj[v].pop())
        else:
            path.append(stack.pop())
    path.reverse()
    out = "".join(path)
    if len(out) != len(win):
        raise RuntimeError("dinuc_shuffle produced incorrect length")
    return out

def window_scan_summary_local(win: str, pwm, bg, threshold: float):
    """Scan a window string directly, return hit count + best score."""
    hits = scan_all_pwm_hits(win, pwm, bg, start=0, end=len(win), threshold=threshold)
    if not hits:
        return 0, None
    best = max(h[2] for h in hits)
    return len(hits), best

def shuffle_control_for_window_fast(seq: str, pwm, bg,
                                    start: int, end: int,
                                    threshold: float,
                                    n_shuffles: int,
                                    rng,
                                    mode: str = "mono"):
    start = max(0, start)
    end = min(len(seq), end)
    if end <= start:
        return None

    win = seq[start:end]
    # real
    real_n, real_best = window_scan_summary_local(win, pwm, bg, threshold)

    # streaming mean/variance (Welford), ignoring "no hit" shuffles for score stats
    n_finite = 0
    mu = 0.0
    m2 = 0.0
    n_hit_shuf = 0

    for _ in range(n_shuffles):
        if mode == "mono":
            win_shuf = mono_shuffle(win, rng)
        elif mode in ("di", "dinuc"):
            win_shuf = dinuc_shuffle(win, rng) if len(win) >= 2 else win
        else:
            raise ValueError("mode must be 'mono' or 'di'/'dinuc'")

        sh_n, sh_best = window_scan_summary_local(win_shuf, pwm, bg, threshold)
        if sh_n > 0:
            n_hit_shuf += 1
        if sh_best is None:
            continue

        n_finite += 1
        delta = sh_best - mu
        mu += delta / n_finite
        m2 += delta * (sh_best - mu)

    sd = math.sqrt(m2 / n_finite) if n_finite > 1 else 0.0
    z = None
    if real_best is not None and n_finite > 1 and sd > 0:
        z = (real_best - mu) / sd

    return {
        "start": start, "end": end,
        "window_len": len(win),
        "window_AT": at_fraction(win),
        "window_TA": ta_fraction(win),
        "real_n_hits": real_n,
        "real_best_score": real_best,
        "shuffle_hit_rate": (n_hit_shuf / n_shuffles) if n_shuffles else 0.0,
        "shuffle_best_mean": (mu if n_finite else None),
        "shuffle_best_sd": (sd if n_finite else None),
        "best_z": z,
        "mode": mode,
    }
    
def spearman_corr(xs, ys):
    try:
        from scipy.stats import spearmanr
        r, p = spearmanr(xs, ys, nan_policy="omit")
        return r, p
    except ImportError:
        # lightweight fallback: rank then Pearson (no p-value)
        def rankdata(a):
            tmp = sorted((v, i) for i, v in enumerate(a))
            ranks = [0.0] * len(a)
            i = 0
            while i < len(tmp):
                j = i
                while j < len(tmp) and tmp[j][0] == tmp[i][0]:
                    j += 1
                avg = (i + 1 + j) / 2.0
                for k in range(i, j):
                    ranks[tmp[k][1]] = avg
                i = j
            return ranks

        def pearson(x, y):
            n = len(x)
            mx = sum(x) / n
            my = sum(y) / n
            num = sum((a - mx) * (b - my) for a, b in zip(x, y))
            denx = math.sqrt(sum((a - mx) ** 2 for a in x))
            deny = math.sqrt(sum((b - my) ** 2 for b in y))
            return num / (denx * deny) if denx and deny else float("nan")

        rx = rankdata(xs)
        ry = rankdata(ys)
        return pearson(rx, ry), None


# ------------------------------
# Argument parsing
# ------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Database-enabled motif scanner with flexible post-analysis filtering.\
            Required inputs: LTR FASTA (-l), U3 FASTA (-u3), TSS summary TSV (-t). Outputs SQLite database (-db) of motif hits and metadata.'
    )
    # Input files and output filenames 
    p.add_argument('-l', '--ltr-fasta', required=True, help='FASTA of full LTR sequences')
    p.add_argument('-u3', '--u3-fasta', required=True, help='FASTA of U3 regions')
    p.add_argument('-t', '--tss-summary', required=True, help='TSV: Feature, TSS1 absolute coordinates')
    p.add_argument('-db', '--database', required=True, help='Output SQLite database file')

    # Sliding window parameters
    p.add_argument('--window-size', type=int, default=10, help='Sliding window size for TA-rich detection')
    p.add_argument('--step-size', type=int, default=1, help='Sliding window step size')
    p.add_argument('--ta-threshold', type=float, default=0.75, help='Threshold for TA-rich windows')

    # PWM controls
    p.add_argument('--pwm-json', type=str, default=None, help='Optional JSON file providing PWMs')
    p.add_argument('--bg-A', type=float, default=0.25,
                    help="Background frequency for A (default 0.25). All bg frequencies should sum to 1.")
    p.add_argument('--bg-C', type=float, default=0.25,
                   help="Background frequency for C (default 0.25). All bg frequencies should sum to 1.")
    p.add_argument('--bg-G', type=float, default=0.25,
                   help="Background frequency for G (default 0.25). All bg frequencies should sum to 1.")
    p.add_argument('--bg-T', type=float, default=0.25,
                   help="Background frequency for T (default 0.25). All bg frequencies should sum to 1.")
    p.add_argument("--store-rel", type=float, default=0.85,
                    help="Relative score fallback for storage (0..1).")
    p.add_argument("--store-slack", type=float, default=1.5,
                    help="Lower storage thresholds by this many LLR units.")
    p.add_argument("--store-floor", type=float, default=0.0,
                    help="Absolute LLR floor for storage (e.g., 0.0 keeps only positive evidence).")
    p.add_argument('--sig-alpha', type=float, default=1e-4,
                help='Per-motif significance level for hits (p-value), FIMO-like default 1e-4.')
    p.add_argument('--bg-n', type=int, default=200000,
                help='Background kmers per motif to estimate null LLR CDF (larger = smoother p-values).')
    p.add_argument('--random-seed', type=int, default=13,
                help='Seed for background sampling reproducibility.')
    p.add_argument('--percentile-cut', type=float, default=95.0,
                help='Percentile of significant scores recorded as a downstream cutoff for STATS.')

    # Mismatch allowances → thresholds
    p.add_argument('--tata_mismatch', type=int, default=0)
    p.add_argument('--ccaat_mismatch', type=int, default=0)
    p.add_argument('--ypatch_mismatch', type=int, default=0)
    p.add_argument('--inr_mismatch', type=int, default=0)
    p.add_argument('--sec_tata_mismatch', type=int, default=0)
    p.add_argument('--dpe_mismatch', type=int, default=0)

    # CA run parameters
    p.add_argument('--ca-min-length', type=int, default=10, help='Minimum length for CA runs')
    p.add_argument('--ca-alpha', type=float, default=0.05, help='Significance level for CA runs')

    # Search windows
    p.add_argument('--tata-max-dist', type=int, default=100, help='Max distance upstream for primary TATA')
    p.add_argument('--inr-max-dist', type=int, default=10, help='Max offset around TSS for Inr')

    return p.parse_args()


# ------------------------------
# Main
# ------------------------------

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Background distribution
    bg = {"A": args.bg_A, "C": args.bg_C, "G": args.bg_G, "T": args.bg_T}
    if abs(sum(bg.values()) - 1.0) > 1e-6:
        s = sum(bg.values())
        bg = {k: v / s for k, v in bg.items()}

    EPS = 1e-9

    # ---- existing PWM loading code stays the same ----
    pwms = DEFAULT_PWMS.copy()
    if args.pwm_json:
        with open(args.pwm_json) as jf:
            loaded = json.load(jf)
        pwms.update(loaded)

    def pwm_len(pwm):
        if isinstance(pwm, list):
            return len(pwm)
        if isinstance(pwm, dict) and all(b in pwm for b in "ACGT"):
            return len(pwm["A"])
        raise TypeError("Unrecognized PWM format")

    for name, pwm in pwms.items():
        L = pwm_len(pwm)
        print(f"{name}\tlen={L}")
    
    # ---- STRICT thresholds from mismatches  ----
    thr = {
        'TATA':   threshold_from_mismatches_anywhere(pwms['TATA'],   bg, args.tata_mismatch),
        'CCAAT':  threshold_from_mismatches_anywhere(pwms['CCAAT'],  bg, args.ccaat_mismatch),
        'YPATCH': threshold_from_mismatches_anywhere(pwms['YPATCH'], bg, args.ypatch_mismatch),
        'INR7':   threshold_from_mismatches_anywhere(pwms['INR7'],   bg, args.inr_mismatch),
        'SEC_TATA':threshold_from_mismatches_anywhere(pwms['TATA'],  bg, args.sec_tata_mismatch),
        'DPE7':   threshold_from_mismatches_anywhere(pwms['DPE7'],   bg, args.dpe_mismatch),
    }

    # ---- PERMISSIVE store-time thresholds:
    # min( mismatch_anywhere, relative ) then subtract slack, floor at store-floor, minus epsilon
    def make_store_threshold(name, pwm, mismatch_count):
        t_mm  = threshold_from_mismatches_anywhere(pwm, bg, mismatch_count)         # (1)
        t_rel = relative_threshold(pwm, bg, args.store_rel)                         # (2)
        t     = min(t_mm, t_rel) - args.store_slack                                 # (3)
        t     = max(t, args.store_floor) - EPS                                      # (4,5)
        return t

    store_thr = {
        'TATA':    make_store_threshold('TATA',    pwms['TATA'],   args.tata_mismatch),
        'CCAAT':   make_store_threshold('CCAAT',   pwms['CCAAT'],  args.ccaat_mismatch),
        'YPATCH':  make_store_threshold('YPATCH',  pwms['YPATCH'], args.ypatch_mismatch),
        'INR7':    make_store_threshold('INR7',    pwms['INR7'],   args.inr_mismatch),
        'SEC_TATA':make_store_threshold('SEC_TATA',pwms['TATA'],   args.sec_tata_mismatch),
        'DPE7':    make_store_threshold('DPE7',    pwms['DPE7'],   args.dpe_mismatch),
    }
    
    logging.info(f"Thresholds: {thr}")
    
    # Initialize database
    conn = init_database(args.database)
    c = conn.cursor()
    
    random.seed(args.random_seed)

    # Build empirical nulls and per-motif significance LLR cutoffs (iid background)
    bg_scores = {}
    sig_cutoff = {}
    for name, pwm in pwms.items():
        # deterministic per-motif seed
        seed = (args.random_seed * 1315423911 + sum(ord(ch) for ch in name)) & 0xFFFFFFFF
        scores = sample_bg_scores_iid(pwm, bg, args.bg_n, seed)
        bg_scores[name] = scores
        sig_cutoff[name] = cutoff_from_bg(scores, args.sig_alpha)

    def put_thr(c, motif, method, param, cutoff):
        c.execute("INSERT OR REPLACE INTO thresholds VALUES (?,?,?,?)",
                (motif, method, param, float(cutoff)))

    for name in store_thr:
        # strict (mismatch), relative, and final storage thresholds
        put_thr(c, name, 'mismatch_anywhere',
                f"m={getattr(args, name.lower()+'_mismatch', 'NA')}",
                thr[name])
        rel_val = relative_threshold(pwms['TATA' if name=='SEC_TATA' else name], bg, args.store_rel)
        put_thr(c, name, 'relative', f"rel={args.store_rel}", rel_val)
        put_thr(c, name, 'store', f"rel={args.store_rel};slack={args.store_slack};floor={args.store_floor}",
                store_thr[name])

        # significance cutoff (SEC_TATA uses TATA PWM)
        base = 'TATA' if name == 'SEC_TATA' else name
        put_thr(c, name, 'empirical_p', f"alpha={args.sig_alpha};bg=iid;n={args.bg_n}",
                sig_cutoff[base])

    # Scan-time threshold: Permissive for wide store
    # Change if needed, keep scan and store separate
    scan_thr = {name: store_thr[name] for name in store_thr}


    # Load U3 coordinates
    u3_map = {}
    for rec in SeqIO.parse(args.u3_fasta, 'fasta'):
        feat, _chrom, s, e, strand = parse_header(rec.id)
        u3_map[feat] = (s, e, strand)
    logging.info(f"Loaded U3 coords for {len(u3_map)} features")

    # Load TSS
    tss_map = {}
    with open(args.tss_summary) as tf:
        rdr = csv.DictReader(tf, delimiter='\t')
        for row in rdr:
            f = row.get('Feature', row.get('Gene'))
            try:
                tss_map[f] = int(row.get('TSS1', ''))
            except (TypeError, ValueError):
                continue
    logging.info(f"Loaded {len(tss_map)} TSS entries")


    total_tests = 0
    prim_best_scores = []
    sec_best_scores  = []

    prim_z_scores = []
    sec_z_scores  = []

    prim_at = []
    prim_ta = []
    prim_best = []
    sec_at  = []
    sec_ta = []
    sec_best = []

    n_prim_total = 0
    n_sec_total  = 0
    n_prim_hit   = 0
    n_sec_hit    = 0
    top_sec = []   # list of tuples (best_score, feat, AT, Z)
    TOP_N = 10
    
    rng = random.Random(12345)
    
    for rec in SeqIO.parse(args.ltr_fasta, 'fasta'):
        feat, chrom, ltr_s, ltr_e, strand = parse_header(rec.id)
        if feat not in u3_map or feat not in tss_map:
            continue
        
        u3_s, u3_e, _u3_strand = u3_map[feat]
        tss_abs = tss_map[feat]

        # Orient sequence
        seq_raw = str(rec.seq).upper()
        if strand == '+':
            seq = seq_raw
            tss_rel = tss_abs - ltr_s
        else:
            seq = str(Seq(seq_raw).reverse_complement())
            tss_rel = ltr_e - tss_abs
        
        mapper = OrientedCoords(ltr_s, ltr_e, strand)

        # U3 interval in oriented coordinates
        u3_r1 = mapper.abs_to_rel(u3_s)
        u3_r2 = mapper.abs_to_rel(u3_e)
        u3_lo, u3_hi = (u3_r1, u3_r2) if u3_r1 <= u3_r2 else (u3_r2, u3_r1)

        # Insert element record
        c.execute('''
            INSERT OR REPLACE INTO elements 
            (feature, chrom, ltr_start, ltr_end, strand, u3_start, u3_end, tss_abs, tss_rel, sequence, sequence_length)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (feat, chrom, ltr_s, ltr_e, strand, u3_s, u3_e, tss_abs, tss_rel, seq, len(seq)))

        # TA-rich regions (restricted to U3)
        ta_windows = []
        for s, e, sub in sliding_windows(seq, args.window_size, args.step_size, start=u3_lo, end=u3_hi):
            if 'N' in sub:
                continue
            if frac_bases(sub, {'T','A'}) >= args.ta_threshold:
                ta_windows.append((s, e))
        ta_regions = merge_intervals(ta_windows)
        
        for ta_s, ta_e in ta_regions:
            ta_s_abs, ta_e_abs = mapper.rel_to_abs_interval(ta_s, ta_e - ta_s)
            c.execute('''
                INSERT INTO ta_regions (feature, start_rel, end_rel, start_abs, end_abs)
                VALUES (?, ?, ?, ?, ?)
            ''', (feat, ta_s, ta_e, ta_s_abs, ta_e_abs))

        # Helper to check if position is in TA region
        def in_ta_region(s_rel: int, e_rel: int) -> bool:
            for a, b in ta_regions:
                if s_rel >= a and e_rel <= b:
                    return True
            return False

        # CA runs across entire LTR
        ca_runs = find_alternating_ca_runs(seq, min_length=args.ca_min_length)
        for run_start, run_end in ca_runs:
            run_length = run_end - run_start
            p_value = binomial_test_ca_run(run_length, bg['C'], bg['A'])
            total_tests += 1
            
            bonferroni_alpha = args.ca_alpha / max(1, total_tests)
            is_significant = 1 if p_value <= bonferroni_alpha else 0
            
            ca_s_abs, ca_e_abs = mapper.rel_to_abs_interval(run_start, run_length)
            c.execute('''
                INSERT INTO ca_runs (feature, start_rel, end_rel, start_abs, end_abs, length, p_value, is_significant)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, run_start, run_end, ca_s_abs, ca_e_abs, run_length, p_value, is_significant))

        # === MOTIF SCANNING - Store ALL hits ===
        
        # 1) Primary TATA upstream
        up_start = max(0, tss_rel - args.tata_max_dist)
        #up_end = max(0, tss_rel)
        up_end = max(0, tss_rel + 8)
        tata_hits = scan_all_pwm_hits(seq, pwms['TATA'], bg, up_start, up_end, scan_thr['TATA'])
        for s_rel, e_rel, score, subseq in tata_hits:
            # p-value under background (TATA PWM)
            pval = score_to_pvalue(bg_scores['TATA'], score)
            
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['TATA']))
            dist = tss_rel - s_rel
            in_ta = 1 if in_ta_region(s_rel, e_rel) else 0
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits 
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs, 
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'TATA', pwm_len(pwms['TATA']), s_rel, e_rel, s_abs, e_abs, 
                  score, subseq, dist, in_ta, yc, up_start, up_end, pval))

        # 2) CCAAT-box upstream 460-140 bp
        cca_band_start = max(0, tss_rel - 460)
        cca_band_end = max(0, tss_rel - 140)
        ccaat_hits = scan_all_pwm_hits(seq, pwms['CCAAT'], bg, cca_band_start, cca_band_end, scan_thr['CCAAT'])
        for s_rel, e_rel, score, subseq in ccaat_hits:
            # p-value under background (CCAAT PWM)
            pval = score_to_pvalue(bg_scores['CCAAT'], score)
            
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['CCAAT']))
            dist = tss_rel - s_rel
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits 
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs, 
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'CCAAT', pwm_len(pwms['CCAAT']), s_rel, e_rel, s_abs, e_abs, 
                  score, subseq, dist, 0, yc, cca_band_start, cca_band_end, pval))

        # 3) Y-PATCH upstream 100–1 bp, 8-mer
        yp_start = max(0, tss_rel - 100)
        yp_end   = max(0, tss_rel)
        ypatch_hits = scan_all_pwm_hits(seq, pwms['YPATCH'], bg, yp_start, yp_end, scan_thr['YPATCH'])
        for s_rel, e_rel, score, subseq in ypatch_hits:
            pval = score_to_pvalue(bg_scores['YPATCH'], score)
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['YPATCH']))
            dist = tss_rel - s_rel
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs,
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'YPATCH', 8, s_rel, e_rel, s_abs, e_abs,
                  score, subseq, dist, 0, yc, yp_start, yp_end, pval))

        # 4) Inr motif near TSS
        inr_band_start = max(0, tss_rel - args.inr_max_dist)
        inr_band_end = min(len(seq), tss_rel + args.inr_max_dist + 1)
        inr_hits = scan_all_pwm_hits(seq, pwms['INR7'], bg, inr_band_start, inr_band_end, scan_thr['INR7'])
        for s_rel, e_rel, score, subseq in inr_hits:
            # p-value under background (INR PWM)
            pval = score_to_pvalue(bg_scores['INR7'], score)
            
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['INR7']))
            dist = abs(s_rel - tss_rel)
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits 
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs, 
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'INR7', pwm_len(pwms['INR7']), s_rel, e_rel, s_abs, e_abs, 
                  score, subseq, dist, 0, yc, inr_band_start, inr_band_end, pval))

        # 5) Secondary TATA -120 bp downstream
        sec_band_start = tss_rel + 1
        sec_band_end = min(len(seq), tss_rel + 120)
        sec_tata_hits = scan_all_pwm_hits(seq, pwms['TATA'], bg, sec_band_start, sec_band_end, scan_thr['SEC_TATA'])
        for s_rel, e_rel, score, subseq in sec_tata_hits:
            # p-value under background (TATA PWM)
            pval = score_to_pvalue(bg_scores['TATA'], score)
            
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['TATA']))
            dist = s_rel - tss_rel
            in_ta = 1 if in_ta_region(s_rel, e_rel) else 0
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits 
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs, 
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'SEC_TATA', pwm_len(pwms['TATA']), s_rel, e_rel, s_abs, e_abs, 
                  score, subseq, dist, in_ta, yc, sec_band_start, sec_band_end, pval))

        # 6) DPE motif 0-50 bp downstream
        dpe_band_start = tss_rel
        dpe_band_end = min(len(seq), tss_rel + 50)
        dpe_hits = scan_all_pwm_hits(seq, pwms['DPE7'], bg, dpe_band_start, dpe_band_end, scan_thr['DPE7'])
        for s_rel, e_rel, score, subseq in dpe_hits:
            # p-value under background (DPE PWM)
            pval = score_to_pvalue(bg_scores['DPE7'], score)
            
            s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, pwm_len(pwms['DPE7']))
            dist = s_rel - tss_rel
            yc = sum(1 for b in subseq if b in 'CT')
            c.execute('''
                INSERT INTO motif_hits 
                (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs, 
                 score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (feat, 'DPE7', pwm_len(pwms['DPE7']), s_rel, e_rel, s_abs, e_abs, 
                  score, subseq, dist, 0, yc, dpe_band_start, dpe_band_end, pval))
            
        # Trouble shooting TATA and SEC TATA 

        prim_mono = shuffle_control_for_window_fast(
            seq, pwms["TATA"], bg, up_start, up_end, scan_thr["TATA"],
            n_shuffles=100, rng=rng, mode="mono"
        )
        prim_di = shuffle_control_for_window_fast(
            seq, pwms["TATA"], bg, up_start, up_end, scan_thr["TATA"],
            n_shuffles=50, rng=rng, mode="di"
        )

        sec_mono = shuffle_control_for_window_fast(
            seq, pwms["TATA"], bg, sec_band_start, sec_band_end, scan_thr["SEC_TATA"],
            n_shuffles=100, rng=rng, mode="mono"
        )
        sec_di = shuffle_control_for_window_fast(
            seq, pwms["TATA"], bg, sec_band_start, sec_band_end, scan_thr["SEC_TATA"],
            n_shuffles=50, rng=rng, mode="di"
        )

        # after you compute prim_mono and sec_mono (or prim_di/sec_di, choose one consistently)
        n_prim_total += 1
        n_sec_total  += 1

        if prim_mono and prim_mono["real_best_score"] is not None:
            prim_at.append(prim_mono["window_AT"])
            prim_ta.append(prim_mono["window_TA"])
            prim_best.append(prim_mono["real_best_score"])
            n_prim_hit += 1
            prim_best_scores.append(prim_mono["real_best_score"])
            if prim_mono["best_z"] is not None:
                prim_z_scores.append(prim_mono["best_z"])

        if sec_mono and sec_mono["real_best_score"] is not None:
            sec_at.append(sec_mono["window_AT"])
            sec_ta.append(sec_mono["window_TA"])
            sec_best.append(sec_mono["real_best_score"])
            n_sec_hit += 1
            sec_best_scores.append(sec_mono["real_best_score"])
            # keep top N secondary hits by score for inspection
            top_sec.append((sec_mono["real_best_score"], feat, sec_mono["window_AT"], sec_mono["best_z"]))
            if sec_mono["best_z"] is not None:
                sec_z_scores.append(sec_mono["best_z"])


    def print_spearman(label, xs, ys):
        # Need at least 3 points to be meaningful
        if len(xs) < 3 or len(ys) < 3:
            print(f"{label} Spearman: not enough points (n={min(len(xs), len(ys))})")
            return
        r, p = spearman_corr(xs, ys)
        # If you’re using the fallback spearman without scipy, p will be None
        if p is None:
            print(f"{label} Spearman rho={r:.3f} (p-value unavailable; install scipy for p)")
        else:
            print(f"{label} Spearman rho={r:.3f}, p={p:.3g} (n={min(len(xs), len(ys))})")

    print("\n=== Correlations (best score vs composition) ===")
    print_spearman("Primary: AT% vs best-score", prim_at, prim_best)
    print_spearman("Primary: TA% vs best-score", prim_ta, prim_best)
    print_spearman("Secondary: AT% vs best-score", sec_at, sec_best)
    print_spearman("Secondary: TA% vs best-score", sec_ta, sec_best)

    def summarize(label, scores, zs, ats, n_total, n_hit):
        print(f"\n=== {label} summary ===")
        print(f"Loci scanned: {n_total}")
        print(f"Loci with >=1 hit: {n_hit} ({(n_hit/n_total*100 if n_total else 0):.2f}%)")

        if scores:
            print(f"Best-score mean:   {mean(scores):.3f}")
            print(f"Best-score median: {median(scores):.3f}")
            print(f"Best-score sd:     {pstdev(scores):.3f}")
        else:
            print("No hits -> no score summary")

        if zs:
            print(f"Z mean:   {mean(zs):.3f}")
            print(f"Z median: {median(zs):.3f}")
            print(f"Z sd:     {pstdev(zs):.3f}")
        else:
            print("No finite Z values to summarize")

        if ats:
            print(f"AT% mean (window): {mean(ats):.3f}")
        else:
            print("No AT% values collected")

    summarize("Primary (mono)", prim_best_scores, prim_z_scores, prim_at, n_prim_total, n_prim_hit)
    summarize("Secondary (mono)", sec_best_scores, sec_z_scores, sec_at, n_sec_total, n_sec_hit)

    # optional: show top N secondary hits
    top_sec.sort(reverse=True, key=lambda x: x[0])
    print(f"\nTop {min(TOP_N, len(top_sec))} secondary hits (by best score):")
    for sc, feat, atv, z in top_sec[:TOP_N]:
        print(f"  {feat}\tscore={sc:.3f}\tAT%={atv:.3f}\tZ={z}")

    def compute_and_store_percentile_cutoffs(conn, perc: float):
        cur = conn.cursor()
        for name in ['TATA','CCAAT','YPATCH','INR7','SEC_TATA','DPE7']:
            rows = cur.execute("SELECT score FROM motif_hits WHERE motif_type=?", (name,)).fetchall()
            if not rows:
                continue
            scores = sorted(r[0] for r in rows)
            # percentile index, inclusive
            idx = max(0, min(len(scores)-1, int(math.ceil(perc/100.0 * len(scores)) - 1)))
            cut = float(scores[idx])
            cur.execute("INSERT OR REPLACE INTO thresholds VALUES (?,?,?,?)",
                        (name, 'empirical_percentile', f'perc={perc}', cut))
        conn.commit()

    compute_and_store_percentile_cutoffs(conn, args.percentile_cut)
    conn.close()
    
    logging.info(f"Database written to {args.database}")
    logging.info(f"Total statistical tests for CA runs: {total_tests}")
    logging.info("Use Query_WSDB.py to analyze results with flexible filtering")


if __name__ == '__main__':
    main()