#!/usr/bin/env python3
"""
Database-enabled motif scanner: stores ALL motif hits for flexible post-analysis filtering.

Every element is anchored on its IsoClassifier TSS and scanned with the same TSS-relative motif
bands. Two schemas decide the scanned sequence and the TA-rich search zone:
  U3          LTR_structural with TSS1 inside the TSS-bearing LTR (LTR FASTA supplied).
              Sequence = that LTR; TA-rich zone = U3, derived as upstream LTR edge .. TSS1.
  TSS_window  Everything else (genes, TIR, Helitron, LINE, SINE, LTR_fragment, and structural
              LTRs without a usable U3). Sequence = genome window TSS-flank_up .. TSS+flank_down;
              TA-rich zone = fixed core band TSS-core_up .. TSS+core_down. No U3 is imputed.

Sequences are always pulled from the genome in the transcript's frame; the LTR FASTA only
supplies coordinates, so its orientation convention cannot flip a scan.

Database Schema:
  elements:   Feature metadata, class, anchor mode, scanned span, TSS, TA-search zone
  ta_regions: TA-rich regions per element
  ca_runs:    All significant CA runs with p-values
  motif_hits: All motif matches above threshold
  thresholds: Per-motif score cutoffs
  run_params: Parameters used to build the database

Inputs:
    - TSS summary TSV from IsoClassifier (*_sense.tss_summary.tsv): Feature, Class, Orientation,
      Strand, Chrom, TSS1 (0-based)
    - Genome FASTA (faidx-indexed)
    - Optional LTR FASTA from IsoClassifier (*_ltr_seqs_<orient>.fa) for the U3 schema
    - Optional IsoClassifier isoform table, to supply Chrom for older TSS summaries lacking it

Usage Example:
  python3 WindowScrubber.py \
    -t iso_sense.tss_summary.tsv \
    -g B73.PLATINUM.pseudomolecules-v1.fasta \
    -l iso_ltr_seqs_sense.fa \
    -db motif_hits.db \
    --tata_mismatch 1 --ccaat_mismatch 0

Query examples (see Query_WSDB.py companion script):
  - Best TATA per element by score
  - All TATA hits within TA-rich regions
  - Motifs within distance ranges, per element class or anchor mode
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


import random
import bisect
from dataclasses import dataclass

import pysam
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
    
    # Run parameters (flanks, core band, inputs) so each DB documents its own schema choices
    c.execute('''
        CREATE TABLE IF NOT EXISTS run_params (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    ''')

    # Elements table: seq_start/seq_end = scanned span (LTR in U3 mode, genome window otherwise);
    # zone_start/zone_end = TA-rich search zone; u3_* is NULL outside U3 mode. All 1-based inclusive.
    c.execute('''
        CREATE TABLE IF NOT EXISTS elements (
            feature TEXT PRIMARY KEY,
            element_class TEXT,
            orientation TEXT,
            anchor_mode TEXT,
            chrom TEXT,
            seq_start INTEGER,
            seq_end INTEGER,
            strand TEXT,
            u3_start INTEGER,
            u3_end INTEGER,
            zone_start INTEGER,
            zone_end INTEGER,
            tss_abs INTEGER,
            tss_rel INTEGER,
            total_reads INTEGER,
            tss_called INTEGER,
            seq_truncated INTEGER,
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

        -- On elements
        CREATE INDEX IF NOT EXISTS idx_el_class_mode ON elements(element_class, anchor_mode);

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
# Element construction (U3 vs TSS_window schemas)
# ------------------------------

ALL_CLASSES = ('LTR_structural', 'LTR_fragment', 'TIR', 'Helitron', 'LINE', 'SINE', 'Gene')
DEFAULT_GENOME = '/home/caleb/data/genome_and_annotations/B73.PLATINUM.pseudomolecules-v1.fasta'


@dataclass
class Element:
    """One TSS-anchored scan unit. seq is oriented 5'->3' in the transcript's frame."""
    feature: str
    cls: str
    orientation: str
    chrom: str
    strand: str
    anchor_mode: str          # 'U3' or 'TSS_window'
    seq_start: int            # 1-based inclusive span of seq
    seq_end: int
    seq: str
    tss_abs: int              # 1-based
    tss_rel: int              # 0-based index into seq
    zone_lo: int              # oriented half-open TA-search zone [zone_lo, zone_hi)
    zone_hi: int
    truncated: int
    total_reads: Optional[int]
    tss_called: Optional[int]


def _int_or_none(x) -> Optional[int]:
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def load_tss_summary(path: str, orientation: str, classes: set,
                     chrom_fallback: Dict[str, str]) -> List[dict]:
    """
    Read an IsoClassifier tss_summary. TSS1 is 0-based there and is converted to 1-based here.
    Rows are filtered to one orientation so Feature is unique; duplicates keep the first row.
    """
    rows, seen = [], set()
    skipped = {}

    def skip(reason):
        skipped[reason] = skipped.get(reason, 0) + 1

    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        missing = {'Feature', 'Class', 'Strand', 'TSS1'} - set(rdr.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks columns {sorted(missing)}; "
                             "use an IsoClassifier *_{sense,antisense}.tss_summary.tsv")
        for r in rdr:
            orient = r.get('Orientation', orientation)
            if orient != orientation:
                continue
            if r['Class'] not in classes:
                skip(f"class {r['Class']} not selected")
                continue
            tss0 = _int_or_none(r['TSS1'])
            if tss0 is None:
                skip('no TSS1')
                continue
            if r['Strand'] not in ('+', '-'):
                skip('unresolved strand')
                continue
            feat = r['Feature']
            if feat in seen:
                skip('duplicate feature')
                continue
            seen.add(feat)
            rows.append({
                'feature': feat, 'cls': r['Class'], 'orientation': orient,
                'strand': r['Strand'],
                'chrom': r.get('Chrom') or chrom_fallback.get(feat),
                'tss_abs': tss0 + 1,
                'total_reads': _int_or_none(r.get('Total_Reads')),
                'tss_called': _int_or_none(r.get('TSS_Called')),
            })
    for k, v in skipped.items():
        logging.info(f"TSS summary: skipped {v} rows ({k})")
    return rows


def load_chrom_map(isoforms_path: Optional[str]) -> Dict[str, str]:
    """Feature -> Chrom from an IsoClassifier isoform table (for TSS summaries without Chrom)."""
    if not isoforms_path:
        return {}
    out = {}
    with open(isoforms_path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out.setdefault(r['Feature'], r['Chrom'])
    return out


def load_ltr_coords(path: Optional[str], orientation: str) -> Dict[str, Tuple[str, int, int, str]]:
    """
    Feature -> (chrom, start1, end1, strand) of the TSS-bearing LTR from IsoClassifier's
    *_ltr_seqs_<orient>.fa. Only header coordinates are used; sequence comes from the genome.
    """
    if not path:
        return {}
    out = {}
    for rec in SeqIO.parse(path, 'fasta'):
        parts = rec.id.split('|')
        if len(parts) >= 3 and parts[2] != orientation:
            continue
        feat, chrom, s, e, strand = parse_header(rec.id)
        out[feat] = (chrom, s, e, strand)
    return out


def fetch_oriented(fa: pysam.FastaFile, chrom: str, s1: int, e1: int, strand: str) -> str:
    """Fetch 1-based inclusive [s1, e1] and orient to the given strand."""
    seq = fa.fetch(chrom, s1 - 1, e1).upper()
    return str(Seq(seq).reverse_complement()) if strand == '-' else seq


def build_element(row: dict, fa: pysam.FastaFile, chrom_len: Dict[str, int],
                  ltr_coords: Dict[str, Tuple[str, int, int, str]], args) -> Tuple[Optional[Element], str]:
    """
    Choose the schema for one TSS-summary row and build its scan unit.
    Returns (element, reason); element is None when the row cannot be scanned.
    """
    feat, strand, tss_abs = row['feature'], row['strand'], row['tss_abs']
    chrom = row['chrom']
    common = dict(feature=feat, cls=row['cls'], orientation=row['orientation'], strand=strand,
                  total_reads=row['total_reads'], tss_called=row['tss_called'])

    # U3 schema: structural LTR whose TSS1 sits inside the TSS-bearing LTR
    if row['cls'] == 'LTR_structural' and feat in ltr_coords:
        l_chrom, l_s, l_e, l_strand = ltr_coords[feat]
        chrom = chrom or l_chrom
        if l_chrom == chrom and l_strand == strand and l_s <= tss_abs <= l_e and chrom in chrom_len:
            mapper = OrientedCoords(l_s, l_e, strand)
            tss_rel = mapper.abs_to_rel(tss_abs)
            # Oriented index 0 is the upstream LTR edge, so U3 = [0, tss_rel)
            return Element(chrom=chrom, anchor_mode='U3', seq_start=l_s, seq_end=l_e,
                           seq=fetch_oriented(fa, chrom, l_s, l_e, strand),
                           tss_abs=tss_abs, tss_rel=tss_rel, zone_lo=0, zone_hi=tss_rel,
                           truncated=0, **common), 'U3'

    # TSS_window schema: fixed window from the genome, fixed core band as TA-search zone
    if not chrom:
        return None, 'no chrom'
    if chrom not in chrom_len:
        return None, 'contig missing from genome'
    if strand == '+':
        want_s, want_e = tss_abs - args.flank_up, tss_abs + args.flank_down
    else:
        want_s, want_e = tss_abs - args.flank_down, tss_abs + args.flank_up
    s1, e1 = max(1, want_s), min(chrom_len[chrom], want_e)
    mapper = OrientedCoords(s1, e1, strand)
    seq = fetch_oriented(fa, chrom, s1, e1, strand)
    tss_rel = mapper.abs_to_rel(tss_abs)
    zone_lo = max(0, tss_rel - args.core_up)
    zone_hi = min(len(seq), tss_rel + args.core_down)
    return Element(chrom=chrom, anchor_mode='TSS_window', seq_start=s1, seq_end=e1, seq=seq,
                   tss_abs=tss_abs, tss_rel=tss_rel, zone_lo=zone_lo, zone_hi=zone_hi,
                   truncated=int((s1, e1) != (want_s, want_e)), **common), 'TSS_window'


# ------------------------------
# Argument parsing
# ------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='TSS-anchored promoter motif scanner. LTR_structural elements with a TSS inside '
                    'an LTR use the U3 schema; all other classes use a fixed genome window around '
                    'the TSS. Outputs an SQLite database (-db) of motif hits and metadata.'
    )
    # Input files and output filenames
    p.add_argument('-t', '--tss-summary', required=True,
                   help='IsoClassifier *_{sense,antisense}.tss_summary.tsv '
                        '(Feature, Class, Orientation, Strand, Chrom, TSS1 0-based)')
    p.add_argument('-g', '--genome-fasta', default=DEFAULT_GENOME,
                   help='faidx-indexed genome FASTA (default: B73v5 PLATINUM)')
    p.add_argument('-l', '--ltr-fasta', default=None,
                   help='IsoClassifier *_ltr_seqs_<orient>.fa; supplies LTR bounds for the U3 schema. '
                        'Without it every element uses the TSS_window schema')
    p.add_argument('-u3', '--u3-fasta', default=None,
                   help='Deprecated, ignored: U3 is derived as upstream LTR edge .. TSS1 so it always '
                        'agrees with the scan anchor')
    p.add_argument('-i', '--isoforms', default=None,
                   help='IsoClassifier *_{sense,antisense}.isoforms.tsv; supplies Chrom when the TSS '
                        'summary predates the Chrom column')
    p.add_argument('-db', '--database', required=True, help='Output SQLite database file')

    # Element selection and TSS_window geometry
    p.add_argument('--orientation', choices=['sense', 'antisense'], default='sense',
                   help='Orientation rows to scan (default sense)')
    p.add_argument('--classes', default=','.join(ALL_CLASSES),
                   help=f'Comma-separated element classes to scan (default: {",".join(ALL_CLASSES)})')
    p.add_argument('--flank-up', type=int, default=500,
                   help='TSS_window: bp upstream of TSS to scan (default 500; CCAAT band needs >=460)')
    p.add_argument('--flank-down', type=int, default=200,
                   help='TSS_window: bp downstream of TSS to scan (default 200; SEC_TATA needs >=120)')
    p.add_argument('--core-up', type=int, default=100,
                   help='TSS_window: TA-rich zone starts this many bp upstream of TSS (default 100)')
    p.add_argument('--core-down', type=int, default=0,
                   help='TSS_window: TA-rich zone ends this many bp downstream of TSS (default 0)')

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
    el_cols = {r[1] for r in c.execute("PRAGMA table_info(elements)")}
    if 'anchor_mode' not in el_cols:
        raise SystemExit(f"{args.database} has the pre-schema-v2 elements table; write to a new DB path")
    
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


    # Record run parameters so each DB documents its schema geometry
    for k in ('tss_summary', 'genome_fasta', 'ltr_fasta', 'orientation', 'classes',
              'flank_up', 'flank_down', 'core_up', 'core_down', 'window_size', 'step_size',
              'ta_threshold', 'tata_max_dist', 'inr_max_dist'):
        c.execute("INSERT OR REPLACE INTO run_params VALUES (?, ?)", (k, str(getattr(args, k))))

    if args.u3_fasta:
        logging.warning("--u3-fasta is ignored; U3 is derived from the LTR bounds and TSS1")

    # Load inputs
    classes = {x.strip() for x in args.classes.split(',') if x.strip()}
    rows = load_tss_summary(args.tss_summary, args.orientation, classes,
                            load_chrom_map(args.isoforms))
    logging.info(f"Loaded {len(rows)} {args.orientation} TSS entries")
    ltr_coords = load_ltr_coords(args.ltr_fasta, args.orientation)
    logging.info(f"Loaded LTR coords for {len(ltr_coords)} features")
    fa = pysam.FastaFile(args.genome_fasta)
    chrom_len = dict(zip(fa.references, fa.lengths))

    total_tests = 0
    outcome = {}

    for row in rows:
        el, reason = build_element(row, fa, chrom_len, ltr_coords, args)
        key = (row['cls'], reason)
        outcome[key] = outcome.get(key, 0) + 1
        if el is None:
            continue

        feat, seq, tss_rel = el.feature, el.seq, el.tss_rel
        mapper = OrientedCoords(el.seq_start, el.seq_end, el.strand)

        if el.zone_hi > el.zone_lo:
            zone_s_abs, zone_e_abs = mapper.rel_to_abs_interval(el.zone_lo, el.zone_hi - el.zone_lo)
        else:
            zone_s_abs = zone_e_abs = None
        u3_s_abs, u3_e_abs = (zone_s_abs, zone_e_abs) if el.anchor_mode == 'U3' else (None, None)

        c.execute('''
            INSERT OR REPLACE INTO elements
            (feature, element_class, orientation, anchor_mode, chrom, seq_start, seq_end, strand,
             u3_start, u3_end, zone_start, zone_end, tss_abs, tss_rel, total_reads, tss_called,
             seq_truncated, sequence, sequence_length)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (feat, el.cls, el.orientation, el.anchor_mode, el.chrom, el.seq_start, el.seq_end,
              el.strand, u3_s_abs, u3_e_abs, zone_s_abs, zone_e_abs, el.tss_abs, tss_rel,
              el.total_reads, el.tss_called, el.truncated, seq, len(seq)))

        # TA-rich regions (restricted to U3 or the core band)
        ta_windows = []
        for s, e, sub in sliding_windows(seq, args.window_size, args.step_size,
                                         start=el.zone_lo, end=el.zone_hi):
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

        def in_ta_region(s_rel: int, e_rel: int) -> bool:
            for a, b in ta_regions:
                if s_rel >= a and e_rel <= b:
                    return True
            return False

        # CA runs across the whole scanned sequence
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
        def scan_band(motif, pwm_key, band_start, band_end, dist_fn, use_ta):
            """Scan one TSS-relative band and store every hit above the scan threshold."""
            pwm = pwms[pwm_key]
            L = pwm_len(pwm)
            for s_rel, e_rel, score, subseq in scan_all_pwm_hits(seq, pwm, bg, band_start,
                                                                 band_end, scan_thr[motif]):
                pval = score_to_pvalue(bg_scores[pwm_key], score)
                s_abs, e_abs = mapper.rel_to_abs_interval(s_rel, L)
                in_ta = 1 if use_ta and in_ta_region(s_rel, e_rel) else 0
                yc = sum(1 for b in subseq if b in 'CT')
                c.execute('''
                    INSERT INTO motif_hits
                    (feature, motif_type, motif_length, start_rel, end_rel, start_abs, end_abs,
                     score, sequence, dist_to_tss, in_ta_region, y_count, search_window_start, search_window_end, p_value)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (feat, motif, L, s_rel, e_rel, s_abs, e_abs,
                      score, subseq, dist_fn(s_rel), in_ta, yc, band_start, band_end, pval))

        upstream = lambda s_rel: tss_rel - s_rel
        downstream = lambda s_rel: s_rel - tss_rel

        # 1) Primary TATA upstream (band extends 8 bp past TSS)
        scan_band('TATA', 'TATA', max(0, tss_rel - args.tata_max_dist), max(0, tss_rel + 8),
                  upstream, True)
        # 2) CCAAT-box upstream 460-140 bp
        scan_band('CCAAT', 'CCAAT', max(0, tss_rel - 460), max(0, tss_rel - 140), upstream, False)
        # 3) Y-PATCH upstream 100-1 bp, 8-mer
        scan_band('YPATCH', 'YPATCH', max(0, tss_rel - 100), max(0, tss_rel), upstream, False)
        # 4) Inr motif near TSS
        scan_band('INR7', 'INR7', max(0, tss_rel - args.inr_max_dist),
                  min(len(seq), tss_rel + args.inr_max_dist + 1),
                  lambda s_rel: abs(s_rel - tss_rel), False)
        # 5) Secondary TATA up to 120 bp downstream
        scan_band('SEC_TATA', 'TATA', tss_rel + 1, min(len(seq), tss_rel + 120), downstream, True)
        # 6) DPE motif 0-50 bp downstream
        scan_band('DPE7', 'DPE7', tss_rel, min(len(seq), tss_rel + 50), downstream, False)

    fa.close()
    for (cls, reason), n in sorted(outcome.items()):
        logging.info(f"  {cls:15s} {reason:28s} {n}")

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
