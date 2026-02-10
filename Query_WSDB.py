#!/usr/bin/env python3
"""
Query interface for motif hit database.
Provides flexible filtering and selection of "best" motifs per element.
The -db should be the output of WindowScrubber.py, which contains all motif hits and metadata.

Usage examples:

1. Export best TATA per element (highest score):
   python Query_WSDB.py -db motif_hits.db --best-per-element TATA --order-by score --output best_tata.tsv

2. Export all TATA hits in TA-rich regions:
   python Query_WSDB.py -db motif_hits.db --motif-type TATA --in-ta-region --output tata_in_ta.tsv

3. Export TATA hits within specific distance range:
   python Query_WSDB.py -db motif_hits.db --motif-type TATA --min-dist 20 --max-dist 100 --output tata_filtered.tsv

4. Get statistics on motif distribution:
   python Query_WSDB.py -db motif_hits.db --stats

5. Export best of each motif type per element:
   python Query_WSDB.py -db motif_hits.db --best-all-motifs --output best_all.tsv

6. Export elements with complete motif set:
   python Query_WSDB.py -db motif_hits.db --require-motifs TATA,CCAAT,INR7 --output complete_promoters.tsv

7. Get CA run summary:
   python Query_WSDB.py -db motif_hits.db --ca-summary --output ca_summary.tsv
   
"""

import argparse
import csv
import sqlite3
import sys
from typing import List, Optional, Dict, Any


def get_connection(db_path: str) -> sqlite3.Connection:
    """Connect to database and set row factory."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def query_best_per_element(conn: sqlite3.Connection, motif_type: str, 
                           order_by: str = 'score', 
                           ascending: bool = False,
                           filters: Optional[Dict[str, Any]] = None,
                           range_pair: Optional[tuple] = None,
                           p_max: Optional[float] = None,
                           thr_method: Optional[str] = None,
                           thr_param: Optional[str] = None,
                           apply_threshold: bool = False) -> List[sqlite3.Row]:    
    """
    Get the best motif hit per element for a given motif type.
    
    Args:
        motif_type: Type of motif (e.g., 'TATA', 'CCAAT', 'INR7', 'YPATCH')
        order_by: Column to use for ranking ('score', 'dist_to_tss', 'y_count')
        ascending: If True, pick lowest value; if False, pick highest
        filters: Optional dict of additional filters (e.g., {'in_ta_region': 1})
        range_pair: Optional (lo, hi) relative distance to TSS
        p_max: Optional max p-value threshold
        thr_method: Threshold method from database
        thr_param: Threshold param from database
        apply_threshold: Whether to apply stored threshold
    """
    c = conn.cursor()

    where_clauses = ["m.motif_type = ?"]
    params = [motif_type]

    if filters:
        for key, value in filters.items():
            where_clauses.append(f"m.{key} = ?")
            params.append(value)

    if range_pair:
        lo, hi = range_pair
        where_clauses.append("m.dist_to_tss BETWEEN ? AND ?")
        params.extend([lo, hi])

    if p_max is not None:
        where_clauses.append("m.p_value <= ?")
        params.append(p_max)

    if apply_threshold:
        cutoff = fetch_threshold(conn, motif_type, thr_method, thr_param)
        if cutoff is None:
            raise ValueError(f"No threshold found for motif={motif_type}, method={thr_method}, param={thr_param}")
        where_clauses.append("m.score >= ?")
        params.append(cutoff)

    where_str = " AND ".join(where_clauses)
    order_dir = "ASC" if ascending else "DESC"

    query = f"""
        WITH ranked AS (
            SELECT m.*,
                   ROW_NUMBER() OVER (
                       PARTITION BY m.feature 
                       ORDER BY m.{order_by} {order_dir}, 
                                ABS(m.dist_to_tss) ASC,
                                m.start_rel ASC
                   ) as rn
            FROM motif_hits m
            WHERE {where_str}
        )
        SELECT * FROM ranked WHERE rn = 1
        ORDER BY feature
    """
    return c.execute(query, params).fetchall()


def query_all_hits(conn: sqlite3.Connection, 
                   motif_type: Optional[str] = None,
                   min_score: Optional[float] = None,
                   max_score: Optional[float] = None,
                   min_dist: Optional[int] = None,
                   max_dist: Optional[int] = None,
                   in_ta_region: Optional[bool] = None,
                   features: Optional[List[str]] = None,
                   p_max: Optional[float] = None,
                   thr_method: Optional[str] = None,
                   thr_param: Optional[str] = None,
                   apply_threshold: bool = False) -> List[sqlite3.Row]:
    """
    Query all motif hits with flexible filtering.
    """
    c = conn.cursor()
    
    where_clauses = []
    params = []
    
    if p_max is not None:
        where_clauses.append("m.p_value <= ?")
        params.append(p_max)

    if apply_threshold and motif_type:
        cutoff = fetch_threshold(conn, motif_type, thr_method, thr_param)
        if cutoff is None:
            raise ValueError(f"No threshold found for motif={motif_type}, method={thr_method}, param={thr_param}")
        where_clauses.append("m.score >= ?")
        params.append(cutoff)
    
    if motif_type:
        where_clauses.append("m.motif_type = ?")
        params.append(motif_type)
    
    if min_score is not None:
        where_clauses.append("m.score >= ?")
        params.append(min_score)
    
    if max_score is not None:
        where_clauses.append("m.score <= ?")
        params.append(max_score)
    
    if min_dist is not None:
        where_clauses.append("ABS(m.dist_to_tss) >= ?")
        params.append(min_dist)
    
    if max_dist is not None:
        where_clauses.append("ABS(m.dist_to_tss) <= ?")
        params.append(max_dist)
    
    if in_ta_region is not None:
        where_clauses.append("m.in_ta_region = ?")
        params.append(1 if in_ta_region else 0)
    
    if features:
        placeholders = ",".join("?" * len(features))
        where_clauses.append(f"m.feature IN ({placeholders})")
        params.extend(features)
    
    where_str = " AND ".join(where_clauses) if where_clauses else "1=1"
    
    query = f"""
        SELECT m.*, e.chrom, e.strand, e.tss_abs
        FROM motif_hits m
        JOIN elements e ON m.feature = e.feature
        WHERE {where_str}
        ORDER BY m.feature, m.motif_type, m.score DESC
    """
    
    return c.execute(query, params).fetchall()


def query_best_all_motifs(conn,
                          order_by: str = 'score',
                          p_max: Optional[float] = None,
                          apply_threshold: bool = False,
                          thr_method: Optional[str] = None,
                          thr_param: Optional[str] = None,
                          ranges_by_motif: Optional[Dict[str, tuple]] = None
                          ) -> List[sqlite3.Row]:
    """
    Query to get best motif hit per element for each motif type,
    with optional motif-specific distance ranges and thresholds.
    YPATCH is now a single 8-mer motif type.
    """
    c = conn.cursor()
    c.execute("PRAGMA temp_store = MEMORY")

    # Discover motif types from the database
    motif_types = [row[0] for row in c.execute(
        "SELECT DISTINCT motif_type FROM motif_hits"
    ).fetchall()]

    # Create temp table for ranges if provided
    if ranges_by_motif:
        c.execute("CREATE TEMP TABLE IF NOT EXISTS _ranges (motif_type TEXT PRIMARY KEY, lo INT, hi INT)")
        c.execute("DELETE FROM _ranges")
        c.executemany("INSERT INTO _ranges VALUES (?,?,?)",
                      [(m, int(lo), int(hi)) for m, (lo, hi) in ranges_by_motif.items()])

    # Build threshold join if needed
    thr_join = ""
    if apply_threshold:
        c.execute("CREATE TEMP TABLE IF NOT EXISTS _thr (motif_type TEXT PRIMARY KEY, cutoff REAL)")
        c.execute("DELETE FROM _thr")
        if thr_param:
            c.execute("""INSERT INTO _thr
                         SELECT motif_type, cutoff
                         FROM thresholds
                         WHERE method=? AND param=?""", (thr_method, thr_param))
        else:
            c.execute("""INSERT INTO _thr
                         SELECT motif_type, MIN(cutoff) AS cutoff
                         FROM thresholds
                         WHERE method=?
                         GROUP BY motif_type""", (thr_method,))
        thr_join = "JOIN _thr t ON t.motif_type = m.motif_type"

    # Build range join if needed
    range_join = "JOIN _ranges r ON r.motif_type = m.motif_type" if ranges_by_motif else ""

    # Build WHERE clause for candidates
    where_bits, params = [], []
    if p_max is not None:
        where_bits.append("m.p_value <= ?")
        params.append(p_max)
    if apply_threshold:
        where_bits.append("m.score >= t.cutoff")
    if ranges_by_motif:
        where_bits.append("m.dist_to_tss BETWEEN r.lo AND r.hi")
    where_str = " AND ".join(where_bits) if where_bits else "1=1"

    # Build CTE with best hit per element per motif
    cte = f"""
    WITH candidates AS (
        SELECT m.*
        FROM motif_hits m
        {thr_join}
        {range_join}
        WHERE {where_str}
    ),
    best AS (
        SELECT *,
               ROW_NUMBER() OVER (
                   PARTITION BY feature, motif_type
                   ORDER BY score DESC, ABS(dist_to_tss) ASC, start_rel ASC
               ) AS rn
        FROM candidates
    ),
    picked AS (
        SELECT * FROM best WHERE rn = 1
    )
    """

    # Build joins and select columns for each motif type
    join_blocks = []
    select_blocks = []
    for m in motif_types:
        alias = f"p_{m}"
        join_blocks.append(
            f"LEFT JOIN picked {alias} "
            f"  ON {alias}.feature = e.feature AND {alias}.motif_type = '{m}'"
        )
        select_blocks.append(
            f"COALESCE({alias}.start_abs, -1) AS {m}_start_abs, "
            f"COALESCE({alias}.end_abs,   -1) AS {m}_end_abs, "
            f"COALESCE({alias}.dist_to_tss, -1) AS Dist_TSS_to_{m}"
        )

    final_sql = f"""
    {cte}
    SELECT
        e.feature,
        e.tss_abs,
        CASE WHEN tr.feature IS NULL THEN 0 ELSE 1 END AS TA_rich_present,
        {", ".join(select_blocks)}
    FROM elements e
    LEFT JOIN ta_regions tr ON tr.feature = e.feature
    {' '.join(join_blocks)}
    GROUP BY e.feature, e.tss_abs
    ORDER BY e.feature
    """

    return c.execute(final_sql, params).fetchall()


def query_elements_with_motifs(conn: sqlite3.Connection, 
                               required_motifs: List[str]) -> List[sqlite3.Row]:
    """
    Find elements that have at least one hit for each required motif type.
    """
    c = conn.cursor()
    
    placeholders = ",".join("?" * len(required_motifs))
    
    query = f"""
        SELECT e.*, 
               GROUP_CONCAT(DISTINCT m.motif_type) as found_motifs
        FROM elements e
        JOIN motif_hits m ON e.feature = m.feature
        WHERE m.motif_type IN ({placeholders})
        GROUP BY e.feature
        HAVING COUNT(DISTINCT m.motif_type) = ?
        ORDER BY e.feature
    """
    
    return c.execute(query, required_motifs + [len(required_motifs)]).fetchall()


def get_statistics(conn: sqlite3.Connection) -> Dict[str, Any]:
    """
    Get summary statistics about the database.
    """
    c = conn.cursor()
    
    stats = {}
    
    # Element count
    stats['total_elements'] = c.execute("SELECT COUNT(*) FROM elements").fetchone()[0]
    
    # Motif hit counts by type
    stats['hits_by_motif'] = {}
    for row in c.execute("""
        SELECT motif_type, COUNT(*) as count, 
               AVG(score) as avg_score, 
               MIN(score) as min_score, 
               MAX(score) as max_score
        FROM motif_hits
        GROUP BY motif_type
    """).fetchall():
        stats['hits_by_motif'][row['motif_type']] = {
            'count': row['count'],
            'avg_score': row['avg_score'],
            'min_score': row['min_score'],
            'max_score': row['max_score']
        }
    
    # Elements with each motif type
    stats['elements_with_motif'] = {}
    for row in c.execute("""
        SELECT motif_type, COUNT(DISTINCT feature) as element_count
        FROM motif_hits
        GROUP BY motif_type
    """).fetchall():
        stats['elements_with_motif'][row['motif_type']] = row['element_count']
    
    # TA region stats
    stats['ta_regions'] = c.execute("SELECT COUNT(*) FROM ta_regions").fetchone()[0]
    stats['elements_with_ta'] = c.execute("SELECT COUNT(DISTINCT feature) FROM ta_regions").fetchone()[0]
    
    # CA run stats
    stats['total_ca_runs'] = c.execute("SELECT COUNT(*) FROM ca_runs").fetchone()[0]
    stats['significant_ca_runs'] = c.execute("SELECT COUNT(*) FROM ca_runs WHERE is_significant = 1").fetchone()[0]
    
    # Motifs in TA regions
    stats['motifs_in_ta'] = {}
    for row in c.execute("""
        SELECT motif_type, COUNT(*) as count
        FROM motif_hits
        WHERE in_ta_region = 1
        GROUP BY motif_type
    """).fetchall():
        stats['motifs_in_ta'][row['motif_type']] = row['count']
    
    return stats


def query_ca_runs(conn: sqlite3.Connection, 
                  significant_only: bool = True,
                  min_length: Optional[int] = None) -> List[sqlite3.Row]:
    """
    Query CA runs with optional filtering.
    """
    c = conn.cursor()
    
    where_clauses = []
    params = []
    
    if significant_only:
        where_clauses.append("is_significant = 1")
    
    if min_length:
        where_clauses.append("length >= ?")
        params.append(min_length)
    
    where_str = " AND ".join(where_clauses) if where_clauses else "1=1"
    
    query = f"""
        SELECT c.*, e.chrom, e.strand
        FROM ca_runs c
        JOIN elements e ON c.feature = e.feature
        WHERE {where_str}
        ORDER BY c.length DESC, c.p_value ASC
    """
    
    return c.execute(query, params).fetchall()


def export_to_tsv(rows: List[sqlite3.Row], output_path: str, columns: Optional[List[str]] = None):
    """
    Export query results to TSV file.
    """
    if not rows:
        print("No results to export.", file=sys.stderr)
        return
    
    if columns is None:
        columns = rows[0].keys()
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))
    
    print(f"Exported {len(rows)} rows to {output_path}")


def print_statistics(stats: Dict[str, Any]):
    """
    Print statistics in readable format.
    """
    print("\n=== Database Statistics ===\n")
    print(f"Total elements: {stats['total_elements']}")
    print(f"TA-rich regions: {stats['ta_regions']} (in {stats['elements_with_ta']} elements)")
    print(f"CA runs: {stats['total_ca_runs']} ({stats['significant_ca_runs']} significant)")
    
    print("\n--- Motif Hits by Type ---")
    for motif_type, info in stats['hits_by_motif'].items():
        elem_count = stats['elements_with_motif'].get(motif_type, 0)
        in_ta = stats['motifs_in_ta'].get(motif_type, 0)
        print(f"\n{motif_type}:")
        print(f"  Total hits: {info['count']}")
        print(f"  Elements with motif: {elem_count}")
        print(f"  Hits in TA regions: {in_ta}")
        print(f"  Score range: {info['min_score']:.2f} - {info['max_score']:.2f}")
        print(f"  Average score: {info['avg_score']:.2f}")


def parse_ranges_arg(ranges_str):
    """
    Parse "TATA:-100:0,CCAAT:-460:-140" into dict {motif: (lo, hi)}
    """
    out = {}
    if not ranges_str:
        return out
    for spec in ranges_str.split(','):
        spec = spec.strip()
        if not spec: 
            continue
        parts = spec.split(':')
        if len(parts) != 3:
            raise ValueError(f"Invalid range spec: {spec}. Expected format: MOTIF:lo:hi")
        name, lo, hi = parts
        out[name] = (int(lo), int(hi))
    return out


def fetch_threshold(conn, motif_type: str, method: str, param: Optional[str]) -> Optional[float]:
    """
    Return cutoff for (motif_type, method, param) from thresholds table.
    If param is None, returns the most permissive (MIN) cutoff for that method.
    """
    c = conn.cursor()
    if param:
        row = c.execute("""
            SELECT cutoff FROM thresholds
            WHERE motif_type=? AND method=? AND param=?
        """, (motif_type, method, param)).fetchone()
        return float(row['cutoff']) if row else None
    else:
        row = c.execute("""
            SELECT MIN(cutoff) AS cutoff FROM thresholds
            WHERE motif_type=? AND method=?
        """, (motif_type, method)).fetchone()
        return float(row['cutoff']) if row and row['cutoff'] is not None else None


def parse_args():
    p = argparse.ArgumentParser(description='Query motif hit database')
    p.add_argument('-db', '--database', required=True, help='SQLite database file')
    p.add_argument('-o', '--output', help='Output TSV file')
    
    # Query modes
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--best-per-element', metavar='MOTIF', help='Get best hit per element for motif type')
    mode.add_argument('--best-all-motifs', action='store_true', help='Get best hit for each motif type per element')
    mode.add_argument('--all-hits', action='store_true', help='Query all hits with filters')
    mode.add_argument('--require-motifs', metavar='MOTIF1,MOTIF2', help='Elements with all specified motifs')
    mode.add_argument('--stats', action='store_true', help='Print database statistics')
    mode.add_argument('--ca-summary', action='store_true', help='Export CA run summary')

    # Filters
    p.add_argument('--motif-type', help='Filter by motif type (e.g., TATA, CCAAT, YPATCH, INR7, DPE7, SEC_TATA)')
    p.add_argument('--min-score', type=float, help='Minimum PWM score')
    p.add_argument('--max-score', type=float, help='Maximum PWM score')
    p.add_argument('--min-dist', type=int, help='Minimum absolute distance to TSS')
    p.add_argument('--max-dist', type=int, help='Maximum absolute distance to TSS')
    p.add_argument('--in-ta-region', action='store_true', help='Only hits in TA-rich regions')
    p.add_argument('--features', help='Comma-separated list of features to include')

    # Range/threshold controls
    p.add_argument('--range', help='Single motif relative lo:hi for distance to TSS (e.g. "-100:0")')
    p.add_argument('--p-max', type=float, help='Keep hits with p_value <= this (e.g. 1e-4)')
    p.add_argument('--apply-threshold', action='store_true',
                help='If set, require score >= threshold(method,param) from the thresholds table')
    p.add_argument('--ranges', help='Comma-separated motif:lo:hi specs (e.g. "TATA:-100:0,CCAAT:-460:-140,INR7:-10:10,DPE7:0:50,YPATCH:-100:0")')
    p.add_argument('--threshold-method', default='empirical_p',
                help='Which stored threshold to apply (e.g., empirical_p, empirical_percentile, mismatch_anywhere, relative, store)')
    p.add_argument('--threshold-param', required=False, default=None,
                help='Exact param string (e.g., "alpha=1e-4;bg=iid;n=200000" or "perc=95.0"). If omitted, use most permissive for method')
    
    # Ordering for "best" queries
    p.add_argument('--order-by', default='score', choices=['score', 'dist_to_tss', 'y_count'],
                   help='Column to rank by for best-per-element')
    p.add_argument('--ascending', action='store_true', help='Pick lowest value instead of highest')
    
    # CA run filters
    p.add_argument('--ca-min-length', type=int, help='Minimum CA run length')
    p.add_argument('--ca-all', action='store_true', help='Include non-significant CA runs')
    
    return p.parse_args()


def main():
    args = parse_args()
    
    conn = get_connection(args.database)
    
    try:
        if args.stats:
            stats = get_statistics(conn)
            print_statistics(stats)
            return
        
        results = None
        
        if args.best_per_element:
            filters = {}
            if args.in_ta_region:
                filters['in_ta_region'] = 1

            range_pair = None
            if args.range:
                parts = args.range.split(':')
                if len(parts) != 2:
                    print("Error: --range must be in format 'lo:hi'", file=sys.stderr)
                    return
                range_pair = (int(parts[0]), int(parts[1]))

            results = query_best_per_element(
                conn,
                args.best_per_element,
                order_by=args.order_by,
                ascending=args.ascending,
                filters=filters,
                range_pair=range_pair,
                p_max=args.p_max,
                thr_method=args.threshold_method,
                thr_param=args.threshold_param,
                apply_threshold=args.apply_threshold
            )

        elif args.best_all_motifs:
            ranges_dict = parse_ranges_arg(args.ranges) if args.ranges else None
            results = query_best_all_motifs(
                conn,
                order_by=args.order_by,
                p_max=args.p_max,
                apply_threshold=args.apply_threshold,
                thr_method=args.threshold_method,
                thr_param=args.threshold_param,
                ranges_by_motif=ranges_dict
            )
        
        elif args.all_hits:
            features = args.features.split(',') if args.features else None
            results = query_all_hits(
                conn,
                motif_type=args.motif_type,
                min_score=args.min_score,
                max_score=args.max_score,
                min_dist=args.min_dist,
                max_dist=args.max_dist,
                in_ta_region=args.in_ta_region if args.in_ta_region else None,
                features=features,
                p_max=args.p_max,
                thr_method=args.threshold_method,
                thr_param=args.threshold_param,
                apply_threshold=args.apply_threshold
            )
        
        elif args.require_motifs:
            required = args.require_motifs.split(',')
            results = query_elements_with_motifs(conn, required)
        
        elif args.ca_summary:
            results = query_ca_runs(conn, 
                                   significant_only=not args.ca_all,
                                   min_length=args.ca_min_length)
        
        else:
            print("Please specify a query mode. Use --help for options.", file=sys.stderr)
            return
        
        if results is not None:
            if args.output:
                export_to_tsv(results, args.output)
            else:
                # Print to stdout
                print(f"\nFound {len(results)} results\n")
                if results:
                    for i, row in enumerate(results[:10]):  # Show first 10
                        print(dict(row))
                    if len(results) > 10:
                        print(f"\n... and {len(results) - 10} more. Use -o to export all.")
    
    finally:
        conn.close()


if __name__ == '__main__':
    main()