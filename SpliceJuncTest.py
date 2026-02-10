#!/usr/bin/env python3
"""
check_splice_junctions_from_gff.py

Given:
  - a splice-aware BAM (e.g. minimap2 -ax splice),
  - a reference FASTA,
  - a GFF with loci (genes, transcripts, or generic regions),

this script:
  * randomly samples up to N loci from the GFF,
  * finds introns (CIGAR N) within those loci,
  * checks their splice motifs on the reference,
  * reports what fraction of those junctions are canonical.

Usage example:
  python check_splice_junctions_from_gff.py \
    -b aln.sorted.bam \
    -f genome.fa \
    -g loci.gff \
    -n 100 \
    -m 3 \
    > junctions_per_locus.tsv

The summary with overall % canonical is printed at the end (prefixed by '#').
"""

import argparse
import random
from collections import Counter
import numpy as np
import scipy.stats as stats
import pysam


def parse_args():
    p = argparse.ArgumentParser(
        description="Check splice motifs for introns within loci from a GFF."
    )
    p.add_argument("-b", "--bam", required=True,
                   help="Coordinate-sorted, indexed BAM file")
    p.add_argument("-f", "--fasta", required=True,
                   help="Reference genome FASTA (indexed with samtools faidx)")
    p.add_argument("-g", "--gff", required=True,
                   help="GFF file with loci (genes/transcripts/regions)")
    p.add_argument("-n", "--num-loci", type=int, default=100,
                   help="Number of loci to sample (default: 100)")
    p.add_argument("-m", "--min-support", type=int, default=1,
                   help="Minimum read support per intron (default: 1)")
    p.add_argument("--seed", type=int, default=1,
                   help="Random seed for sampling loci (default: 1)")
    return p.parse_args()


# CIGAR code for N (ref skip; intron)
CIGAR_N = 3

# Canonical motifs (on the genomic strand)
CANONICAL_MOTIFS = {
    "GT-AG", "GC-AG", "AT-AC",    # forward
    "CT-AC", "CT-GC", "GT-AT"     # reverse complements (for genes on - strand)
}


def parse_gff(gff_path):
    """
    Parse GFF and return a list of loci:
    (chrom, start0, end0, locus_id)

    GFF is 1-based inclusive; convert to 0-based half-open for comparison.
    """
    loci = []
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature_type, start, end, _, _, _, attrs = fields
            start = int(start)
            end = int(end)

            # 0-based half-open
            start0 = start - 1
            end0 = end

            # Extract ID/Name for helpful labeling
            locus_id = f"{feature_type}:{chrom}:{start}-{end}"
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    locus_id = attr.split("=", 1)[1]
                    break
                if attr.startswith("Name="):
                    locus_id = attr.split("=", 1)[1]
                    # don't break, ID takes priority if present
            loci.append((chrom, start0, end0, locus_id))
    return loci


def get_motif(fasta, chrom, start, end):
    """
    Get 2bp donor and acceptor motifs on the reference (0-based, end-exclusive).

    Returns motif like 'GT-AG' or None if out-of-bounds / chrom not found.
    """
    if start < 0 or end - 2 < 0:
        return None
    try:
        donor = fasta.fetch(chrom, start, start + 2).upper()
        acceptor = fasta.fetch(chrom, end - 2, end).upper()
    except KeyError:
        return None
    if len(donor) != 2 or len(acceptor) != 2:
        return None
    return f"{donor}-{acceptor}"


def collect_introns_in_locus(bam, chrom, locus_start, locus_end, min_support):
    """
    For a single locus, collect introns fully contained in [locus_start, locus_end)
    and return a dict: (start, end) -> read_support.
    """
    intron_counts = Counter()

    # Fetch only alignments overlapping the locus
    for aln in bam.fetch(chrom, locus_start, locus_end):
        if aln.is_unmapped or aln.cigartuples is None:
            continue

        ref_pos = aln.reference_start  # 0-based
        for op, length in aln.cigartuples:
            if op in (0, 2, 7, 8):  # M, D, =, X (consume reference)
                ref_pos += length
            elif op == CIGAR_N:     # intron
                intron_start = ref_pos
                intron_end = ref_pos + length

                # require intron fully within locus
                if intron_start >= locus_start and intron_end <= locus_end:
                    intron_counts[(intron_start, intron_end)] += 1

                ref_pos += length
            else:
                # I, S, H, P etc. do not consume reference
                continue

    # Filter by support
    return {k: c for k, c in intron_counts.items() if c >= min_support}


def main():
    args = parse_args()
    random.seed(args.seed)

    loci = parse_gff(args.gff)
    if not loci:
        raise SystemExit("No loci parsed from GFF.")

    if len(loci) > args.num_loci:
        loci = random.sample(loci, args.num_loci)

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)

    print("\t".join([
        "locus_id", "chrom",
        "locus_start", "locus_end",
        "intron_start", "intron_end",
        "intron_length", "read_support",
        "motif", "is_canonical"
    ]))

    total_introns = 0
    canonical_introns = 0
    canonical_supports = []
    noncanonical_supports = []

    for chrom, locus_start, locus_end, locus_id in loci:
        introns = collect_introns_in_locus(
            bam, chrom, locus_start, locus_end, args.min_support
        )

        for (istart, iend), support in sorted(introns.items()):
            motif = get_motif(fasta, chrom, istart, iend)
            if motif is None:
                is_canonical = "NA"
            else:
                is_canonical = "yes" if motif in CANONICAL_MOTIFS else "no"

            # Update global stats (ignore NA in stats)
            total_introns += 1
            if is_canonical == "yes":
                canonical_introns += 1
                canonical_supports.append(support)
            elif is_canonical == "no":
                noncanonical_supports.append(support)

            print("\t".join(map(str, [
                locus_id,
                chrom,
                locus_start + 1,  # back to 1-based
                locus_end,
                istart,
                iend,
                iend - istart,
                support,
                motif if motif is not None else "NA",
                is_canonical
            ])))

    bam.close()
    fasta.close()


    # Summary
    if total_introns > 0:
        pct = 100.0 * canonical_introns / total_introns
        print(
            f"# Summary: {canonical_introns}/{total_introns} "
            f"({pct:.2f}%) introns have canonical splice motifs "
            f"(min_support >= {args.min_support})"
        )
        # Basic summaries
        print("# Canonical median:", np.median(canonical_supports))
        print("# Noncanonical median:", np.median(noncanonical_supports))

        # Wilcoxon rank-sum test (non-parametric)
        stat, pvalue = stats.ranksums(canonical_supports, noncanonical_supports)
        print(f"# Wilcoxon test p={pvalue}, test statistic={stat}")
    else:
        print("# Summary: No introns found in sampled loci with the given settings.")


if __name__ == "__main__":
    main()
