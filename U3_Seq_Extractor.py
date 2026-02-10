#!/usr/bin/env python3
"""
Unified extractor for U3 regions and gene promoters.

Modes
=====
1) LTR mode (mode=LTR)
   - Extracts U3 regions upstream of primary TSS from LTR elements
     and saves both U3 and full LTR sequences.

   Inputs:
     --summary-tsv : TSV with columns: Feature, TSS1, Top_Isoform_Strand
                     (Top_Isoform_Strand is used as fallback strand)
     --ltr-gff     : GFF with lLTR/rLTR annotations and Parent=...
     --genome-fasta: Reference genome FASTA
     --u3-output   : Output FASTA for U3 sequences
     --ltr-output  : Output FASTA for full LTR sequences

2) Gene mode (mode=Gene)
   - Generates promoter regions (±1000 bp around TSS) for genes and
     outputs:
       * BED file for ±1000 bp around TSS (2 kb window)
       * FASTA of that BED (strand-aware, like bedtools getfasta -s)
       * Dummy U3 FASTA entries for the upstream 1 kb region only

   Inputs:
     --summary-tsv : TSV of TSS summary. Must contain:
                     - Gene or Feature column (identifier)
                     - TSS1 or TSS_abs column (TSS coordinate, 1-based)
     --gtf         : GTF/GFF with 'gene' features and gene_id/ID attribute
     --genome-fasta: Reference genome FASTA
     --promoter-bed-output   : BED file of ±1000 bp around TSS
     --promoter-fasta-output : FASTA for that BED (promoter sequences)
     --u3-output   : Dummy U3 FASTA for each promoter
                     (sequence is 'NNNNN', header encodes upstream 1 kb)

Usage examples:

# LTR mode
python3 U3_seq_extractor.py \
  --mode LTR \
  --summary-tsv summary.tsv \
  --ltr-gff annotations.gff \
  --genome-fasta genome.fa \
  --u3-output LTR_U3.fa \
  --ltr-output full_LTRs.fa

# Gene mode
python3 U3_seq_extractor.py \
  --mode Gene \
  --summary-tsv tss_summary.tsv \
  --gtf genes.gtf \
  --genome-fasta genome.fa \
  --promoter-bed-output Gene_2kb_proms.bed \
  --promoter-fasta-output Gene_2kb_proms.fa \
  --u3-output Gene_dummy_U3.fa
"""

import argparse
import logging
from typing import Dict, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ---------- Argument parsing ----------

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract U3 regions / promoters in LTR or Gene mode."
    )
    p.add_argument(
        "--mode",
        required=True,
        choices=["LTR", "Gene"],
        help="Operating mode: 'LTR' for LTR U3 extraction, 'Gene' for gene promoters."
    )
    # common summary TSV (used by both modes)
    p.add_argument(
        "--summary-tsv", "--tss",
        dest="summary_tsv",
        required=True,
        help="Summary TSV. "
             "LTR mode: Feature,TSS1,Top_Isoform_Strand. "
             "Gene mode: must have Gene/Feature and TSS1/TSS_abs."
    )

    # LTR-mode specific
    p.add_argument("--ltr-gff", help="LTR GFF with lLTR/rLTR and Parent=… (LTR mode)")
    p.add_argument("--genome-fasta", help="Reference genome FASTA (LTR or Gene mode)")
    p.add_argument("--ltr-output", help="Output FASTA for full LTR sequences (LTR mode)")

    # Gene-mode specific
    p.add_argument("--gtf", help="Gene annotation GTF/GFF (Gene mode)")
    p.add_argument(
        "--promoter-bed-output",
        help="Output BED file for ±1000 bp around TSS (Gene mode)"
    )
    p.add_argument(
        "--promoter-fasta-output",
        help="Output FASTA for promoter regions (Gene mode)"
    )

    # U3 FASTA (used in both modes, real in LTR, dummy in Gene)
    p.add_argument(
        "--u3-output",
        required=True,
        help="Output FASTA for U3 sequences "
             "(real for LTR mode, dummy for upstream 1 kb in Gene mode)."
    )

    args = p.parse_args()

    # Mode-specific sanity checks
    if args.mode == "LTR":
        missing = []
        if not args.ltr_gff:
            missing.append("--ltr-gff")
        if not args.genome_fasta:
            missing.append("--genome-fasta")
        if not args.ltr_output:
            missing.append("--ltr-output")
        if missing:
            raise SystemExit(
                f"ERROR: In LTR mode, the following arguments are required: {', '.join(missing)}"
            )
    elif args.mode == "Gene":
        missing = []
        if not args.gtf:
            missing.append("--gtf")
        if not args.genome_fasta:
            missing.append("--genome-fasta")
        if not args.promoter_bed_output:
            missing.append("--promoter-bed-output")
        if not args.promoter_fasta_output:
            missing.append("--promoter-fasta-output")
        if missing:
            raise SystemExit(
                f"ERROR: In Gene mode, the following arguments are required: {', '.join(missing)}"
            )

    return args


# ---------- Shared helpers ----------

def extract_region(chrom_seq: Seq, start: int, end: int, strand: str) -> Seq:
    """
    Extract 1-based inclusive [start, end] from chrom_seq and
    reverse-complement if strand == '-'.
    """
    seq = chrom_seq[start - 1:end]
    return seq.reverse_complement() if strand == '-' else seq


# ---------- LTR MODE IMPLEMENTATION ----------

def load_ltr_summary(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        tsv_path,
        sep='\t',
        usecols=["Feature", "TSS1", "Top_Isoform_Strand"],
        dtype={"Feature": str, "TSS1": int, "Top_Isoform_Strand": str}
    )
    df = df.rename(columns={"Top_Isoform_Strand": "Strand"})
    logging.info(f"[LTR] Loaded {len(df)} TSS entries from {tsv_path}")
    return df


def load_ltr_gff(gff_path: str) -> Dict[str, dict]:
    ltr_map: Dict[str, dict] = {}
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split('\t')
            if len(cols) < 9:
                continue
            chrom, start, end, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]
            attrs = dict(
                p.split("=", 1) for p in cols[8].split(";") if "=" in p
            )
            parent = attrs.get("Parent")
            if not parent:
                continue
            entry = ltr_map.setdefault(parent, {
                "chrom": chrom,
                "strand": None,
                "lLTR": None,
                "rLTR": None
            })
            eid = attrs.get("ID", "")
            # full-element LTRRT_ carries strand
            if eid.startswith("LTRRT_"):
                entry["strand"] = strand
            # identify left/right LTRs
            if "lLTR" in eid or "lLTR" in cols[8]:
                entry['lLTR'] = (start, end)
            elif "rLTR" in eid or "rLTR" in cols[8]:
                entry['rLTR'] = (start, end)
    logging.info(f"[LTR] Parsed {len(ltr_map)} LTR parent entries from GFF")
    return ltr_map


def load_genome(fasta_path: str):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))
    logging.info(f"Loaded {len(genome)} sequences from genome fasta")
    return genome


def run_ltr_mode(args):
    df = load_ltr_summary(args.summary_tsv)
    ltr_map = load_ltr_gff(args.ltr_gff)
    genome = load_genome(args.genome_fasta)

    u3_records = []
    ltr_records = []

    for _, row in df.iterrows():
        feat, tss, fallback_strand = row.Feature, int(row.TSS1), row.Strand
        info = ltr_map.get(feat)
        if not info:
            logging.warning(f"[LTR] No annotation for {feat}")
            continue
        chrom = info['chrom']
        strand = info['strand']
        if strand not in ('+', '-'):
            if fallback_strand in ('+', '-'):
                strand = fallback_strand
                logging.info(f"[LTR] Using fallback strand {strand} for {feat}")
            else:
                logging.warning(f"[LTR] No full-element strand for {feat} and no valid fallback")
                continue

        chosen = None
        side = None
        for key in ('lLTR', 'rLTR'):
            iv = info[key]
            if iv and iv[0] <= tss <= iv[1]:
                chosen = iv
                side = key
                break
        if not chosen:
            logging.warning(f"[LTR] TSS {tss} not in lLTR/rLTR of {feat}")
            continue

        boundary = chosen[0] if strand == '+' else chosen[1]
        u3_start, u3_end = sorted((tss, boundary))

        chrom_seq_record = genome.get(chrom)
        if not chrom_seq_record:
            logging.error(f"[LTR] Chrom {chrom} not in genome")
            return

        # U3
        u3_seq = extract_region(chrom_seq_record.seq, u3_start, u3_end, strand)
        u3_header = f"{feat}|{chrom}:{u3_start}-{u3_end}({strand})"
        u3_records.append(SeqRecord(u3_seq, id=u3_header, description=""))

        # full LTR (same side as TSS)
        full_ltr_start, full_ltr_end = info[side]
        full_seq = extract_region(chrom_seq_record.seq, full_ltr_start, full_ltr_end, strand)
        ltr_header = f"{feat}|{chrom}:{full_ltr_start}-{full_ltr_end}({strand})"
        ltr_records.append(SeqRecord(full_seq, id=ltr_header, description=""))

    SeqIO.write(u3_records, args.u3_output, 'fasta')
    logging.info(f"[LTR] Wrote {len(u3_records)} U3 sequences to {args.u3_output}")
    SeqIO.write(ltr_records, args.ltr_output, 'fasta')
    logging.info(f"[LTR] Wrote {len(ltr_records)} full LTR sequences to {args.ltr_output}")


# ---------- GENE MODE IMPLEMENTATION ----------

def load_gene_annotation(gtf_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse GTF/GFF and return mapping:
      gene_id -> (chrom, strand)

    Expects 'gene' features with attributes containing gene_id or ID.
    """
    gene_info: Dict[str, Tuple[str, str]] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, feature_type, strand = cols[0], cols[2], cols[6]
            if feature_type != "gene":
                continue
            attrs_field = cols[8].rstrip(";")
            attrs = dict(
                kv.split("=", 1) for kv in attrs_field.split(";") if "=" in kv
            )
            gene_id = attrs.get("gene_id") or attrs.get("ID")
            if gene_id:
                gene_info[gene_id] = (chrom, strand)
    logging.info(f"[Gene] Parsed {len(gene_info)} gene entries from {gtf_path}")
    return gene_info


def run_gene_mode(args):
    # Load TSS summary
    tss_df = pd.read_csv(args.summary_tsv, sep='\t')
    # Determine identifier column
    id_col = 'Gene' if 'Gene' in tss_df.columns else 'Feature' if 'Feature' in tss_df.columns else None
    if id_col is None:
        raise KeyError("[Gene] TSS summary must contain a 'Gene' or 'Feature' column")

    # Determine TSS coordinate column
    tss_col = 'TSS1' if 'TSS1' in tss_df.columns else 'TSS_abs' if 'TSS_abs' in tss_df.columns else None
    if tss_col is None:
        raise KeyError("[Gene] TSS summary must contain 'TSS1' or 'TSS_abs' column")

    gene_info = load_gene_annotation(args.gtf)
    genome = load_genome(args.genome_fasta)

    u3_records = []        # dummy U3s (upstream 1 kb)
    promoter_records = []  # actual promoter sequences (±1000 bp)
    n_written = 0

    with open(args.promoter_bed_output, 'w') as bed_out:
        for _, row in tss_df.iterrows():
            gene = str(row[id_col])
            if gene not in gene_info:
                logging.warning(f"[Gene] {gene} not found in gene annotation; skipping")
                continue
            tss = int(row[tss_col])  # 1-based
            chrom, strand = gene_info[gene]
            if strand not in ('+', '-'):
                logging.warning(f"[Gene] {gene} has invalid strand '{strand}'; skipping")
                continue

            chrom_seq_record = genome.get(chrom)
            if not chrom_seq_record:
                logging.warning(f"[Gene] Chrom {chrom} not in genome; skipping {gene}")
                continue
            chrom_len = len(chrom_seq_record.seq)

            # Promoter region: ±1000 bp around TSS (in genomic coords)
            promoter_start0 = max(0, tss - 1000)  # 0-based start
            promoter_end1 = min(tss + 1000, chrom_len)  # 1-based end, clamp to chrom len

            # BED line (0-based start, 1-based end, half-open [start0, end1))
            bed_out.write(f"{chrom}\t{promoter_start0}\t{promoter_end1}\t{gene}\t0\t{strand}\n")

            # Promoter FASTA sequence (strand-aware)
            promoter_start1 = promoter_start0 + 1       # convert to 1-based inclusive
            promoter_end1_1based = promoter_end1        # already 1-based end
            promoter_seq = extract_region(chrom_seq_record.seq,
                                          promoter_start1, promoter_end1_1based, strand)
            promoter_header = f"{gene}|{chrom}:{promoter_start0}-{promoter_end1_1based}({strand})"
            promoter_records.append(SeqRecord(promoter_seq, id=promoter_header, description=""))

            # Upstream 1 kb region for dummy U3 header
            if strand == '+':
                u3_start0 = max(0, tss - 1000)  # upstream [TSS-1000, TSS)
                u3_end1 = tss
            else:  # strand == '-'
                u3_start0 = tss                 # upstream [TSS, TSS+1000)
                u3_end1 = min(tss + 1000, chrom_len)

            u3_header = f"{gene}|{chrom}:{u3_start0}-{u3_end1}({strand})"
            dummy_seq = Seq("NNNNN")
            u3_records.append(SeqRecord(dummy_seq, id=u3_header, description=""))

            n_written += 1

    # Write outputs
    SeqIO.write(promoter_records, args.promoter_fasta_output, "fasta")
    logging.info(f"[Gene] Wrote {n_written} promoter entries to {args.promoter_bed_output} "
                 f"and {args.promoter_fasta_output}")
    SeqIO.write(u3_records, args.u3_output, "fasta")
    logging.info(f"[Gene] Wrote {n_written} dummy U3 sequences to {args.u3_output}")


# ---------- MAIN ----------

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    if args.mode == "LTR":
        run_ltr_mode(args)
    else:  # Gene
        run_gene_mode(args)


if __name__ == '__main__':
    main()
