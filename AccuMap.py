#!/usr/bin/env python3
"""
Wrapper for processes that directionalize, clean, and map long reads to genome of interest for purposes of 
LTR study. 
The pipeline can be started at any point depending on the input files. For instance, if aligning reads
from Illumina or PacBio, omit -run_pyc and start with Cutadapt using -run_cut and -run_map. 
Check required file types in the function definitions below "main" and adjsut inputs (--bam vs --fq) as necessary.
Running PyChopper assumes ONT reads are raw and were not cleaned prior to input no 
trimming or demultiplexxing). Running from Cutadapt assumes reads are demultiplexxed but not cleaned 
or aligned.

Minimum inputs: 
    --fq, Input FASTQ
    --pyc, PyChopper output FASTQ filename
    --sample, Sample name prefix
    --run_pyc, Run PyChopper
    --run_cut, Run Cutadapt
    --run_map, Run Minimap2
    --ref, Reference genome
    --kit, for ONT reads, default="PCB114"

Outputs (using provided sample name prefix):
    - cutadapt.fastq : FASTQ after Cutadapt trimming
    - pychopped.fastq : FASTQ after PyChopper primer removal
    - minimap2.sam : SAM file from Minimap2 alignment
    - strandtags.tsv : TSV file of read names and PyChopper strand orientations
    - STtagged.bam : BAM file with PyChopper strand tags (ST) added
    - STtagged.bed : BED file conversion of BAM with strand information
    - STtagged.sorted.bam : Sorted BAM file with PyChopper strand tags
    - pychopper.log : Log file from PyChopper run
    - pychopper.report.pdf : PDF report from PyChopper
    - pychopper.rescued.fastq : FASTQ of rescued reads from PyChopper
    - cutadapt.log : Log file from Cutadapt run
    - minimap2.log : Log file from Minimap2 run

Dependencies: subprocess, pysam, argparse

Usage (Default parameters):
    python3 AccuMap.py \
        --fq sample.fastq \
        --pyc sample.pychopped.fastq \
        --sample sample_name \
        --ref reference_genome.fasta \
        --run_pyc \
        --run_cut \
        --run_map  
"""
import argparse
import subprocess
import pysam
import os

## Run a shell command, log the output and print progress to terminal.
def run_command(command, log_file, step_name=None, sample_name=None):
    if step_name and sample_name:
        print(f"[INFO] Starting {step_name} for sample {sample_name}...")
    with open(log_file, "a") as log:
        log.write(f"Running: {command}\n")
        process = subprocess.run(command, shell=True, stdout=log, stderr=log, text=True)
        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {command}")
    if step_name and sample_name:
        print(f"[INFO] Finished {step_name} for sample {sample_name}.")

## Remove pychopper orientations
def extract_pychopper_tags(pychopped, out):
    if not pychopped:
        print("[INFO] No PyChopper FASTQ provided; skipping tag extraction.")
        return False
    command = f"grep \"^@\" {pychopped} | awk -F '[ =]' '{{print $1 \"\\t\" $3}}' | sed 's/^@//' > {out}"
    run_command(command, f"{out}.log", step_name="Extract PyChopper Tags", sample_name=os.path.basename(pychopped))
    return True

def load_pychopper_tags(path):
    tag_dict = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) == 2:
                qname, strand = fields
                tag_dict[qname] = strand.replace("strand=", "")
    return tag_dict

## Use pychopper strand directions to annotate alignment
def annotate_bam_with_strand(in_bam, out_bam, tag_dict, overwrite=False):
    print(f"[INFO] Annotating BAM with ST tags for sample {os.path.basename(in_bam)}...")
    n_seen = n_written = n_preserved = 0

    with pysam.AlignmentFile(in_bam, "rb") as infile, \
         pysam.AlignmentFile(out_bam, "wb", template=infile) as outfile:
        for read in infile:
            if tag_dict:
                st = tag_dict.get(read.query_name)
                if st in {"+", "-"}:
                    n_seen += 1
                    if overwrite or not read.has_tag("ST"):
                        read.set_tag("ST", st)
                        n_written += 1
                    else:
                        n_preserved += 1
            outfile.write(read)

    print(f"[INFO] Finished annotating BAM: {out_bam} "
          f"(tagged: {n_written}, preserved existing ST: {n_preserved}, seen in dict: {n_seen}).")


def sort_bam(in_bam, out_bam):
    command = f"samtools sort -o {out_bam} {in_bam}"
    run_command(command, f"{out_bam}.log", step_name="Sort BAM", sample_name=os.path.basename(in_bam))

def infer_strand(read):
    """
    Determine strand for a read with the following priority:
      1) ST tag  (PyChopper-derived tag)
      2) minimap2 'ts' tag (transcript strand)
      3) generic 'XS' tag (used by some aligners)
      4) minimap2 intron motif 'jM' tag (if spliced)
      5) alignment orientation (is_reverse)
      6) '.' as last resort
    Returns '+', '-', or '.'
    """
    # 1) ST (if already annotated)
    if read.has_tag("ST"):
        st = read.get_tag("ST")
        if st in {"+", "-"}:
            return st

    # 2) minimap2 transcript strand
    if read.has_tag("ts"):
        ts = read.get_tag("ts")
        if ts in {"+", "-"}:
            return ts

    # 3) XS (common in STAR/TopHat style)
    if read.has_tag("XS"):
        xs = read.get_tag("XS")
        if xs in {"+", "-"}:
            return xs

    # 4) minimap2 intron motif-based inference
    #    jM = list of motif codes per intron. We'll use the first informative one.
    #    Map: 1=GT-AG(+), 2=CT-AC(-), 3=GC-AG(+), 4=CT-GC(-)
    if read.has_tag("jM"):
        try:
            jm = read.get_tag("jM")
            if isinstance(jm, (list, tuple)) and len(jm) > 0:
                # take the first motif that matches our known set
                for code in jm:
                    # jM can be stored as ints or chars; coerce to int when possible
                    try:
                        c = int(code)
                    except Exception:
                        continue
                    if c == 1 or c == 3:
                        return "+"
                    if c == 2 or c == 4:
                        return "-"
        except Exception:
            pass

    # 5) fall back to alignment orientation (not transcript-aware, but better than a blanket '+')
    try:
        return "-" if read.is_reverse else "+"
    except Exception:
        pass

    # 6) unknown
    return "."


## Make custom bed file using bam alignments and pychopper orientation
def bam_to_bed_with_strand(bam, bed):
    print(f"[INFO] Converting BAM to BED for sample {os.path.basename(bam)}...")
    with pysam.AlignmentFile(bam, "rb") as infile, open(bed, "w") as out:
        for read in infile:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            name = read.query_name
            score = read.mapping_quality
            strand = infer_strand(read)
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    print(f"[INFO] Finished BED output: {bed}")

## Cutadapt with defualt settings unless specified by user
def run_cutadapt(fq_in, fq_out, log, a="A{15}", g="T{15}", t=16):
    info_file = f"{fq_out}.info"
    command = f"cutadapt -j {t} -a \"{a}\" -g \"{g}\" -o {fq_out} --info-file {info_file} {fq_in}"
    run_command(command, log, step_name="Cutadapt", sample_name=os.path.basename(fq_in))

## Minimap2 alignemnt with default settings unless specified by the user
def run_minimap2(ref, fq, sam_out, log, t=24, sec="no", gap=20000):
    command = (
        f"minimap2 -ax splice -uf -k14 --secondary={sec} -G {gap} -t {t} {ref} {fq}"
        f" | samtools view -b -o {bam_out} -"
    )
    run_command(command, log, step_name="Minimap2", sample_name=os.path.basename(fq))

## Pychopper primer removal with default (multiplex) settings unless specified by the user
def run_pychopper(fq_in, fq_out, unc, resc, report, log, kit="PCB114", t=8):
    command = f"pychopper -k {kit} -t {t} -r {report} -u {unc} -w {resc} {fq_in} {fq_out}"
    run_command(command, log, step_name="PyChopper", sample_name=os.path.basename(fq_in))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq", help="Input FASTQ")
    parser.add_argument("--bam", help="Input BAM file")
    parser.add_argument("--pyc", help="PyChopper output FASTQ")
    parser.add_argument("--sample", required=True, help="Sample name prefix")
    parser.add_argument("--run_pyc", action="store_true")
    parser.add_argument("--run_cut", action="store_true")
    parser.add_argument("--run_map", action="store_true")
    parser.add_argument("--ref", help="Reference genome")
    parser.add_argument("--kit", default="PCB114")
    parser.add_argument("--pyc_threads", type=int, default=8)
    parser.add_argument("--cut_threads", type=int, default=16)
    parser.add_argument("--map_threads", type=int, default=24)
    parser.add_argument("--map_sec", default="no")
    parser.add_argument("--map_gap", type=int, default=5000)
    args = parser.parse_args()

    current_fq = args.fq

    if args.run_pyc:
        pyc_out = f"{args.sample}.pychopped.fastq"
        unc = f"{args.sample}.unclassified.fastq"
        rpt = f"{args.sample}.pychopper.report.pdf"
        log = f"{args.sample}.pychopper.log"
        resc = f"{args.sample}pychopper.rescued.fastq"
        run_pychopper(current_fq, pyc_out, unc, resc, rpt, log, kit=args.kit, t=args.pyc_threads)
        args.pyc = pyc_out
        current_fq = pyc_out

    if args.run_cut:
        trim_out = f"{args.sample}.cutadapt.fastq"
        log = f"{args.sample}.cutadapt.log"
        run_cutadapt(current_fq, trim_out, log, t=args.cut_threads)
        current_fq = trim_out

    if args.run_map:
        bam_out = f"{args.sample}.minimap2.bam"
        log = f"{args.sample}.minimap2.log"
        run_minimap2(args.ref, current_fq, bam_out, log, t=args.map_threads, sec=args.map_sec, gap=args.map_gap)
        args.bam = bam_out

    tag_file = f"{args.sample}.strandtags.tsv"
    have_tags = extract_pychopper_tags(args.pyc, tag_file)
    tags = load_pychopper_tags(tag_file) if have_tags else {}

    tagged = f"{args.sample}.STtagged.bam"
    annotate_bam_with_strand(args.bam, tagged, tags)

    sorted_bam = f"{args.sample}.STtagged.sorted.bam"
    sort_bam(tagged, sorted_bam)

    bed = f"{args.sample}.STtagged.bed"
    bam_to_bed_with_strand(sorted_bam, bed)
