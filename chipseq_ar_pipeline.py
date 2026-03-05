#!/usr/bin/env python3
"""
ChIP-seq AR Motif Discovery Pipeline
=====================================
A complete bioinformatics pipeline for identifying Androgen Receptor (AR)
binding motifs from public ChIP-seq data (NCBI GEO).

Usage:
    python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984
    python chipseq_ar_pipeline.py -c SRR1615985,SRR1615986 -i SRR1615984 -p PE -t 16
    python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 --resume

Author: MOUNEM HOUF
License: MIT
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import textwrap
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional, Dict, Tuple

# ============================================================================
# CONSTANTS
# ============================================================================

VERSION = "1.1.0"
HOMER_URL = "http://homer.ucsd.edu/homer/configureHomer.pl"

# Genome configurations: URL, MACS2 genome size, HOMER name, R TxDb package
GENOME_CONFIG = {
    "hg38": {
        "name": "GRCh38 (hg38)",
        "description": "Human genome assembly GRCh38 — UCSC hg38 (2013, most widely used)",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "macs2_genome": "hs",
        "homer_genome": "hg38",
        "homer_install": "hg38",
        "txdb_package": "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "txdb_r": "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "orgdb": "org.Hs.eg.db",
        "species": "hsa",
    },
    "hg19": {
        "name": "GRCh37 (hg19)",
        "description": "Human genome assembly GRCh37 — UCSC hg19 (2009, legacy, many older datasets)",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "macs2_genome": "hs",
        "homer_genome": "hg19",
        "homer_install": "hg19",
        "txdb_package": "TxDb.Hsapiens.UCSC.hg19.knownGene",
        "txdb_r": "TxDb.Hsapiens.UCSC.hg19.knownGene",
        "orgdb": "org.Hs.eg.db",
        "species": "hsa",
    },
    "t2t": {
        "name": "T2T-CHM13v2.0",
        "description": "Telomere-to-Telomere complete human genome (2022, gapless, most complete)",
        "url": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
        "macs2_genome": "2.9e9",
        "homer_genome": "chm13v2.0",
        "homer_install": None,  # HOMER does not have T2T pre-built — use custom FASTA
        "txdb_package": None,  # No pre-built TxDb for T2T
        "txdb_r": None,
        "orgdb": "org.Hs.eg.db",
        "species": "hsa",
    },
}

STEP_NAMES = {
    0: "Download data",
    1: "Quality control (FastQC)",
    2: "Adapter trimming (Trim Galore)",
    3: "Alignment (Bowtie2)",
    4: "Post-alignment processing",
    5: "ChIP-seq QC (deepTools)",
    6: "Peak calling (MACS2)",
    7: "Peak annotation (HOMER + ChIPseeker)",
    8: "Motif discovery (HOMER + MEME-ChIP)",
    9: "Motif analysis & visualization (R)",
}

# ============================================================================
# LOGGING SETUP
# ============================================================================

class ColorFormatter(logging.Formatter):
    """Custom formatter with colors for terminal output."""
    COLORS = {
        "DEBUG": "\033[0;37m",
        "INFO": "\033[0;32m",
        "WARNING": "\033[0;33m",
        "ERROR": "\033[0;31m",
        "CRITICAL": "\033[1;31m",
    }
    RESET = "\033[0m"
    ICONS = {
        "DEBUG": "    ",
        "INFO": " ✓  ",
        "WARNING": " ⚠  ",
        "ERROR": " ✗  ",
        "CRITICAL": " ✗✗ ",
    }

    def format(self, record):
        color = self.COLORS.get(record.levelname, self.RESET)
        icon = self.ICONS.get(record.levelname, "")
        record.msg = f"{color}[{icon}]{self.RESET} {record.msg}"
        return super().format(record)


def setup_logging(log_dir: Path) -> logging.Logger:
    """Configure logging to both file and console."""
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"pipeline_{timestamp}.log"

    logger = logging.getLogger("chipseq")
    logger.setLevel(logging.DEBUG)

    # File handler (detailed)
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logger.addHandler(fh)

    # Console handler (colored)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(ColorFormatter("%(message)s"))
    logger.addHandler(ch)

    logger.info(f"Log file: {log_file}")
    return logger


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def run_cmd(cmd: str, log: logging.Logger, description: str = "",
            log_file: Optional[Path] = None, check: bool = True,
            env: Optional[dict] = None) -> subprocess.CompletedProcess:
    """Execute a shell command with logging and error handling."""
    if description:
        log.info(description)
    log.debug(f"Command: {cmd}")

    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    stderr_dest = subprocess.PIPE
    stdout_dest = subprocess.PIPE

    try:
        result = subprocess.run(
            cmd, shell=True, stdout=stdout_dest, stderr=stderr_dest,
            text=True, env=merged_env
        )

        # Write to log file if specified
        if log_file:
            with open(log_file, "a") as f:
                if result.stdout:
                    f.write(result.stdout)
                if result.stderr:
                    f.write(result.stderr)

        if check and result.returncode != 0:
            log.error(f"Command failed (exit code {result.returncode}): {cmd}")
            if result.stderr:
                log.error(f"STDERR: {result.stderr[:500]}")
            raise subprocess.CalledProcessError(result.returncode, cmd)

        return result

    except FileNotFoundError:
        log.error(f"Command not found: {cmd.split()[0]}")
        raise


def check_tool(tool: str, log: logging.Logger) -> bool:
    """Check if a command-line tool is available."""
    result = shutil.which(tool)
    if result:
        log.info(f"{tool}")
        return True
    else:
        log.warning(f"{tool} — NOT FOUND")
        return False


def file_exists_and_nonempty(path: Path) -> bool:
    """Check if a file exists and is not empty."""
    return path.exists() and path.stat().st_size > 0


def format_duration(seconds: float) -> str:
    """Format seconds into human-readable duration."""
    td = timedelta(seconds=int(seconds))
    hours, remainder = divmod(td.seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if td.days:
        return f"{td.days}d {hours}h {minutes}m {secs}s"
    elif hours:
        return f"{hours}h {minutes}m {secs}s"
    elif minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def print_banner(text: str):
    """Print a formatted step banner."""
    width = 70
    print(f"\n{'═' * width}")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  │  {text}")
    print(f"{'═' * width}\n")


# ============================================================================
# STATE MANAGEMENT (for resume capability)
# ============================================================================

class PipelineState:
    """Track pipeline progress for resume capability."""

    def __init__(self, state_file: Path):
        self.state_file = state_file
        self.data = self._load()

    def _load(self) -> dict:
        if self.state_file.exists():
            with open(self.state_file) as f:
                return json.load(f)
        return {"completed_steps": [], "current_step": 0, "params": {}}

    def save(self):
        with open(self.state_file, "w") as f:
            json.dump(self.data, f, indent=2)

    def mark_complete(self, step: int):
        if step not in self.data["completed_steps"]:
            self.data["completed_steps"].append(step)
        self.data["current_step"] = step + 1
        self.save()

    def is_complete(self, step: int) -> bool:
        return step in self.data["completed_steps"]

    def get_resume_step(self) -> int:
        return self.data.get("current_step", 0)

    def save_params(self, params: dict):
        self.data["params"] = params
        self.save()

    def reset(self):
        self.data = {"completed_steps": [], "current_step": 0, "params": {}}
        self.save()


# ============================================================================
# PIPELINE CLASS
# ============================================================================

class ChIPseqPipeline:
    """Main pipeline class for AR ChIP-seq motif discovery."""

    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.chip_srrs = args.chip.split(",")
        self.input_srr = args.input
        self.num_chip = len(self.chip_srrs)
        self.threads = args.threads
        self.read_type = args.read_type
        self.genome_build = args.genome
        self.macs2_qvalue = args.qvalue
        self.fe_cutoff = args.fold_enrichment
        self.summit_window = args.window

        # Genome configuration
        if self.genome_build not in GENOME_CONFIG:
            print(f"ERROR: Unknown genome '{self.genome_build}'.")
            print(f"Available genomes: {', '.join(GENOME_CONFIG.keys())}")
            for key, cfg in GENOME_CONFIG.items():
                print(f"  {key:6s} — {cfg['description']}")
            sys.exit(1)

        self.genome_cfg = GENOME_CONFIG[self.genome_build]
        self.macs2_genome_size = self.genome_cfg["macs2_genome"]
        self.homer_genome = self.genome_cfg["homer_genome"]
        self.genome_url = self.genome_cfg["url"]

        # Genome FASTA filename
        if self.genome_build == "t2t":
            self.genome_fa_name = "chm13v2.0.fa"
        else:
            self.genome_fa_name = f"{self.genome_build}.fa"

        # Directories
        self.project_dir = Path(args.output).resolve()
        self.raw_dir = self.project_dir / "data" / "raw"
        self.genome_dir = self.project_dir / "genome"
        self.fastqc_dir = self.project_dir / "results" / "fastqc"
        self.trim_dir = self.project_dir / "results" / "trimmed"
        self.align_dir = self.project_dir / "results" / "alignment"
        self.chipqc_dir = self.project_dir / "results" / "chipqc"
        self.peaks_dir = self.project_dir / "results" / "peaks"
        self.annot_dir = self.project_dir / "results" / "annotation"
        self.motifs_dir = self.project_dir / "results" / "motifs"
        self.logs_dir = self.project_dir / "logs"

        # Key file paths
        self.genome_fa = self.genome_dir / self.genome_fa_name
        self.genome_index = self.genome_dir / f"{self.genome_build}_bt2"

        # Create directories
        for d in [self.raw_dir, self.genome_dir, self.fastqc_dir, self.trim_dir,
                  self.align_dir, self.chipqc_dir, self.peaks_dir, self.annot_dir,
                  self.motifs_dir / "homer", self.motifs_dir / "meme",
                  self.motifs_dir / "tomtom", self.motifs_dir / "fimo",
                  self.logs_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Logging
        self.log = setup_logging(self.logs_dir)

        # State management
        self.state = PipelineState(self.project_dir / ".pipeline_state.json")

        # Sample names
        self.chip_names = []
        for idx in range(self.num_chip):
            if self.num_chip == 1:
                self.chip_names.append("ChIP")
            else:
                self.chip_names.append(f"ChIP_rep{idx + 1}")
        self.all_names = self.chip_names + ["Input"]
        self.all_srrs = self.chip_srrs + [self.input_srr]

        # HOMER path
        homer_bin = Path.home() / "homer" / "bin"
        if homer_bin.exists():
            os.environ["PATH"] = str(homer_bin) + ":" + os.environ.get("PATH", "")

        # Step timing
        self.step_times: Dict[int, float] = {}

    # ========================================================================
    # STEP 0: DOWNLOAD DATA
    # ========================================================================

    def step0_download(self):
        """Download SRA data and reference genome."""
        print_banner("Step 0: Download SRA data + reference genome")

        # Download SRA files
        for srr, name in zip(self.all_srrs, self.all_names):
            if self.read_type == "SE":
                target = self.raw_dir / f"{name}.fastq.gz"
            else:
                target = self.raw_dir / f"{name}_R1.fastq.gz"

            if file_exists_and_nonempty(target):
                self.log.info(f"SKIP {name} ({srr}) — already exists")
                continue

            run_cmd(
                f"prefetch {srr} -O {self.raw_dir}/",
                self.log, f"Downloading {srr} → {name}...",
                log_file=self.logs_dir / f"fasterq_{name}.log"
            )

            if self.read_type == "SE":
                run_cmd(
                    f"fasterq-dump {srr} --outdir {self.raw_dir} --threads {self.threads}",
                    self.log, f"Converting {srr} to FASTQ...",
                    log_file=self.logs_dir / f"fasterq_{name}.log"
                )
                src = self.raw_dir / f"{srr}.fastq"
                dst = self.raw_dir / f"{name}.fastq"
                if src.exists():
                    src.rename(dst)
                run_cmd(f"gzip {dst}", self.log)
            else:
                run_cmd(
                    f"fasterq-dump {srr} --split-files --outdir {self.raw_dir} --threads {self.threads}",
                    self.log, f"Converting {srr} to FASTQ (paired-end)...",
                    log_file=self.logs_dir / f"fasterq_{name}.log"
                )
                (self.raw_dir / f"{srr}_1.fastq").rename(self.raw_dir / f"{name}_R1.fastq")
                (self.raw_dir / f"{srr}_2.fastq").rename(self.raw_dir / f"{name}_R2.fastq")
                run_cmd(f"gzip {self.raw_dir / f'{name}_R1.fastq'} {self.raw_dir / f'{name}_R2.fastq'}", self.log)

            self.log.info(f"{name} downloaded")

        # Download genome
        if not file_exists_and_nonempty(self.genome_fa):
            self.log.info(f"Downloading {self.genome_cfg['name']} genome...")
            self.log.info(f"  {self.genome_cfg['description']}")
            run_cmd(
                f"wget -q -O {self.genome_fa}.gz {self.genome_url}",
                self.log, f"Downloading {self.genome_build} genome..."
            )
            run_cmd(f"gunzip {self.genome_fa}.gz", self.log)
            self.log.info("Genome downloaded")
        else:
            self.log.info(f"SKIP genome ({self.genome_cfg['name']}) — already exists")

        # Build Bowtie2 index
        bt2_file = Path(f"{self.genome_index}.1.bt2")
        if not bt2_file.exists():
            run_cmd(
                f"bowtie2-build --threads {self.threads} {self.genome_fa} {self.genome_index}",
                self.log, "Building Bowtie2 index (~1 hour)..."
            )
            self.log.info("Bowtie2 index built")
        else:
            self.log.info("SKIP Bowtie2 index — already exists")

        # Chromosome sizes
        chrom_sizes = self.genome_dir / f"{self.genome_build}.chrom.sizes"
        if not chrom_sizes.exists():
            run_cmd(f"samtools faidx {self.genome_fa}", self.log)
            run_cmd(
                f"awk '{{print $1\"\\t\"$2}}' {self.genome_fa}.fai > {chrom_sizes}",
                self.log
            )

    # ========================================================================
    # STEP 1: FASTQC
    # ========================================================================

    def step1_fastqc(self):
        """Run FastQC quality control on raw reads."""
        print_banner("Step 1: Quality control (FastQC)")

        run_cmd(
            f"fastqc -o {self.fastqc_dir} -t {self.threads} {self.raw_dir}/*.fastq.gz",
            self.log, "Running FastQC...",
            log_file=self.logs_dir / "fastqc.log"
        )

        multiqc_out = self.fastqc_dir / "multiqc"
        run_cmd(
            f"multiqc {self.fastqc_dir} -o {multiqc_out} --force",
            self.log, "Generating MultiQC report...",
            log_file=self.logs_dir / "multiqc_raw.log"
        )

        self.log.info(f"Reports: {multiqc_out}/multiqc_report.html")

    # ========================================================================
    # STEP 2: TRIMMING
    # ========================================================================

    def step2_trimming(self):
        """Trim adapters and low-quality bases."""
        print_banner("Step 2: Adapter trimming (Trim Galore)")

        for name in self.all_names:
            self.log.info(f"Trimming {name}...")

            if self.read_type == "SE":
                cmd = (
                    f"trim_galore --quality 20 --length 20 --fastqc --cores 4 "
                    f"--output_dir {self.trim_dir} "
                    f"{self.raw_dir}/{name}.fastq.gz"
                )
            else:
                cmd = (
                    f"trim_galore --quality 20 --length 20 --paired --fastqc --cores 4 "
                    f"--output_dir {self.trim_dir} "
                    f"{self.raw_dir}/{name}_R1.fastq.gz {self.raw_dir}/{name}_R2.fastq.gz"
                )

            run_cmd(cmd, self.log, log_file=self.logs_dir / f"trim_{name}.log")
            self.log.info(f"{name} trimmed")

        run_cmd(
            f"multiqc {self.trim_dir} -o {self.trim_dir}/multiqc --force",
            self.log, log_file=self.logs_dir / "multiqc_trim.log"
        )

    # ========================================================================
    # STEP 3: ALIGNMENT
    # ========================================================================

    def step3_alignment(self):
        """Align reads to reference genome with Bowtie2."""
        print_banner(f"Step 3: Alignment (Bowtie2 → {self.genome_cfg['name']})")

        for name in self.all_names:
            self.log.info(f"Aligning {name}...")
            bam_out = self.align_dir / f"{name}_raw.bam"
            bt2_log = self.logs_dir / f"bowtie2_{name}.log"

            if self.read_type == "SE":
                trimmed = self.trim_dir / f"{name}_trimmed.fq.gz"
                cmd = (
                    f"bowtie2 -x {self.genome_index} -U {trimmed} "
                    f"--very-sensitive --threads {self.threads} "
                    f"2> {bt2_log} "
                    f"| samtools view -bS -@ {self.threads} - > {bam_out}"
                )
            else:
                r1 = self.trim_dir / f"{name}_R1_val_1.fq.gz"
                r2 = self.trim_dir / f"{name}_R2_val_2.fq.gz"
                cmd = (
                    f"bowtie2 -x {self.genome_index} -1 {r1} -2 {r2} "
                    f"--very-sensitive --no-mixed --no-discordant --maxins 500 "
                    f"--threads {self.threads} "
                    f"2> {bt2_log} "
                    f"| samtools view -bS -@ {self.threads} - > {bam_out}"
                )

            run_cmd(cmd, self.log)

            # Parse alignment rate
            try:
                with open(bt2_log) as f:
                    for line in f:
                        if "overall alignment rate" in line:
                            rate = line.strip().split()[0]
                            self.log.info(f"{name} — {rate} alignment rate")
                            break
            except FileNotFoundError:
                pass

    # ========================================================================
    # STEP 4: POST-ALIGNMENT
    # ========================================================================

    def step4_post_alignment(self):
        """Sort, filter, deduplicate, and index BAM files."""
        print_banner("Step 4: Post-alignment processing")

        for name in self.all_names:
            self.log.info(f"Processing {name}...")

            raw_bam = self.align_dir / f"{name}_raw.bam"
            sorted_bam = self.align_dir / f"{name}_sorted.bam"
            filt_bam = self.align_dir / f"{name}_filt.bam"
            dedup_bam = self.align_dir / f"{name}_dedup.bam"
            final_bam = self.align_dir / f"{name}_final.bam"
            dedup_metrics = self.align_dir / f"{name}_dedup_metrics.txt"

            # Sort
            run_cmd(
                f"samtools sort -@ {self.threads} -o {sorted_bam} {raw_bam}",
                self.log
            )

            # Filter by MAPQ
            run_cmd(
                f"samtools view -b -q 10 -F 1804 -@ {self.threads} {sorted_bam} > {filt_bam}",
                self.log
            )

            # Remove PCR duplicates
            run_cmd(
                f"picard MarkDuplicates INPUT={filt_bam} OUTPUT={dedup_bam} "
                f"METRICS_FILE={dedup_metrics} REMOVE_DUPLICATES=true "
                f"VALIDATION_STRINGENCY=LENIENT",
                self.log,
                log_file=self.logs_dir / f"picard_{name}.log"
            )

            # Remove non-standard chromosomes
            run_cmd(
                f"samtools view -h {dedup_bam} | "
                f"grep -v -E 'chrUn|chrM|_random|_alt|_hap|_fix' | "
                f"samtools view -b -@ {self.threads} - > {final_bam}",
                self.log
            )

            # Index
            run_cmd(f"samtools index {final_bam}", self.log)

            # Stats
            run_cmd(
                f"samtools flagstat {final_bam} > {self.align_dir / f'{name}_flagstat.txt'}",
                self.log
            )

            # Count final reads
            result = run_cmd(f"samtools view -c {final_bam}", self.log)
            reads = result.stdout.strip()
            self.log.info(f"{name} — {reads} reads after filtering")

            # Cleanup intermediates
            for f in [raw_bam, sorted_bam, filt_bam, dedup_bam]:
                if f.exists():
                    f.unlink()

    # ========================================================================
    # STEP 5: ChIP-seq QC
    # ========================================================================

    def step5_chipqc(self):
        """ChIP-seq specific quality control with deepTools."""
        print_banner("Step 5: ChIP-seq quality control")

        chip_bam = self.align_dir / f"{self.chip_names[0]}_final.bam"
        input_bam = self.align_dir / "Input_final.bam"

        # Fingerprint plot
        run_cmd(
            f"plotFingerprint -b {chip_bam} {input_bam} "
            f"--labels {self.chip_names[0]} Input "
            f"--plotFile {self.chipqc_dir}/fingerprint.png "
            f"--outQualityMetrics {self.chipqc_dir}/fingerprint_metrics.txt "
            f"--numberOfProcessors {self.threads}",
            self.log, "Generating fingerprint plot...",
            log_file=self.logs_dir / "deeptools.log"
        )

        # BigWig tracks
        for bam, label in [(chip_bam, self.chip_names[0]), (input_bam, "Input")]:
            run_cmd(
                f"bamCoverage -b {bam} -o {self.chipqc_dir}/{label}_RPKM.bw "
                f"--normalizeUsing RPKM --binSize 10 "
                f"--numberOfProcessors {self.threads}",
                self.log, f"Generating bigWig for {label}...",
                log_file=self.logs_dir / "deeptools.log"
            )

        # Replicate correlation
        if self.num_chip > 1:
            bams = " ".join(str(self.align_dir / f"{n}_final.bam") for n in self.chip_names)
            run_cmd(
                f"multiBamSummary bins --bamfiles {bams} "
                f"--outFileName {self.chipqc_dir}/multibam_summary.npz "
                f"--numberOfProcessors {self.threads}",
                self.log, "Computing replicate correlation...",
                log_file=self.logs_dir / "deeptools.log"
            )
            run_cmd(
                f"plotCorrelation -in {self.chipqc_dir}/multibam_summary.npz "
                f"--corMethod pearson --whatToPlot heatmap "
                f"--plotFile {self.chipqc_dir}/replicate_correlation.png",
                self.log,
                log_file=self.logs_dir / "deeptools.log"
            )

        self.log.info("ChIP-seq QC complete")

    # ========================================================================
    # STEP 6: PEAK CALLING
    # ========================================================================

    def step6_peak_calling(self):
        """Call peaks with MACS2."""
        print_banner("Step 6: Peak calling (MACS2)")

        chip_bams = " ".join(str(self.align_dir / f"{n}_final.bam") for n in self.chip_names)
        input_bam = self.align_dir / "Input_final.bam"

        run_cmd(
            f"macs2 callpeak -t {chip_bams} -c {input_bam} "
            f"-f BAM -g {self.macs2_genome_size} --outdir {self.peaks_dir} -n AR_ChIP "
            f"-q {self.macs2_qvalue} --keep-dup 1 --call-summits --bdg",
            self.log, "Running MACS2...",
            log_file=self.logs_dir / "macs2.log"
        )

        # Count peaks
        peaks_file = self.peaks_dir / "AR_ChIP_peaks.narrowPeak"
        total = sum(1 for _ in open(peaks_file))
        self.log.info(f"Total peaks: {total}")

        # Filter by fold enrichment
        filtered = self.peaks_dir / "AR_ChIP_peaks_filtered.narrowPeak"
        with open(peaks_file) as fin, open(filtered, "w") as fout:
            for line in fin:
                cols = line.strip().split("\t")
                if float(cols[6]) > self.fe_cutoff:
                    fout.write(line)

        filt_count = sum(1 for _ in open(filtered))
        self.log.info(f"Filtered peaks (FE > {self.fe_cutoff}): {filt_count}")

        # Summit windows
        summits = self.peaks_dir / "AR_ChIP_summits.bed"
        summit_window = self.peaks_dir / "AR_summits_window.bed"
        w = self.summit_window

        with open(summits) as fin, open(summit_window, "w") as fout:
            for line in fin:
                cols = line.strip().split("\t")
                chrom = cols[0]
                pos = int(cols[1])
                start = max(0, pos - w)
                end = pos + w
                fout.write(f"{chrom}\t{start}\t{end}\t{cols[3]}\t{cols[4]}\n")

        # Merge overlapping regions
        run_cmd(
            f"sort -k1,1 -k2,2n {summit_window} | bedtools merge -i - > {self.peaks_dir / 'AR_summits_merged.bed'}",
            self.log
        )

        self.log.info("Peak calling complete")

    # ========================================================================
    # STEP 7: PEAK ANNOTATION
    # ========================================================================

    def step7_annotation(self):
        """Annotate peaks with HOMER and ChIPseeker."""
        print_banner("Step 7: Peak annotation")

        peaks_filtered = self.peaks_dir / "AR_ChIP_peaks_filtered.narrowPeak"

        # HOMER annotation
        run_cmd(
            f"annotatePeaks.pl {peaks_filtered} {self.homer_genome} "
            f"-annStats {self.annot_dir}/annotation_stats.txt "
            f"> {self.annot_dir}/AR_peaks_annotated.txt",
            self.log, "Running HOMER annotatePeaks...",
            log_file=self.logs_dir / "homer_annotate.log"
        )
        self.log.info("HOMER annotation done")

        # ChIPseeker R script — adapt TxDb to genome build
        r_script = self.logs_dir / "_annotation.R"
        txdb_r = self.genome_cfg.get("txdb_r")
        txdb_pkg = self.genome_cfg.get("txdb_package")
        orgdb = self.genome_cfg.get("orgdb", "org.Hs.eg.db")
        species = self.genome_cfg.get("species", "hsa")

        if txdb_r is None:
            self.log.warning(
                f"No pre-built TxDb available for {self.genome_build} (T2T). "
                "ChIPseeker annotation will be skipped. "
                "HOMER annotation is still available above."
            )
            return

        r_code = textwrap.dedent(f"""\
        args <- commandArgs(trailingOnly = TRUE)
        peaks_file <- args[1]; output_dir <- args[2]
        suppressPackageStartupMessages({{
          library(ChIPseeker); library({txdb_pkg})
          library({orgdb}); library(clusterProfiler); library(ggplot2)
        }})
        cat("Loading peaks...\\n")
        peaks <- readPeakFile(peaks_file)
        txdb <- {txdb_r}
        peakAnno <- annotatePeak(peaks, TxDb=txdb, annoDb="{orgdb}", level="gene")
        anno_df <- as.data.frame(peakAnno)
        write.csv(anno_df, file.path(output_dir,"AR_peaks_chipseeker.csv"), row.names=FALSE)
        pdf(file.path(output_dir,"annotation_pie.pdf"), w=8, h=8); plotAnnoPie(peakAnno); dev.off()
        pdf(file.path(output_dir,"distance_to_TSS.pdf"), w=10, h=6)
        plotDistToTSS(peakAnno, title="AR peaks: Distance to TSS"); dev.off()
        promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
        tagMatrix <- getTagMatrix(peaks, windows=promoter)
        pdf(file.path(output_dir,"TSS_profile.pdf"), w=10, h=6)
        plotAvgProf(tagMatrix, xlim=c(-3000,3000), xlab="Distance to TSS (bp)",
                    ylab="Read count frequency"); dev.off()
        genes <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
        ego <- enrichGO(gene=genes, OrgDb={orgdb}, ont="BP",
                        pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
        if(!is.null(ego) && nrow(as.data.frame(ego))>0) {{
          write.csv(as.data.frame(ego), file.path(output_dir,"GO_BP_enrichment.csv"), row.names=FALSE)
          pdf(file.path(output_dir,"GO_BP_dotplot.pdf"), w=10, h=8)
          print(dotplot(ego, showCategory=20, title="GO — AR targets")); dev.off()
        }}
        ekegg <- enrichKEGG(gene=genes, organism="{species}", pAdjustMethod="BH", qvalueCutoff=0.05)
        if(!is.null(ekegg) && nrow(as.data.frame(ekegg))>0) {{
          write.csv(as.data.frame(ekegg), file.path(output_dir,"KEGG_enrichment.csv"), row.names=FALSE)
          pdf(file.path(output_dir,"KEGG_dotplot.pdf"), w=10, h=8)
          print(dotplot(ekegg, showCategory=15, title="KEGG — AR targets")); dev.off()
        }}
        cat("Annotation complete.\\n")
        """)

        with open(r_script, "w") as f:
            f.write(r_code)

        # Resolve Rscript path explicitly
        rscript_path = shutil.which("Rscript") or "Rscript"
        run_cmd(
            f"{rscript_path} {r_script} {peaks_filtered} {self.annot_dir}",
            self.log, "Running ChIPseeker (R)...",
            log_file=self.logs_dir / "chipseq.log",
            check=False
        )
        self.log.info("Annotation complete")

    # ========================================================================
    # STEP 8: MOTIF DISCOVERY
    # ========================================================================

    def step8_motif_discovery(self):
        """Discover motifs with HOMER and MEME-ChIP."""
        print_banner("Step 8: Motif discovery (HOMER + MEME-ChIP)")

        summit_bed = self.peaks_dir / "AR_summits_window.bed"
        summit_fa = self.motifs_dir / "AR_summit_sequences.fa"

        # Extract FASTA sequences
        run_cmd(
            f"bedtools getfasta -fi {self.genome_fa} -bed {summit_bed} -fo {summit_fa}",
            self.log, "Extracting summit sequences..."
        )

        result = run_cmd(f"grep -c '^>' {summit_fa}", self.log)
        num_seq = int(result.stdout.strip())
        self.log.info(f"{num_seq} sequences extracted")

        # Subsample if too many
        if num_seq > 10000:
            self.log.info("Subsampling to top 10,000 peaks...")
            top_bed = self.motifs_dir / "top_summits.bed"
            run_cmd(f"sort -k5,5rn {summit_bed} | head -10000 > {top_bed}", self.log)
            summit_fa_top = self.motifs_dir / "AR_summit_top10k.fa"
            run_cmd(
                f"bedtools getfasta -fi {self.genome_fa} -bed {top_bed} -fo {summit_fa_top}",
                self.log
            )
            summit_fa = summit_fa_top
            summit_bed = top_bed

        # HOMER
        self.log.info("Running HOMER findMotifsGenome.pl...")
        homer_out = self.motifs_dir / "homer"

        # For T2T genome, HOMER needs the FASTA path instead of genome name
        if self.genome_cfg["homer_install"] is None:
            homer_genome_arg = str(self.genome_fa)
            self.log.info(f"Using custom FASTA for HOMER (T2T): {self.genome_fa}")
        else:
            homer_genome_arg = self.homer_genome

        run_cmd(
            f"findMotifsGenome.pl {summit_bed} {homer_genome_arg} {homer_out}/ "
            f"-size 200 -mask -p {self.threads}",
            self.log,
            log_file=self.logs_dir / "homer_motifs.log",
            check=False
        )
        homer_html = homer_out / "homerResults.html"
        if homer_html.exists():
            self.log.info(f"HOMER results: {homer_html}")
        else:
            self.log.warning("HOMER did not produce results — check logs")

        # MEME-ChIP (uses meme_env)
        self.log.info("Running MEME-ChIP...")
        meme_out = self.motifs_dir / "meme"

        # Try to find motif database
        meme_db_arg = ""
        db_search_paths = [
            self.project_dir / "motif_databases" / "JASPAR",
        ]
        for db_dir in db_search_paths:
            if db_dir.exists():
                for f in db_dir.glob("*.meme"):
                    meme_db_arg = f"-db {f}"
                    break
                for f in db_dir.glob("*.txt"):
                    meme_db_arg = f"-db {f}"
                    break

        meme_cmd = (
            f"meme-chip -oc {meme_out} {meme_db_arg} "
            f"-maxw 20 -minw 6 "
            f"-meme-nmotifs 10 -meme-mod zoops "
            f"-centrimo-local "
            f"{summit_fa}"
        )

        # Try running in meme_env
        # We use a bash -c wrapper to properly activate the conda env
        conda_base = os.environ.get("CONDA_EXE", "").replace("/bin/conda", "")
        if not conda_base:
            conda_base = str(Path.home() / "miniforge3")
            if not Path(conda_base).exists():
                conda_base = str(Path.home() / "anaconda3")
            if not Path(conda_base).exists():
                conda_base = str(Path("/opt/anaconda3"))

        meme_activate = (
            f"source {conda_base}/etc/profile.d/conda.sh && "
            f"conda activate meme_env && "
        )

        run_cmd(
            f"bash -c '{meme_activate}{meme_cmd}'",
            self.log,
            log_file=self.logs_dir / "memechip.log",
            check=False
        )

        meme_html = meme_out / "meme-chip.html"
        if meme_html.exists():
            self.log.info(f"MEME-ChIP results: {meme_html}")
        else:
            self.log.warning(
                "MEME-ChIP did not produce results. Run manually:\n"
                f"  conda activate meme_env\n"
                f"  {meme_cmd}"
            )

        # TOMTOM
        meme_txt = meme_out / "meme_out" / "meme.txt"
        if meme_txt.exists() and meme_db_arg:
            db_path = meme_db_arg.replace("-db ", "")
            run_cmd(
                f"bash -c '{meme_activate}tomtom -oc {self.motifs_dir / 'tomtom'} {meme_txt} {db_path}'",
                self.log, "Running TOMTOM...",
                log_file=self.logs_dir / "tomtom.log",
                check=False
            )

        # FIMO
        if meme_txt.exists():
            run_cmd(
                f"bash -c '{meme_activate}fimo --oc {self.motifs_dir / 'fimo'} --thresh 1e-4 {meme_txt} {self.genome_fa}'",
                self.log, "Running FIMO...",
                log_file=self.logs_dir / "fimo.log",
                check=False
            )

    # ========================================================================
    # STEP 9: MOTIF ANALYSIS
    # ========================================================================

    def step9_motif_analysis(self):
        """Analyze and visualize discovered motifs with R."""
        print_banner("Step 9: Motif analysis & visualization")

        r_script = self.logs_dir / "_motif_analysis.R"
        r_code = textwrap.dedent(f"""\
        motifs_dir <- "{self.motifs_dir}"
        suppressPackageStartupMessages({{
          library(universalmotif); library(ggplot2); library(ggseqlogo)
        }})
        cat("=== AR Motif Analysis ===\\n\\n")
        meme_file <- file.path(motifs_dir,"meme","meme_out","meme.txt")
        if(!file.exists(meme_file)){{cat("MEME output not found\\n"); quit(status=0)}}
        motifs <- read_meme(meme_file)
        if(!is.list(motifs)) motifs <- list(motifs)
        cat(sprintf("%d motifs discovered\\n\\n", length(motifs)))
        summary_df <- data.frame(
          Motif=paste0("Motif_",seq_along(motifs)),
          Consensus=sapply(motifs, function(m) m@consensus),
          Width=sapply(motifs, function(m) ncol(m@motif)),
          stringsAsFactors=FALSE)
        print(summary_df)
        write.csv(summary_df, file.path(motifs_dir,"motif_summary.csv"), row.names=FALSE)
        for(i in seq_along(motifs)){{
          m <- motifs[[i]]
          pdf(file.path(motifs_dir,sprintf("motif_%d_logo.pdf",i)),w=max(6,ncol(m@motif)*0.5),h=3)
          print(ggseqlogo(m@motif,method="bits")+ggtitle(sprintf("Motif %d: %s",i,m@consensus)))
          dev.off()
        }}
        if(length(motifs)>1){{
          pwm_list <- lapply(motifs,function(m) m@motif)
          names(pwm_list) <- paste0("Motif ",seq_along(motifs),": ",sapply(motifs,function(m) m@consensus))
          pdf(file.path(motifs_dir,"all_motifs.pdf"),w=12,h=3*min(length(motifs),10))
          print(ggseqlogo(pwm_list,method="bits",ncol=1)); dev.off()
        }}
        are_pwm <- matrix(c(
          .85,.05,.05,.05,.05,.05,.85,.05,.85,.05,.05,.05,
          .85,.05,.05,.05,.05,.85,.05,.05,.85,.05,.05,.05,
          .25,.25,.25,.25,.25,.25,.25,.25,.25,.25,.25,.25,
          .05,.05,.05,.85,.05,.05,.85,.05,.05,.05,.05,.85,
          .05,.05,.05,.85,.05,.85,.05,.05,.05,.05,.05,.85),ncol=4,byrow=TRUE)
        colnames(are_pwm)<-c("A","C","G","T")
        pdf(file.path(motifs_dir,"ARE_reference_vs_discovered.pdf"),w=12,h=6)
        comp_list <- list("ARE canonical"=t(are_pwm))
        comp_list[[paste0("Discovered: ",motifs[[1]]@consensus)]] <- motifs[[1]]@motif
        print(ggseqlogo(comp_list,method="bits",ncol=1)+
              theme(strip.text=element_text(size=11,face="bold")))
        dev.off()
        homer_known <- file.path(motifs_dir,"homer","knownResults.txt")
        if(file.exists(homer_known)){{
          cat("\\n--- Top HOMER known motifs ---\\n")
          hk <- read.delim(homer_known,header=TRUE,check.names=FALSE)
          for(i in 1:min(10,nrow(hk))) cat(sprintf("  %2d. %-35s p=%s\\n",i,hk[i,1],hk[i,3]))
          write.csv(hk[1:min(10,nrow(hk)),], file.path(motifs_dir,"homer_top_known.csv"),row.names=FALSE)
        }}
        cat("\\n[DONE]\\n")
        """)

        with open(r_script, "w") as f:
            f.write(r_code)

        # Find Rscript path explicitly (may be lost after conda env switches)
        rscript_path = shutil.which("Rscript")
        if not rscript_path:
            # Try common conda paths
            conda_base_result = run_cmd("conda info --base", self.log, check=False)
            cbase = conda_base_result.stdout.strip() if conda_base_result.stdout else ""
            for candidate in [
                f"{cbase}/envs/chipseq/bin/Rscript",
                f"{cbase}/bin/Rscript",
                "/usr/local/bin/Rscript",
                "/usr/bin/Rscript",
            ]:
                if Path(candidate).exists():
                    rscript_path = candidate
                    break

        if rscript_path:
            run_cmd(
                f"{rscript_path} {r_script}",
                self.log, "Running motif analysis (R)...",
                log_file=self.logs_dir / "motif_analysis.log",
                check=False
            )
            self.log.info("Motif analysis complete")
        else:
            self.log.warning(
                "Rscript not found. Run manually:\n"
                f"  conda activate chipseq\n"
                f"  Rscript {r_script}"
            )

    # ========================================================================
    # HTML REPORT
    # ========================================================================

    def generate_report(self):
        """Generate an HTML summary report."""
        print_banner("Generating HTML report")

        report_file = self.project_dir / "report.html"

        # Gather stats
        stats = {}

        # Alignment stats
        for name in self.all_names:
            flagstat = self.align_dir / f"{name}_flagstat.txt"
            if flagstat.exists():
                with open(flagstat) as f:
                    for line in f:
                        if "mapped (" in line:
                            stats[f"{name}_mapped"] = line.strip()
                            break

        # Peak counts
        peaks_file = self.peaks_dir / "AR_ChIP_peaks.narrowPeak"
        filtered_file = self.peaks_dir / "AR_ChIP_peaks_filtered.narrowPeak"
        stats["total_peaks"] = sum(1 for _ in open(peaks_file)) if peaks_file.exists() else "N/A"
        stats["filtered_peaks"] = sum(1 for _ in open(filtered_file)) if filtered_file.exists() else "N/A"

        # Motif summary
        motif_csv = self.motifs_dir / "motif_summary.csv"
        motif_rows = ""
        if motif_csv.exists():
            with open(motif_csv) as f:
                next(f)  # skip header
                for line in f:
                    cols = line.strip().split(",")
                    motif_rows += f"<tr><td>{cols[0]}</td><td><code>{cols[1]}</code></td><td>{cols[2]}</td></tr>\n"

        # Step timings
        timing_rows = ""
        for step_num, elapsed in sorted(self.step_times.items()):
            timing_rows += (
                f"<tr><td>Step {step_num}</td>"
                f"<td>{STEP_NAMES.get(step_num, '')}</td>"
                f"<td>{format_duration(elapsed)}</td></tr>\n"
            )

        # Image paths (relative)
        def img_tag(path: Path, alt: str) -> str:
            if path.exists():
                rel = os.path.relpath(path, self.project_dir)
                return f'<img src="{rel}" alt="{alt}" style="max-width:600px;">'
            return f"<p><em>{alt} — not found</em></p>"

        html = textwrap.dedent(f"""\
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>ChIP-seq AR Pipeline Report</title>
            <style>
                body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; background: #f8f9fa; color: #2c3e50; }}
                h1 {{ color: #1b4f72; border-bottom: 3px solid #2e86c1; padding-bottom: 10px; }}
                h2 {{ color: #2e86c1; margin-top: 30px; }}
                table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
                th {{ background: #1b4f72; color: white; padding: 10px 15px; text-align: left; }}
                td {{ padding: 8px 15px; border-bottom: 1px solid #ddd; }}
                tr:nth-child(even) {{ background: #eaf2f8; }}
                .card {{ background: white; border-radius: 8px; padding: 20px; margin: 15px 0;
                         box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
                code {{ background: #e8e8e8; padding: 2px 6px; border-radius: 3px; font-size: 0.95em; }}
                .config {{ display: grid; grid-template-columns: 200px 1fr; gap: 5px 15px; }}
                .config dt {{ font-weight: bold; }}
                .success {{ color: #27ae60; font-weight: bold; }}
                .warning {{ color: #f39c12; font-weight: bold; }}
                img {{ border-radius: 4px; margin: 10px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.15); }}
                .motif-box {{ background: #1b4f72; color: white; font-family: monospace;
                              font-size: 1.3em; padding: 15px; text-align: center;
                              border-radius: 6px; letter-spacing: 3px; margin: 10px 0; }}
            </style>
        </head>
        <body>
            <h1>ChIP-seq AR Motif Discovery — Pipeline Report</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

            <div class="card">
                <h2>Configuration</h2>
                <dl class="config">
                    <dt>ChIP SRR(s)</dt><dd>{', '.join(self.chip_srrs)}</dd>
                    <dt>Input SRR</dt><dd>{self.input_srr}</dd>
                    <dt>Read type</dt><dd>{self.read_type}</dd>
                    <dt>Genome</dt><dd>{self.genome_build}</dd>
                    <dt>Threads</dt><dd>{self.threads}</dd>
                    <dt>MACS2 q-value</dt><dd>{self.macs2_qvalue}</dd>
                    <dt>FE cutoff</dt><dd>{self.fe_cutoff}</dd>
                    <dt>Summit window</dt><dd>±{self.summit_window} bp</dd>
                    <dt>Project directory</dt><dd><code>{self.project_dir}</code></dd>
                </dl>
            </div>

            <div class="card">
                <h2>Peak Calling Results</h2>
                <table>
                    <tr><th>Metric</th><th>Value</th></tr>
                    <tr><td>Total peaks</td><td>{stats.get('total_peaks', 'N/A')}</td></tr>
                    <tr><td>Filtered peaks (FE &gt; {self.fe_cutoff})</td><td>{stats.get('filtered_peaks', 'N/A')}</td></tr>
                </table>
            </div>

            <div class="card">
                <h2>Expected Motif: ARE (Androgen Response Element)</h2>
                <div class="motif-box">A G A A C A n n n T G T T C T</div>
                <p>The canonical ARE is a 15 bp inverted palindrome composed of two hexameric half-sites.</p>
            </div>

            <div class="card">
                <h2>Discovered Motifs</h2>
                <table>
                    <tr><th>Motif</th><th>Consensus</th><th>Width (bp)</th></tr>
                    {motif_rows if motif_rows else '<tr><td colspan="3">Run Step 9 to generate motif summary</td></tr>'}
                </table>
            </div>

            <div class="card">
                <h2>Pipeline Timing</h2>
                <table>
                    <tr><th>Step</th><th>Description</th><th>Duration</th></tr>
                    {timing_rows if timing_rows else '<tr><td colspan="3">No timing data available</td></tr>'}
                </table>
            </div>

            <div class="card">
                <h2>Key Output Files</h2>
                <table>
                    <tr><th>File</th><th>Description</th></tr>
                    <tr><td><code>results/motifs/homer/homerResults.html</code></td><td>HOMER motif report</td></tr>
                    <tr><td><code>results/motifs/meme/meme-chip.html</code></td><td>MEME-ChIP motif report</td></tr>
                    <tr><td><code>results/motifs/all_motifs.pdf</code></td><td>All motif logos</td></tr>
                    <tr><td><code>results/motifs/ARE_reference_vs_discovered.pdf</code></td><td>ARE comparison</td></tr>
                    <tr><td><code>results/annotation/annotation_pie.pdf</code></td><td>Peak genomic distribution</td></tr>
                    <tr><td><code>results/annotation/GO_BP_dotplot.pdf</code></td><td>GO enrichment</td></tr>
                    <tr><td><code>results/fastqc/multiqc/multiqc_report.html</code></td><td>QC report</td></tr>
                    <tr><td><code>results/chipqc/fingerprint.png</code></td><td>ChIP enrichment QC</td></tr>
                </table>
            </div>

            <div class="card">
                <h2>ChIP-seq QC</h2>
                {img_tag(self.chipqc_dir / "fingerprint.png", "Fingerprint plot")}
            </div>

            <hr>
            <p style="color: #888; font-size: 0.9em;">
                Pipeline v{VERSION} | Generated by chipseq_ar_pipeline.py |
                <a href="https://github.com/mounemhf/chipseq-ar-motifs-python">GitHub</a>
            </p>
        </body>
        </html>
        """)

        with open(report_file, "w") as f:
            f.write(html)

        self.log.info(f"HTML report: {report_file}")

    # ========================================================================
    # RUN PIPELINE
    # ========================================================================

    def run(self):
        """Execute the full pipeline with resume support."""
        pipeline_start = time.time()

        # Print configuration
        print_banner("PIPELINE ChIP-seq AR — Motif Discovery (Python)")
        print(f"  Version       : {VERSION}")
        print(f"  ChIP SRR(s)   : {', '.join(self.chip_srrs)}")
        print(f"  Input SRR     : {self.input_srr}")
        print(f"  Read type     : {self.read_type}")
        print(f"  Genome        : {self.genome_cfg['name']}")
        print(f"                  {self.genome_cfg['description']}")
        print(f"  Threads       : {self.threads}")
        print(f"  MACS2 q-value : {self.macs2_qvalue}")
        print(f"  FE cutoff     : {self.fe_cutoff}")
        print(f"  Summit window : ±{self.summit_window} bp")
        print(f"  Project dir   : {self.project_dir}")
        print()

        # Save parameters for resume
        self.state.save_params({
            "chip": self.args.chip,
            "input": self.args.input,
            "read_type": self.read_type,
            "genome": self.genome_build,
        })

        # Check dependencies
        print_banner("Checking dependencies")
        missing = []
        for tool in ["prefetch", "fasterq-dump", "fastqc", "multiqc",
                      "trim_galore", "bowtie2", "samtools", "picard",
                      "macs2", "bedtools", "plotFingerprint", "bamCoverage",
                      "Rscript", "wget"]:
            if not check_tool(tool, self.log):
                missing.append(tool)

        # HOMER
        if not check_tool("findMotifsGenome.pl", self.log):
            missing.append("homer")

        # MEME (in meme_env — check via bash activation)
        conda_exe = shutil.which("conda")
        if conda_exe:
            result = run_cmd("conda info --base", self.log, check=False)
            cbase = result.stdout.strip() if result.stdout else ""
            if cbase:
                meme_check = run_cmd(
                    f"bash -c 'source {cbase}/etc/profile.d/conda.sh && "
                    f"conda activate meme_env && which meme-chip'",
                    self.log, check=False
                )
                if meme_check.returncode == 0:
                    self.log.info("MEME Suite (in meme_env)")
                else:
                    self.log.warning("MEME Suite — NOT FOUND in meme_env")
                    missing.append("meme")
            else:
                self.log.warning("MEME Suite — could not detect conda base")
                missing.append("meme")
        else:
            self.log.warning("MEME Suite — conda not found")
            missing.append("meme")

        if missing:
            self.log.warning(f"Missing tools: {', '.join(missing)}")
            self.log.warning("Continuing — some steps may fail")
        else:
            self.log.info("All dependencies found")

        # Determine start step
        start_step = self.args.skip_to
        if self.args.resume:
            start_step = self.state.get_resume_step()
            self.log.info(f"Resuming from Step {start_step}")

        # Steps mapping
        steps = {
            0: self.step0_download,
            1: self.step1_fastqc,
            2: self.step2_trimming,
            3: self.step3_alignment,
            4: self.step4_post_alignment,
            5: self.step5_chipqc,
            6: self.step6_peak_calling,
            7: self.step7_annotation,
            8: self.step8_motif_discovery,
            9: self.step9_motif_analysis,
        }

        # Execute steps
        for step_num, step_func in steps.items():
            if step_num < start_step:
                continue

            step_start = time.time()
            step_name = STEP_NAMES[step_num]

            try:
                step_func()
                elapsed = time.time() - step_start
                self.step_times[step_num] = elapsed
                self.state.mark_complete(step_num)
                self.log.info(f"Step {step_num} completed in {format_duration(elapsed)}")
            except Exception as e:
                elapsed = time.time() - step_start
                self.step_times[step_num] = elapsed
                self.log.error(f"Step {step_num} ({step_name}) failed: {e}")
                self.log.error(f"Resume with: python {sys.argv[0]} --resume -o {self.project_dir}")
                self.state.save()

                # Generate partial report
                try:
                    self.generate_report()
                except Exception:
                    pass

                if not self.args.force:
                    sys.exit(1)
                else:
                    self.log.warning("--force specified, continuing to next step...")

        # Generate final report
        self.generate_report()

        # Final summary
        total_time = time.time() - pipeline_start
        print_banner("PIPELINE COMPLETE")
        print(f"  Total runtime: {format_duration(total_time)}")
        print()
        print("  ┌──────────────────────────────────────────────────────────┐")
        print("  │  ★ MAIN RESULTS                                         │")
        print("  ├──────────────────────────────────────────────────────────┤")
        print(f"  │  Report    : {self.project_dir / 'report.html'}")
        print(f"  │  Peaks     : {self.peaks_dir / 'AR_ChIP_peaks_filtered.narrowPeak'}")
        print(f"  │  HOMER     : {self.motifs_dir / 'homer' / 'homerResults.html'}")
        print(f"  │  MEME-ChIP : {self.motifs_dir / 'meme' / 'meme-chip.html'}")
        print(f"  │  Logos     : {self.motifs_dir / 'all_motifs.pdf'}")
        print(f"  │  ARE comp  : {self.motifs_dir / 'ARE_reference_vs_discovered.pdf'}")
        print(f"  │  Annotation: {self.annot_dir / 'annotation_pie.pdf'}")
        print("  └──────────────────────────────────────────────────────────┘")
        print()
        print("  Expected top motif: ARE (AGAACAnnnTGTTCT)")
        print("  Co-enriched: FOXA1, HOXB13, GATA2/3, ETS")
        print()
        print("  Done! 🧬")


# ============================================================================
# CLI
# ============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="chipseq_ar_pipeline.py",
        description="ChIP-seq AR Motif Discovery Pipeline — "
                    "Identify Androgen Receptor binding motifs from public GEO data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        Examples:
          python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984
          python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -g hg19
          python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -g t2t
          python chipseq_ar_pipeline.py -c SRR1615985,SRR1615986 -i SRR1615984 -p PE -t 16
          python chipseq_ar_pipeline.py --resume -o AR_ChIP_project

        Available genomes:
          hg38  — GRCh38 (2013) : Most widely used, recommended for new analyses
          hg19  — GRCh37 (2009) : Legacy, use for older datasets aligned to hg19
          t2t   — T2T-CHM13v2.0 (2022) : Gapless telomere-to-telomere assembly

        Suggested GEO datasets:
          GSE65066 (LNCaP): -c SRR1615985 -i SRR1615984
          GSE28126 (LNCaP): -c SRR039291  -i SRR039292
        """)
    )

    # Required
    required = parser.add_argument_group("required arguments")
    required.add_argument("-c", "--chip", type=str,
                          help="ChIP SRR accession(s), comma-separated for replicates")
    required.add_argument("-i", "--input", type=str,
                          help="Input/control SRR accession")

    # Optional
    parser.add_argument("-o", "--output", type=str, default="AR_ChIP_project",
                        help="Output project directory (default: AR_ChIP_project)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("-p", "--read-type", type=str, default="SE",
                        choices=["SE", "PE"],
                        help="Read type: SE (single-end) or PE (paired-end) (default: SE)")
    parser.add_argument("-g", "--genome", type=str, default="hg38",
                        choices=["hg38", "hg19", "t2t"],
                        help="Genome build: hg38 (GRCh38, default), hg19 (GRCh37, legacy), "
                             "t2t (T2T-CHM13v2.0, telomere-to-telomere)")
    parser.add_argument("-q", "--qvalue", type=float, default=0.01,
                        help="MACS2 q-value threshold (default: 0.01)")
    parser.add_argument("-f", "--fold-enrichment", type=float, default=4,
                        help="Fold enrichment cutoff (default: 4)")
    parser.add_argument("-w", "--window", type=int, default=150,
                        help="Window around summit in bp (default: 150)")
    parser.add_argument("-s", "--skip-to", type=int, default=0,
                        help="Skip to step number 0-9 (default: 0)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last completed step")
    parser.add_argument("--force", action="store_true",
                        help="Continue pipeline even if a step fails")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {VERSION}")

    args = parser.parse_args()

    # Handle resume mode
    if args.resume:
        state_file = Path(args.output).resolve() / ".pipeline_state.json"
        if state_file.exists():
            with open(state_file) as f:
                saved = json.load(f)
            params = saved.get("params", {})
            if not args.chip:
                args.chip = params.get("chip")
            if not args.input:
                args.input = params.get("input")

    if not args.chip or not args.input:
        parser.error("Both -c/--chip and -i/--input are required (unless using --resume)")

    return args


# ============================================================================
# MAIN
# ============================================================================

def main():
    args = parse_args()
    pipeline = ChIPseqPipeline(args)
    pipeline.run()


if __name__ == "__main__":
    main()
