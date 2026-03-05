# ChIP-seq AR Motif Discovery Pipeline (Python)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)]()
[![Platform](https://img.shields.io/badge/Platform-macOS%20|%20Linux%20|%20WSL-lightgrey.svg)]()

A Python-based bioinformatics pipeline for identifying **Androgen Receptor (AR) binding motifs** from public ChIP-seq data (NCBI GEO).

> This is the **Python version** of the pipeline. For the Bash version, see [chipseq-ar-motifs](https://github.com/mounemhf/chipseq-ar-motifs).

---

## Features

- **Single script** — One Python file runs the entire 10-step pipeline
- **CLI with argparse** — Clean argument parsing with help messages
- **Advanced logging** — Color-coded console output + detailed log files
- **Auto-resume** — Pipeline saves state and can resume after interruption with `--resume`
- **HTML report** — Automatic summary report generated at the end
- **Error handling** — Graceful failure with `--force` to skip broken steps

---

## Quick Start

```bash
# Install dependencies (see Installation section below)
conda activate chipseq

# Run the pipeline
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984

# Resume after interruption
python chipseq_ar_pipeline.py --resume -o AR_ChIP_project

# Skip to a specific step
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -s 6

# Continue even if a step fails
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 --force
```

---

## Installation

### macOS (Apple Silicon M1/M2/M3/M4)

```bash
# 1. Prerequisites
xcode-select --install
softwareupdate --install-rosetta --agree-to-license

# 2. Install Miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3
~/miniforge3/bin/conda init zsh
# Close and reopen Terminal

# 3. Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 4. Main environment
conda create -n chipseq -y python=3.10 \
    sra-tools fastqc multiqc trim-galore bowtie2 samtools \
    picard-slim deeptools macs2 bedtools wget

# 5. MEME environment (separate due to dependency conflicts)
conda create -n meme_env -y meme

# 6. HOMER (manual install)
mkdir -p ~/homer && cd ~/homer
curl -L -O http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
perl configureHomer.pl -install hg38
echo 'export PATH=$HOME/homer/bin:$PATH' >> ~/.zshrc
source ~/.zshrc

# 7. R packages
conda activate chipseq
conda install -y r-base r-ggplot2 r-ggrepel
Rscript -e '
install.packages("BiocManager", repos="https://cloud.r-project.org")
install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("YuLab-SMU/ggtree", force=TRUE)
BiocManager::install(c("ChIPseeker","TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db","clusterProfiler"), ask=FALSE, update=FALSE)
install.packages("ggseqlogo", repos="https://cloud.r-project.org")
BiocManager::install("universalmotif", ask=FALSE, update=FALSE)
'

# 8. Verify
conda activate chipseq
python verify_install.py
```

### Linux (Ubuntu/Debian)

```bash
# 1. Prerequisites
sudo apt update
sudo apt install -y build-essential curl wget git perl zip unzip \
    default-jre libxml2-dev libcurl4-openssl-dev libssl-dev

# 2. Install Miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
~/miniforge3/bin/conda init bash
# Close and reopen Terminal

# 3-7. Same as macOS (replace ~/.zshrc with ~/.bashrc for HOMER PATH)

# 8. Verify
conda activate chipseq
python verify_install.py
```

### Linux (CentOS/RHEL/Fedora)

```bash
# 1. Prerequisites
sudo yum groupinstall -y "Development Tools"  # or: sudo dnf groupinstall
sudo yum install -y curl wget git perl java-11-openjdk \
    libxml2-devel libcurl-devel openssl-devel

# 2-8. Same as Ubuntu
```

### Windows (WSL2)

```powershell
# In PowerShell as Administrator:
wsl --install -d Ubuntu-22.04
```

Then open Ubuntu terminal and follow the **Linux (Ubuntu)** instructions above.

---

## Usage

### All Options

```
python chipseq_ar_pipeline.py [options]

Required:
  -c, --chip          ChIP SRR accession(s), comma-separated
  -i, --input         Input/control SRR accession

Optional:
  -o, --output        Output directory          [default: AR_ChIP_project]
  -t, --threads       Number of threads         [default: 8]
  -p, --read-type     SE or PE                  [default: SE]
  -g, --genome        Genome build              [default: hg38]
  -q, --qvalue        MACS2 q-value threshold   [default: 0.01]
  -f, --fold-enrichment  FE cutoff              [default: 4]
  -w, --window        Summit window (bp)        [default: 150]
  -s, --skip-to       Skip to step (0-9)        [default: 0]
  --resume            Resume from last step
  --force             Continue on errors
  --version           Show version
```

### Examples
# hg38 (défaut) — le plus utilisé
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -g hg38

# hg19 — legacy, pour les anciennes données
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -g hg19

# T2T — assemblage complet télomère-à-télomère (2022)
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -g t2t
```bash
# Basic
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984

# Multiple replicates + paired-end
python chipseq_ar_pipeline.py -c SRR1615985,SRR1615986 -i SRR1615984 -p PE -t 16

# Resume interrupted pipeline
python chipseq_ar_pipeline.py --resume -o AR_ChIP_project

# Skip to motif discovery
python chipseq_ar_pipeline.py -c SRR1615985 -i SRR1615984 -s 8 -o AR_ChIP_project
```

### Running HOMER and MEME Manually

If the automatic MEME step fails (due to separate conda environments), run them manually:

```bash
# HOMER
conda activate chipseq
export PATH="$HOME/homer/bin:$PATH"
findMotifsGenome.pl <project>/results/motifs/top_summits.bed \
    hg38 <project>/results/motifs/homer/ -size 200 -mask -p 8

# MEME-ChIP
conda activate meme_env
meme-chip -oc <project>/results/motifs/meme \
    -maxw 20 -minw 6 -meme-nmotifs 10 -meme-mod zoops \
    -centrimo-local <project>/results/motifs/AR_summit_sequences.fa
```

---

## Pipeline Steps

| Step | Description | Tools |
|------|-------------|-------|
| 0 | Download SRA data + genome | SRA Toolkit, wget, Bowtie2 |
| 1 | Quality control | FastQC, MultiQC |
| 2 | Adapter trimming | Trim Galore, Cutadapt |
| 3 | Alignment to hg38 | Bowtie2, Samtools |
| 4 | Post-alignment filtering | Samtools, Picard |
| 5 | ChIP-seq QC | deepTools |
| 6 | Peak calling | MACS2 |
| 7 | Peak annotation | HOMER, ChIPseeker (R) |
| 8 | Motif discovery | HOMER, MEME-ChIP |
| 9 | Motif visualization | R (ggseqlogo, universalmotif) |

---

## Output

The pipeline generates an **HTML report** (`report.html`) at the end, plus all individual result files in the project directory. Key outputs:

| File | Description |
|------|-------------|
| `report.html` | Interactive HTML summary |
| `results/motifs/homer/homerResults.html` | HOMER motif report |
| `results/motifs/meme/meme-chip.html` | MEME-ChIP report |
| `results/motifs/all_motifs.pdf` | All motif logos |
| `results/motifs/ARE_reference_vs_discovered.pdf` | ARE comparison |
| `results/annotation/annotation_pie.pdf` | Peak distribution |
| `results/annotation/GO_BP_dotplot.pdf` | GO enrichment |

---

## Expected Results

The top enriched motif should be the **ARE (Androgen Response Element)**:

```
5'-AGAACAnnnTGTTCT-3'
```

Co-enriched motifs typically include FOXA1, HOXB13, GATA2/3, and ETS factors.

---

## Troubleshooting

See the detailed [troubleshooting section](https://github.com/mounemhf/chipseq-ar-motifs#troubleshooting) in the Bash version README — the same solutions apply.

Additional Python-specific tips:

- **Resume after crash**: `python chipseq_ar_pipeline.py --resume -o AR_ChIP_project`
- **Force continue on error**: add `--force` flag
- **Check state**: `cat AR_ChIP_project/.pipeline_state.json`
- **Reset state**: delete `.pipeline_state.json` to start fresh

---

## Citation

If you use this pipeline, please cite the tools listed in the [Bash version README](https://github.com/mounemhf/chipseq-ar-motifs#citation).

---

## License

MIT License. See [LICENSE](LICENSE) for details.
