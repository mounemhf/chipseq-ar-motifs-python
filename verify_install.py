#!/usr/bin/env python3
"""
Verify all dependencies for the ChIP-seq AR pipeline.
Usage: python verify_install.py
"""

import shutil
import subprocess
import sys

GREEN = "\033[0;32m"
RED = "\033[0;31m"
YELLOW = "\033[0;33m"
NC = "\033[0m"

passed = 0
failed = 0


def check(name: str, cmd: str, use_shell: bool = False) -> bool:
    global passed, failed
    try:
        if use_shell:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
            ok = result.returncode == 0
        else:
            ok = shutil.which(cmd) is not None

        if ok:
            print(f"  {name:<30} : {GREEN}OK{NC}")
            passed += 1
        else:
            print(f"  {name:<30} : {RED}MISSING{NC}")
            failed += 1
        return ok
    except Exception:
        print(f"  {name:<30} : {RED}ERROR{NC}")
        failed += 1
        return False


print("\n=== Checking chipseq environment ===")
for tool in ["prefetch", "fasterq-dump", "fastqc", "multiqc", "trim_galore",
             "bowtie2", "samtools", "picard", "macs2", "bedtools",
             "plotFingerprint", "bamCoverage", "Rscript", "wget"]:
    check(tool, tool)

print("\n=== Checking HOMER ===")
check("findMotifsGenome.pl", "findMotifsGenome.pl")
check("annotatePeaks.pl", "annotatePeaks.pl")

print("\n=== Checking MEME Suite (meme_env) ===")
check("meme-chip", "bash -c 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate meme_env && which meme-chip'", use_shell=True)

print("\n=== Checking R packages ===")
for pkg in ["ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
            "clusterProfiler", "ggplot2", "ggseqlogo", "universalmotif"]:
    check(pkg, f'Rscript -e \'if(!requireNamespace("{pkg}",quietly=TRUE)) quit(status=1)\'', use_shell=True)

print(f"\n{'─' * 45}")
if failed == 0:
    print(f"{GREEN}All {passed} checks passed!{NC}")
else:
    print(f"{YELLOW}Passed: {passed}{NC} | {RED}Failed: {failed}{NC}")
    print("See README.md for installation instructions.")
print()

sys.exit(0 if failed == 0 else 1)
