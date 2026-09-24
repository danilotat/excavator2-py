#!/usr/bin/env bash
# Run the unmodified source snapshot; logs and artifacts stay in /work.
set -euo pipefail
export LC_ALL=C TZ=UTC
cd /work
mkdir -p runtime
micromamba list -n excavator2 --explicit > runtime/conda-explicit.txt
Rscript -e 'sessionInfo(); print(installed.packages()[,c("Package","Version")]); print(RNGkind())' > runtime/R.txt
samtools --version > runtime/samtools.txt
bedtools --version > runtime/bedtools.txt
perl -V > runtime/perl.txt
R CMD config FC > runtime/fortran-compiler.txt
R CMD config FFLAGS > runtime/fortran-flags.txt
uname -a > runtime/platform.txt
bash --version > runtime/bash.txt
(command -v awk; awk -W version) > runtime/awk.txt 2>&1 || true

# The pinned image has working original binaries but an incomplete compiler.
# Use its binaries only when accompanying Fortran sources exactly match ours.
for name in F4R FastJointSLMLibraryI; do
    cmp "/opt/excavator2/lib/F77/$name.f" "/work/source/excavator2/lib/F77/$name.f"
    cp "/opt/excavator2/lib/F77/$name.so" "/work/source/excavator2/lib/F77/$name.so"
done
sha256sum /work/source/excavator2/lib/F77/*.so > runtime/native-binaries.txt
printf '%s\n' 'Original binaries from pinned image; source identity verified with cmp.' > runtime/native-policy.txt
program=/work/source/excavator2
perl "$program/TargetPerla.pl" -v -s /work/config.yaml -o /work/target > target.log 2>&1
# Target path is generated from the fixture settings by the Python runner.
target=$(cat /work/target-path.txt)
perl "$program/EXCAVATORDataPrepare.pl" -v -s /work/samples.yaml -t "$target" -o /work/prepared -@ "$1" > prepare.log 2>&1
perl "$program/EXCAVATORDataAnalysis.pl" -v -s /work/design.yaml -i /work/prepared -t "$target" -o /work/results -@ "$1" -e paired -p "$program/parameters.yaml" > analyze.log 2>&1
