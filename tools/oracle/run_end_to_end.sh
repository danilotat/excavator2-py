#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C TZ=UTC
cd /work
program=/work/source/excavator2
for name in F4R FastJointSLMLibraryI; do
    cmp "/opt/excavator2/lib/F77/$name.f" "$program/lib/F77/$name.f"
    cp "/opt/excavator2/lib/F77/$name.so" "$program/lib/F77/$name.so"
done
sha256sum "$program"/lib/F77/*.so > native-binaries.sha256
Rscript -e 'sessionInfo()' > legacy-runtime.txt
perl "$program/TargetPerla.pl" -s /work/legacy-config.yaml -o /work/legacy-target > target.log 2>&1
target=/work/legacy-target/synthetic/panel/w_100
perl "$program/EXCAVATORDataPrepare.pl" -s /work/legacy-samples.yaml -t "$target" -o /work/legacy-prepared -@ 1 > prepare.log 2>&1
for mode in paired pooling; do
    perl "$program/EXCAVATORDataAnalysis.pl" -s "/work/$mode.yaml" -i /work/legacy-prepared -t "$target" -o "/work/legacy-$mode" -@ 1 -e "$mode" -p "$program/parameters.yaml" > "$mode.log" 2>&1
done
Rscript /repo/tools/oracle/export.R /work/legacy-prepared /work/exports/prepared
Rscript /repo/tools/oracle/export.R "$target" /work/exports/target
