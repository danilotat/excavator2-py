#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C TZ=UTC
cd /work
program=/work/source/excavator2
cmp /opt/excavator2/lib/F77/F4R.f "$program/lib/F77/F4R.f"
cp /opt/excavator2/lib/F77/F4R.so "$program/lib/F77/F4R.so"
sha256sum "$program/lib/F77/F4R.so" > preparation-native-binary.sha256
Rscript /repo/tools/oracle/characterize_preparation_pipeline.R /work "$program"
Rscript /repo/tools/oracle/export.R /work/legacy-prepared /work/exports/legacy-prepared
Rscript /repo/tools/oracle/export.R /work/legacy/target /work/exports/preparation-target
