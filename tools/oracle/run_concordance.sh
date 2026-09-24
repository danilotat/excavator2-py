#!/usr/bin/env bash
# Fresh old-version execution from the pinned source, not the container's R wrapper.
set -euo pipefail
export LC_ALL=C TZ=UTC
cd /work
program=/work/source/excavator2
cmp /opt/excavator2/lib/F77/FastJointSLMLibraryI.f "$program/lib/F77/FastJointSLMLibraryI.f"
cp /opt/excavator2/lib/F77/FastJointSLMLibraryI.so "$program/lib/F77/FastJointSLMLibraryI.so"
sha256sum "$program/lib/F77/FastJointSLMLibraryI.so" > native-binary.sha256
Rscript -e 'sessionInfo()' > legacy-runtime.txt
Rscript /repo/tools/oracle/characterize_analysis.R /work/legacy "$program"
Rscript /repo/tools/oracle/export.R /work/legacy/prepared /work/exports/prepared
Rscript /repo/tools/oracle/export.R /work/legacy/target /work/exports/target
