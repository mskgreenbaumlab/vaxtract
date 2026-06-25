#!/usr/bin/env bash
# Build the antVac extractor SIF on MSKCC HPC (RHEL 8, no sudo), per the singularity-build skill.
#
# The HPC has no /etc/subuid entry for the user, so apptainer builds in a ROOT-MAPPED namespace.
# That mode CANNOT run apt-get (setgroups blocked) and the miniforge base ships a fakeroot binary
# built against a newer glibc than RHEL 8 has — so `--ignore-fakeroot-command` is MANDATORY.
#
# Usage:
#   bash scale/singularity/build_sif.sh [PACKAGE_DIR] [SIF_OUT] [LOG]
# Defaults assume the iris layout; override positionally as needed. Run it from (or pass) the
# package dir that contains scale/ — do NOT rely on $BASH_SOURCE (SLURM copies scripts to spool).
set -uo pipefail   # deliberately NOT -e: the smoke-test pipes can SIGPIPE (exit 141) on success

PACKAGE_DIR="${1:-$PWD}"
SIF_OUT="${2:-/path/to/containers/antvac_extractor.sif}"
LOG="${3:-$HOME/sif_build.log}"
DEF="$PACKAGE_DIR/scale/singularity/antvac_extractor.def"

# prefer apptainer; the conda `singularity` on this HPC is apptainer underneath (accepts the flags)
BUILDER="$(command -v apptainer || command -v singularity || true)"
if [[ -z "$BUILDER" ]]; then echo "ERROR: no apptainer/singularity on PATH" | tee -a "$LOG"; exit 1; fi
if [[ ! -f "$DEF" ]]; then echo "ERROR: def not found: $DEF" | tee -a "$LOG"; exit 1; fi

export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$HOME/apptainer_cache}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$(dirname "$SIF_OUT")" "$APPTAINER_CACHEDIR"
unset APPTAINER_BIND SINGULARITY_BIND   # a stale bind source breaks %post with a fatal mount error

echo "[build] $(date) builder=$BUILDER def=$DEF sif=$SIF_OUT cache=$APPTAINER_CACHEDIR" | tee -a "$LOG"
"$BUILDER" build --fakeroot --ignore-fakeroot-command "$SIF_OUT" "$DEF" 2>&1 | tee -a "$LOG"

if [[ -f "$SIF_OUT" ]]; then
    echo "[smoke] $(date)" | tee -a "$LOG"
    "$BUILDER" exec "$SIF_OUT" claude --version 2>&1 | tee -a "$LOG" || true
    "$BUILDER" exec "$SIF_OUT" python -c \
        "import pydantic, openpyxl, pypdf, PIL, docx, claude_agent_sdk; print('deps OK')" 2>&1 | tee -a "$LOG" || true
    echo "[build] DONE -> $SIF_OUT" | tee -a "$LOG"
else
    echo "[build] ERROR: SIF not created (read the lines above the FATAL line)" | tee -a "$LOG"
    exit 1
fi
