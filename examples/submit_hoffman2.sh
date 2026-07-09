#!/bin/bash
# Batch submission (UGE/SGE, e.g. UCLA Hoffman2) for a scVarSim run.
#   submit:   qsub examples/submit_hoffman2.sh
#   monitor:  qstat -u $USER            (state: qw=queued, r=running)
#
# Cores : `-pe shared 16` -> $NSLOTS=16, auto-used by the runner's NCORES.
#         Most steps parallelize ~linearly; raise to `-pe shared 32` for more speed.
# Memory: h_data is PER slot -> 16 x 4G = 64G total (headroom for scDesign2/R).
# Wall  : raise h_rt if you expect a long run (large chromosomes / high depth).
#
#$ -cwd
#$ -N scVarSim
#$ -o scVarSim.$JOB_ID.log
#$ -j y
#$ -pe shared 16
#$ -l h_data=4G,h_rt=24:00:00,highp
# #$ -M you@example.edu
# #$ -m bea

set -uo pipefail

# ---- EDIT THESE for your environment ----------------------------------------
CONDA_ENV="scIsoSim"                                   # conda env from environment.yml
CONFIG="examples/config.chr19_GM12878.yaml"            # your scVarSim config
# -----------------------------------------------------------------------------

# Activate the conda env. scReadSim imports rpy2 -> R must be discoverable.
. /u/local/Modules/default/init/modules.sh 2>/dev/null || true
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# rpy2 needs R_HOME + R libs on LD_LIBRARY_PATH (adjust prefix if your env differs).
export R_HOME="${CONDA_PREFIX}/lib/R"
export LD_LIBRARY_PATH="${R_HOME}/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

echo "[submit] host=$(hostname) NSLOTS=${NSLOTS:-NA} start=$(date)"
python examples/run_simulation.py --config "$CONFIG"
status=$?
echo "[submit] python exit=$status end=$(date)"
exit $status
