#!/bin/bash
#SBATCH --job-name=stickleback_illumina
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=OUT/%x_%j.out
# Paired-end insertion-site mapping for one sample.
# USAGE: sbatch stickleback_illumina.sbatch R1.fastq.gz R2.fastq.gz QUERY template.fasta out/prefix [extra options]
#   extra options are passed through, e.g. --circular, --max-dist 2, --flank 16
# Many samples: one job each, e.g.
#   for r1 in fastq/*_R1_001.fastq.gz; do
#     s=$(basename "$r1" _R1_001.fastq.gz)
#     sbatch stickleback_illumina.sbatch "$r1" "${r1/_R1_/_R2_}" QUERY template.fasta "out/$s" --circular
#   done
set -euo pipefail
if [ $# -lt 5 ]; then
    echo "USAGE: sbatch $0 R1.fastq.gz R2.fastq.gz QUERY template.fasta out/prefix [extra options]" >&2
    exit 1
fi
R1=$1; R2=$2; QUERY=$3; TEMPLATE=$4; OUT=$5; shift 5

module purge
module load python
module load pigz 2>/dev/null || true   # optional: faster decompression

mkdir -p "$(dirname "$OUT")"
SCRIPT_DIR=${STICKLEBACK_DIR:-$SLURM_SUBMIT_DIR}
# R1 and R2 are each decompressed by pigz in parallel with the Python process.
python "$SCRIPT_DIR/stickleback_illumina.py" -1 "$R1" -2 "$R2" -q "$QUERY" -t "$TEMPLATE" -o "$OUT" \
    --threads $(( ${SLURM_CPUS_PER_TASK:-4} >= 2 ? ${SLURM_CPUS_PER_TASK:-4} / 2 : 1 )) "$@"
