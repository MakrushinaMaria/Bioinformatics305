#!/bin/bash
#SBATCH --job-name=soil_core
#SBATCH --partition=AMD9554
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_core_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_core_%j.err

set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/soil
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw
image=$student/qiime2-amplicon-2024.10.sif
metadata=$shared/soil_metadata_full.tsv
sampling_depth=4500

cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=20 OPENBLAS_NUM_THREADS=20 MKL_NUM_THREADS=20
q2() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
    "$image" qiime "$@"
}

test -s qza/soil_rooted_tree.qza; test -s qza/soil_ASV_table.qza
test ! -e soil_core_metrics_results
printf 'START %s\n' "$(date -Is)"
q2 diversity core-metrics-phylogenetic \
  --i-phylogeny qza/soil_rooted_tree.qza --i-table qza/soil_ASV_table.qza \
  --p-sampling-depth "$sampling_depth" --m-metadata-file "$metadata" \
  --p-n-jobs-or-threads 20 --output-dir soil_core_metrics_results
q2 diversity alpha-group-significance \
  --i-alpha-diversity soil_core_metrics_results/shannon_vector.qza \
  --m-metadata-file "$metadata" \
  --o-visualization soil_core_metrics_results/soil_shannon_significance.qzv
for file in soil_core_metrics_results/shannon_vector.qza \
  soil_core_metrics_results/weighted_unifrac_distance_matrix.qza \
  soil_core_metrics_results/weighted_unifrac_emperor.qzv \
  soil_core_metrics_results/soil_shannon_significance.qzv; do
  q2 tools validate "$file" --level max
done
printf 'SAMPLING_DEPTH %s (all 36 samples retained; minimum ASV count 4809)\n' "$sampling_depth"
printf 'SOIL_CORE_COMPLETE %s\n' "$(date -Is)"
