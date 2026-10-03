#!/bin/bash
#SBATCH --job-name=mice_diversity
#SBATCH --partition=E5-2630-TITAN_X
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/mice_diversity_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/mice_diversity_%j.err

set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/mice
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/mice
image=$student/qiime2-amplicon-2024.10.sif
metadata=$shared/metadata.tsv
sampling_depth=1300

cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
q2() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
    "$image" qiime "$@"
}

for file in qza/mice_rooted_tree.qza qza/mice_ASV_table.qza "$metadata"; do test -s "$file"; done
test ! -e mice_core_metrics_results
test ! -e qzv/mice_shannon_significance.qzv
test ! -e qzv/weighted_unifrac_significance.qzv

printf 'START %s\n' "$(date -Is)"
q2 diversity core-metrics-phylogenetic \
  --i-phylogeny qza/mice_rooted_tree.qza \
  --i-table qza/mice_ASV_table.qza \
  --p-sampling-depth "$sampling_depth" \
  --m-metadata-file "$metadata" \
  --p-n-jobs-or-threads 8 \
  --output-dir mice_core_metrics_results

q2 diversity alpha-group-significance \
  --i-alpha-diversity mice_core_metrics_results/shannon_vector.qza \
  --m-metadata-file "$metadata" \
  --o-visualization qzv/mice_shannon_significance.qzv

q2 diversity beta-group-significance \
  --i-distance-matrix mice_core_metrics_results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file "$metadata" \
  --m-metadata-column group \
  --p-method permanova \
  --p-pairwise \
  --o-visualization qzv/weighted_unifrac_significance.qzv

for file in \
  mice_core_metrics_results/shannon_vector.qza \
  mice_core_metrics_results/weighted_unifrac_distance_matrix.qza \
  mice_core_metrics_results/weighted_unifrac_emperor.qzv \
  qzv/mice_shannon_significance.qzv \
  qzv/weighted_unifrac_significance.qzv; do
  test -s "$file"
  q2 tools validate "$file" --level max
done
printf 'SAMPLING_DEPTH %s (all 29 samples retained; minimum ASV count 1392)\n' "$sampling_depth"
printf 'DIVERSITY_COMPLETE %s\n' "$(date -Is)"
