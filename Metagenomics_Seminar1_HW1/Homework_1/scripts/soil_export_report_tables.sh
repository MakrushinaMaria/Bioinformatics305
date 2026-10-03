#!/bin/bash
# Export machine-readable tables used for the Homework 1 report.
set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/soil
image=$student/qiime2-amplicon-2024.10.sif

cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1
q2() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
    "$image" qiime "$@"
}
raw() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 "$image" "$@"
}

for dir in tables/export_feature_table tables/export_taxonomy tables/export_dada2 \
  tables/export_shannon tables/export_weighted_pcoa; do test ! -e "$dir"; done
q2 tools export --input-path qza/soil_ASV_table.qza --output-path tables/export_feature_table
q2 tools export --input-path qza/soil_taxonomy.qza --output-path tables/export_taxonomy
q2 tools export --input-path qza/soil_reads.dada2.stats.qza --output-path tables/export_dada2
q2 tools export --input-path soil_core_metrics_results/shannon_vector.qza \
  --output-path tables/export_shannon
q2 tools export --input-path soil_core_metrics_results/weighted_unifrac_pcoa_results.qza \
  --output-path tables/export_weighted_pcoa
raw biom convert -i tables/export_feature_table/feature-table.biom \
  -o tables/feature-table.tsv --to-tsv
cp tables/export_taxonomy/taxonomy.tsv tables/taxonomy.tsv
cp tables/export_dada2/stats.tsv tables/dada2_stats.tsv
cp tables/export_shannon/alpha-diversity.tsv tables/shannon.tsv
cp tables/export_weighted_pcoa/ordination.txt tables/weighted_unifrac_pcoa.txt
for file in tables/feature-table.tsv tables/taxonomy.tsv tables/dada2_stats.tsv \
  tables/shannon.tsv tables/weighted_unifrac_pcoa.txt; do test -s "$file"; done
printf 'REPORT_TABLE_EXPORT_COMPLETE %s\n' "$(date -Is)"
