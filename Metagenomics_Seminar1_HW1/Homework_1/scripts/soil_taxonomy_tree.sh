#!/bin/bash
#SBATCH --job-name=soil_tax_tree
#SBATCH --partition=AMD9554
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_tax_tree_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_tax_tree_%j.err

set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/soil
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw
image=$student/qiime2-amplicon-2024.10.sif
classifier=$shared/2022.10.backbone.v4.nb.sklearn-1.4.2.qza
metadata=$shared/soil_metadata_full.tsv

cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
q2() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
    "$image" qiime "$@"
}

for file in qza/soil_taxonomy.qza qzv/soil_taxonomy.qzv qzv/soil_taxonomy_barplot.qzv \
  qza/aligned_soil_rep_seq.qza qza/masked_aligned_soil_rep_seq.qza \
  qza/soil_unrooted_tree.qza qza/soil_rooted_tree.qza; do test ! -e "$file"; done
test -s qza/soil_rep_seq.qza; test -s qza/soil_ASV_table.qza

printf 'START %s\n' "$(date -Is)"
q2 feature-classifier classify-sklearn --i-classifier "$classifier" \
  --i-reads qza/soil_rep_seq.qza --o-classification qza/soil_taxonomy.qza \
  --p-n-jobs 16
q2 metadata tabulate --m-input-file qza/soil_taxonomy.qza \
  --o-visualization qzv/soil_taxonomy.qzv
q2 taxa barplot --i-table qza/soil_ASV_table.qza \
  --i-taxonomy qza/soil_taxonomy.qza --m-metadata-file "$metadata" \
  --o-visualization qzv/soil_taxonomy_barplot.qzv
q2 phylogeny align-to-tree-mafft-fasttree --i-sequences qza/soil_rep_seq.qza \
  --o-alignment qza/aligned_soil_rep_seq.qza \
  --o-masked-alignment qza/masked_aligned_soil_rep_seq.qza \
  --o-tree qza/soil_unrooted_tree.qza \
  --o-rooted-tree qza/soil_rooted_tree.qza --p-n-threads 16
for file in qza/soil_taxonomy.qza qzv/soil_taxonomy.qzv qzv/soil_taxonomy_barplot.qzv \
  qza/aligned_soil_rep_seq.qza qza/masked_aligned_soil_rep_seq.qza \
  qza/soil_unrooted_tree.qza qza/soil_rooted_tree.qza; do
  q2 tools validate "$file" --level max
done
printf 'SOIL_TAXONOMY_TREE_COMPLETE %s\n' "$(date -Is)"
