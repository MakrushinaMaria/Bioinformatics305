#!/bin/bash
# Reproducibility script for the successful server-side taxonomy and tree stage.
# Run inside the Seminar 1 mice work directory with the existing QIIME2 image.
set -euo pipefail

student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/mice
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/mice
image=$student/qiime2-amplicon-2024.10.sif
classifier=$shared/2022.10.backbone.v4.nb.sklearn-1.4.2.qza
metadata=$shared/metadata.tsv

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

q2 feature-classifier classify-sklearn --i-classifier "$classifier" \
  --i-reads qza/mice_rep_seq.qza --p-n-jobs 8 \
  --o-classification qza/mice_taxonomy.qza
q2 metadata tabulate --m-input-file qza/mice_taxonomy.qza \
  --o-visualization qzv/mice_taxonomy.qzv
q2 taxa barplot --i-table qza/mice_ASV_table.qza \
  --i-taxonomy qza/mice_taxonomy.qza --m-metadata-file "$metadata" \
  --o-visualization qzv/mice_taxonomy_barplot.qzv
q2 phylogeny align-to-tree-mafft-fasttree --i-sequences qza/mice_rep_seq.qza \
  --o-alignment qza/aligned_mice_rep_seq.qza \
  --o-masked-alignment qza/masked_aligned_mice_rep_seq.qza \
  --o-tree qza/mice_unrooted_tree.qza \
  --o-rooted-tree qza/mice_rooted_tree.qza --p-n-threads 8
