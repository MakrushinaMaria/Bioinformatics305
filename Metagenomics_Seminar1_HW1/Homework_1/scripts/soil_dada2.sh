#!/bin/bash
#SBATCH --job-name=soil_dada2
#SBATCH --partition=AMD9554
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_dada2_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_dada2_%j.err

set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/soil
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw
image=$student/qiime2-amplicon-2024.10.sif
metadata=$shared/soil_metadata_full.tsv

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

test -s qza/soil_reads.qza
for file in qza/soil_ASV_table.qza qza/soil_rep_seq.qza qza/soil_reads.dada2.stats.qza \
  qzv/soil_reads.dada2.stats.qzv qzv/soil_ASV_table.qzv qzv/soil_rep_seqs.qzv; do
  test ! -e "$file"
done

printf 'START %s\n' "$(date -Is)"
q2 dada2 denoise-single --i-demultiplexed-seqs qza/soil_reads.qza \
  --p-trim-left 25 --p-trunc-len 200 --p-max-ee 3 --p-n-threads 20 \
  --p-pooling-method pseudo --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 4 \
  --o-table qza/soil_ASV_table.qza \
  --o-representative-sequences qza/soil_rep_seq.qza \
  --o-denoising-stats qza/soil_reads.dada2.stats.qza --verbose
q2 metadata tabulate --m-input-file qza/soil_reads.dada2.stats.qza \
  --o-visualization qzv/soil_reads.dada2.stats.qzv
q2 feature-table summarize --i-table qza/soil_ASV_table.qza \
  --m-sample-metadata-file "$metadata" --o-visualization qzv/soil_ASV_table.qzv
q2 feature-table tabulate-seqs --i-data qza/soil_rep_seq.qza \
  --o-visualization qzv/soil_rep_seqs.qzv
for file in qza/soil_ASV_table.qza qza/soil_rep_seq.qza qza/soil_reads.dada2.stats.qza \
  qzv/soil_reads.dada2.stats.qzv qzv/soil_ASV_table.qzv qzv/soil_rep_seqs.qzv; do
  q2 tools validate "$file" --level max
done
printf 'SOIL_DADA2_COMPLETE %s\n' "$(date -Is)"
