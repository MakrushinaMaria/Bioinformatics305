#!/bin/bash
#SBATCH --job-name=mice_dada2
#SBATCH --partition=E5-2630-TITAN_X
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/dada2_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/dada2_%j.err
set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/mice
for dir in "$work" "$work/qza" "$work/qzv" "$work/logs" "$student/tmp" "$student/q2home"; do
 case "$(realpath "$dir")" in "$student"/*) ;; *) exit 90 ;; esac
done
cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
q2() {
 singularity exec -H "$student/q2home:/home/qiime2" -B "$student:$student" \
 -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
 --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
 "$student/qiime2-amplicon-2024.10.sif" qiime "$@"
}
printf 'START %s\n' "$(date -Is)"
hostname
for file in qza/mice_ASV_table.qza qza/mice_rep_seq.qza qza/mice_reads.dada2.stats.qza qzv/mice_reads.dada2.stats.qzv qzv/mice_ASV_table.qzv qzv/rep-seqs.qzv; do test ! -e "$file"; done
q2 dada2 denoise-paired --i-demultiplexed-seqs qza/mice_reads.qza \
 --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 220 --p-trunc-len-r 200 \
 --p-n-threads 2 --o-table qza/mice_ASV_table.qza \
 --o-representative-sequences qza/mice_rep_seq.qza \
 --o-denoising-stats qza/mice_reads.dada2.stats.qza --verbose
for file in qza/mice_ASV_table.qza qza/mice_rep_seq.qza qza/mice_reads.dada2.stats.qza; do
 test -s "$file"
 q2 tools peek "$file"
 q2 tools validate "$file" --level max
done
q2 metadata tabulate --m-input-file qza/mice_reads.dada2.stats.qza --o-visualization qzv/mice_reads.dada2.stats.qzv
q2 feature-table summarize --i-table qza/mice_ASV_table.qza --o-visualization qzv/mice_ASV_table.qzv --m-sample-metadata-file /home/STUDY/FBMF/bioinformatics/metagenomes/mice/metadata.tsv
q2 feature-table tabulate-seqs --i-data qza/mice_rep_seq.qza --o-visualization qzv/rep-seqs.qzv
for file in qzv/mice_reads.dada2.stats.qzv qzv/mice_ASV_table.qzv qzv/rep-seqs.qzv; do
 test -s "$file"
 q2 tools validate "$file" --level max
done
printf 'DADA2_AND_SUMMARIES_COMPLETE %s\n' "$(date -Is)"
