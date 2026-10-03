#!/bin/bash
#SBATCH --job-name=mice_qc_import
#SBATCH --partition=E5-2630-TITAN_X
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/checkpoint1_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/mice/logs/checkpoint1_%j.err
set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/mice
source_dir=/home/STUDY/FBMF/bioinformatics/metagenomes/mice
for dir in "$work" "$work/qza" "$work/qzv" "$work/fastqc" "$work/multiqc" "$work/logs" "$student/tmp" "$student/q2home"; do
    case "$(realpath "$dir")" in "$student"/*) ;; *) exit 90 ;; esac
done
cd "$work"
export TMPDIR=$student/tmp
export MPLCONFIGDIR=$student/tmp
export XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export PATH=/home/STUDY/FBMF/bioinformatics/anaconda3/envs/qc/bin:$PATH
printf 'START %s\n' "$(date -Is)"
hostname
fastqc --version
multiqc --version
shopt -s nullglob
fq=("$source_dir"/*.fq)
test "${#fq[@]}" -eq 58
for file in "${fq[@]}"; do
    stem=$(basename "$file" .fq)
    test ! -e "$work/fastqc/${stem}_fastqc.html"
    test ! -e "$work/fastqc/${stem}_fastqc.zip"
done
fastqc --threads 2 --outdir "$work/fastqc" "${fq[@]}"
html=("$work"/fastqc/*_fastqc.html)
zips=("$work"/fastqc/*_fastqc.zip)
test "${#html[@]}" -eq 58
test "${#zips[@]}" -eq 58
printf 'FASTQC_COMPLETE: 58 HTML and 58 ZIP reports\n'
test ! -e "$work/multiqc/multiqc_report.html"
test ! -e "$work/multiqc/multiqc_data"
multiqc "$work/fastqc" --outdir "$work/multiqc" --no-version-check
test -s "$work/multiqc/multiqc_report.html"
printf 'MULTIQC_COMPLETE\n'
q2() {
    singularity exec -H "$student/q2home:/home/qiime2" \
      -B "$student:$student" \
      -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
      --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
      "$student/qiime2-amplicon-2024.10.sif" qiime "$@"
}
test ! -e "$work/qza/mice_reads.qza"
q2 tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$source_dir/metadata_for_input.tsv" \
    --output-path "$work/qza/mice_reads.qza" \
    --input-format PairedEndFastqManifestPhred33V2
test -s "$work/qza/mice_reads.qza"
q2 tools peek "$work/qza/mice_reads.qza"
q2 tools validate "$work/qza/mice_reads.qza" --level max
printf 'IMPORT_COMPLETE\n'
test ! -e "$work/qzv/mice_reads.qzv"
q2 demux summarize --i-data "$work/qza/mice_reads.qza" \
    --o-visualization "$work/qzv/mice_reads.qzv"
test -s "$work/qzv/mice_reads.qzv"
q2 tools peek "$work/qzv/mice_reads.qzv"
q2 tools validate "$work/qzv/mice_reads.qzv" --level max
printf 'CHECKPOINT_1_STOP %s\n' "$(date -Is)"
