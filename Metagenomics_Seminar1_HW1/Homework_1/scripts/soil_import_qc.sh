#!/bin/bash
#SBATCH --job-name=soil_import_qc
#SBATCH --partition=E5-2630-TITAN_X
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_import_qc_%j.out
#SBATCH --error=/home/STUDY/FBMF/studfbmf01_18/bioinformatics/metagenomics_hw/soil/logs/soil_import_qc_%j.err

set -euo pipefail
student=/home/STUDY/FBMF/studfbmf01_18
work=$student/bioinformatics/metagenomics_hw/soil
shared=/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw
image=$student/qiime2-amplicon-2024.10.sif
manifest=$shared/soil_metadata_for_input.tsv
metadata=$shared/soil_metadata_full.tsv
sample=$shared/raw_reads/SRR17307316_1.fastq

for dir in "$work" "$work/qza" "$work/qzv" "$work/logs" "$work/fastqc" "$student/tmp" "$student/q2home"; do
  case "$(realpath "$dir")" in "$student"/*) ;; *) exit 90 ;; esac
done
test -s "$manifest"; test -s "$metadata"; test -s "$sample"
test ! -e "$work/qza/soil_reads.qza"
test ! -e "$work/qzv/soil_reads.qzv"

cd "$work"
export TMPDIR=$student/tmp MPLCONFIGDIR=$student/tmp XDG_CACHE_HOME=$student/tmp
export PYTHONDONTWRITEBYTECODE=1
export PATH=/home/STUDY/FBMF/bioinformatics/anaconda3/envs/qc/bin:$PATH
q2() {
  singularity exec -H "$student/q2home:/home/qiime2" \
    -B "$student:$student" \
    -B /home/STUDY/FBMF/bioinformatics:/home/STUDY/FBMF/bioinformatics:ro \
    --env MPLCONFIGDIR="$student/tmp",TMPDIR="$student/tmp",PYTHONDONTWRITEBYTECODE=1 \
    "$image" qiime "$@"
}

printf 'START %s\n' "$(date -Is)"
fastqc --threads 2 --outdir fastqc "$sample"
q2 tools import --type 'SampleData[SequencesWithQuality]' \
  --input-path "$manifest" --input-format SingleEndFastqManifestPhred33V2 \
  --output-path qza/soil_reads.qza
q2 demux summarize --i-data qza/soil_reads.qza --o-visualization qzv/soil_reads.qzv
q2 tools validate qza/soil_reads.qza --level max
q2 tools validate qzv/soil_reads.qzv --level max
printf 'SOIL_IMPORT_QC_COMPLETE %s\n' "$(date -Is)"
