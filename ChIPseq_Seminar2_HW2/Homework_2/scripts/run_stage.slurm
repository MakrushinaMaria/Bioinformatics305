#!/usr/bin/env bash
#SBATCH --job-name=chipseq_stage
#SBATCH --partition=AMD9554-common
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch run_stage.slurm <75ep|1-2ss|14ss>" >&2
    exit 2
fi

stage="$1"
case "$stage" in
    75ep|1-2ss|14ss) ;;
    *) echo "Unsupported stage: $stage" >&2; exit 2 ;;
esac

work="/home/STUDY/FBMF/studfbmf01_18/chipseq_seminar"
env_dir="/home/STUDY/FBMF/bioinformatics/chipseq/software/seminar1-chipseq"
export PATH="$env_dir/bin:$PATH"

chip="$work/data/$stage/foxd3_${stage}_chip.nuclear.bam"
input="$work/data/$stage/foxd3_${stage}_input.nuclear.bam"
signal_dir="$work/results/signal"
peak_dir="$work/results/peaks/$stage"
figure_dir="$work/results/figures"
table_dir="$work/results/tables"

mkdir -p "$signal_dir" "$peak_dir" "$figure_dir" "$table_dir" "$work/logs"

chip_bw="$signal_dir/${stage}_chip.CPM.bw"
input_bw="$signal_dir/${stage}_input.CPM.bw"
log2_bw="$signal_dir/${stage}.log2_IP_Input.bw"
peaks="$peak_dir/${stage}_peaks.narrowPeak"
summits="$peak_dir/${stage}_summits.bed"
matrix="$figure_dir/${stage}.matrix.gz"

samtools quickcheck -v "$chip" "$input"
samtools flagstat -@ 4 "$chip" > "$table_dir/${stage}_chip.flagstat.txt"
samtools flagstat -@ 4 "$input" > "$table_dir/${stage}_input.flagstat.txt"

bamCoverage \
    -b "$chip" \
    -o "$chip_bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 4

bamCoverage \
    -b "$input" \
    -o "$input_bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 4

bigwigCompare \
    -b1 "$chip_bw" \
    -b2 "$input_bw" \
    --operation log2 \
    --pseudocount 1 \
    --binSize 10 \
    --numberOfProcessors 4 \
    -o "$log2_bw"

macs3 callpeak \
    -t "$chip" \
    -c "$input" \
    -f BAMPE \
    -g 1360000000 \
    -n "$stage" \
    -q 0.01 \
    --outdir "$peak_dir"

sort -k9,9nr "$peaks" | sed -n '1,10p' > "$table_dir/${stage}_top10_peaks.tsv"

total_alignments=$(samtools view -@ 4 -c "$chip")
reads_in_peaks=$(samtools view -@ 4 -c -L "$peaks" "$chip")
frip=$(awk -v n="$reads_in_peaks" -v d="$total_alignments" 'BEGIN { if (d > 0) printf "%.8f", n/d; else print "NA" }')
{
    printf 'stage\treads_in_peaks\ttotal_chip_alignments\tFRiP\n'
    printf '%s\t%s\t%s\t%s\n' "$stage" "$reads_in_peaks" "$total_alignments" "$frip"
} > "$table_dir/${stage}_frip.tsv"

computeMatrix reference-point \
    --referencePoint center \
    -R "$summits" \
    -S "$chip_bw" "$input_bw" \
    -b 2000 \
    -a 2000 \
    --binSize 20 \
    --numberOfProcessors 4 \
    -o "$matrix"

plotHeatmap \
    -m "$matrix" \
    -out "$figure_dir/${stage}.heatmap.pdf" \
    --samplesLabel "Foxd3 ChIP" "Input" \
    --refPointLabel "summit"

plotProfile \
    -m "$matrix" \
    -out "$figure_dir/${stage}.profile.pdf" \
    --samplesLabel "Foxd3 ChIP" "Input" \
    --refPointLabel "summit"

for f in \
    "$chip_bw" "$input_bw" "$log2_bw" "$peaks" "$summits" \
    "$peak_dir/${stage}_peaks.xls" "$matrix" \
    "$figure_dir/${stage}.heatmap.pdf" "$figure_dir/${stage}.profile.pdf"; do
    test -s "$f"
done

wc -l "$peaks" "$summits"
head -n 3 "$peaks"
cat "$table_dir/${stage}_frip.tsv"
echo "STAGE_PIPELINE_COMPLETE $stage"
