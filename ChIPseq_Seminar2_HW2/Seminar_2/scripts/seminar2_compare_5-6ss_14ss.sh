#!/usr/bin/env bash
# Exact Seminar 2 comparison of teacher 5-6ss peaks with completed 14ss peaks.
# Run on MIPT; all outputs stay inside the student's workspace.

set -euo pipefail

work="/home/STUDY/FBMF/studfbmf01_18/chipseq_seminar"
env_dir="/home/STUDY/FBMF/bioinformatics/chipseq/software/seminar1-chipseq"
export PATH="$env_dir/bin:$PATH"

peaks_5_6ss="$work/teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.narrowPeak"
peaks_14ss="$work/results/peaks/14ss/14ss_peaks.narrowPeak"
out_dir="$work/results/stage_comparison/seminar2"

mkdir -p "$out_dir"

bedtools intersect \
    -a "$peaks_5_6ss" \
    -b "$peaks_14ss" \
    -u > "$out_dir/shared_5-6ss_14ss.bed"

bedtools intersect \
    -a "$peaks_14ss" \
    -b "$peaks_5_6ss" \
    -v > "$out_dir/only_14ss.bed"

{
    printf 'category\tcount\n'
    printf 'all_5-6ss_peaks\t%s\n' "$(wc -l < "$peaks_5_6ss")"
    printf 'all_14ss_peaks\t%s\n' "$(wc -l < "$peaks_14ss")"
    printf 'shared_5-6ss_peaks_overlapping_14ss\t%s\n' "$(wc -l < "$out_dir/shared_5-6ss_14ss.bed")"
    printf 'only_14ss_peaks_vs_5-6ss\t%s\n' "$(wc -l < "$out_dir/only_14ss.bed")"
} > "$out_dir/5-6ss_vs_14ss_counts.tsv"

cat "$out_dir/5-6ss_vs_14ss_counts.tsv"
