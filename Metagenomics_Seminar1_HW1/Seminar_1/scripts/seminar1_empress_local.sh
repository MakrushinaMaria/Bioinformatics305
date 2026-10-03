#!/bin/bash
# Run Empress in the existing isolated Mac environment; no server installation.
set -euo pipefail

project='/Users/masamakrusina/Documents/УНИВЕРСИТЕТ/Биоинформатика рабочий'
env_dir="$project/.tools/qiime2-empress"
work="$project/assignments/seminar1_metagenomics"
mkdir -p "$project/.tools/tmp" "$work/qzv"

CONDA_PREFIX="$env_dir" TMPDIR="$project/.tools/tmp" \
  "$env_dir/bin/python" "$env_dir/bin/qiime" empress tree-plot \
  --i-tree "$work/qza/mice_rooted_tree.qza" \
  --m-feature-metadata-file "$work/qza/mice_taxonomy.qza" \
  --o-visualization "$work/qzv/mice_tree.qzv"
CONDA_PREFIX="$env_dir" TMPDIR="$project/.tools/tmp" \
  "$env_dir/bin/python" "$env_dir/bin/qiime" tools validate \
  "$work/qzv/mice_tree.qzv" --level max
