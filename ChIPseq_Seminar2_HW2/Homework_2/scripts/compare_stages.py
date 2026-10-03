#!/usr/bin/env python3
"""Build four-stage Foxd3 ChIP-seq comparison tables and IGV candidates."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys


WORK = Path(os.environ.get("CHIPSEQ_WORK", "/home/STUDY/FBMF/studfbmf01_18/chipseq_seminar"))
OUT = WORK / "results" / "stage_comparison"
TABLES = WORK / "results" / "tables"
STAGES = ["75ep", "1-2ss", "5-6ss", "14ss"]
PEAKS = {
    "75ep": WORK / "results/peaks/75ep/75ep_peaks.narrowPeak",
    "1-2ss": WORK / "results/peaks/1-2ss/1-2ss_peaks.narrowPeak",
    "5-6ss": WORK / "teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.narrowPeak",
    "14ss": WORK / "results/peaks/14ss/14ss_peaks.narrowPeak",
}


def run(args: list[str], *, stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        return subprocess.run(args, check=True, text=True, capture_output=True).stdout
    with stdout_path.open("w") as handle:
        subprocess.run(args, check=True, text=True, stdout=handle)
    return ""


def rows(path: Path) -> list[list[str]]:
    result = []
    with path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                result.append(line.rstrip("\n").split("\t"))
    return result


def count_lines(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def flagstat_values(path: Path) -> dict[str, str]:
    if not path.exists():
        return {"total": "NA", "mapped": "NA", "properly_paired": "NA"}
    result = {"total": "NA", "mapped": "NA", "properly_paired": "NA"}
    for line in path.read_text().splitlines():
        value = line.split(" + ", 1)[0]
        if " in total " in line:
            result["total"] = value
        elif " mapped (" in line and "primary mapped" not in line:
            result["mapped"] = value
        elif " properly paired (" in line:
            result["properly_paired"] = value
    return result


def scientific_q(neg_log10_q: float) -> str:
    exponent = math.floor(neg_log10_q)
    mantissa = 10 ** (exponent - neg_log10_q)
    return f"{mantissa:.4g}e-{exponent}"


def write_tsv(path: Path, header: list[str], data: list[list[object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(data)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    sorted_dir = OUT / "sorted_peaks"
    sorted_dir.mkdir(exist_ok=True)

    missing = [str(path) for path in PEAKS.values() if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing peak files: " + ", ".join(missing))

    sorted_peaks: dict[str, Path] = {}
    for stage in STAGES:
        path = sorted_dir / f"{stage}.sorted.narrowPeak"
        run(["bedtools", "sort", "-i", str(PEAKS[stage])], stdout_path=path)
        sorted_peaks[stage] = path

    peak_counts = {stage: count_lines(sorted_peaks[stage]) for stage in STAGES}
    comparison_rows: list[list[object]] = []
    for i, a in enumerate(STAGES):
        for b in STAGES[i + 1 :]:
            a_overlap = int(run(["bedtools", "intersect", "-a", str(sorted_peaks[a]), "-b", str(sorted_peaks[b]), "-u"]).count("\n"))
            b_overlap = int(run(["bedtools", "intersect", "-a", str(sorted_peaks[b]), "-b", str(sorted_peaks[a]), "-u"]).count("\n"))
            jaccard_lines = run(["bedtools", "jaccard", "-a", str(sorted_peaks[a]), "-b", str(sorted_peaks[b])]).strip().splitlines()
            keys = jaccard_lines[0].split("\t")
            vals = jaccard_lines[1].split("\t")
            jac = dict(zip(keys, vals))
            comparison_rows.append([
                a, b, peak_counts[a], peak_counts[b],
                a_overlap, f"{a_overlap / peak_counts[a]:.6f}",
                b_overlap, f"{b_overlap / peak_counts[b]:.6f}",
                peak_counts[a] - a_overlap, peak_counts[b] - b_overlap,
                jac["intersection"], jac["union"], jac["jaccard"], jac["n_intersections"],
            ])
    write_tsv(
        OUT / "pairwise_peak_comparison.tsv",
        ["stage_A", "stage_B", "peaks_A", "peaks_B", "A_peaks_overlapping_B", "A_overlap_fraction",
         "B_peaks_overlapping_A", "B_overlap_fraction", "A_only_vs_B", "B_only_vs_A",
         "intersection_bp", "union_bp", "jaccard", "n_intersections"],
        comparison_rows,
    )

    global_specific_rows: list[list[object]] = []
    global_specific_paths: dict[str, Path] = {}
    for stage in STAGES:
        combined = OUT / f"{stage}.other_stages.concat.bed"
        with combined.open("w") as handle:
            for other in STAGES:
                if other != stage:
                    handle.write(sorted_peaks[other].read_text())
        combined_sorted = OUT / f"{stage}.other_stages.sorted.bed"
        other_union = OUT / f"{stage}.other_stages.union.bed"
        run(["bedtools", "sort", "-i", str(combined)], stdout_path=combined_sorted)
        run(["bedtools", "merge", "-i", str(combined_sorted)], stdout_path=other_union)
        specific = OUT / f"{stage}.global_stage_specific.narrowPeak"
        run(["bedtools", "intersect", "-a", str(sorted_peaks[stage]), "-b", str(other_union), "-v"], stdout_path=specific)
        global_specific_paths[stage] = specific
        n_specific = count_lines(specific)
        global_specific_rows.append([stage, peak_counts[stage], n_specific, f"{n_specific / peak_counts[stage]:.6f}"])
        combined.unlink()
        combined_sorted.unlink()
        other_union.unlink()
    write_tsv(OUT / "global_stage_specific_counts.tsv", ["stage", "total_peaks", "global_stage_specific_peaks", "fraction"], global_specific_rows)

    summary_rows: list[list[object]] = []
    top_rows: list[list[object]] = []
    for stage in STAGES:
        peak_rows = rows(PEAKS[stage])
        signals = [float(r[6]) for r in peak_rows]
        q_scores = [float(r[8]) for r in peak_rows]
        top = max(peak_rows, key=lambda r: float(r[8]))
        if stage == "5-6ss":
            chip_stats = {"total": "18853396", "mapped": "NA", "properly_paired": "NA"}
            input_stats = {"total": "57051062", "mapped": "NA", "properly_paired": "NA"}
            frip = "0.00235952"
            depth_note = "Total alignments: ChIP from teacher FRiP denominator; Input inferred as 2 x 28,525,531 paired fragments in MACS3 peaks.xls."
        else:
            chip_stats = flagstat_values(TABLES / f"{stage}_chip.flagstat.txt")
            input_stats = flagstat_values(TABLES / f"{stage}_input.flagstat.txt")
            frip_rows = rows(TABLES / f"{stage}_frip.tsv")
            frip = frip_rows[-1][3]
            depth_note = "samtools flagstat total alignments"
        summary_rows.append([
            stage, chip_stats["total"], input_stats["total"], chip_stats["mapped"], chip_stats["properly_paired"],
            peak_counts[stage], frip, f"{statistics.median(signals):.6f}", f"{statistics.mean(signals):.6f}",
            f"{max(signals):.6f}", f"{statistics.median(q_scores):.6f}", f"{max(q_scores):.6f}", depth_note,
        ])
        q_score = float(top[8])
        top_rows.append([stage, top[0], top[1], top[2], top[3], top[6], top[7], top[8], scientific_q(q_score)])
    write_tsv(
        OUT / "four_stage_summary.tsv",
        ["stage", "total_chip_alignments", "total_input_alignments", "mapped_chip", "properly_paired_chip",
         "number_of_peaks", "FRiP", "median_signalValue", "mean_signalValue", "max_signalValue",
         "median_neglog10q", "max_neglog10q", "depth_source_note"],
        summary_rows,
    )
    write_tsv(
        OUT / "top_peak_by_stage.tsv",
        ["stage", "chrom", "start", "end", "peak_name", "signalValue", "neglog10p", "neglog10q", "q_value"],
        top_rows,
    )

    # IGV candidate 1: strongest 5-6ss peak by signalValue.
    teacher_rows = rows(PEAKS["5-6ss"])

    def nearest_genes(peak_rows: list[list[str]], label: str) -> dict[tuple[str, str, str, str], tuple[str, int, int, int]]:
        temp = OUT / f".{label}.candidate_pool.bed"
        with temp.open("w") as handle:
            for r in peak_rows:
                handle.write("\t".join(r) + "\n")
        text = run(["bedtools", "closest", "-a", str(temp), "-b", str(WORK / "reference/genes.bed"), "-d"])
        temp.unlink()
        result = {}
        for line in text.splitlines():
            r = line.split("\t")
            result[(r[0], r[1], r[2], r[3])] = (r[13], int(r[11]), int(r[12]), int(r[16]))
        return result

    teacher_genes = nearest_genes(teacher_rows, "teacher")
    teacher_near_gene = [r for r in teacher_rows if teacher_genes[(r[0], r[1], r[2], r[3])][3] <= 20000]
    strongest_teacher = max(teacher_near_gene or teacher_rows, key=lambda r: float(r[6]))

    # Candidate 2: strongest 5-6ss peak shared with the most other stages (at least one).
    overlap_membership: dict[tuple[str, str, str, str], set[str]] = {}
    for other in STAGES:
        if other == "5-6ss":
            continue
        text = run(["bedtools", "intersect", "-a", str(sorted_peaks["5-6ss"]), "-b", str(sorted_peaks[other]), "-u"])
        for line in text.splitlines():
            r = line.split("\t")
            overlap_membership.setdefault((r[0], r[1], r[2], r[3]), set()).add(other)
    shared_pool = [r for r in teacher_rows if (r[0], r[1], r[2], r[3]) in overlap_membership]
    if not shared_pool:
        raise RuntimeError("No shared 5-6ss peaks found")
    shared_near_gene = [r for r in shared_pool if teacher_genes[(r[0], r[1], r[2], r[3])][3] <= 20000]
    shared = max(
        shared_near_gene or shared_pool,
        key=lambda r: (len(overlap_membership[(r[0], r[1], r[2], r[3])]), float(r[6])),
    )
    shared_stages = sorted(overlap_membership[(shared[0], shared[1], shared[2], shared[3])])

    # Candidate 3: strongest globally stage-specific peak by signalValue.
    specific_candidates: list[tuple[str, list[str], tuple[str, int, int, int]]] = []
    for stage, path in global_specific_paths.items():
        stage_rows = rows(path)
        stage_genes = nearest_genes(stage_rows, f"specific_{stage}") if stage_rows else {}
        specific_candidates.extend((stage, r, stage_genes[(r[0], r[1], r[2], r[3])]) for r in stage_rows)
    near_gene_specific = [item for item in specific_candidates if item[2][3] <= 20000]
    specific_stage, specific, specific_gene = max(near_gene_specific or specific_candidates, key=lambda item: float(item[1][6]))

    def candidate_row(label: str, peak: list[str], kind: str, rationale: str, stage: str,
                      gene: tuple[str, int, int, int] | None = None) -> list[object]:
        if gene is None:
            gene = teacher_genes[(peak[0], peak[1], peak[2], peak[3])]
        gene_name, gene_start, gene_end, distance = gene
        start = max(0, min(int(peak[1]) - 2500, gene_start - 500))
        end = max(int(peak[2]) + 2500, min(gene_end, gene_start + 5000) + 500)
        return [label, peak[0], start, end, int(peak[1]), int(peak[2]), stage, kind,
                peak[6], peak[8], gene_name, distance, rationale]

    candidates = [
        candidate_row("IGV_1", strongest_teacher, "strong_5-6ss", "Highest 5-6ss signalValue; shows a clear teacher-example peak.", "5-6ss"),
        candidate_row("IGV_2", shared, "multi_stage_shared", "5-6ss peak overlapping peaks in: " + ", ".join(shared_stages) + ".", "5-6ss"),
        candidate_row("IGV_3", specific, "global_stage_specific", f"Highest-signal peak unique to {specific_stage} among all four called-peak sets.", specific_stage, specific_gene),
    ]
    write_tsv(
        OUT / "IGV_candidates.tsv",
        ["candidate", "chrom", "view_start", "view_end", "peak_start", "peak_end", "anchor_stage", "type",
         "anchor_signalValue", "anchor_neglog10q", "nearest_gene", "gene_distance_bp", "rationale"],
        candidates,
    )

    print((OUT / "four_stage_summary.tsv").read_text())
    print((OUT / "pairwise_peak_comparison.tsv").read_text())
    print((OUT / "global_stage_specific_counts.tsv").read_text())
    print((OUT / "IGV_candidates.tsv").read_text())


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
