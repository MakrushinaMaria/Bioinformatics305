#!/usr/bin/env python3
"""Create compact, lossless BigWig slices and annotations for selected IGV loci."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import shutil

import pyBigWig


WORK = Path(os.environ.get("CHIPSEQ_WORK", "/home/STUDY/FBMF/studfbmf01_18/chipseq_seminar"))
OUT = WORK / "results" / "igv_assets"
CANDIDATES = WORK / "results" / "stage_comparison" / "IGV_candidates.tsv"

TRACKS = {
    "75ep": {
        "chip": WORK / "results/signal/75ep_chip.CPM.bw",
        "input": WORK / "results/signal/75ep_input.CPM.bw",
        "log2": WORK / "results/signal/75ep.log2_IP_Input.bw",
        "peaks": WORK / "results/peaks/75ep/75ep_peaks.narrowPeak",
    },
    "1-2ss": {
        "chip": WORK / "results/signal/1-2ss_chip.CPM.bw",
        "input": WORK / "results/signal/1-2ss_input.CPM.bw",
        "log2": WORK / "results/signal/1-2ss.log2_IP_Input.bw",
        "peaks": WORK / "results/peaks/1-2ss/1-2ss_peaks.narrowPeak",
    },
    "5-6ss": {
        "chip": WORK / "teacher/5-6ss/signal/foxd3_5-6ss_chip.nuclear.CPM.bw",
        "input": WORK / "teacher/5-6ss/signal/foxd3_5-6ss_input.nuclear.CPM.bw",
        "log2": WORK / "teacher/5-6ss/signal/foxd3_5-6ss_chip.log2_IP_Input.nuclear.bw",
        "peaks": WORK / "teacher/5-6ss/peaks/foxd3_5-6ss_chip_peaks.narrowPeak",
    },
    "14ss": {
        "chip": WORK / "results/signal/14ss_chip.CPM.bw",
        "input": WORK / "results/signal/14ss_input.CPM.bw",
        "log2": WORK / "results/signal/14ss.log2_IP_Input.bw",
        "peaks": WORK / "results/peaks/14ss/14ss_peaks.narrowPeak",
    },
}


def read_candidates() -> list[dict[str, str]]:
    with CANDIDATES.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def subset_bigwig(source: Path, destination: Path, candidates: list[dict[str, str]]) -> None:
    with pyBigWig.open(str(source)) as src:
        chroms = src.chroms()
        wanted: set[tuple[str, int, int, float]] = set()
        for candidate in candidates:
            chrom = candidate["chrom"]
            start = int(candidate["view_start"])
            end = min(int(candidate["view_end"]), chroms[chrom])
            for interval in src.intervals(chrom, start, end) or []:
                wanted.add((chrom, int(interval[0]), int(interval[1]), float(interval[2])))
        with pyBigWig.open(str(destination), "w") as dst:
            dst.addHeader(list(chroms.items()))
            ordered = sorted(wanted, key=lambda x: (x[0], x[1], x[2]))
            if ordered:
                dst.addEntries(
                    [x[0] for x in ordered],
                    [x[1] for x in ordered],
                    ends=[x[2] for x in ordered],
                    values=[x[3] for x in ordered],
                )


def overlaps(chrom: str, start: int, end: int, candidates: list[dict[str, str]]) -> bool:
    return any(
        chrom == c["chrom"] and start < int(c["view_end"]) and end > int(c["view_start"])
        for c in candidates
    )


def subset_gtf(candidates: list[dict[str, str]]) -> None:
    source = WORK / "reference/GRCz11.gtf"
    destination = OUT / "candidate_gene_annotations.gtf"
    with source.open() as src, destination.open("w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and overlaps(fields[0], int(fields[3]) - 1, int(fields[4]), candidates):
                dst.write(line)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = read_candidates()
    shutil.copy2(CANDIDATES, OUT / "IGV_candidates.tsv")
    manifest = []
    for stage, tracks in TRACKS.items():
        for kind in ("chip", "input", "log2"):
            source = tracks[kind]
            destination = OUT / f"{stage}.{kind}.candidate_regions.bw"
            subset_bigwig(source, destination, candidates)
            if not destination.exists() or destination.stat().st_size == 0:
                raise RuntimeError(f"Failed to create {destination}")
            manifest.append([stage, kind, destination.name, source.as_posix()])
        peak_destination = OUT / f"{stage}.peaks.narrowPeak"
        shutil.copy2(tracks["peaks"], peak_destination)
        manifest.append([stage, "peaks", peak_destination.name, tracks["peaks"].as_posix()])
    subset_gtf(candidates)
    with (OUT / "track_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "track_type", "local_asset", "source_file"])
        writer.writerows(manifest)
    for path in sorted(OUT.iterdir()):
        print(path.name, path.stat().st_size)


if __name__ == "__main__":
    main()
