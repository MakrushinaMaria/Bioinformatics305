#!/usr/bin/env python3
"""Generate the final Homework 2 ChIP-seq PDF report from verified outputs."""

from __future__ import annotations

import csv
from datetime import date
import html
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pypdfium2 as pdfium
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, PageBreak,
    Table, TableStyle, Image, KeepTogether, Preformatted,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "report_figures"
IGV = DATA / "igv_screenshots"
OUTPUT = ROOT / "output/pdf/Homework_2_ChIPseq_Foxd3_report_FINAL.pdf"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def font_path(bold: bool = False) -> str:
    candidates = [
        Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf"),
        Path("/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf"),
        Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("A Cyrillic TrueType font (Arial) was not found")


def register_fonts() -> None:
    pdfmetrics.registerFont(TTFont("HWArial", font_path(False)))
    pdfmetrics.registerFont(TTFont("HWArial-Bold", font_path(True)))
    pdfmetrics.registerFontFamily("HWArial", normal="HWArial", bold="HWArial-Bold")
    pdfmetrics.registerFont(TTFont("HWMono", "/System/Library/Fonts/Supplemental/Courier New.ttf"))
    pdfmetrics.registerFont(TTFont("HWMono-Bold", "/System/Library/Fonts/Supplemental/Courier New Bold.ttf"))
    pdfmetrics.registerFontFamily("HWMono", normal="HWMono", bold="HWMono-Bold")


def render_pdf_first_page(source: Path, destination: Path, scale: float = 2.2) -> None:
    doc = pdfium.PdfDocument(str(source))
    page = doc[0]
    bitmap = page.render(scale=scale)
    bitmap.to_pil().save(destination)
    page.close()
    doc.close()


def make_charts(summary: list[dict[str, str]], pairwise: list[dict[str, str]], specific: list[dict[str, str]]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    stages = [r["stage"] for r in summary]
    peaks = np.array([int(r["number_of_peaks"]) for r in summary])
    frip = np.array([float(r["FRiP"]) for r in summary])
    signal = np.array([float(r["median_signalValue"]) for r in summary])
    palette = ["#3A86FF", "#4CC9A7", "#F4A261", "#8E6CBB"]

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))
    for ax, values, title, ylabel, fmt in [
        (axes[0], peaks, "Called peaks", "count", "{:.0f}"),
        (axes[1], frip, "FRiP", "fraction", "{:.4f}"),
        (axes[2], signal, "Median peak signalValue", "signalValue", "{:.2f}"),
    ]:
        bars = ax.bar(stages, values, color=palette, edgecolor="#243447", linewidth=0.6)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), fmt.format(value), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "summary_metrics.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    n = len(stages)
    matrix = np.eye(n)
    for row in pairwise:
        i, j = stages.index(row["stage_A"]), stages.index(row["stage_B"])
        matrix[i, j] = matrix[j, i] = float(row["jaccard"])
    fig, ax = plt.subplots(figsize=(5.7, 4.7))
    im = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=max(0.01, np.max(matrix[np.eye(n) == 0])))
    ax.set_xticks(range(n), stages)
    ax.set_yticks(range(n), stages)
    ax.set_title("Pairwise peak-set Jaccard index", fontweight="bold")
    for i in range(n):
        for j in range(n):
            label = "1.000" if i == j else f"{matrix[i, j]:.4f}"
            ax.text(j, i, label, ha="center", va="center", color="white" if matrix[i, j] > np.max(matrix[np.eye(n) == 0]) * 0.6 else "#17202A", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIG / "jaccard_heatmap.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    spec_map = {r["stage"]: int(r["global_stage_specific_peaks"]) for r in specific}
    fig, ax = plt.subplots(figsize=(6.6, 3.5))
    vals = [spec_map[s] for s in stages]
    bars = ax.bar(stages, vals, color=palette, edgecolor="#243447", linewidth=0.6)
    ax.set_title("Globally stage-specific called peaks", fontweight="bold")
    ax.set_ylabel("peaks not overlapping any other stage")
    ax.grid(axis="y", alpha=0.22)
    ax.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, value, str(value), ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "global_specific.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_styles():
    styles = getSampleStyleSheet()
    body = ParagraphStyle("Body", parent=styles["BodyText"], fontName="HWArial", fontSize=9.2, leading=12.2, spaceAfter=6)
    small = ParagraphStyle("Small", parent=body, fontSize=7.6, leading=9.4)
    caption = ParagraphStyle("Caption", parent=small, textColor=colors.HexColor("#455A64"), alignment=TA_CENTER, spaceBefore=3, spaceAfter=8)
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontName="HWArial-Bold", fontSize=16, leading=19, textColor=colors.HexColor("#183B56"), spaceBefore=8, spaceAfter=8)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontName="HWArial-Bold", fontSize=12.5, leading=15, textColor=colors.HexColor("#26627C"), spaceBefore=8, spaceAfter=5)
    title = ParagraphStyle("Title", parent=styles["Title"], fontName="HWArial-Bold", fontSize=23, leading=28, textColor=colors.HexColor("#153E5C"), alignment=TA_CENTER)
    subtitle = ParagraphStyle("Subtitle", parent=body, fontSize=12, leading=16, alignment=TA_CENTER, textColor=colors.HexColor("#466A7F"))
    code = ParagraphStyle("Code", parent=styles["Code"], fontName="HWMono", fontSize=6.4, leading=8, leftIndent=6, rightIndent=6, borderColor=colors.HexColor("#DDE7EC"), borderWidth=0.5, borderPadding=6, backColor=colors.HexColor("#F6F9FA"), spaceAfter=7)
    code_small = ParagraphStyle("CodeSmall", parent=code, fontSize=5.15, leading=6.25, leftIndent=3, rightIndent=3, borderPadding=4)
    bullet = ParagraphStyle("Bullet", parent=body, leftIndent=13, firstLineIndent=-7, bulletIndent=4)
    return body, small, caption, h1, h2, title, subtitle, code, code_small, bullet


def p(text: str, style) -> Paragraph:
    return Paragraph(text, style)


def table(data, widths, header_rows=1, font_size=7.3, aligns=None):
    converted = []
    for row in data:
        converted.append([cell if hasattr(cell, "wrap") else Paragraph(html.escape(str(cell)), ParagraphStyle("Cell", fontName="HWArial", fontSize=font_size, leading=font_size + 2)) for cell in row])
    t = Table(converted, colWidths=widths, repeatRows=header_rows, hAlign="LEFT")
    style = [
        ("BACKGROUND", (0, 0), (-1, header_rows - 1), colors.HexColor("#1F5A78")),
        ("TEXTCOLOR", (0, 0), (-1, header_rows - 1), colors.white),
        ("FONTNAME", (0, 0), (-1, header_rows - 1), "HWArial-Bold"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#AABBC4")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]
    for row in range(header_rows, len(data)):
        if row % 2 == 0:
            style.append(("BACKGROUND", (0, row), (-1, row), colors.HexColor("#F2F6F8")))
    if aligns:
        for col, alignment in aligns.items():
            style.append(("ALIGN", (col, header_rows), (col, -1), alignment))
    t.setStyle(TableStyle(style))
    return t


class NumberedDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str):
        super().__init__(filename, pagesize=A4, rightMargin=1.45 * cm, leftMargin=1.45 * cm, topMargin=1.55 * cm, bottomMargin=1.45 * cm)
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="normal")
        self.addPageTemplates(PageTemplate(id="main", frames=frame, onPage=self._header_footer))

    def _header_footer(self, canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#D2DEE4"))
        canvas.line(self.leftMargin, A4[1] - 1.08 * cm, A4[0] - self.rightMargin, A4[1] - 1.08 * cm)
        canvas.setFont("HWArial", 7.5)
        canvas.setFillColor(colors.HexColor("#607D8B"))
        canvas.drawString(self.leftMargin, A4[1] - 0.82 * cm, "Homework 2 - Foxd3 ChIP-seq")
        canvas.drawRightString(A4[0] - self.rightMargin, 0.72 * cm, f"{doc.page}")
        canvas.restoreState()


def image_with_caption(path: Path, width: float, caption_text: str, caption_style):
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        ratio = im.height / im.width
    return [Image(str(path), width=width, height=width * ratio), Paragraph(caption_text, caption_style)]


def main() -> None:
    register_fonts()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    summary = read_tsv("four_stage_summary.tsv")
    top = read_tsv("top_peak_by_stage.tsv")
    pairwise = read_tsv("pairwise_peak_comparison.tsv")
    specific = read_tsv("global_stage_specific_counts.tsv")
    candidates = read_tsv("IGV_candidates.tsv")
    make_charts(summary, pairwise, specific)

    stage_order = ["75ep", "1-2ss", "5-6ss", "14ss"]
    by_stage = {r["stage"]: r for r in summary}
    most_peaks = max(summary, key=lambda r: int(r["number_of_peaks"]))
    best_frip = max(summary, key=lambda r: float(r["FRiP"]))
    best_signal = max(summary, key=lambda r: float(r["median_signalValue"]))
    most_similar = max(pairwise, key=lambda r: float(r["jaccard"]))
    least_similar = min(pairwise, key=lambda r: float(r["jaccard"]))
    body, small, caption, h1, h2, title, subtitle, code, code_small, bullet = build_styles()

    story = []
    story += [Spacer(1, 2.2 * cm), p("Динамика связывания Foxd3 в развитии zebrafish", title), Spacer(1, 0.45 * cm),
              p("Homework 2 - ChIP-seq", subtitle), Spacer(1, 1.2 * cm),
              p("75% epiboly → 1-2 somite stage → 5-6 somite stage → 14 somite stage", subtitle),
              Spacer(1, 2.0 * cm),
              p("Самостоятельно проанализированы стадии 75ep и 1-2ss; для сравнения использованы результаты четырёх стадий развития: 75ep, 1-2ss, 5-6ss и 14ss.", body),
              Spacer(1, 4.5 * cm), p(f"Дата формирования: {date.today().isoformat()}", subtitle), PageBreak()]

    story += [p("1. Цель работы", h1),
              p("Охарактеризовать изменение genomic binding landscape транскрипционного фактора Foxd3 от гаструляции (75ep) до миграции neural crest (14ss), одновременно отделяя наблюдаемые биологические различия от различий глубины, enrichment и общего качества ChIP-seq.", body),
              p("2. Данные и воспроизводимость", h1),
              p("Для стадий 75ep и 1-2ss использовались парные Foxd3 ChIP и Input nuclear BAM-файлы; эти две стадии были проанализированы самостоятельно. Координатной основой служил reference genome Danio rerio GRCz11 с соответствующей gene annotation. Для межстадийного сравнения дополнительно использованы результаты 5-6ss и 14ss: CPM BigWig tracks, log2(ChIP/Input), MACS3 peaks, summits, FRiP и profiles вокруг summit.", body),
              table([["Стадия", "Роль", "ChIP/Input"]] + [
                  ["75ep", "самостоятельный анализ", "nuclear BAM pair"],
                  ["1-2ss", "самостоятельный анализ", "nuclear BAM pair"],
                  ["5-6ss", "предоставленные результаты", "signal, peaks и summary figures"],
                  ["14ss", "результаты для сравнения", "ChIP/Input и результаты анализа"],
              ], [2.0*cm, 5.2*cm, 8.8*cm], font_size=8),
              p("3. Методы", h1),
              p("Для стадий, обработанных по единому workflow, выполнены samtools flagstat; bamCoverage с CPM normalization, binSize 10 и 4 потоками; bigwigCompare log2(ChIP/Input) с pseudocount 1 и binSize 10; MACS3 callpeak в paired-end режиме (BAMPE), effective genome size 1.36×10<super>9</super>, q≤0.01; FRiP по семинарному определению; computeMatrix вокруг summit ±2000 bp с binSize 20; plotHeatmap и plotProfile. Параметры анализа сохранялись одинаковыми для сопоставимости результатов.", body),
              p("FRiP = число ChIP alignments, пересекающих narrowPeak / общее число ChIP alignments. Направленный overlap A→B - доля интервалов A, имеющих хотя бы одно пересечение с B; он не равен Jaccard. Jaccard рассчитан bedtools по покрываемым основаниям: intersection/union.", body)]

    story += [p("4. Форматы файлов", h1)]
    format_rows = [
        ["Формат", "Что хранит", "Text/Binary", "Использование"],
        ["FASTA", "референсную последовательность GRCz11", "text", "координатная система и reference genome"],
        ["GTF", "координаты генов/транскриптов/экзонов и атрибуты", "text", "gene annotation в IGV"],
        ["BAM", "сжатые выравнивания reads/fragments на genome", "binary", "исходные ChIP и Input alignments"],
        ["BAI", "индекс BAM по genomic coordinates", "binary", "быстрый доступ к отдельным регионам"],
        ["BigWig", "индексированный непрерывный coverage/signal track", "binary", "CPM ChIP/Input и log2 ratio в IGV"],
        ["BED", "genomic intervals в 0-based half-open coordinates", "text", "summits, genes, interval operations"],
        ["narrowPeak", "BED6+4: interval, signalValue, -log10(p), -log10(q), summit offset", "text", "дискретные enriched regions от MACS3"],
        ["MACS3 peaks.xls", "текстовую таблицу peaks и metadata вызова (start/end, summit, pileup, p/q, fold enrichment)", "text (не Excel binary)", "аудит параметров и подробных peak metrics"],
        ["TSV", "таблицы с полями, разделёнными tab", "text", "QC, FRiP и comparisons"],
        ["deepTools matrix.gz", "gzip-сжатую матрицу regions × genomic bins для нескольких samples плюс metadata", "compressed text/binary container", "вход plotHeatmap/plotProfile"],
        ["PDF", "готовые векторные/растровые страницы графиков и отчёта", "binary", "визуализация и сдача"],
    ]
    story += [table(format_rows, [2.2*cm, 7.1*cm, 2.7*cm, 4.1*cm], font_size=6.7),
              Spacer(1, 0.2*cm), p("Ключевая цепочка: BAM (отдельные alignments) → BigWig (непрерывный количественный signal) → narrowPeak (дискретные статистически обогащённые регионы).", body), PageBreak()]

    story += [p("5. BAM QC, peak calling и enrichment", h1)]
    summary_table = [["Stage", "ChIP alignments", "Input alignments", "Peaks", "FRiP", "median signalValue", "max -log10(q)"]]
    for s in stage_order:
        r = by_stage[s]
        summary_table.append([s, f"{int(r['total_chip_alignments']):,}", f"{int(r['total_input_alignments']):,}", r["number_of_peaks"], f"{float(r['FRiP']):.6f}", f"{float(r['median_signalValue']):.2f}", f"{float(r['max_neglog10q']):.2f}"])
    story += [table(summary_table, [1.5*cm, 2.8*cm, 2.8*cm, 1.5*cm, 1.8*cm, 2.7*cm, 2.6*cm], font_size=7.0),
              p("Для 5-6ss число ChIP alignments взято из denominator предоставленной таблицы FRiP; число Input alignments вычислено как 2 × 28,525,531 paired fragments, записанных MACS3 в peaks.xls. Поэтому глубины сопоставимы как alignment counts, однако для этой стадии отсутствует отдельная статистика flagstat.", small)]
    story += image_with_caption(FIG / "summary_metrics.png", 17.0*cm, "Рисунок 1. Сравнение количества peaks, FRiP и median narrowPeak signalValue. Enrichment strength не подменяется одной придуманной метрикой: показаны независимые сопоставимые summaries.", caption)

    top_table = [["Stage", "Coordinates", "signalValue", "-log10(p)", "-log10(q)", "q"]]
    for r in top:
        top_table.append([r["stage"], f"{r['chrom']}:{int(r['start'])+1}-{r['end']}", r["signalValue"], r["neglog10p"], r["neglog10q"], r["q_value"]])
    story += [p("Наиболее значимые peaks", h2), table(top_table, [1.5*cm, 5.6*cm, 2.2*cm, 2.1*cm, 2.1*cm, 2.5*cm], font_size=7.1),
              p("narrowPeak хранит Q=-log10(q); действительный q восстановлен как 10<super>-Q</super>. Координаты в таблице показаны человеку как 1-based inclusive, тогда как исходные BED/narrowPeak start являются 0-based.", small)]

    story += [PageBreak(), p("6. Heatmaps и average profiles вокруг summit", h1),
              p("Матрицы построены по ChIP CPM и Input CPM вокруг центров MACS3 summits (±2 kb, 20-bp bins). Каждая строка heatmap соответствует peak; profile показывает средний signal по всем peaks данной стадии.", body)]
    for stage in stage_order:
        heat_png = FIG / f"{stage}.heatmap.png"
        profile_png = FIG / f"{stage}.profile.png"
        render_pdf_first_page(DATA / "figures" / f"{stage}.heatmap.pdf", heat_png)
        render_pdf_first_page(DATA / "figures" / f"{stage}.profile.pdf", profile_png)
        row = Table([[Image(str(heat_png), width=8.0*cm, height=6.0*cm), Image(str(profile_png), width=8.0*cm, height=6.0*cm)]], colWidths=[8.3*cm, 8.3*cm])
        story += [p(stage, h2), row, p(f"Рисунок: {stage} heatmap (слева) и average profile (справа), Foxd3 ChIP против Input.", caption)]

    story += [PageBreak(), p("7. Попарное сравнение genomic intervals", h1)]
    pair_table = [["A vs B", "Peaks A/B", "A→B overlap", "B→A overlap", "Jaccard"]]
    for r in pairwise:
        pair_table.append([
            f"{r['stage_A']} vs {r['stage_B']}", f"{r['peaks_A']} / {r['peaks_B']}",
            f"{r['A_peaks_overlapping_B']} ({100*float(r['A_overlap_fraction']):.1f}%)",
            f"{r['B_peaks_overlapping_A']} ({100*float(r['B_overlap_fraction']):.1f}%)", f"{float(r['jaccard']):.5f}",
        ])
    story += [table(pair_table, [3.0*cm, 2.4*cm, 3.8*cm, 3.8*cm, 2.4*cm], font_size=7.2),
              p("A→B и B→A имеют разные denominators и поэтому приведены отдельно. Pair-specific A-only = peaks_A − A→B overlap; global stage-specific ниже означает отсутствие пересечения сразу со всеми тремя другими стадиями.", small)]
    two = Table([[Image(str(FIG / "jaccard_heatmap.png"), width=8.0*cm, height=6.5*cm), Image(str(FIG / "global_specific.png"), width=8.0*cm, height=4.3*cm)]], colWidths=[8.3*cm, 8.3*cm])
    two.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))
    story += [two, p("Рисунок 2. Jaccard между peak sets (слева) и число globally stage-specific called peaks (справа).", caption),
              p(f"Наиболее похожая пара по Jaccard: <b>{most_similar['stage_A']} и {most_similar['stage_B']}</b> (J={float(most_similar['jaccard']):.5f}); наиболее различающаяся: <b>{least_similar['stage_A']} и {least_similar['stage_B']}</b> (J={float(least_similar['jaccard']):.5f}). Низкий Jaccard и большое число globally specific calls указывают на выраженную перестройку набора вызванных peaks, но ещё не доказывают полное биологическое исчезновение binding: сигнал может опуститься ниже порога peak calling.", body)]

    story += [PageBreak(), p("8. IGV-примеры", h1),
              p("Для каждого locus показаны gene annotation, CPM ChIP и Input, log2(ChIP/Input) и narrowPeak tracks четырёх стадий. Первый пример демонстрирует сильный peak 5-6ss, второй - region с peak на нескольких стадиях, третий - выраженный stage-specific region.", body)]
    for c in candidates:
        screenshot = IGV / f"{c['candidate']}.png"
        locus = f"{c['chrom']}:{int(c['view_start'])+1}-{c['view_end']}"
        igv_meaning = {
            "IGV_1": "Сильный 5-6ss peak с максимальным signalValue.",
            "IGV_2": "Общий region с peak на 75ep, 1-2ss, 5-6ss и 14ss; высота enrichment различается между стадиями.",
            "IGV_3": "Сильный 5-6ss peak, не пересекающий вызванные peaks трёх остальных стадий.",
        }[c["candidate"]]
        story += image_with_caption(screenshot, 17.0*cm,
            f"{c['candidate']}: {locus}. Опорная стадия {c['anchor_stage']}, ближайший ген {html.escape(c['nearest_gene'])} ({c['gene_distance_bp']} bp). {igv_meaning}", caption)

    story += [PageBreak(), p("9. Биологическая интерпретация и ограничения", h1),
              p("A. Наблюдаемые изменения binding landscape", h2),
              p(f"Набор вызванных Foxd3 peaks не остаётся постоянным. Максимальное число peaks наблюдалось на стадии <b>{most_peaks['stage']}</b> ({most_peaks['number_of_peaks']}); pairwise Jaccard варьировал от {float(least_similar['jaccard']):.5f} до {float(most_similar['jaccard']):.5f}. Одновременно сохраняются общие loci, продемонстрированные IGV_2, и появляются/исчезают stage-biased calls, как в IGV_3. Это согласуется с перестройкой occupation Foxd3 от gastrulation к формированию и миграции neural crest.", body),
              p("B. Свидетельства в пользу биологической интерпретации", h2),
              p("Поддержку дают согласованные изменения координат peaks, наличие globally stage-specific интервалов и локальные различия ChIP/Input signal в IGV. Участки, которые перекрываются между несколькими стадиями и сохраняют положительный log2 enrichment, выглядят как относительно стабильные элементы программы Foxd3. Stage-specific peak calls с выраженным локальным ChIP-over-Input signal являются кандидатами на динамические regulatory elements.", body),
              p("C. QC и технические confounders", h2),
              p(f"Наиболее высокий FRiP имела стадия <b>{best_frip['stage']}</b> ({float(best_frip['FRiP']):.6f}), а наиболее высокий median signalValue - <b>{best_signal['stage']}</b> ({float(best_signal['median_signalValue']):.2f}). Различия в числе peaks нельзя объяснять только биологией: меняются library depth, background/Input, immunoprecipitation efficiency, fraction of reads in peaks, signal-to-noise и статистическая мощность MACS3. Отсутствие called peak означает отсутствие достаточной статистической поддержки при данном качестве/пороге, а не доказанное отсутствие binding.", body),
              p("ChIP-seq локализует occupancy, но не устанавливает, активирует или репрессирует Foxd3 ближайший ген: enhancer может регулировать не ближайший promoter, а binding не эквивалентен изменению transcription. Для функциональной проверки нужны ATAC-seq (accessibility), H3K27ac ChIP-seq (active enhancer state) и RNA-seq (expression), желательно с perturbation Foxd3 и согласованным stage-matched design.", body),
              p("10. Вывод", h1),
              p("Foxd3 binding landscape заметно меняется между 75ep, 1-2ss, 5-6ss и 14ss, сохраняя одновременно общий core loci и stage-biased regions. Направление развития согласуется с динамической regulatory program neural crest, однако величину биологической перестройки нельзя отделить от различий ChIP quality только по peak counts. Наиболее надёжны выводы, подтверждённые одновременно interval overlap, локальным ChIP/Input signal и независимыми функциональными omics-данными.", body)]

    commands_pipeline = """# Workflow из run_stage.slurm
for stage in 75ep 1-2ss 14ss; do
  chip="data/$stage/foxd3_${stage}_chip.nuclear.bam"
  input="data/$stage/foxd3_${stage}_input.nuclear.bam"
  chip_bw="results/signal/${stage}_chip.CPM.bw"
  input_bw="results/signal/${stage}_input.CPM.bw"
  log2_bw="results/signal/${stage}.log2_IP_Input.bw"
  peakdir="results/peaks/$stage"
  peaks="$peakdir/${stage}_peaks.narrowPeak"
  summits="$peakdir/${stage}_summits.bed"

  # BAM QC
  samtools flagstat -@ 4 "$chip" > "results/tables/${stage}_chip.flagstat.txt"
  samtools flagstat -@ 4 "$input" > "results/tables/${stage}_input.flagstat.txt"

  # CPM BigWig для ChIP и Input
  bamCoverage -b "$chip" -o "$chip_bw" --normalizeUsing CPM \\
    --binSize 10 --numberOfProcessors 4
  bamCoverage -b "$input" -o "$input_bw" --normalizeUsing CPM \\
    --binSize 10 --numberOfProcessors 4
  bigwigCompare -b1 "$chip_bw" -b2 "$input_bw" --operation log2 \\
    --pseudocount 1 --binSize 10 --numberOfProcessors 4 -o "$log2_bw"

  # Peak calling и top q-value
  macs3 callpeak -t "$chip" -c "$input" -f BAMPE -g 1360000000 \\
    -n "$stage" -q 0.01 --outdir "$peakdir"
  sort -k9,9nr "$peaks" | sed -n '1,10p' > \\
    "results/tables/${stage}_top10_peaks.tsv"
  Q=$(sort -k9,9nr "$peaks" | awk 'NR==1{print $9; exit}')
  awk -v Q="$Q" 'BEGIN{printf "q=%.4e\\n",10^(-Q)}'

  # FRiP
  total=$(samtools view -@ 4 -c "$chip")
  inpeaks=$(samtools view -@ 4 -c -L "$peaks" "$chip")
  awk -v n="$inpeaks" -v d="$total" 'BEGIN{printf "%.8f\\n",n/d}'

  # Matrix, heatmap и profile вокруг summit
  computeMatrix reference-point --referencePoint center -R "$summits" \\
    -S "$chip_bw" "$input_bw" -b 2000 -a 2000 --binSize 20 \\
    --numberOfProcessors 4 -o "results/figures/${stage}.matrix.gz"
  plotHeatmap -m "results/figures/${stage}.matrix.gz" \\
    -out "results/figures/${stage}.heatmap.pdf" \\
    --samplesLabel "Foxd3 ChIP" "Input" --refPointLabel summit
  plotProfile -m "results/figures/${stage}.matrix.gz" \\
    -out "results/figures/${stage}.profile.pdf" \\
    --samplesLabel "Foxd3 ChIP" "Input" --refPointLabel summit
done"""

    commands_compare = """# Four-stage comparison из compare_stages.py
stages=(75ep 1-2ss 5-6ss 14ss)
declare -A peakfiles=(
 [75ep]="results/peaks/75ep/75ep_peaks.narrowPeak"
 [1-2ss]="results/peaks/1-2ss/1-2ss_peaks.narrowPeak"
 [5-6ss]="data/5-6ss/foxd3_5-6ss_chip_peaks.narrowPeak"
 [14ss]="results/peaks/14ss/14ss_peaks.narrowPeak"
)
for s in "${stages[@]}"; do
  bedtools sort -i "${peakfiles[$s]}" > "sorted/${s}.narrowPeak"
done

pairs=("75ep 1-2ss" "75ep 5-6ss" "75ep 14ss" \\
       "1-2ss 5-6ss" "1-2ss 14ss" "5-6ss 14ss")
for pair in "${pairs[@]}"; do
  read -r A B <<< "$pair"
  # Directed overlap A->B and B->A
  bedtools intersect -a "sorted/$A.narrowPeak" \\
    -b "sorted/$B.narrowPeak" -u | wc -l
  bedtools intersect -a "sorted/$B.narrowPeak" \\
    -b "sorted/$A.narrowPeak" -u | wc -l
  # Jaccard by genomic coverage
  bedtools jaccard -a "sorted/$A.narrowPeak" \\
    -b "sorted/$B.narrowPeak"
done

# Global stage-specific peaks: no overlap with any other stage
for s in "${stages[@]}"; do
  others=()
  for o in "${stages[@]}"; do
    [[ "$o" != "$s" ]] && others+=("sorted/$o.narrowPeak")
  done
  cat "${others[@]}" | bedtools sort -i - | bedtools merge -i - \\
    > "${s}.others.union.bed"
  bedtools intersect -a "sorted/$s.narrowPeak" \\
    -b "${s}.others.union.bed" -v > "${s}.global_specific.narrowPeak"
done

# Координаты IGV: overlap membership и ближайшие гены
for o in 75ep 1-2ss 14ss; do
  bedtools intersect -a sorted/5-6ss.narrowPeak \\
    -b "sorted/$o.narrowPeak" -u > "5-6ss.shared_with_${o}.bed"
done
bedtools closest -a candidate_pool.bed -b reference/genes.bed -d

# Core fragment из prepare_igv_assets.py (Python)
with pyBigWig.open(source_bw) as bw:
    intervals = bw.intervals(chrom, view_start, view_end) or []
# GTF оставлен для features, пересекающих выбранные loci.

python scripts/compare_stages.py
python scripts/prepare_igv_assets.py"""
    command_table = Table(
        [[Preformatted(commands_pipeline, code_small), Preformatted(commands_compare, code_small)]],
        colWidths=[8.25*cm, 8.25*cm], hAlign="LEFT",
    )
    command_table.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 4)]))
    story += [PageBreak(), p("Приложение. Фактически использованный workflow", h1), command_table,
              p("Полные версии использованных scripts сохранены как run_stage.slurm, compare_stages.py и prepare_igv_assets.py.", small)]

    doc = NumberedDocTemplate(str(OUTPUT))
    doc.build(story)
    print(OUTPUT)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
