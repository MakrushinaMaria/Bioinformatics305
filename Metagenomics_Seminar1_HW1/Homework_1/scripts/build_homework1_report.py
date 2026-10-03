#!/usr/bin/env python3
"""Собрать финальный русскоязычный отчёт по Homework 1 Metagenomics."""
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer,
    Table, TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
OUT = ROOT / "output" / "Homework_1_Metagenomics_report.pdf"
OUT.parent.mkdir(parents=True, exist_ok=True)

pdfmetrics.registerFont(TTFont("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", "/System/Library/Fonts/Supplemental/Arial Bold.ttf"))

PAGE_W, PAGE_H = A4
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name="TitleArial", parent=styles["Title"], fontName="Arial-Bold", fontSize=22,
    leading=27, textColor=colors.HexColor("#17324D"), alignment=TA_CENTER,
    spaceAfter=16,
))
styles.add(ParagraphStyle(
    name="Subtitle", parent=styles["Normal"], fontName="Arial", fontSize=11,
    leading=16, textColor=colors.HexColor("#455A64"), alignment=TA_CENTER,
))
styles.add(ParagraphStyle(
    name="H1Arial", parent=styles["Heading1"], fontName="Arial-Bold", fontSize=15,
    leading=19, textColor=colors.HexColor("#17324D"), spaceBefore=8, spaceAfter=8,
))
styles.add(ParagraphStyle(
    name="H2Arial", parent=styles["Heading2"], fontName="Arial-Bold", fontSize=12,
    leading=15, textColor=colors.HexColor("#244A67"), spaceBefore=7, spaceAfter=5,
))
styles.add(ParagraphStyle(
    name="BodyArial", parent=styles["BodyText"], fontName="Arial", fontSize=9.7,
    leading=14.0, spaceAfter=7, alignment=TA_LEFT,
))
styles.add(ParagraphStyle(
    name="SmallArial", parent=styles["BodyText"], fontName="Arial", fontSize=8.2,
    leading=11.2, textColor=colors.HexColor("#37474F"), spaceAfter=5,
))
styles.add(ParagraphStyle(
    name="Caption", parent=styles["BodyText"], fontName="Arial", fontSize=8.3,
    leading=11, textColor=colors.HexColor("#455A64"), alignment=TA_CENTER,
    spaceBefore=3, spaceAfter=9,
))


def P(text, style="BodyArial"):
    return Paragraph(text, styles[style])


def figure(path, caption, width_cm, height_cm):
    return KeepTogether([
        Image(str(FIG / path), width=width_cm * cm, height=height_cm * cm),
        P(caption, "Caption"),
    ])


def table(data, widths=None, font=8.5):
    result = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#DCEAF5")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#17324D")),
        ("FONTNAME", (0, 0), (-1, 0), "Arial-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Arial"),
        ("FONTSIZE", (0, 0), (-1, -1), font),
        ("LEADING", (0, 0), (-1, -1), font + 3),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#90A4AE")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return result


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#CFD8DC"))
    canvas.line(1.7 * cm, 1.25 * cm, PAGE_W - 1.7 * cm, 1.25 * cm)
    canvas.setFont("Arial", 8)
    canvas.setFillColor(colors.HexColor("#607D8B"))
    canvas.drawString(1.7 * cm, 0.85 * cm, "Домашняя работа 1 - Metagenomics")
    canvas.drawRightString(PAGE_W - 1.7 * cm, 0.85 * cm, f"Страница {doc.page}")
    canvas.restoreState()


doc = SimpleDocTemplate(
    str(OUT), pagesize=A4, rightMargin=1.7 * cm, leftMargin=1.7 * cm,
    topMargin=1.55 * cm, bottomMargin=1.55 * cm,
    title="Домашняя работа 1 - Metagenomics: ответ почвенного микробиома на керосин",
    author="Учебный отчёт",
)

story = []
story += [Spacer(1, 2.0 * cm), P("Домашняя работа 1: Metagenomics", "TitleArial")]
story += [P("Изменение почвенного микробиома после загрязнения керосином", "Subtitle"), Spacer(1, 0.8 * cm)]
story += [P(
    "В эксперименте сравнивали контрольную почву без керосина (k0), слабое загрязнение "
    "(k5, 5 г керосина на кг почвы) и сильное загрязнение (k25, 25 г/кг). Образцы собирали "
    "на 3, 90, 180 и 360 день, по три биологических повтора для каждой комбинации условий "
    "(всего 36 образцов). Ампликоны 16S rRNA V4-V5 секвенировали на Illumina MiSeq PE250, "
    "но в анализе использовали только forward R1 reads.", "BodyArial")]
story += [Spacer(1, 0.4 * cm), table([
    ["Этап", "Использованные параметры"],
    ["DADA2", "trim-left 25; trunc-len 200; max-ee 3; pseudo-pooling; consensus chimeras; parent fold 4"],
    ["Taxonomy", "Учебный GreenGenes2 V4 naive-Bayes classifier"],
    ["Core metrics", "Phylogenetic metrics при sampling depth 4500; сохранены все 36 образцов"],
], widths=[4.0 * cm, 12.8 * cm], font=8.5)]
story += [Spacer(1, 0.6 * cm), P("Результаты и ответы", "H1Arial")]

story += [P("Задача 1. FASTQC одного R1 образца", "H1Arial")]
story += [P(
    "Для проверки выбрали образец <b>SRR17307316</b> (k0, день 3, повтор 1). В нём было "
    "10 833 reads длиной 251 nt, GC-состав 57%, reads низкого качества не отмечены. Модуль "
    "per-base sequence quality получил статус pass. В первых пяти позициях среднее качество "
    "составляло примерно Q32-Q33, в центральной части - в основном Q34-Q38, а в последнем "
    "интервале 250-251 - Q32,67. К 3'-концу разброс качества увеличивался (10-й перцентиль "
    "Q19,5), но медиана оставалась Q36. В целом качество R1 высокое, хотя на конце read заметно "
    "умеренное снижение.")]
story += [figure(
    "task1_fastqc_per_base_quality.png",
    "Рисунок 1. Стандартный FASTQC plot качества по позициям для SRR17307316 R1.",
    16.8, 11.6,
)]

story += [P("Задача 2. Почему слева удаляли 25 нуклеотидов?", "H1Arial")]
story += [P(
    "В начале read находится forward primer и прилегающая техническая область (515F: "
    "GTGCCAGCMGCCGCGGTAA, 19 nt). Параметр <b>--p-trim-left 25</b> удаляет primer и "
    "несколько следующих начальных оснований, как задано в workflow. Благодаря этому primer "
    "и возможные несовпадения в нём не воспринимаются как биологические различия между ASV. "
    "Все reads обрезаются одинаково перед построением error model.")]

story += [P("Задача 3. Статистика DADA2 и качество данных", "H1Arial")]
story += [table([
    ["Этап", "Число reads", "% от input"],
    ["Input", "641 362", "100,00%"],
    ["Filtered", "638 869", "99,61%"],
    ["Denoised", "622 542", "97,07%"],
    ["Non-chimeric", "608 362", "94,85%"],
], widths=[6.2 * cm, 5.0 * cm, 4.0 * cm])]
story += [Spacer(1, 0.25 * cm), P(
    "После полного workflow DADA2 сохранилось <b>94,85%</b> исходных reads. Максимальное "
    "число non-chimeric reads было у SRR17307475 - <b>43 518</b>, минимальное у SRR17307380 - "
    "<b>4 809</b>. Высокая доля сохранённых reads согласуется с результатом FASTQC. Глубина "
    "между образцами различается примерно в девять раз, поэтому sampling depth 4500 выбрали "
    "ниже наблюдаемого минимума; при этом в анализе сохранились все образцы.")]
story += [figure(
    "task3_dada2_retention.png",
    "Рисунок 2. Сохранение reads на этапах DADA2 для всех 36 образцов.",
    16.2, 9.45,
)]

story += [P("Задача 4. Почему подходит V4 classifier?", "H1Arial")]
story += [P(
    "PCR primers охватывали область V4-V5, но мы анализировали только forward R1 reads. После "
    "удаления 25 оснований слева и truncation на позиции 200 сохранённая последовательность "
    "начинается со стороны 515F и в основном перекрывает V4. Поэтому classifier, обученный на "
    "соответствующем участке V4, может классифицировать эти частичные reads. Для taxonomy не "
    "обязательно иметь полный V4-V5 amplicon, но read должен перекрываться с областью, на которой "
    "обучен classifier.")]

story += [PageBreak(), P("Задача 5. Таксономический состав", "H1Arial")]
story += [P("Level 3 (class)", "H2Arial")]
story += [P(
    "На 3 день состав трёх вариантов был похож. Acidobacteriae составляли 18,1-19,5%, "
    "Alphaproteobacteria - 16,1-17,0%, Actinomycetia - 10,5-11,2%, Gammaproteobacteria - "
    "9,3-10,4%, Thermoleophilia - 7,4-8,7%. Затем траектории разошлись. На 360 день в k0 "
    "и k5 преобладали Acidobacteriae (47,9% и 46,1%). В k25 доля Gammaproteobacteria достигла "
    "<b>87,5%</b>, а Acidobacteriae снизилась до 1,2%. Значит, сильное загрязнение сопровождалось "
    "выраженным отбором и доминированием одного класса.")]
story += [figure(
    "task5_taxonomy_level3.png",
    "Рисунок 3. Средний состав на Level 3 для каждой группы и временной точки (n=3).",
    17.0, 8.3,
)]
story += [P("Level 6 (genus)", "H2Arial")]
story += [P(
    "Среди классифицированных родов в загрязнённых образцах наиболее обильным был "
    "<b>Paraburkholderia_580243</b>: в среднем 25,84% для k5 и k25. Особенно заметно его доля "
    "увеличилась в k25 - 26,7% на 90 день, 50,7% на 180 день и 83,6% на 360 день. По нашим "
    "данным этот genus явно связан с сильным загрязнением керосином и поздними временными "
    "точками. Однако 16S taxonomy показывает только состав сообщества и ассоциацию с условием; "
    "она не доказывает функцию этих ASV и не позволяет напрямую утверждать, что именно они "
    "разрушают керосин.")]
story += [figure(
    "task5_taxonomy_level6.png",
    "Рисунок 4. Средний состав на Level 6 для каждой группы и временной точки (n=3).",
    17.0, 8.3,
)]

story += [PageBreak(), P("Задача 6. Weighted UniFrac", "H1Arial")]
story += [P(
    "Первые три PCoA axes объясняют <b>55,7%, 25,3% и 7,0%</b> вариации соответственно, "
    "вместе - 88,0%. Поэтому основное разделение видно по Axis 1, затем по Axis 2 и Axis 3. "
    "На 3 день k0, k5 и k25 расположены близко: трёхмерное расстояние от центроидов k5 и k25 "
    "до соответствующего центроида k0 составляет только 0,026 и 0,047. В следующие сроки обе "
    "загрязнённые группы смещаются в положительную сторону Axis 1. Для k25 это смещение "
    "последовательно усиливается до 360 дня.")]
story += [P(
    "Для k5 расстояние до k0 увеличивается до 0,146 на 90 день и 0,206 на 180 день, а затем "
    "уменьшается до 0,159 на 360 день. Это похоже на <b>частичное восстановление</b>, главным "
    "образом за счёт возврата по Axis 1. При этом поздний центроид k5 всё ещё отличается по "
    "Axis 2 и Axis 3, поэтому восстановление нельзя считать полным. Для k25 расстояние растёт "
    "от 0,047 до 0,155, 0,275 и 0,379. Признаков восстановления после сильного загрязнения за "
    "360 дней не видно.")]
story += [figure(
    "task6_weighted_unifrac_pcoa.png",
    "Рисунок 5. Weighted UniFrac PCoA. Полупрозрачные символы показывают отдельные образцы. "
    "Сплошные линии соединяют центроиды одной группы строго по времени: d3 -> d90 -> d180 -> d360; "
    "это временные траектории, а не fitted curves.",
    16.2, 12.5,
)]

story += [P("Задача 7. Alpha diversity по Shannon", "H1Arial")]
story += [P(
    "На 3 день Shannon diversity была почти одинаковой во всех вариантах: средние значения "
    "8,98 для k0, 9,05 для k5 и 8,91 для k25. В контроле показатель снижался постепенно: "
    "8,87 на 90 день, 8,65 на 180 день и 8,28 на 360 день. В k5 diversity снизилась до 7,79 "
    "и 6,58, а затем немного выросла до 6,80 на 360 день. Это согласуется с частичным, но "
    "неполным восстановлением. В k25 снижение было непрерывным и намного сильнее: 7,61, 5,91 "
    "и 2,95. Значит, сильное загрязнение керосином сопровождалось большой потерей alpha "
    "diversity без восстановления к 360 дню.")]
story += [figure(
    "task7_shannon.png",
    "Рисунок 6. Shannon entropy по вариантам и времени: среднее +/- SD, n=3.",
    16.5, 10.2,
)]

story += [PageBreak(), P("Воспроизводимость анализа", "H1Arial")]
story += [P(
    "В QIIME2 импортировали 36 forward R1 FASTQ как SampleData[SequencesWithQuality] в формате "
    "SingleEndFastqManifestPhred33V2. Основные команды и параметры:", "BodyArial")]
commands = [
    "qiime dada2 denoise-single --p-trim-left 25 --p-trunc-len 200 --p-max-ee 3 --p-n-threads 20 --p-pooling-method pseudo --p-chimera-method consensus --p-min-fold-parent-over-abundance 4",
    "qiime feature-classifier classify-sklearn --p-n-jobs 16 (учебный GreenGenes2 V4 classifier)",
    "qiime phylogeny align-to-tree-mafft-fasttree",
    "qiime diversity core-metrics-phylogenetic --p-sampling-depth 4500",
    "qiime diversity alpha-group-significance (Shannon vector и полная metadata)",
]
for command in commands:
    story.append(P("- " + command.replace("&", "&amp;"), "SmallArial"))
story += [Spacer(1, 0.25 * cm), P(
    "Sampling depth 4500 выбран ниже минимального наблюдаемого значения 4809 reads, поэтому "
    "сохранились все 36 образцов. На рисунках 3 и 4 показаны средние sample-level relative "
    "abundances. Для оценки восстановления по Weighted UniFrac использовали расстояния между "
    "центроидами, рассчитанные по PCoA Axes 1-3.", "SmallArial")]

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUT)
