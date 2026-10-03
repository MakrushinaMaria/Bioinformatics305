# Семинар 2: Foxd3 ChIP-seq на стадии 14ss

## Что мы делали

На семинаре мы анализировали ChIP-seq Foxd3 у Danio rerio на стадии 14 somites. ChIP показывает участки генома, связанные с Foxd3, а Input нужен для оценки фонового сигнала. Мы проверили BAM, построили signal tracks, вызвали peaks с помощью MACS3, рассчитали FRiP, посмотрели signal вокруг summit и сравнили 14ss с готовыми данными для 5-6ss.

## BAM и исходные библиотеки

Оба BAM содержали paired-end nuclear alignments. Все reads были mapped и properly paired; secondary, supplementary, duplicate и singleton alignments не было.

| Образец | Всего alignments | Пар reads | Mapped | Properly paired |
|---|---:|---:|---:|---:|
| Foxd3 ChIP | 41 553 880 | 20 776 940 | 41 553 880 (100%) | 41 553 880 (100%) |
| Input | 27 794 668 | 13 897 334 | 27 794 668 (100%) | 27 794 668 (100%) |

ChIP-библиотека была примерно в 1,5 раза глубже Input. Поэтому для просмотра signal мы использовали CPM normalization, а при peak calling сравнивали ChIP с Input как с background.

## Signal tracks и peaks

Для ChIP и Input мы построили CPM-normalized BigWig с bin size 10 bp. Затем рассчитали отдельный track `log2((ChIP + 1)/(Input + 1))`: положительное значение означает, что ChIP signal выше Input.

Peaks вызывали в MACS3 в режиме `BAMPE`, с effective genome size 1 360 000 000 и порогом `q <= 0.01`. Для 14ss получилось **170 peaks**.

Число peaks нельзя напрямую считать числом отдельных сайтов связывания Foxd3. Peak - это обогащённый genomic interval. Внутри него может быть несколько событий связывания, а слабые настоящие сайты могут не пройти выбранный порог.

## Как читать narrowPeak

Формат narrowPeak содержит 10 колонок:

1. `chrom` - chromosome или reference sequence;
2. `chromStart` - начало interval в zero-based coordinates;
3. `chromEnd` - конец half-open BED interval;
4. `name` - имя peak;
5. `score` - score для отображения;
6. `strand` - для ChIP-seq обычно `.`;
7. `signalValue` - сила enrichment;
8. `pValue` - `-log10(p-value)`;
9. `qValue` - `-log10(q-value)`;
10. `peak` - положение summit относительно начала interval.

Поэтому команда `sort -k9,9nr` сортирует peaks по девятой колонке, то есть по `-log10(q-value)`, начиная с самого значимого.

Самым значимым был `14ss_peak_74`:

- координаты: `NC_007126.7:35 583 506-35 583 991`;
- `signalValue = 24,9995`;
- `-log10(p) = 330,35`;
- `-log10(q) = 321,453`;
- summit offset: 230 bp, genomic position 35 583 736;
- q-value примерно **3,52 x 10^-322**.

Полный top-10 сохранён в `results/14ss_top10_peaks.tsv`.

## FRiP

FRiP мы считали по логике семинара как долю ChIP alignments, которые попали в peaks:

- alignments внутри peaks: **67 346**;
- всего ChIP alignments: **41 553 880**;
- FRiP = 67 346 / 41 553 880 = **0,00162069**, или **0,1621%**.

FRiP получился низким. Это значит, что только небольшая часть всей ChIP-библиотеки сосредоточена в вызванных peaks. Отдельные peaks при этом могут быть очень значимыми, но общую силу enrichment и качество эксперимента всё равно нужно учитывать в выводах.

## Heatmap и average profile

С помощью `computeMatrix reference-point` мы собрали CPM signal в окне 2 kb слева и справа от каждого summit с bin size 20 bp. Каждая строка matrix соответствует одному peak, а столбцы показывают положение относительно summit.

На heatmap виден signal отдельных peaks. Average profile усредняет все 170 regions. У ChIP около summit есть выраженный центральный максимум, а у Input он слабее. Значит, вызванные peaks действительно центрированы на локальном обогащении Foxd3 signal. При этом такой aggregate plot показывает общую форму сигнала, но не доказывает функцию каждого отдельного peak.

## Сравнение 5-6ss и 14ss

Затем мы сравнили peaks 14ss с готовым набором 5-6ss. Команда `bedtools intersect -u` оставляла peaks 5-6ss, которые пересекаются с 14ss, а `bedtools intersect -v` - peaks 14ss без пересечения с 5-6ss.

| Результат | Число peaks |
|---|---:|
| Все peaks 5-6ss | 839 |
| Все peaks 14ss | 170 |
| Peaks 5-6ss, пересекающиеся с 14ss | 40 |
| Peaks 14ss без пересечения с 5-6ss | 130 |

С 14ss пересекаются 40 из 839 peaks 5-6ss, то есть **4,77%**. Если смотреть со стороны набора 14ss, общими являются 40 из 170 peaks, или **23,53%**. Остальные **130 peaks 14ss** не пересекаются с 5-6ss и в этом pairwise comparison считаются специфичными для 14ss.

Видно, что binding landscape Foxd3 между стадиями сильно меняется. Однако всю разницу нельзя автоматически считать биологической: наборы также различаются по глубине, FRiP, enrichment и чувствительности peak calling. Поэтому stage-specific peak в этом анализе означает отсутствие пересечения между двумя полученными наборами, а не доказанное полное отсутствие связывания Foxd3 на другой стадии.

## Итог

Для 14ss мы получили CPM и log2 BigWig tracks, 170 MACS3 peaks, heatmap и average profile. Самый значимый peak имел q-value около 3,52 x 10^-322, но общий FRiP был низким - 0,1621%. При сравнении с 5-6ss нашли 40 общих peaks и 130 peaks, специфичных для 14ss в рамках pairwise overlap. Результаты показывают заметную перестройку связывания Foxd3 во времени, но интерпретировать её нужно вместе с различиями качества ChIP-seq между стадиями.
