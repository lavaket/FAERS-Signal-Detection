# FAERS Drug Safety Signal Detection Pipeline

FAERS Drug Safety Signal Detection Pipeline

Using FDA Adverse Event Reporting System (FAERS) Data
Author: Tim LaVake
Quarter Analyzed: 2025 Q3

Project Badges

These are optional but standard on high-quality GitHub repos:

![R](https://img.shields.io/badge/language-R-blue.svg)
![FAERS](https://img.shields.io/badge/data-FAERS_2025Q3-green.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/status-active-lightgrey.svg)


If you add these to your README, they will render as clean badges at the top of your repo.

Project Overview

This project implements a reproducible pharmacoepidemiologic signal-detection pipeline using real FAERS (FDA Adverse Event Reporting System) ASCII data. It performs:

Modular data ingestion

Case-level aggregation

Drug–event 2×2 table construction

ROR and PRR estimation with confidence intervals

Export of tables in CSV, Markdown, and HTML

Generation of visual outputs

The goal is a practical, transparent example of data engineering + pharmacovigilance analytics suitable for industry and academic review.

Repository Structure
faers-signal-detection/
│
├── faers_signal_detection.R
│
├── data_raw/
│     ├── DRUG25Q3.txt
│     ├── REAC25Q3.txt
│     ├── DEMO25Q3.txt
│     ├── OUTC25Q3.txt
│     ├── INDI25Q3.txt
│     └── THER25Q3.txt
│
├── results/
│     ├── 2x2_pretty_table.csv
│     ├── 2x2_matrix.csv
│     ├── 2x2_table.md
│     ├── 2x2_table.html
│     ├── signal_results.csv
│
├── figures/
│     └── 2x2_counts_plot.png
│
└── README.md

Methods
1. Data Sources

FAERS ASCII data for 2025 Q3 was downloaded from:

FDA FAERS Quarterly Data Extracts:
https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html

Files used:

DRUG25Q3.txt

REAC25Q3.txt

Supporting tables (DEMO, OUTC, INDI, THER)

FAERS version metadata (from PDFs):

Extract built: 20-OCT-2025 14:49:34

Start Date: 01-JUL-2025

End Date: 30-SEP-2025

2. Case-Level Aggregation

The DRUG + REAC files are joined using primaryid.

Each case is represented as:

primaryid | vector of drug names | vector of MedDRA PTs

3. Exposure & Outcome Flagging

For a user-specified:

DRUG_PATTERN (e.g., “OZEMPIC”)

EVENT_PATTERN (e.g., “PANCREATITIS”)

We flag each case:

has_drug

has_event

4. Epidemiologic 2×2 Table
Actual 2×2 Table (Markdown Export)

Below is the real output from your run, inserted directly into README:

| Exposure       | Event_Present | Event_Absent |
|----------------|---------------|--------------|
| Drug Present   | 220           | 11927        |
| Drug Absent    | 1684          | 424681       |


This table is also exported as:

results/2x2_table.md

results/2x2_pretty_table.csv

results/2x2_matrix.csv

results/2x2_table.html

5. Disproportionality Analysis Results

Computed using continuity correction for zeros.

Your ROR / PRR Output
| Measure | Estimate | Lower | Upper |
|---------|----------|-------|-------|
| ROR     | 4.65     | 4.04  | 5.36  |
| PRR     | 4.59     | 3.99  | 5.27  |


Interpretation:
There is a statistically significant disproportionality signal for Ozempic and pancreatitis in FAERS Q3 2025.
As with all spontaneous reporting analyses, this does not imply causation.

Visualization

The following figure is automatically saved:

figures/2x2_counts_plot.png


It displays the a, b, c, and d cell counts from the 2×2 table.

How to Run

Place all FAERS .TXT files into:
data_raw/

Open RStudio and run:

source("faers_signal_detection.R")


All outputs will be generated in results/ and figures/.

Pharmacoepidemiology Notes

This project illustrates key real-world evidence considerations:

Spontaneous reporting systems have no denominator

Underreporting and reporting bias

Duplicate case submissions

Confounding by indication

Signals require follow-up studies

This system generates hypotheses, not causal conclusions.

Planned Enhancements

Multi-quarter FAERS ingestion

EBGM (Empirical Bayes Geometric Mean)

Bayesian shrinkage for sparse cells

Time-trend signal visualizations

Shiny dashboard

Python version of the pipeline