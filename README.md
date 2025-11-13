FAERS Drug Safety Signal Detection Pipeline

Using FDA Adverse Event Reporting System (FAERS) Data
Author: Tim LaVake
Quarter Analyzed: 2025 Q3

Project Badges








Project Overview

This project implements a reproducible pharmacoepidemiologic signal-detection pipeline using real FAERS (FDA Adverse Event Reporting System) ASCII data. It includes:

Modular data ingestion

Case-level aggregation

Drug–event 2×2 table construction

Disproportionality analysis (ROR & PRR)

CSV, Markdown, and HTML output exports

Automated visualizations

The goal is a transparent, industry-ready workflow that demonstrates data engineering, R programming, and pharmacovigilance analytics.

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
https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html

Files used:

DRUG25Q3.txt

REAC25Q3.txt

Support tables: DEMO, OUTC, INDI, THER

Metadata:

Extract built: 20-OCT-2025 14:49:34

Start Date: 01-JUL-2025

End Date: 30-SEP-2025

2. Case-Level Aggregation

Drug and reaction files are merged using primaryid.

Each case is represented as:

primaryid | list of drug names | list of MedDRA PTs

3. Exposure & Outcome Flagging

User-defined:

DRUG_PATTERN (e.g., "OZEMPIC")

EVENT_PATTERN (e.g., "PANCREATITIS")

Flags created per case:

has_drug

has_event

4. Epidemiologic 2×2 Table
Actual Output (From Your Run)
| Exposure       | Event_Present | Event_Absent |
|----------------|---------------|--------------|
| Drug Present   | 220           | 11927        |
| Drug Absent    | 1684          | 424681       |


Exports:

results/2x2_table.md

results/2x2_pretty_table.csv

results/2x2_matrix.csv

results/2x2_table.html

5. Disproportionality Analysis Results
| Measure | Estimate | Lower | Upper |
|---------|----------|-------|-------|
| ROR     | 4.65     | 4.04  | 5.36  |
| PRR     | 4.59     | 3.99  | 5.27  |


Interpretation:
There is a statistically significant disproportionality signal for Ozempic and pancreatitis in FAERS Q3 2025.
As always, spontaneous reporting signals are hypothesis-generating, not causal.

Visualization

Figure generated:

figures/2x2_counts_plot.png


This barplot displays the a, b, c, d cell counts from the 2×2 table.

How to Run

Place all FAERS .TXT files into:

data_raw/


In RStudio, run:

source("faers_signal_detection.R")


Outputs appear in:

results/

figures/

Pharmacoepidemiology Notes

This project demonstrates real-world pharmacovigilance considerations:

FAERS lacks a denominator

Reporting is voluntary (underreporting bias)

Duplicate reports occur

Confounding by indication

Signals require follow-up epidemiologic studies

FAERS signal detection should be used for surveillance and hypothesis generation.

Planned Enhancements

Add multi-quarter ingestion

Implement EBGM / MGPS shrinkage (Empirical Bayes)

Add Bayesian hierarchical disproportionality

Time-trend visualization module

Build a Shiny dashboard

Add a Python version of the pipeline
