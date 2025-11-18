\# Genomic Signature Explorer

Unlike simple command-line genome scripts, this project provides an interactive Streamlit UI with six-frame ORF detection, ORF density heatmaps, and direct comparison of codon usage between two organisms. Designed for quick exploratory analysis and reproducible reporting: export ORF tables (CSV), codon frequency tables, and publication-quality PNGs directly from the UI. 

\## features

A) ORF Analysis
- Detection of ORFs across all six reading frames
- ORF length distribution visualization
- ORF density heatmap across genomic windows
- Export ORF tables and figures

B) Codon Usage
- Codon frequency and per-1000 usage metrics
- Relative Synonymous Codon Usage (RSCU)
- Codon usage comparison between two groups using log2 fold change
- Global chi-square statistics for codon usage differences

C) GC3 Content
- GC content at third codon positions for all detected ORFs
- GC3 distribution plotting

### Installation

```bash
pip install -r requirements.txt


## Run locally

```bash

streamlit run app.py

