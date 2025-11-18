\# Genomic Signature Explorer

A Streamlit-based bioinformatics tool for comprehensive DNA sequence analysis including:

\- 6-frame ORF detection  
\- ORF density heatmap  
\- ORF length distribution  
\- Codon usage analysis (counts, frequencies, RSCU)  
\- Codon usage comparison (log2 FC)  
\- GC3 content analysis  

\## features

\A) ORF Analysis
\- Detection of ORFs across all six reading frames
\- ORF length distribution visualization
\- ORF density heatmap across genomic windows
\- Export ORF tables and figures
\B) Codon Usage
\- Codon frequency and per-1000 usage metrics
\- Relative Synonymous Codon Usage (RSCU)
\- Codon usage comparison between two groups using log2 fold change
\- Global chi-square statistics for codon usage differences
\C) GC3 Content
\- GC content at third codon positions for all detected ORFs
\- GC3 distribution plotting

\### Installation

```bash
pip install -r requirements.txt


\## Run locally

```bash

streamlit run app.py

