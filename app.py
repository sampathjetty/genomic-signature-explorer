"""
orfheatmap.py

ORF Heatmap Generator + Codon Bias & GC3 Analyzer
Cleaned, unified and optimized single-file Streamlit app.

Features:
- 6-frame ORF finder (forward + reverse) with robust scanning (no premature frame breaks)
- ORF table export (CSV)
- ORF length distribution (hist + KDE)
- ORF density heatmap (windows, configurable)
- Codon counts from ORFs, frequencies, RSCU, log2 fold-change between two groups
- GC3 per-ORF distribution
- Streamlit UI with tabs, caching, and responsive tables/plots
- Small utility functions and defensive coding

Author: ChatGPT (GPT-5 Thinking mini)
Date: 2025-11
"""

from io import StringIO, BytesIO
import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import textwrap
import matplotlib
matplotlib.use("Agg")

# ---------------------------
# Constants / Genetic code
# ---------------------------
START_CODONS = ("ATG",)  # default; UI can expose alternatives later
STOP_CODONS = ("TAA", "TAG", "TGA")

# canonical 64 codons in TCAG order (keeps heatmap layout deterministic)
ALL_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]

# Basic genetic code mapping (DNA codon -> AA single letter)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# ---------------------------
# Helper functions
# ---------------------------

def safe_seq_upper(seq):
    """Return uppercase DNA sequence with U->T and only standard letters kept."""
    return str(seq).upper().replace('U', 'T')

def translate_nt_seq(nt_seq):
    """Translate a DNA sequence (Biopython Seq) to AA sequence without trailing '*' (use to_stop=True)."""
    # Using Biopython Seq.translate with to_stop=True to avoid terminal '*' in result
    s = Seq(nt_seq)
    try:
        aa = s.translate(table=1, to_stop=True)  # table=1 standard
        return str(aa)
    except Exception:
        # fallback: naive codon translation using GENETIC_CODE (not including termination)
        aa_chars = []
        for i in range(0, (len(nt_seq)//3)*3, 3):
            cod = nt_seq[i:i+3]
            aa_chars.append(GENETIC_CODE.get(cod, 'X'))
        # trim at first '*' if present
        aa_str = ''.join(aa_chars)
        if '*' in aa_str:
            aa_str = aa_str.split('*')[0]
        return aa_str

# ---------------------------
# ORF finder (robust)
# ---------------------------

def find_orfs_in_seq(seq_str, min_length_nt=150, starts=START_CODONS, stops=STOP_CODONS):
    """
    Scan sequence in all 6 frames and return ORFs list.
    Each ORF: dict with keys: frame (1..6), strand (1 or -1), orf_start, orf_end, length_nt, length_aa, aa_seq
    Coordinates: 1-based inclusive relative to forward strand.
    Robust scanning: if a start has no stop, continue scanning frame (do not break frame early).
    """
    seq = safe_seq_upper(seq_str)
    seq_len = len(seq)
    results = []

    # Forward strand frames (offset 0,1,2)
    for f in range(3):
        i = f
        while i + 3 <= seq_len:
            cod = seq[i:i+3]
            if cod in starts:
                # search for next stop codon downstream
                j = i + 3
                found_stop = False
                while j + 3 <= seq_len:
                    codj = seq[j:j+3]
                    if codj in stops:
                        found_stop = True
                        orf_start = i + 1
                        orf_end = j + 3
                        length_nt = orf_end - orf_start + 1
                        if length_nt >= min_length_nt:
                            aa_seq = translate_nt_seq(seq[i:orf_end])
                            results.append({
                                'frame': f+1,
                                'strand': 1,
                                'orf_start': orf_start,
                                'orf_end': orf_end,
                                'length_nt': length_nt,
                                'length_aa': len(aa_seq),
                                'aa_seq': aa_seq
                            })
                        # advance scan to just after this stop to allow downstream ORFs
                        i = j + 3
                        break
                    else:
                        j += 3
                if not found_stop:
                    # No stop found for this start; move one codon forward and keep scanning frame
                    i += 3
            else:
                i += 3

    # Reverse complement frames
    rc_seq = str(Seq(seq).reverse_complement())
    for f in range(3):
        i = f
        while i + 3 <= len(rc_seq):
            cod = rc_seq[i:i+3]
            if cod in starts:
                j = i + 3
                found_stop = False
                while j + 3 <= len(rc_seq):
                    codj = rc_seq[j:j+3]
                    if codj in stops:
                        found_stop = True
                        rc_start = i + 1
                        rc_end = j + 3
                        # map rc coords to forward coords (1-based)
                        fwd_end = seq_len - rc_start + 1
                        fwd_start = seq_len - rc_end + 1
                        length_nt = fwd_end - fwd_start + 1
                        if length_nt >= min_length_nt:
                            aa_seq = translate_nt_seq(rc_seq[i:rc_end])
                            # frame numbering: 4,5,6 for reverse frames (mapped)
                            results.append({
                                'frame': 3 + f + 1,
                                'strand': -1,
                                'orf_start': fwd_start,
                                'orf_end': fwd_end,
                                'length_nt': length_nt,
                                'length_aa': len(aa_seq),
                                'aa_seq': aa_seq
                            })
                        i = j + 3
                        break
                    else:
                        j += 3
                if not found_stop:
                    i += 3
            else:
                i += 3

    # Sort ORFs by start coordinate (useful)
    results_sorted = sorted(results, key=lambda r: (r['orf_start'], -r['length_nt']))
    return results_sorted

# ---------------------------
# Conversion to DataFrame
# ---------------------------

def orfs_to_dataframe(orfs, seq_id):
    """Convert ORF dict list to pandas DataFrame"""
    if not orfs:
        cols = ['seq_id', 'frame', 'strand_sym', 'orf_start', 'orf_end', 'length_nt', 'length_aa', 'aa_seq']
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(orfs)
    df['seq_id'] = seq_id
    df['strand_sym'] = df['strand'].map({1: '+', -1: '-'})
    # reorder columns
    df = df[['seq_id', 'frame', 'strand_sym', 'orf_start', 'orf_end', 'length_nt', 'length_aa', 'aa_seq']]
    return df

# ---------------------------
# Window density computation
# ---------------------------

def compute_window_density(orfs, seq_len, window_size=500, step=None):
    """Return windows list and counts matrix (6 x n_windows)"""
    if step is None:
        step = window_size
    starts = list(range(1, seq_len + 1, step))
    windows = [(s, min(s + window_size - 1, seq_len)) for s in starts]
    counts = np.zeros((6, len(windows)), dtype=int)
    for orf in orfs:
        orf_s = orf.get('orf_start', 0)
        orf_e = orf.get('orf_end', 0)
        frame_idx = int(orf.get('frame', 1)) - 1
        for w_idx, (ws, we) in enumerate(windows):
            if not (orf_e < ws or orf_s > we):
                counts[frame_idx, w_idx] += 1
    return windows, counts

# ---------------------------
# Codon counts / RSCU / stats
# ---------------------------

def codon_counts_from_orfs(seq_str, orfs):
    """Count codons for a sequence string across given ORF coordinate list (1-based inclusive)."""
    seq = safe_seq_upper(seq_str)
    counts = {c: 0 for c in ALL_CODONS}
    for orf in orfs:
        s = max(0, orf.get('orf_start', 1) - 1)
        e = min(len(seq), orf.get('orf_end', len(seq)))
        seg = seq[s:e]
        L = (len(seg) // 3) * 3
        seg = seg[:L]
        for i in range(0, L, 3):
            cod = seg[i:i+3]
            if len(cod) == 3 and all(ch in "TCAG" for ch in cod):
                counts[cod] += 1
    return counts

def counts_to_dataframe(counts_dict):
    """Return DataFrame with count, aa, freq, per1000."""
    df = pd.DataFrame.from_dict(counts_dict, orient='index', columns=['count'])
    df.index.name = 'codon'
    df = df.reindex(ALL_CODONS, fill_value=0)
    df['aa'] = df.index.map(lambda c: GENETIC_CODE.get(c, 'X'))
    total = df['count'].sum()
    df['freq'] = df['count'] / total if total > 0 else 0.0
    df['per1000'] = df['freq'] * 1000.0
    return df

def compute_rscu(df_counts):
    """
    Compute RSCU (Relative Synonymous Codon Usage).
    For each amino acid, RSCU = observed_count / (expected_count if all synonymous codons used equally).
    """
    res = []
    for aa, group in df_counts.groupby('aa'):
        if aa == '*':
            res.append(group.assign(RSCU=np.nan))
            continue
        total = group['count'].sum()
        ncod = len(group)
        if total == 0:
            res.append(group.assign(RSCU=np.nan))
            continue
        expected = total / ncod
        rscu = group['count'] / expected
        res.append(group.assign(RSCU=rscu))
    return pd.concat(res).sort_index()

def log2_fold_change(dfA, dfB, pseudocount=1e-8):
    """Return Series of log2(freqA / freqB) in codon index order."""
    fa = dfA['freq'] + pseudocount
    fb = dfB['freq'] + pseudocount
    return np.log2(fa / fb)

def chi_square_test(dfA, dfB):
    """Global chi-square test across codons (2 x 64 table). Returns chi2, p"""
    a = dfA['count'].values
    b = dfB['count'].values
    table = np.vstack([a, b])
    try:
        chi2, p, dof, exp = chi2_contingency(table)
        return chi2, p
    except Exception:
        return None, None

# ---------------------------
# GC3 computation
# ---------------------------

def gc3_from_orfs(seq_str, orfs):
    """Compute GC3 per ORF; returns DataFrame with orf_start, orf_end, gc3, length."""
    seq = safe_seq_upper(seq_str)
    rows = []
    for orf in orfs:
        s = max(0, orf.get('orf_start', 1) - 1)
        e = min(len(seq), orf.get('orf_end', len(seq)))
        seg = seq[s:e]
        L = (len(seg) // 3) * 3
        seg = seg[:L]
        if L < 3:
            rows.append({'orf_start': s + 1, 'orf_end': e, 'gc3': np.nan, 'length': len(seg)})
            continue
        third_pos = seg[2::3]
        gc3 = sum(1 for ch in third_pos if ch in 'GC') / len(third_pos)
        rows.append({'orf_start': s + 1, 'orf_end': e, 'gc3': gc3, 'length': len(seg)})
    return pd.DataFrame(rows)

# ---------------------------
# Plotting helpers
# ---------------------------

def plot_length_histogram(df, bins=50):
    fig, ax = plt.subplots(figsize=(8, 4))
    if df.empty:
        ax.text(0.5, 0.5, 'No ORFs found', ha='center', va='center')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        return fig
    sns.histplot(df['length_aa'], bins=bins, kde=True, ax=ax)
    ax.set_xlabel('ORF length (aa)'); ax.set_ylabel('Count'); ax.set_title('ORF length distribution (aa)')
    plt.tight_layout()
    return fig

def plot_heatmap(windows, counts):
    fig, ax = plt.subplots(figsize=(12, 4))
    if counts.sum() == 0:
        ax.text(0.5, 0.5, 'No ORFs found', ha='center', va='center')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        return fig
    sns.heatmap(counts, cmap='viridis', ax=ax, cbar_kws={'label': 'ORF count'})
    ax.set_yticks(np.arange(6) + 0.5)
    ax.set_yticklabels(['Frame 1 (+)', 'Frame 2 (+)', 'Frame 3 (+)', 'Frame 4 (-)', 'Frame 5 (-)', 'Frame 6 (-)'])
    n_windows = len(windows)
    # create readable x-tick labels (s-e)
    xlabels = [f"{s}-{e}" for (s, e) in windows]
    if n_windows > 30:
        step = max(1, n_windows // 20)
        xlabels = [lbl if (i % step == 0) else '' for i, lbl in enumerate(xlabels)]
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    ax.set_xlabel('Genomic window (bp)')
    ax.set_title('ORF density heatmap (counts per window by frame)')
    plt.tight_layout()
    return fig

def plot_codon_log2fc_heatmap(log2fc_series):
    """Arrange 64 codons as 8x8 grid following ALL_CODONS order for visualization."""
    # Clip to reasonable range for color scale
    arr = np.array(log2fc_series).reshape(8, 8)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(arr, cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=ax, cbar_kws={'label': 'log2(freq A / freq B)'})
    ax.set_title('Codon usage log2 fold-change (A vs B)')
    # annotate codons in grid
    codon_grid = np.array(ALL_CODONS).reshape(8, 8)
    for (i, j), cod in np.ndenumerate(codon_grid):
        ax.text(j + 0.5, i + 0.5, cod, ha='center', va='center', fontsize=8, color='black')
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    return fig

def plot_gc3_distribution(gc3_df, title='GC3 distribution'):
    fig, ax = plt.subplots(figsize=(8, 3))
    if gc3_df.empty or gc3_df['gc3'].isna().all():
        ax.text(0.5, 0.5, 'No GC3 data available', ha='center', va='center')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        return fig
    sns.histplot(gc3_df['gc3'].dropna(), bins=30, kde=True, ax=ax)
    ax.set_xlabel('GC3'); ax.set_ylabel('Count'); ax.set_title(title)
    plt.tight_layout()
    return fig

# ---------------------------
# Caching heavy functions
# ---------------------------

@st.cache_data(show_spinner=False)
def cached_find_orfs(seq_str, min_length_nt):
    return find_orfs_in_seq(seq_str, min_length_nt=min_length_nt)

@st.cache_data(show_spinner=False)
def cached_codon_counts(seq_str, orfs):
    return codon_counts_from_orfs(seq_str, orfs)

# ---------------------------
# Streamlit UI
# ---------------------------

st.set_page_config(page_title='Genomic Signature Explorer', layout='wide', page_icon='🧬')

st.title('Genomic Signature Explorer')
st.markdown('Interactive ORF analysis, codon bias visualizer, and GC3 explorer. Upload FASTA(s) or paste sequence.')

tab_orf, tab_codon = st.tabs(['📊 ORF Analysis', '🔬 Codon & GC3'])

# ---------------------------
# ORF Analysis tab
# ---------------------------
with tab_orf:
    st.header('ORF Finder & Density Heatmap')
    with st.expander('How it works'):
        st.markdown(textwrap.dedent("""\
            - Scans both strands in 3 frames each (6-frame ORF finding).
            - Finds start-to-stop ORFs (configurable minimum length).
            - Visualizes length distribution and density heatmap (counts per window).
            - You can download ORF tables and plot images.
            """))

    col_cfg, col_input = st.columns([1, 2])
    with col_cfg:
        min_orf_nt = st.number_input('Minimum ORF length (nt)', min_value=30, max_value=100000, value=150, step=3)
        window_size = st.number_input('Window size (bp) for heatmap', min_value=50, max_value=100000, value=500, step=50)
        step_size = st.number_input('Window step size (bp)', min_value=1, max_value=100000, value=window_size, step=50)
        frames_select = st.multiselect('Frames to include', options=['1','2','3','4','5','6'], default=['1','2','3','4','5','6'])
    with col_input:
        uploaded = st.file_uploader('Upload FASTA file(s)', type=['fa','fasta','fas','txt','fna'], accept_multiple_files=True, key='orf_files')
        paste = st.text_area('Or paste FASTA (single or multiple records)', height=120)

    # Parse sequences
    records = []
    if uploaded:
        for f in uploaded:
            try:
                txt = f.read().decode('utf-8')
            except Exception:
                f.seek(0)
                txt = f.read().decode('latin-1')
            try:
                stream = StringIO(txt)
                for rec in SeqIO.parse(stream, 'fasta'):
                    records.append(rec)
            except Exception as e:
                st.error(f'Error parsing {f.name}: {str(e)}')

    if paste and paste.strip():
        try:
            stream = StringIO(paste)
            for rec in SeqIO.parse(stream, 'fasta'):
                records.append(rec)
        except Exception as e:
            st.error('Failed to parse pasted FASTA: ' + str(e))

    if not records:
        st.info('Upload or paste FASTA to run ORF analysis')
    else:
        # run ORF finding per record (cached)
        all_orf_dfs = []
        seq_summaries = []
        for rec in records:
            orfs = cached_find_orfs(str(rec.seq), min_orf_nt)
            df_orfs = orfs_to_dataframe(orfs, rec.id)
            if not df_orfs.empty:
                all_orf_dfs.append(df_orfs)
            seq_summaries.append({'Sequence ID': rec.id, 'Length (bp)': len(rec.seq), 'ORFs found': len(orfs)})

        if all_orf_dfs:
            df_all = pd.concat(all_orf_dfs, ignore_index=True)
        else:
            df_all = pd.DataFrame(columns=['seq_id','frame','strand_sym','orf_start','orf_end','length_nt','length_aa','aa_seq'])

        # filter by frames selected
        if frames_select:
            frames_int = [int(x) for x in frames_select]
            df_all = df_all[df_all['frame'].astype(int).isin(frames_int)]

        # Sequence summary table
        st.subheader('Sequence summary')
        st.dataframe(pd.DataFrame(seq_summaries), use_container_width=True)

        st.subheader('Detected ORFs')
        if not df_all.empty:
            st.dataframe(df_all, use_container_width=True)
            csv_bytes = df_all.to_csv(index=False).encode('utf-8')
            st.download_button('Download ORF table (CSV)', data=csv_bytes, file_name='orfs_table.csv', mime='text/csv')

            # ORF length distribution
            st.subheader('ORF length distribution')
            fig_len = plot_length_histogram(df_all)
            st.pyplot(fig_len)
            buf_len = BytesIO(); fig_len.savefig(buf_len, format='png', dpi=150, bbox_inches='tight'); buf_len.seek(0)
            st.download_button('Download length distribution (PNG)', data=buf_len, file_name='orf_length_distribution.png', mime='image/png')

            # Heatmap per sequence
            st.subheader('ORF density heatmap (by sequence)')
            seq_choice = st.selectbox('Choose sequence', options=[r.id for r in records])
            seq_rec = next(r for r in records if r.id == seq_choice)
            orfs_seq = df_all[df_all['seq_id'] == seq_choice].to_dict('records')
            if orfs_seq:
                windows, counts = compute_window_density(orfs_seq, len(seq_rec.seq), window_size=window_size, step=step_size)
                fig_hm = plot_heatmap(windows, counts)
                st.pyplot(fig_hm)
                buf_hm = BytesIO(); fig_hm.savefig(buf_hm, format='png', dpi=200, bbox_inches='tight'); buf_hm.seek(0)
                st.download_button('Download heatmap (PNG)', data=buf_hm, file_name='orf_density_heatmap.png', mime='image/png')

                # small metrics
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric('Total ORFs', len(orfs_seq))
                with col2:
                    st.metric('Avg length (aa)', f"{df_all[df_all['seq_id'] == seq_choice]['length_aa'].mean():.1f}")
                with col3:
                    st.metric('Max length (aa)', int(df_all[df_all['seq_id'] == seq_choice]['length_aa'].max()))
            else:
                st.warning('No ORFs in selected sequence with current filters')
        else:
            st.warning('No ORFs found with current settings. Try lowering minimum ORF length.')

# ---------------------------
# Codon & GC3 tab
# ---------------------------
with tab_codon:
    st.header('Codon usage & GC3 analysis')
    st.markdown('Compute codon frequencies from detected ORFs and compare two groups.')

    col_inputs, col_opts = st.columns([2, 1])
    with col_inputs:
        uploaded_A = st.file_uploader('Upload FASTA for Group A (organism A / set A)', type=['fa','fasta','fas','txt','fna'], key='codonA')
        uploaded_B = st.file_uploader('Upload FASTA for Group B (optional)', type=['fa','fasta','fas','txt','fna'], key='codonB')
        paste_A = st.text_area('Or paste FASTA for Group A (optional)', height=80)
        paste_B = st.text_area('Paste FASTA for Group B (optional)', height=80)
    with col_opts:
        use_orfs_for_counts = st.checkbox('Count codons from ORFs (recommended)', value=True)
        min_orf_for_codon = st.number_input('Min ORF length (nt) for codon counts', min_value=30, max_value=10000, value=150, step=3)
        show_rscu = st.checkbox('Show RSCU (Relative Synonymous Codon Usage)', value=True)

    if st.button('Run codon analysis'):
        # parse A
        recsA = []
        if uploaded_A:
            try:
                txt = uploaded_A.read().decode('utf-8')
            except Exception:
                uploaded_A.seek(0)
                txt = uploaded_A.read().decode('latin-1')
            stream = StringIO(txt)
            for rec in SeqIO.parse(stream, 'fasta'):
                recsA.append(rec)
        if paste_A and paste_A.strip():
            stream = StringIO(paste_A)
            for rec in SeqIO.parse(stream, 'fasta'):
                recsA.append(rec)

        # parse B
        recsB = []
        if uploaded_B:
            try:
                txt = uploaded_B.read().decode('utf-8')
            except Exception:
                uploaded_B.seek(0)
                txt = uploaded_B.read().decode('latin-1')
            stream = StringIO(txt)
            for rec in SeqIO.parse(stream, 'fasta'):
                recsB.append(rec)
        if paste_B and paste_B.strip():
            stream = StringIO(paste_B)
            for rec in SeqIO.parse(stream, 'fasta'):
                recsB.append(rec)

        if not recsA:
            st.error('Provide at least Group A sequences (upload or paste).')
            st.stop()

        # For each record in A, detect ORFs (cached) and count codons
        countsA = {c: 0 for c in ALL_CODONS}
        all_orfs_A = []
        for rec in recsA:
            orfs = cached_find_orfs(str(rec.seq), min_orf_for_codon)
            # optional: filter ORFs by frame selection? here we keep all
            all_orfs_A.extend([o for o in orfs])
            if use_orfs_for_counts:
                cts = cached_codon_counts(str(rec.seq), orfs)
            else:
                # count from whole sequence (frame 1)
                fake_orf = [{'orf_start':1, 'orf_end':len(rec.seq)}]
                cts = cached_codon_counts(str(rec.seq), fake_orf)
            for k, v in cts.items():
                countsA[k] += v

        dfA = counts_to_dataframe(countsA)
        st.subheader('Group A codon usage')
        display_cols = ['count', 'per1000', 'aa']
        st.dataframe(dfA[display_cols], use_container_width=True)
        st.download_button('Download Group A codon table (CSV)', data=dfA.to_csv().encode('utf-8'), file_name='groupA_codon_table.csv', mime='text/csv')

        if show_rscu:
            rA = compute_rscu(dfA)
            st.subheader('Group A RSCU')
            st.dataframe(rA[['count', 'RSCU', 'aa']], use_container_width=True)

        # GC3 for Group A
        if all_orfs_A:
            gc3_A = gc3_from_orfs(''.join(str(r.seq) for r in recsA) if len(recsA)>1 else str(recsA[0].seq), all_orfs_A)
            if not gc3_A.empty:
                st.subheader('Group A GC3 distribution (per ORF)')
                fig_gc3_A = plot_gc3_distribution(gc3_A, 'GC3 - Group A')
                st.pyplot(fig_gc3_A)
                st.download_button('Download Group A GC3 (PNG)', data=(BytesIO(fig_gc3_A.savefig(BytesIO(), format='png'),) if False else b''), file_name='groupA_gc3.png')  # placeholder no-op; below we also provide proper buffer

                # proper buffer and metrics
                buf_gc3A = BytesIO(); fig_gc3_A.savefig(buf_gc3A, format='png', dpi=150, bbox_inches='tight'); buf_gc3A.seek(0)
                st.download_button('Download Group A GC3 (PNG)', data=buf_gc3A, file_name='groupA_gc3.png', mime='image/png')
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric('Mean GC3', f"{gc3_A['gc3'].mean():.3f}")
                with col2:
                    st.metric('Median GC3', f"{gc3_A['gc3'].median():.3f}")
                with col3:
                    st.metric('Std Dev', f"{gc3_A['gc3'].std():.3f}")
        else:
            st.info('No ORFs found in Group A to compute GC3.')

        # If Group B provided - compute counts and comparison
        if recsB:
            countsB = {c: 0 for c in ALL_CODONS}
            all_orfs_B = []
            for rec in recsB:
                orfs = cached_find_orfs(str(rec.seq), min_orf_for_codon)
                all_orfs_B.extend([o for o in orfs])
                if use_orfs_for_counts:
                    cts = cached_codon_counts(str(rec.seq), orfs)
                else:
                    fake_orf = [{'orf_start':1, 'orf_end':len(rec.seq)}]
                    cts = cached_codon_counts(str(rec.seq), fake_orf)
                for k, v in cts.items():
                    countsB[k] += v

            dfB = counts_to_dataframe(countsB)
            st.subheader('Group B codon usage')
            st.dataframe(dfB[display_cols], use_container_width=True)
            st.download_button('Download Group B codon table (CSV)', data=dfB.to_csv().encode('utf-8'), file_name='groupB_codon_table.csv', mime='text/csv')

            if show_rscu:
                rB = compute_rscu(dfB)
                st.subheader('Group B RSCU')
                st.dataframe(rB[['count', 'RSCU', 'aa']], use_container_width=True)

            # Comparison heatmap (log2 FC)
            st.subheader('Codon usage comparison (log2 fold-change)')
            l2 = log2_fold_change(dfA, dfB)
            fig_l2 = plot_codon_log2fc_heatmap(l2.values)
            st.pyplot(fig_l2)
            buf_l2 = BytesIO(); fig_l2.savefig(buf_l2, format='png', dpi=200, bbox_inches='tight'); buf_l2.seek(0)
            st.download_button('Download codon log2FC heatmap (PNG)', data=buf_l2, file_name='codon_log2fc.png', mime='image/png')

            # chi-square
            chi2, p = chi_square_test(dfA, dfB)
            if chi2 is not None:
                st.markdown(f"**Chi-square test**: χ² = {chi2:.2f}, p = {p:.3e}")
                if p < 0.001:
                    st.success('Codon usage patterns differ significantly between groups (p < 0.001)')
                else:
                    st.info('No strong evidence of overall codon usage difference (p >= 0.001)')

            # GC3 for B
            if all_orfs_B:
                gc3_B = gc3_from_orfs(''.join(str(r.seq) for r in recsB) if len(recsB)>1 else str(recsB[0].seq), all_orfs_B)
                if not gc3_B.empty:
                    st.subheader('Group B GC3 distribution (per ORF)')
                    fig_gc3_B = plot_gc3_distribution(gc3_B, 'GC3 - Group B')
                    st.pyplot(fig_gc3_B)
                    buf_gc3B = BytesIO(); fig_gc3_B.savefig(buf_gc3B, format='png', dpi=150, bbox_inches='tight'); buf_gc3B.seek(0)
                    st.download_button('Download Group B GC3 (PNG)', data=buf_gc3B, file_name='groupB_gc3.png', mime='image/png')

        st.success('Codon analysis complete.')

# Footer
st.markdown('---')
st.caption('Genomic Signature Explorer — ORF & Codon visual analysis (Streamlit + Biopython)')

