######## step 4 How to visualize the alignment using dotplotly in wga_env
import pandas as pd
import plotly.express as px

def create_wga_dotplot(paf_file, output_html):
    """
    Create an interactive whole-genome dotplot using a minimap2 PAF file.
    Optimized for speed, large datasets, and correct chromosome ordering.
    """

    # Expected PAF columns (minimap2 may output 12–>20 columns)
    base_cols = [
        'q_name', 'q_len', 'q_start', 'q_end',
        'strand', 'r_name', 'r_len', 'r_start', 'r_end',
        'n_match', 'n_bp', 'mapq'
    ]
    
    # 1. Read the PAF file
    try:
        df = pd.read_csv(paf_file, sep='\t', header=None, comment='#')
    except FileNotFoundError:
        print(f"File not found: {paf_file}")
        return
    except Exception as e:
        print("Error reading PAF:", e)
        return

    # Assign only the columns present
    df.columns = base_cols[:len(df.columns)] + list(df.columns[len(base_cols):])

    # Filter out very small or noisy alignments (recommended; optional)
    if 'n_bp' in df.columns:
        df = df[df['n_bp'] > 2000]

    # -------------------------------
    # 2. Genome-order contig sorting
    # -------------------------------
    # Order = first appearance in PAF (reflects actual genome)
    ref_order = df['r_name'].drop_duplicates().tolist()
    qry_order = df['q_name'].drop_duplicates().tolist()

    # Build reference contig layout
    ref_len_df = (
        df[['r_name', 'r_len']]
        .drop_duplicates()
        .set_index('r_name')
        .loc[ref_order]                # enforce correct genome order
        .reset_index()
    )
    ref_len_df['r_end_cum'] = ref_len_df['r_len'].cumsum()
    ref_len_df['r_start_cum'] = ref_len_df['r_end_cum'].shift(1, fill_value=0)

    # Build query contig layout
    qry_len_df = (
        df[['q_name', 'q_len']]
        .drop_duplicates()
        .set_index('q_name')
        .loc[qry_order]
        .reset_index()
    )
    qry_len_df['q_end_cum'] = qry_len_df['q_len'].cumsum()
    qry_len_df['q_start_cum'] = qry_len_df['q_end_cum'].shift(1, fill_value=0)

    # Convert to lookup dictionaries
    ref_lookup = dict(zip(ref_len_df['r_name'], ref_len_df['r_start_cum']))
    qry_lookup = dict(zip(qry_len_df['q_name'], qry_len_df['q_start_cum']))

    # -------------------------------------------------------
    # 3. Compute global coordinates (vectorized, fast)
    # -------------------------------------------------------
    df['r_global'] = df['r_start'] + df['r_name'].map(ref_lookup)
    df['q_global'] = df['q_start'] + df['q_name'].map(qry_lookup)

    # Color by strand (synteny vs inversion)
    if 'strand' in df.columns:
        df['orientation'] = df['strand'].map({
            '+': 'Forward (Synteny)',
            '-': 'Reverse (Inversion)'
        }).fillna('Unknown')
    else:
        df['orientation'] = 'Unknown'

    # -------------------------------------------------------
    # 4. Plotly dotplot (WebGL for large datasets)
    # -------------------------------------------------------
    fig = px.scatter(
        df,
        x='r_global',
        y='q_global',
        color='orientation',
        render_mode='webgl',      # critical for speed
        opacity=0.6,
        hover_data={
            'r_name': True,
            'q_name': True,
            'strand': True if 'strand' in df.columns else False,
            'r_global': False,
            'q_global': False
        },
        title="Whole Genome Alignment Dotplot"
    )

    # Axis labels
    fig.update_layout(
        xaxis_title="Reference genome (cumulative position)",
        yaxis_title="Query genome (cumulative position)",
        width=1100,
        height=1100
    )

    # -------------------------------------------------------
    # 5. Save output
    # -------------------------------------------------------
    fig.write_html(output_html)
    print(f"✅ Dotplot written to {output_html}")


# ---------------------------
# Script entry point
# ---------------------------
if __name__ == "__main__":
    INPUT_PAF = "1_vs_2.paf"
    OUTPUT_HTML = "dotplot_1_vs_2.html"
    create_wga_dotplot(INPUT_PAF, OUTPUT_HTML)
