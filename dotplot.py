 ######## step 4 How to visualize the alignment using dotplotly in wga_env
import pandas as pd
import plotly.express as px

def create_wga_dotplot(paf_file, output_html):
    """
    Reads a Minimap2 PAF file, calculates cumulative coordinates for whole-genome
    visualization, and generates an interactive Plotly scatter dot plot.
    """
    
    # 1. Define PAF column headers
    cols = [
        'q_name', 'q_len', 'q_start', 'q_end', 'strand', 'r_name', 'r_len', 'r_start', 'r_end', 
        'n_match', 'n_bp', 'mapq', 'cigar', 'tp', 'cm', 's1', 's2', 'dv', 'de', 'rl'
    ]
    
    # 2. Read the PAF file
    try:
        df = pd.read_csv(paf_file, sep='\t', header=None, comment='#')
        df.columns = cols[:len(df.columns)]
    except FileNotFoundError:
        print(f"Error: PAF file not found at '{paf_file}'")
        return
    except Exception as e:
        print(f"Error reading PAF file: {e}")
        return

    # Filter out short or low-quality alignments for cleaner visualization (Optional)
    # df = df[df['n_bp'] > 1000] 
    
    # --- 3. Calculate Cumulative Coordinates (CRITICAL STEP for WGA) ---
    
    # Prepare dataframes for calculating offsets
    # Reference Genome Offsets
    ref_len_df = df.drop_duplicates(subset=['r_name']).sort_values(by='r_name')[['r_name', 'r_len']]
    ref_len_df['r_end_cumulative'] = ref_len_df['r_len'].cumsum()
    ref_len_df['r_start_cumulative'] = ref_len_df['r_end_cumulative'].shift(1, fill_value=0)
    ref_len_dict = dict(zip(ref_len_df['r_name'], ref_len_df['r_start_cumulative']))
    
    # Apply offsets to generate global X-coordinates
    df['r_cumulative_start'] = df.apply(lambda row: row['r_start'] + ref_len_dict[row['r_name']], axis=1)

    # Query Genome Offsets
    query_len_df = df.drop_duplicates(subset=['q_name']).sort_values(by='q_name')[['q_name', 'q_len']]
    query_len_df['q_end_cumulative'] = query_len_df['q_len'].cumsum()
    query_len_df['q_start_cumulative'] = query_len_df['q_end_cumulative'].shift(1, fill_value=0)
    query_len_dict = dict(zip(query_len_df['q_name'], query_len_df['q_start_cumulative']))
    
    # Apply offsets to generate global Y-coordinates
    df['q_cumulative_start'] = df.apply(lambda row: row['q_start'] + query_len_dict[row['q_name']], axis=1)
    
    # Define color based on alignment strand (for inversions)
    df['color'] = df['strand'].apply(lambda x: 'Forward (Synteny)' if x == '+' else 'Reverse (Inversion)')

    # --- 4. Create the Plotly Scatter Plot ---
    fig = px.scatter(
        df,
        x='r_cumulative_start',      # Global X-axis coordinate
        y='q_cumulative_start',      # Global Y-axis coordinate
        color='color',               # Color by strand
        hover_data={
            'r_cumulative_start': False, 
            'q_cumulative_start': False, 
            'r_name': True,          # Show Reference Contig on hover
            'q_name': True,          # Show Query Contig on hover
            'strand': True
        },
        title="Whole Genome Alignment Dot Plot: *M. aguillicaudatus* (Ref) vs *C. carpio* (Query)"
    )

    # --- 5. Configure Layout and Save ---
    fig.update_layout(
        xaxis_title="Reference Genome: *M. aguillicaudatus* (Cumulative Position)",
        yaxis_title="Query Genome: *C. carpio* (Cumulative Position)",
        # Set aspect ratio to 1 for better visual comparison of genome sizes
        autosize=False, 
        width=1000, 
        height=1000 
    )

    # Add chromosome/scaffold separation lines (highly recommended but complex to automate perfectly here)
    # Use ref_len_df['r_end_cumulative'] and query_len_df['q_end_cumulative'] to draw vertical/horizontal lines
    
    fig.write_html(output_html)
    print(f"Dot plot saved to: {output_html}")
    print("Open the HTML file in your web browser for an interactive plot.")

# --- EXECUTION ---
if __name__ == "__main__":
    # Define your input and output files
    INPUT_PAF = "ag_vs_cc.paf"
    OUTPUT_HTML = "dotplot_ag_vs_cc.html"
    
    # Make sure to run the Minimap2 alignment command first to create the PAF file:
    # minimap2 -x asm5 -t 16 M_aguillicaudatus.fa C_carpio.fa > ag_vs_cc.paf
    
    create_wga_dotplot(INPUT_PAF, OUTPUT_HTML)