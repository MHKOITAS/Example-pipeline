Host login.anunna.wur.nl
  HostName login.anunna.wur.nl
  User katra001
   ##### Step 1  1. Load the module (adjust version if needed)
module load 2024 # inside "mamba activate repeatmasker_env"
module load repeatmasker/4.0.7
# 2. Mask the Reference Genome (M. aguilicaudatus)
RepeatMasker aguillicaudatus.fa -species vertebrates -pa 16 -dir ./masked_ref -e rmblast 
# To only mask the first chromosome, use the following command instead:
####samtools faidx aguillicaudatus.fa 
####samtools faidx aguillicaudatus.fa NC_073337.2 > chr1_ref.fa
####RepeatMasker chr1_ref.fa -species vertebrates -pa 16 -dir ./masked_ref_chr1 -e rmblast
# 3. Mask the Query Genome (Cyprinus carpio)
RepeatMasker common_carp.fa -species vertebrates -pa 16 -dir ./masked_query -e rmblast
# To only mask the first chromosome, use the following command instead:
####samtools faidx common_carp.fa 
####samtools faidx common_carp.fa NC_031697.1 > chr1_query.fa
####RepeatMasker chr1_query.fa -species vertebrates -pa 16 -dir ./masked_query_chr1 -e rmblast

######## step 2 How to align the masked genomes using minimap2
module load legacy
module load WUR/RIKILT/minimap2/2.17

##repeatmodller
# 1) PAF output (default/minimap2 lightweight)
 minimap2 -x asm5 -t 16 masked_ref/aguillicaudatus.masked.fa masked_query/common_carp.masked.fa > ag_vs_cc.paf

# 2) SAM -> sorted BAM (for inspection/ downstream tools)
minimap2 -x asm5 -t 16 -a aguillicaudatus.fa common_carp.fa | samtools sort -@ 8 -o ag_vs_cc.sorted.bam
samtools index ag_vs_cc.sorted.bam

####### step 3 Chain and Net Steps
# Define file names
REF_FA="M_aguilicaudatus.fa.masked"
QUERY_FA="Cyprinus_carpio.fa.masked"
INPUT_AXT="ag_vs_cc.axt"
OUTPUT_CHAIN="ag_vs_cc.chain"

# --- 1. axtChain: Chain Local Alignments ---
# This step groups the small, local AXT alignments into larger, ordered blocks (chains).
echo "Starting axtChain..."
axtChain -minScore=5000 -linearGap=loose \
    $INPUT_AXT $REF_FA $QUERY_FA \
    $OUTPUT_CHAIN

# --- 2. chainNet: Net the Chains ---
# This step resolves overlaps between chains and builds a hierarchical alignment structure (net).
# This is where inversions and translocations are structurally defined.
echo "Starting chainNet..."
chainNet $OUTPUT_CHAIN $REF_FA $QUERY_FA \
    ag_vs_cc.net cc_vs_ag.net

echo "Chain and Net steps complete."

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