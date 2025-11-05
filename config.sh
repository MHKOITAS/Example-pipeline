#!/bin/bash
#SBATCH --comment=dotplot_wga
#SBATCH --time=0-03:00:00
#SBATCH --mem=2048 
#SBATCH --cpus-per-task=1
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=wga_dotplot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=apostolos.katrachouras@wur.nl
# Load conda/mamba

source /home/WUR/katra001/miniforge3/bin/conda.sh

conda activate wga_env
##### step 0 repeat modeler
echo "Starting RepeatModeler for custom library generation"

RepeatModeler -database aguillicaudatus_db -pa 1 -engine rmblast aguillicaudatus.fa > repeatmodeler.out
##### Step 1  1. Load the module (adjust version if needed)
# 2. Mask the Reference Genome (M. aguilicaudatus)
echo "strarting RepeatMasker for Reference Genome"
RepeatMasker aguillicaudatus.fa -species vertebrates -pa 1 -dir ./masked_ref -e rmblast 
# To only mask the first chromosome, use the following command instead:
###samtools faidx aguillicaudatus.fa 
###samtools faidx aguillicaudatus.fa NC_073337.2 > chr1_ref.fa
###RepeatMasker chr1_ref.fa -species vertebrates -pa 16 -dir ./masked_ref_chr1 -e rmblast
# 3. Mask the Query Genome (Cyprinus carpio)
echo "strarting RepeatMasker for Query Genome"
RepeatMasker common_carp.fa -species vertebrates -pa 1 -dir ./masked_query -e rmblast
# To only mask the first chromosome, use the following command instead:
###samtools faidx common_carp.fa 
###samtools faidx common_carp.fa NC_031697.1 > chr1_query.fa
###RepeatMasker chr1_query.fa -species vertebrates -pa 16 -dir ./masked_query_chr1 -e rmblast

######## step 2 How to align the masked genomes using minimap2
##repeatmodller
# 1) PAF output (default/minimap2 lightweight)
echo "Starting minimap2 alignment"
minimap2 -x asm5 -t 1 masked_ref/aguillicaudatus.masked.fa masked_query/common_carp.masked.fa > ag_vs_cc.paf

# 2) SAM -> sorted BAM (for inspection/ downstream tools)
minimap2 -x asm5 -t 1 -a aguillicaudatus.fa common_carp.fa | samtools sort -@ 1 -o ag_vs_cc.sorted.bam
samtools index ag_vs_cc.sorted.bam

######## step 3 visualization  python
echo "Starting dotplot visualization"
python dotplot1.py

mamba deactivate 

####### step 3 Chain and Net Steps
# Define file names
#REF_FA="M_aguilicaudatus.fa.masked"
#QUERY_FA="Cyprinus_carpio.fa.masked"
#INPUT_AXT="ag_vs_cc.axt"
#OUTPUT_CHAIN="ag_vs_cc.chain"

# --- 1. axtChain: Chain Local Alignments ---
# This step groups the small, local AXT alignments into larger, ordered blocks (chains).
#echo "Starting axtChain..."
#axtChain -minScore=5000 -linearGap=loose \
# $INPUT_AXT $REF_FA $QUERY_FA \
# $OUTPUT_CHAIN

# --- 2. chainNet: Net the Chains ---
# This step resolves overlaps between chains and builds a hierarchical alignment structure (net).
# This is where inversions and translocations are structurally defined.
#echo "Starting chainNet..."
#chainNet $OUTPUT_CHAIN $REF_FA $QUERY_FA \
    ag_vs_cc.net cc_vs_ag.net

#echo "Chain and Net steps complete."

