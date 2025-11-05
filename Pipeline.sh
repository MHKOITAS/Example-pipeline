#!/bin/bash
#SBATCH --comment=dotplot_wga
#SBATCH --time=0-03:00:00
#SBATCH --mem=2048 
#SBATCH --cpus-per-task=2
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=wga_dotplot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=apostolos.katrachouras@wur.nl
# Load conda/mamba

#!/bin/bash
source ~/.bashrc
mamba activate wga_env

##### Step 1 Build Data bases for Repeat Modeler
echo "Starting BuildDatabases for custom library generation"

# Build databases
BuildDatabase -name my_my my_my.fa
BuildDatabase -name myc_pne myc_pne.fa

mamba deactivate

mamba activate sam_t

##### Step 2 — Run RepeatModeler to identify repeats
echo "Starting RepeatModeler"
RepeatModeler -database my_my -threads 2 -engine rmblast > remomymy.out
RepeatModeler -database myc_pne -threads 2 -engine rmblast > remomycpne.out

mamba deactivate

mamba activate wga_env

##### Step 3 — Mask genomes using RepeatMasker
echo "Starting RepeatMasker for Reference Genome"
RepeatMasker my_my.fa -lib my_my-families.fa -pa 2 -dir masked_ref -e rmblast

echo "Starting RepeatMasker for Query Genome"
RepeatMasker myc_pne.fa -lib myc_pne-families.fa -pa 2 -dir masked_query -e rmblast

mamba deactivate

mamba activate wga_env

##### Step 4 — Align masked genomes with minimap2
echo "Starting minimap2 alignment"
minimap2 -x asm10 -t 2 masked_ref/my_my.fa.masked masked_query/myc_pne.fa.masked > 1_vs_2.paf

mamba deactivate

mamba activate wga_env

##### Step 5 — Visualization
echo "Starting dotplot visualization"

python dotplot1.py

mamba deactivate

# Optional SAM/BAM output for IGV use
#minimap2 -x asm10 -t 2 -a masked_ref/myc_gen.fa.masked masked_query/myc_pne.fa.masked | samtools sort -@ 2 -o 1_vs_2.sorted.bam
#samtools index 1_vs_2.sorted.bam

