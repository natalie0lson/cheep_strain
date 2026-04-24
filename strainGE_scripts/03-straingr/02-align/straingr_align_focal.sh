#!/bin/bash
#SBATCH --job-name="straingr_align"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same as the number of cores
#SBATCH --cpus-per-task=2  # Ensure this matches the `-t` value in bwa mem
#SBATCH --partition=largemem
#SBATCH --time=96:00:00
#SBATCH --output=straingr_align-%j.o
#SBATCH --error=straingr_align-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

# Load the conda environment
source /projects/mnadimlab/conda/bin/activate
conda init bash
source ~/.bashrc  # Ensure the changes made by `conda init` take effect
conda activate /projects/mnadimlab/nmolson/condaenvs/bwa-mem2
 
# Define the directory containing the sample files and the reference genome
SAMPLES_DIR="/projects/mnadimlab/nmolson/data/ecoli_FOCAL/rawdata"
OUTPUT_DIR="/projects/mnadimlab/nmolson/data/straingr/focal2_align2"
REFERENCE="/projects/mnadimlab/nmolson/data/straingr/refs_concat4.fasta"

# Loop through all unique sample prefixes in the directory
for sample in $(ls $SAMPLES_DIR/*_R1.fq.gz | sed 's/_R1.fq.gz//' | xargs -n 1 basename); do
    # Define the input and output files
    READ1="$SAMPLES_DIR/${sample}_R1.fq.gz"
    READ2="$SAMPLES_DIR/${sample}_R2.fq.gz"
    BAM_OUT="$OUTPUT_DIR/${sample}.bam"

    # Run bwa mem and samtools sort
    bwa mem -I 300 -t 2 $REFERENCE $READ1 $READ2 | samtools sort -@ 2 -O BAM -o $BAM_OUT -

    # Create BAM index
    samtools index $BAM_OUT
done

# Deactivate the conda environment
conda deactivate