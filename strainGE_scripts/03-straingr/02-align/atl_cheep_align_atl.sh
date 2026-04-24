#!/bin/bash
#SBATCH --job-name="atl_align_atl"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same as the number of cores
#SBATCH --cpus-per-task=2  # Ensure this matches the `-t` value in bwa mem
#SBATCH --partition=largemem
#SBATCH --time=96:00:00
#SBATCH --output=atl_align_atl-%j.o
#SBATCH --error=atl_align_atl-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

# Load the conda environment
source /projects/mnadimlab/conda/bin/activate
conda init bash
source ~/.bashrc  # Ensure the changes made by `conda init` take effect
conda activate /projects/mnadimlab/nmolson/condaenvs/bwa-mem2

# Define the directory containing the sample files and the reference genome
SAMPLES_DIR="/projects/mnadimlab/nmolson/data/atl_fastq"
OUTPUT_DIR="/projects/mnadimlab/nmolson/data/straingr/atl_cheep_aligned/atl"
REFERENCE="/projects/mnadimlab/nmolson/data/straingr/refs_concat_atl.fasta"
# Loop through all unique sample prefixes in the directory
for sample in $(ls $SAMPLES_DIR/*_1.fastq.gz | sed 's/_1.fastq.gz//' | xargs -n 1 basename); do
    # Define the input and output files
    READ1="$SAMPLES_DIR/${sample}_1.fastq.gz"
    READ2="$SAMPLES_DIR/${sample}_2.fastq.gz"
    BAM_OUT="$OUTPUT_DIR/${sample}.bam"

    # Run bwa mem and samtools sort
    bwa mem -I 300 -t 2 $REFERENCE $READ1 $READ2 | samtools sort -@ 2 -O BAM -o $BAM_OUT -

    # Create BAM index
    samtools index $BAM_OUT
done

# Deactivate the conda environment
conda deactivate