#!/bin/bash
#SBATCH --job-name="straingr_var"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same as the number of cores
#SBATCH --cpus-per-task=2  # Ensure this matches the `-t` value in bwa mem
#SBATCH --partition=week-long-cpu
#SBATCH --time=72:00:00
#SBATCH --output=straingr_var-%j.o
#SBATCH --error=straingr_var-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

# Load the conda environment
source /projects/mnadimlab/conda/bin/activate
conda init bash
source ~/.bashrc  # Ensure the changes made by `conda init` take effect
conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

# Define the directory containing the BAM files and the reference genome
BAM_DIR="/projects/mnadimlab/nmolson/data/straingr/cheep"
REFERENCE="/projects/mnadimlab/nmolson/data/straingr/refs_concat.fasta"
OUTPUT_DIR="/projects/mnadimlab/nmolson/data/straingr/cheep"
# Loop through all BAM files in the directory
for bam_file in $BAM_DIR/*.bam; do
    # Get the base name of the BAM file (without the directory and extension)
    base_name=$(basename "$bam_file" .bam)
    
    # Define the output files
    HDF5_OUT="$OUTPUT_DIR/${base_name}.hdf5"
    SUMMARY_OUT="$OUTPUT_DIR/${base_name}.tsv"
    
    # Check if the output files already exist
    if [ -f "$HDF5_OUT" ] && [ -f "$SUMMARY_OUT" ]; then
        echo "Output files for $base_name already exist, skipping..."
        continue
    fi
    
    # Run the StrainGR variant caller
    straingr call $REFERENCE $bam_file --hdf5-out $HDF5_OUT --summary $SUMMARY_OUT --tracks all
done

# Deactivate the conda environment
conda deactivate

