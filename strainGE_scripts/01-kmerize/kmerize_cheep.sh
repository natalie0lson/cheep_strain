#!/bin/bash
#SBATCH --job-name="create_sample_kmers"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=week-long-cpu
#SBATCH --time=32:00:00
#SBATCH --output=strainge-%j.o
#SBATCH --error=strainge-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

conda init bash

conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

# Loop through every subdirectory in the specified directory
for dir in /projects/mnadimlab/mnadim/cheep/data/usftp21.novogene.com/01.RawData/EUP{215,216,217}*/; do 

   # Navigate to the current subdirectory
    cd "$dir" || continue

    # Loop through every .fq.gz file in the current subdirectory
    for f in *.fq.gz; do
        output_file="/projects/mnadimlab/nmolson/data/cheep2_kmers/$(basename ${f%.fq.gz}.hdf5)"  # Define the output file path in the new directory
        straingst kmerize -o "$output_file" "$f";  # Execute the kmerize command with updated output path
    done;

    # Navigate back to the parent directory
    cd - >/dev/null || exit
done;


# Deactivate the conda environment
conda deactivate


