#!/bin/bash
#SBATCH --job-name="create_sample_kmers"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=week-long-cpu
#SBATCH --time=96:00:00
#SBATCH --output=strainge-%j.o
#SBATCH --error=strainge-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

conda init bash

conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

cd /projects/mnadimlab/nmolson/data/mapsan_fastq2

for f in /projects/mnadimlab/nmolson/data/mapsan_fastq2/*.fastq.gz; do
    output_file="/projects/mnadimlab/nmolson/data/mapsan_kmers/kmerized$(basename ${f%.fastq.gz}.hdf5)"  # Define the output file path in the new directory
    straingst kmerize -o $output_file $f;  # Execute the kmerize command with updated output path
done;

conda deactivate
