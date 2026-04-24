#!/bin/bash
#SBATCH --job-name="download_mapsan"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=week-long-cpu
#SBATCH --time=32:00:00
#SBATCH --output=ncbi-%j.o
#SBATCH --error=ncbi-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

# Create a directory to store downloaded files
mkdir -p /projects/mnadimlab/nmolson/data/mapsan_fastq

# Change to the specified directory
cd /projects/mnadimlab/nmolson/data/mapsan_fastq

# Loop through the range of IDs and download the files
for i in {2760..2936}; do
    curl -L "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/${i:0:3}/SRR152${i}/SRR152${i}_1.fastq.gz" -o "SRR152${i}_shotgun_metagenomes_of_Maputo_Mozambique_children_gut_microbiome_1.fastq.gz"
    curl -L "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/${i:0:3}/SRR152${i}/SRR152${i}_2.fastq.gz" -o "SRR152${i}_shotgun_metagenomes_of_Maputo_Mozambique_children_gut_microbiome_2.fastq.gz"
done
