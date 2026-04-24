#!/bin/bash
#SBATCH --job-name="strain_gst_focal"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=week-long-cpu
#SBATCH --time=96:00:00
#SBATCH --output=straingst-%j.o
#SBATCH --error=straingst-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

conda init bash

conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

#!/bin/bash

# MapSan Directories
output_dir="/projects/mnadimlab/nmolson/data/straingst/focal2"
dir="/projects/mnadimlab/nmolson/data/kmers/focal_kmers"
cd /projects/mnadimlab/nmolson/data/kmers/focal_kmers
# Loop through all files in the directory
for file in "$dir"/*; do
    # Extract filename without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Run straingst command
    straingst run -o "${output_dir}/${filename_no_ext}.tsv" /projects/mnadimlab/mschwab/pan-genome-db.hdf5 "${filename_no_ext}.hdf5"
done

conda deactivate
