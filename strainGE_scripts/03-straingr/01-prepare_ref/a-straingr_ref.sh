#!/bin/bash
#SBATCH --job-name="straingr_ref"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=day-long-cpu
#SBATCH --time=20:00:00
#SBATCH --output=straingst-%j.o
#SBATCH --error=straingst-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

conda init bash

conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

#!/bin/bash

straingr prepare-ref -s /projects/mnadimlab/nmolson/data/straingst/all/*.tsv \
   -p "/projects/mnadimlab/mschwab/strainge_db/{ref}.fa.gz" \
   -S /projects/mnadimlab/mschwab/similarities.tsv \
   -o /projects/mnadimlab/nmolson/data/straingr/refs_concat.fasta


conda deactivate
