#!/bin/bash
#SBATCH --job-name="straingr_ref_atl"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=largemem
#SBATCH --time=96:00:00
#SBATCH --output=straingr-ref-atl-%j.o
#SBATCH --error=straingr-ref-atl-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

source /projects/mnadimlab/conda/bin/activate

conda init bash

conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

export PATH=/projects/mnadimlab/nmolson/condaenvs/mummer/bin:$PATH

#!/bin/bash

straingr prepare-ref -s /projects/mnadimlab/nmolson/data/straingst/atl_cheep/*.tsv \
   -p "/projects/mnadimlab/mschwab/strainge_db/{ref}.fa.gz" \
   -S /projects/mnadimlab/mschwab/similarities.tsv \
   -o /projects/mnadimlab/nmolson/data/straingr/refs_concat_atl.fasta

conda deactivate
