#!/bin/bash
#SBATCH --job-name="strain_gst_output"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same that number of cores
#SBATCH --cpus-per-task=1
#SBATCH --partition=day-long-cpu
#SBATCH --time=06:00:00
#SBATCH --output=straingst-%j.o
#SBATCH --error=straingst-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

# Load required modules
# module load python/3.8

cd /projects/mnadimlab/nmolson/data/straingst/focal2

# Set directory for StrainGST
STRAINGST_DIR="/projects/mnadimlab/nmolson/data/straingst/focal2/"

# Directory to save the output DataFrame
OUTPUT_DIR="/projects/mnadimlab/nmolson/data/straingst/output"

# Activate your Python environment
source /projects/mnadimlab/conda/bin/activate
conda activate /projects/mnadimlab/nmolson/condaenvs/straingst_env

# Run Python script
python << END
import pandas as pd
from pathlib import Path
import os

STRAINGST_DIR = Path("$STRAINGST_DIR".strip('"'))
OUTPUT_DIR = Path("$OUTPUT_DIR".strip('"'))

df_list = []
sample_names = []
for f in STRAINGST_DIR.glob("*.tsv"):
    sample_name = f.stem
    df = pd.read_csv(f, sep='\t', comment='#', skiprows=2, index_col=1)
    df_list.append(df)
    sample_names.append(sample_name)

# Combine all StrainGST results from each sample into a single DataFrame.
straingst_df = pd.concat(df_list, keys=sample_names, names=["sample"])

# Sort sample names
straingst_df = straingst_df.sort_index()

# Save the DataFrame to a specified directory
output_file = OUTPUT_DIR / "straingst_focal_results.tsv"
straingst_df.to_csv(output_file, sep='\t')
END
