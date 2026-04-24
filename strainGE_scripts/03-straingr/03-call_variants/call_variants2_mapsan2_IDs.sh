#!/bin/bash
#SBATCH --job-name="id_map2_var"
#SBATCH --nodes=1
#SBATCH --account=nmolson
#SBATCH --tasks-per-node=1  # The same as the number of cores
#SBATCH --cpus-per-task=2  # Ensure this matches the `-t` value in bwa mem
#SBATCH --partition=largemem
#SBATCH --time=96:00:00
#SBATCH --output=id_map2_var-%j.o
#SBATCH --error=id_map2_var-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.olson@emory.edu

# Load the conda environment
source /projects/mnadimlab/conda/bin/activate
conda init bash
source ~/.bashrc  # Ensure the changes made by `conda init` take effect
conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

# Define the directory containing the BAM files, the reference genome, and the output directory
BAM_DIR="/projects/mnadimlab/nmolson/data/straingr/mapsan_align2"
REFERENCE="/projects/mnadimlab/nmolson/data/straingr/refs_concat4.fasta"
OUTPUT_DIR="/projects/mnadimlab/nmolson/data/straingr/call_var2/mapsan"
# Specify the CSV file paths
CSV_FILES=(
    "/projects/mnadimlab/nmolson/data/straingr/compare_ids/df_13C1079T_IDS_MAP_CHEEP.csv"
    "/projects/mnadimlab/nmolson/data/straingr/compare_ids/df_81_IDS_MAPSAN_CHEEP.csv"
    "/projects/mnadimlab/nmolson/data/straingr/compare_ids/df_APEC_102026_IDS_MAPSAN_CHEEP.csv"
    "/projects/mnadimlab/nmolson/data/straingr/compare_ids/df_NCTC9087_IDS_MAP_CHEEP.csv"
)

# Aggregate all IDs from the specified CSV files into a single array
ids=()
for csv_file in "${CSV_FILES[@]}"; do
    echo "Processing IDs from $csv_file..."
    ids+=( $(awk -F',' '{print $1}' "$csv_file" | tail -n +2) )  # Assuming the first column contains IDs and the CSVs have headers
done

# Remove duplicates from the aggregated ID list
unique_ids=($(echo "${ids[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

# Loop through all BAM files in reverse alphanumerical order
for bam_file in $(ls -r $BAM_DIR/*.bam); do
    # Get the base name of the BAM file (without the directory and extension)
    base_name=$(basename "$bam_file" .bam)
    
    # Check if the base name is in the list of IDs
    if [[ ! " ${unique_ids[*]} " =~ " ${base_name} " ]]; then
        echo "Skipping $bam_file: Not in the ID list."
        continue
    fi
    
    # Define the output files
    HDF5_OUT="$OUTPUT_DIR/${base_name}.hdf5"
    SUMMARY_OUT="$OUTPUT_DIR/${base_name}.tsv"
    
    # Check if the output files already exist
    if [[ -f "$HDF5_OUT" && -f "$SUMMARY_OUT" ]]; then
        echo "Skipping $bam_file: Output files already exist."
        continue
    fi
    
    # Run the StrainGR variant caller
    straingr call $REFERENCE $bam_file --hdf5-out $HDF5_OUT --summary $SUMMARY_OUT --tracks all
done

# Deactivate the conda environment
conda deactivate