#!/bin/bash
#SBATCH --job-name=straingr_array
#SBATCH --nodes=1
#SBATCH --account=mnadimp
#SBATCH --tasks-per-node=1  #The same that number of cores
#SBATCH --cpus-per-task=8          # Increased to 8 for faster alignment
#SBATCH --output=log/straingr.out
#SBATCH --error=log/straingr.err
#SBATCH --array=12-70
#SBATCH --partition=largemem
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maya.l.nadimpalli@emory.edu

source /projects/mnadimlab/conda/bin/activate
conda init bash
export PATH=/projects/mnadimlab/nmolson/condaenvs/mummer/bin:$PATH

# Set variables
WORKDIR=/projects/mnadimlab/mnadim/kukusafe/03_output/cef_ecoli/strainGST_ecoli
TRIM_DIR=/projects/mnadimlab/mnadim/kukusafe/03_output/cef_ecoli/trim
SIM_TSV=/projects/mnadimlab/mnadim/kukusafe/02_scripts/02_strainGE/03-straingr/01-prepare_ref/similarities.tsv
REF_DB_PATH="/projects/mnadimlab/mschwab/strainge_db/{ref}"
OUTDIR=/projects/mnadimlab/mnadim/kukusafe/03_output/cef_ecoli/strainGR_by_sample

work_units=/projects/mnadimlab/mnadim/scratch/slice/${SLURM_ARRAY_TASK_ID}

cd $WORKDIR

cat "${work_units}" | while read sample ; do

# 1. Define specific files for this sample
# We assume .strains.tsv is in $WORKDIR
STRAIN_TSV="$WORKDIR/${sample}.strains.tsv"
REF_FASTA="$OUTDIR/${sample}.fasta"
BAM_OUT="$OUTDIR/${sample}.bam"
R1="$TRIM_DIR/${sample}_1.trim.fq.gz"
R2="$TRIM_DIR/${sample}_2.trim.fq.gz"

# 2. Process
conda activate /projects/mnadimlab/nmolson/condaenvs/strainge
straingr prepare-ref -s $STRAIN_TSV -p "$REF_DB_PATH" -S "$SIM_TSV" -o "$REF_FASTA"
conda deactivate

#Alignment
cd $OUTDIR
conda activate /projects/mnadimlab/nmolson/condaenvs/bwa-mem2
bwa index $REF_FASTA

bwa mem -t $SLURM_CPUS_PER_TASK "$REF_FASTA" $R1 $R2 | \
samtools sort -@ 2 -O BAM -o $BAM_OUT -
samtools index $BAM_OUT

conda deactivate

#Call variants
conda activate /projects/mnadimlab/nmolson/condaenvs/strainge
straingr call "$REF_FASTA" $BAM_OUT --hdf5-out $OUTDIR/${sample}.hdf5 \
--summary $OUTDIR/${sample}.tsv --tracks all
rm -f $BAM_OUT ${BAM_OUT}.bai $REF_FASTA ${sample}*wig ${sample}*bed

conda deactivate

done