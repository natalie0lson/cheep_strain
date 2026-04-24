#!/bin/bash
#SBATCH --job-name=straingr_cheep_atl
#SBATCH --nodes=1
#SBATCH --account=mnadimp
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=log/straingr_cheep_atl_%a.out
#SBATCH --error=log/straingr_cheep_atl_%a.err
#SBATCH --array=1-6
#SBATCH --partition=day-long-cpu
#SBATCH --time=0-02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maya.l.nadimpalli@emory.edu

source /projects/mnadimlab/conda/bin/activate
conda init bash
export PATH=/projects/mnadimlab/nmolson/condaenvs/mummer/bin:$PATH

# ── Paths ──────────────────────────────────────────────────────────────────────
PAIRS_TSV=/projects/mnadimlab/nmolson/scripts/strainGE/03-straingr/05-all_commands/cheep_atl_shared_ref_pairs.tsv
CHEEP_STRAINS_TSV=/projects/mnadimlab/nmolson/data/straingst/cheep/cheep_straingst_strains.tsv
ATL_STRAINS_TSV=/projects/mnadimlab/nmolson/data/straingst/atl/atl_straingst_strains.tsv

CHEEP_STRAINS_DIR=/projects/mnadimlab/nmolson/data/straingst/cheep
ATL_STRAINS_DIR=/projects/mnadimlab/nmolson/data/straingst/atl

TRIMDIR1=/projects/mnadimlab/mnadim/cheep/data/output/trimmed        # cheep fastq.gz files
TRIMDIR2=/projects/mnadimlab/nmolson/data/atl_fastq            # atl fastq.gz files

SIM_TSV=/projects/mnadimlab/mschwab/strainge_db/similarities.tsv
REF_DB_PATH="/projects/mnadimlab/mschwab/strainge_db/{ref}"

SIM_TSV=/projects/mnadimlab/mnadim/kukusafe/02_scripts/02_strainGE/03-straingr/01-prepare_ref/similarities.tsv


OUTDIR=/projects/mnadimlab/nmolson/data/straingr/straingr_cheep_atl
TMPDIR=/projects/mnadimlab/nmolson/data/straingr/tmp

#mkdir -p "$OUTDIR" "$TMPDIR"

# ── Get this task's sample pair ────────────────────────────────────────────────
# Build a deduplicated list of unique pairs (skip header, sort unique on cols 1-2)
PAIR=$(awk -F'\t' 'NR>1 {print $1"\t"$2}' "$PAIRS_TSV" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")

CHEEP_SAMPLE=$(echo "$PAIR" | cut -f1)
ATL_SAMPLE=$(echo "$PAIR"   | cut -f2)

echo "Processing pair: cheep=${CHEEP_SAMPLE}  atl=${ATL_SAMPLE}"

# ── Locate per-sample strains TSV files ───────────────────────────────────────
#CHEEP_SAMPLE_TSV="${CHEEP_STRAINS_DIR}/${CHEEP_SAMPLE}.strains.tsv"
#ATL_SAMPLE_TSV="${ATL_STRAINS_DIR}/${ATL_SAMPLE}.strains.tsv"

#if [[ ! -f "$CHEEP_SAMPLE_TSV" ]]; then
 #   echo "ERROR: strains TSV not found for cheep sample ${CHEEP_SAMPLE}: ${CHEEP_SAMPLE_TSV}" >&2
 #   exit 1
#fi
#if [[ ! -f "$ATL_SAMPLE_TSV" ]]; then
 #   echo "ERROR: strains TSV not found for atl sample ${ATL_SAMPLE}: ${ATL_SAMPLE_TSV}" >&2
  #  exit 1
#fi

# ── Step 1: prepare combined reference FASTA ──────────────────────────────────
#conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

#PAIR_REF_FASTA="$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}.fasta"

#straingr prepare-ref \
#    -s "$CHEEP_SAMPLE_TSV" "$ATL_SAMPLE_TSV" \
#    -p "$REF_DB_PATH" \
#    -S "$SIM_TSV" \
 #   -o "$PAIR_REF_FASTA"

#conda deactivate

# ── Step 2: index the combined reference once ──────────────────────────────────
#conda activate /projects/mnadimlab/nmolson/condaenvs/bwa-mem2
#bwa index "$PAIR_REF_FASTA"
#conda deactivate

# ── Step 3: align + variant-call each sample ──────────────────────────────────
#for SAMPLE in "$CHEEP_SAMPLE" "$ATL_SAMPLE"; do

    # Choose the correct trim directory and file extension
  #  if [[ "$SAMPLE" == "$CHEEP_SAMPLE" ]]; then
    #    R1="${TRIMDIR1}/${SAMPLE}_1.trim.fq.gz"
   #     R2="${TRIMDIR1}/${SAMPLE}_2.trim.fq.gz"
  #  else
      #  R1=$(ls "${TRIMDIR2}/${SAMPLE}"*_1.fastq.gz 2>/dev/null | head -1)
      #  R2=$(ls "${TRIMDIR2}/${SAMPLE}"*_2.fastq.gz 2>/dev/null | head -1)
      #  if [[ -z "$R1" || -z "$R2" ]]; then
      #      echo "ERROR: could not find fastq.gz files for ${SAMPLE} in ${TRIMDIR2}" >&2
     #       exit 1
    #    fi
   # fi
   # BAM_OUT="$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}_${SAMPLE}.bam"

    # Alignment
   # conda activate /projects/mnadimlab/nmolson/condaenvs/bwa-mem2

   # bwa mem -I 300 -t "$SLURM_CPUS_PER_TASK" "$PAIR_REF_FASTA" "$R1" "$R2" | \
    #    samtools sort -@ 2 -O BAM -o "$BAM_OUT" -
   # samtools index "$BAM_OUT"

   # conda deactivate

    # Variant calling
   # conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

   # straingr call "$PAIR_REF_FASTA" "$BAM_OUT" \
   #     --hdf5-out "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}_${SAMPLE}.hdf5" \
  #      --summary  "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}_${SAMPLE}.tsv" \
 #       --tracks all

#    rm -f "$BAM_OUT" "${BAM_OUT}.bai"

  #  conda deactivate

#done

# ── Step 4: compare strains across the two samples ────────────────────────────
conda activate /projects/mnadimlab/nmolson/condaenvs/strainge

straingr compare \
    "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}_${CHEEP_SAMPLE}.hdf5" \
    "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}_${ATL_SAMPLE}.hdf5" \
    -o "$OUTDIR/${CHEEP_SAMPLE}.vs.${ATL_SAMPLE}.summary.tsv" \
    -d "$OUTDIR/${CHEEP_SAMPLE}.vs.${ATL_SAMPLE}.details.tsv"

conda deactivate
# ── Cleanup shared files ───────────────────────────────────────────────────────
#rm -f "$PAIR_REF_FASTA" "${PAIR_REF_FASTA}".*
#rm -f "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}"*wig \
 #     "$OUTDIR/${CHEEP_SAMPLE}_${ATL_SAMPLE}"*bed

