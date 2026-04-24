#!/bin/bash
#SBATCH --job-name="atl-cheep-compile"
#SBATCH --nodes=1
#SBATCH --account=mnadim
#SBATCH --tasks-per-node=1  # The same as the number of cores
#SBATCH --cpus-per-task=1  # 
#SBATCH --partition=largemem
#SBATCH --time=0-01:00:00
#SBATCH --output=log/atl-cheep-compile-%j.o
#SBATCH --error=log/atl-cheep-compile-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maya.l.nadimpalli@emory.edu


#Set paths
COMPARE_DIR=/projects/mnadimlab/nmolson/data/straingr/straingr_cheep_atl
OUTPUT1=$COMPARE_DIR/cheep_atl_sample_comparisons.tsv

cd $COMPARE_DIR

# Compile all summary TSV files and filter only by commonPCT >=0.5
# - Extract last field (split by "_") for sample1 and sample2
# - Filter: commonPct >= 0.5
# - Retain: sample1, sample2, ref, commonPct, singleAgreePct, gapJaccardSim

shopt -s nullglob
tsv_files=( *.summary.tsv )
if [ ${#tsv_files[@]} -eq 0 ]; then
    echo "No summary TSV files found in $COMPARE_DIR"
    exit 1
fi

awk '
BEGIN { FS="\t"; OFS="\t"; header_printed=0 }
FNR==1 {
    if (!header_printed) {
        # Map column names to indices
        for (i=1; i<=NF; i++) {
            col[$i] = i
        }
        print "sample1", "sample2", "ref", "commonPct", "singleAgreePct", "gapJaccardSim"
        header_printed=1
    }
    next
}
{
	
    s1 = $col["sample1"]
    s2 = $col["sample2"]
    commonPct     = $col["commonPct"]
    singleAgreePct = $col["singleAgreePct"]
    gapJaccardSim = $col["gapJaccardSim"]
    ref           = $col["ref"]

    # Extract last field when split by "_"
    n1 = split(s1, a1, "_"); s1_clean = a1[n1]
    n2 = split(s2, a2, "_"); s2_clean = a2[n2]

    if (commonPct+0 >= 0.5) {
        print s1_clean, s2_clean, ref, commonPct, singleAgreePct, gapJaccardSim
    }
}
' "${tsv_files[@]}" > $OUTPUT1

echo "Done. Filtered results written to $OUTPUT1"


