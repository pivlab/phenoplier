#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=genotype_dosage
#SBATCH --output=genotype_dosage.%j.out
#SBATCH --time=18:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milton.pividori@cuanschutz.edu

module load bcftools/1.16

# make sure we use the number of CPUs specified
export PHENOPLIER_N_JOBS=1
export NUMBA_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

GTEX_V8_DIR=${PHENOPLIER_EXTERNAL_GTEX_V8_DIR}
CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/03_gtex_v8

# Function taken and adapted from:
#  https://github.com/hakyimlab/summary-gwas-imputation/wiki/Reference-Data-Set-Compilation
filter_and_convert ()
{
echo -ne "varID\t" | gzip > $3
bcftools view $1 -S $2 --force-samples -Ou |  bcftools query -l | tr '\n' '\t' | sed 's/\t$/\n/' | gzip >> $3

#The first python inline script will check if a variant is blacklisted
NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Starting at $NOW"
bcftools view -S $2 --force-samples $1 -Ou | \
bcftools query -f '%ID[\t%GT]\n' | \
awk '
{
for (i = 1; i <= NF; i++) {
    if (substr($i,0,1) == "c") {
        printf("%s",$i)
    } else if ( substr($i, 0, 1) == ".") {
        printf("\tNA")
    } else if ($i ~ "[0-9]|[0-9]") {
        n = split($i, array, "|")
        printf("\t%d",array[1]+array[2])
    } else {
        #printf("\t%s",$i)
        printf("Unexpected: %s",$i)
        exit 1
    }
}
printf("\n")
}
' | gzip >> $3

NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Ending at $NOW"
}

# Run
filter_and_convert \
  ${GTEX_V8_DIR}/haplotype_phasing/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
  ${GTEX_V8_DIR}/generated/eur_samples.txt \
  ${GTEX_V8_DIR}/generated/gtex_v8_eur_filtered.txt.gz

