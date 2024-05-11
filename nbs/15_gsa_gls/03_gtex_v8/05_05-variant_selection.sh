#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=variant_selection
#SBATCH --output=variant_selection.%j.out
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milton.pividori@cuanschutz.edu

# Code taken and adapted from:
#  https://github.com/hakyimlab/summary-gwas-imputation/wiki/Reference-Data-Set-Compilation

# make sure we use the number of CPUs specified
export PHENOPLIER_N_JOBS=1
export NUMBA_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

GTEX_V8_DIR=${PHENOPLIER_EXTERNAL_GTEX_V8_DIR}
CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/03_gtex_v8

OUTPUT_DIR=${PHENOPLIER_PHENOMEXCAN_LD_BLOCKS_GTEX_V8_GENOTYPE_DIR}
mkdir -p ${OUTPUT_DIR}

PYTHON_EXECUTABLE="${PHENOPLIER_GWAS_IMPUTATION_CONDA_ENV}/bin/python"
if [ ! -f ${PYTHON_EXECUTABLE} ]; then
    >&2 echo "The python executable does not exist: ${PYTHON_EXECUTABLE}"
    exit 1
fi

# The file that this command generates has this name in the summary-gwas-imputation pipeline:
#  gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz
${PYTHON_EXECUTABLE} ${PHENOPLIER_GWAS_IMPUTATION_BASE_DIR}/src/get_reference_metadata.py \
  -genotype ${GTEX_V8_DIR}/generated/gtex_v8_eur_filtered.txt.gz \
  -annotation ${GTEX_V8_DIR}/analysis_supplement/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz \
  -rsid_column rs_id_dbSNP151_GRCh38p7 \
  -filter MAF 0.01 \
  -filter TOP_CHR_POS_BY_FREQ \
  -output ${OUTPUT_DIR}/variant_metadata.txt.gz
