#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=genotype_compilation
#SBATCH --output=genotype_compilation.%j.out
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40GB
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
OUTPUT_DIR=${PHENOPLIER_EXTERNAL_GTEX_V8_DIR}/generated/genotype
mkdir -p ${OUTPUT_DIR}
CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/03_gtex_v8

OUTPUT_DIR=${PHENOPLIER_PHENOMEXCAN_LD_BLOCKS_GTEX_V8_GENOTYPE_DIR}

PYTHON_EXECUTABLE="${PHENOPLIER_GWAS_IMPUTATION_CONDA_ENV}/bin/python"
if [ ! -f ${PYTHON_EXECUTABLE} ]; then
    >&2 echo "The python executable does not exist: ${PYTHON_EXECUTABLE}"
    exit 1
fi

${PYTHON_EXECUTABLE} ${PHENOPLIER_GWAS_IMPUTATION_BASE_DIR}/src/model_training_genotype_to_parquet.py \
  -input_genotype_file ${GTEX_V8_DIR}/generated/gtex_v8_eur_filtered.txt.gz \
  -snp_annotation_file ${OUTPUT_DIR}/variant_metadata.txt.gz METADATA \
  -parsimony 9 \
  --impute_to_mean \
  --split_by_chromosome \
  --only_in_key \
  -output_prefix ${OUTPUT_DIR}/gtex
