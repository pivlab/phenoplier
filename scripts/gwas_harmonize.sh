#!/bin/bash
set -euo pipefail
# IFS=$'\n\t'

# Runs the harmonization step of the pipeline here: https://github.com/hakyimlab/summary-gwas-imputation

# read arguments
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input-gwas-file)
      INPUT_GWAS_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -l|--liftover-chain-file)
      LIFTOVER_CHAIN_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output-dir)
      OUTPUT_DIR="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--sample-size)
      SAMPLE_SIZE="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--sample-n-cases)
      SAMPLES_N_CASES="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

#
# check arguments
#
if [ -z "${INPUT_GWAS_FILE}" ]; then
    >&2 echo "Error, --input-gwas-file <value> not provided"
    exit 1
fi

if [ -z "${OUTPUT_DIR}" ]; then
    >&2 echo "Error, --output-dir <value> not provided"
    exit 1
fi

if [ -z "${SAMPLE_SIZE}" ]; then
    >&2 echo "Error, --samples-size <value> not provided"
    exit 1
fi

ARGS_SAMPLE_N_CASES=""
if [ -z "${SAMPLES_N_CASES}" ]; then
    echo "Parameter --samples-n-cases was not provided"
else
    # if value is NA, then we don't need to pass it
    if [ "${SAMPLES_N_CASES}" = "NA" ]; then
      echo "Parameter --samples-n-cases is NA, skipping"
    else
      ARGS_SAMPLE_N_CASES="--insert_value n_cases ${SAMPLES_N_CASES}"
    fi
fi

# liftover argument (optional)
LIFTOVER_ARG=""
if [ ! -z "${LIFTOVER_CHAIN_FILE}" ]; then
    LIFTOVER_ARG="-liftover ${LIFTOVER_CHAIN_FILE}"
fi

#
# Global PhenoPLIER environmental variables
#

# make sure we have environment variables with configuration
if [ -z "${PHENOPLIER_ROOT_DIR}" ] || [ -z "${PHENOPLIER_GWAS_IMPUTATION_BASE_DIR}" ]; then
    >&2 echo "PhenoPLIER configuration was not loaded"
    exit 1
fi

PYTHON_EXECUTABLE="${PHENOPLIER_GWAS_IMPUTATION_CONDA_ENV}/bin/python"
if [ ! -f ${PYTHON_EXECUTABLE} ]; then
    >&2 echo "The python executable does not exist: ${PYTHON_EXECUTABLE}"
    exit 1
fi

A1000G_VARIANTS_METADATA_FILE="${PHENOPLIER_PHENOMEXCAN_LD_BLOCKS_1000G_GENOTYPE_DIR}/variant_metadata.txt.gz"
if [ ! -f ${A1000G_VARIANTS_METADATA_FILE} ]; then
    >&2 echo "The 1000 Genomes variants metadata file does not exist: ${A1000G_VARIANTS_METADATA_FILE}"
    exit 1
fi

# Build GWAS columns arguments
GWAS_CONFIG_FILE="${INPUT_GWAS_FILE}.config"
if [ ! -f ${GWAS_CONFIG_FILE} ]; then
    >&2 echo "The GWAS config file does not exist: ${GWAS_CONFIG_FILE}"
    exit 1
fi

echo "Loading GWAS config file: ${GWAS_CONFIG_FILE}"
. ${GWAS_CONFIG_FILE}

ARGS_CHROMOSOME_COL="-output_column_map ${PHENOPLIER_GWAS_COL_CHROMOSOME} chromosome"
if [ -z "${PHENOPLIER_GWAS_COL_CHROMOSOME}" ]; then
    ARGS_CHROMOSOME_COL=""
fi

ARGS_VAR_ID_COL="-output_column_map ${PHENOPLIER_GWAS_COL_VARIANT_ID} variant_id"
if [ -z "${PHENOPLIER_GWAS_COL_VARIANT_ID}" ]; then
    ARGS_VAR_ID_COL=""
fi

ARGS_POS_COL="-output_column_map ${PHENOPLIER_GWAS_COL_POSITION} position"
if [ -z "${PHENOPLIER_GWAS_COL_POSITION}" ]; then
    ARGS_POS_COL=""
fi

ARGS_EFFECT_ALLELE_COL="-output_column_map ${PHENOPLIER_GWAS_COL_EFFECT_ALLELE} effect_allele"
if [ -z "${PHENOPLIER_GWAS_COL_EFFECT_ALLELE}" ]; then
    ARGS_EFFECT_ALLELE_COL=""
fi

ARGS_REF_ALLELE_COL="-output_column_map ${PHENOPLIER_GWAS_COL_NON_EFFECT_ALLELE} non_effect_allele"
if [ -z "${PHENOPLIER_GWAS_COL_NON_EFFECT_ALLELE}" ]; then
    ARGS_REF_ALLELE_COL=""
fi

ARGS_OR_COL="-output_column_map ${PHENOPLIER_GWAS_COL_OR} or"
if [ -z "${PHENOPLIER_GWAS_COL_OR}" ]; then
    ARGS_OR_COL=""
fi

ARGS_EFFECT_SIZE_COL="-output_column_map ${PHENOPLIER_GWAS_COL_EFFECT_SIZE} effect_size"
if [ -z "${PHENOPLIER_GWAS_COL_EFFECT_SIZE}" ]; then
    ARGS_EFFECT_SIZE_COL=""
fi

ARGS_SE_COL="-output_column_map ${PHENOPLIER_GWAS_COL_STANDARD_ERROR} standard_error"
if [ -z "${PHENOPLIER_GWAS_COL_STANDARD_ERROR}" ]; then
    ARGS_SE_COL=""
fi

ARGS_PVAL_COL="-output_column_map ${PHENOPLIER_GWAS_COL_PVALUE} pvalue"
if [ -z "${PHENOPLIER_GWAS_COL_PVALUE}" ]; then
    ARGS_PVAL_COL=""
fi

# Create output directory
mkdir -p ${OUTPUT_DIR}

INPUT_GWAS_FILENAME=$(basename ${INPUT_GWAS_FILE})
OUTPUT_FILENAME_BASE="${INPUT_GWAS_FILENAME%.*}-harmonized"

set -x
${PYTHON_EXECUTABLE} ${PHENOPLIER_GWAS_IMPUTATION_BASE_DIR}/src/gwas_parsing.py \
    -gwas_file ${INPUT_GWAS_FILE} \
    -separator $'\t' \
    -snp_reference_metadata ${A1000G_VARIANTS_METADATA_FILE} METADATA \
    --chromosome_format ${ARGS_CHROMOSOME_COL} ${ARGS_VAR_ID_COL} ${ARGS_POS_COL} ${ARGS_EFFECT_ALLELE_COL} ${ARGS_REF_ALLELE_COL} ${ARGS_EFFECT_SIZE_COL} ${ARGS_OR_COL} ${ARGS_SE_COL} ${ARGS_PVAL_COL} \
    --insert_value sample_size ${SAMPLE_SIZE} ${ARGS_SAMPLE_N_CASES} \
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
    ${LIFTOVER_ARG} -output "${OUTPUT_DIR}/${OUTPUT_FILENAME_BASE}.txt" \
>"${OUTPUT_DIR}/${OUTPUT_FILENAME_BASE}.log" 2>&1
set +x
