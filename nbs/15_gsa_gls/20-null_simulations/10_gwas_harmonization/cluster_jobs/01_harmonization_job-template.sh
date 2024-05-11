#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=random_pheno${pheno_id}
#SBATCH --output=_tmp/harmonization/random_pheno${pheno_id}.%j.out
#SBATCH --time=00:15:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milton.pividori@cuanschutz.edu

# make sure we use the number of CPUs specified
export n_jobs=1
export PHENOPLIER_N_JOBS=${n_jobs}
export NUMBA_NUM_THREADS=${n_jobs}
export MKL_NUM_THREADS=${n_jobs}
export OPEN_BLAS_NUM_THREADS=${n_jobs}
export NUMEXPR_NUM_THREADS=${n_jobs}
export OMP_NUM_THREADS=${n_jobs}

CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/20-null_simulations/10_gwas_harmonization
GWAS_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/gwas"
OUTPUT_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/harmonized_gwas"

GWAS_JOBINDEX=${pheno_id}

bash ${CODE_DIR}/01_harmonize.sh \
  --input-gwas-file ${GWAS_DIR}/random.pheno${GWAS_JOBINDEX}.glm.linear.tsv.gz \
  --liftover-chain-file ${PHENOPLIER_GENERAL_LIFTOVER_HG19_TO_HG38} \
  --output-dir ${OUTPUT_DIR}

