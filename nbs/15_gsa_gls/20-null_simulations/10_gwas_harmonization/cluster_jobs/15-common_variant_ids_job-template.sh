#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=random_pheno${pheno_id}
#SBATCH --output=_tmp/common_var_ids/random_pheno${pheno_id}.%j.out
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#DISABLED SBATCH --mail-type=ALL
#DISABLED SBATCH --mail-user=milton.pividori@cuanschutz.edu

# make sure we use the number of CPUs specified
export n_jobs=1
export PHENOPLIER_N_JOBS=${n_jobs}
export NUMBA_NUM_THREADS=${n_jobs}
export MKL_NUM_THREADS=${n_jobs}
export OPEN_BLAS_NUM_THREADS=${n_jobs}
export NUMEXPR_NUM_THREADS=${n_jobs}
export OMP_NUM_THREADS=${n_jobs}

CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/20-null_simulations/10_gwas_harmonization
POST_GWAS_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/post_imputed_gwas"
OUTPUT_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/final_imputed_gwas"
mkdir -p ${OUTPUT_DIR}

GWAS_JOBINDEX=${pheno_id}

python ${CODE_DIR}/15_common_variant_ids.py \
  --input-gwas-file ${POST_GWAS_DIR}/random.pheno${GWAS_JOBINDEX}.glm-imputed.txt.gz \
  --common-variant-ids-file ${POST_GWAS_DIR}/common_variant_ids.pkl \
  --output-dir ${OUTPUT_DIR}

