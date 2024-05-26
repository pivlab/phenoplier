#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=random_pheno${pheno_id}
#SBATCH --output=_tmp/smultixcan/random_pheno${pheno_id}.%j.out
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
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

CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/20-null_simulations/15_spredixcan
FINAL_IMPUTED_GWAS_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/final_imputed_gwas"
SPREDIXCAN_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/spredixcan"
OUTPUT_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/smultixcan"

GWAS_JOBINDEX=${pheno_id}

# +
bash ${CODE_DIR}/05_smultixcan.sh \
  --input-gwas-file ${FINAL_IMPUTED_GWAS_DIR}/random.pheno${GWAS_JOBINDEX}.glm-imputed.txt.gz \
  --spredixcan-folder ${SPREDIXCAN_DIR} \
  --phenotype-name "random.pheno${GWAS_JOBINDEX}" \
  --output-dir ${OUTPUT_DIR}