#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=random_pheno${pheno_id}-${chromosome}-${batch_id}
#SBATCH --output=_tmp/imputation/random_pheno${pheno_id}-${chromosome}-${batch_id}.%j.out
#SBATCH --time=00:45:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
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
#GWAS_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/gwas"
HARMONIZED_GWAS_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/harmonized_gwas"
OUTPUT_DIR="${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/imputed_gwas"

bash ${CODE_DIR}/05_impute.sh \
  --input-gwas-file ${HARMONIZED_GWAS_DIR}/random.pheno${pheno_id}.glm.linear.tsv-harmonized.txt \
  --chromosome ${chromosome} \
  --n-batches 10 \
  --batch-id ${batch_id} \
  --output-dir ${OUTPUT_DIR}

