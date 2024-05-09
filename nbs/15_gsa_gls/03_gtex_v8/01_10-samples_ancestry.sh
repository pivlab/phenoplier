#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=samples_ancestry
#SBATCH --output=samples_ancestry.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milton.pividori@cuanschutz.edu

# make sure we use the number of CPUs specified
export PHENOPLIER_N_JOBS=1
export NUMBA_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPEN_BLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

CODE_DIR=${PHENOPLIER_CODE_DIR}/nbs/15_gsa_gls/03_gtex_v8

module purge
module load mambaforge/23.1.0-1
conda activate phenoplier_light

. ${PHENOPLIER_CODE_DIR}/scripts/alpine/env.sh
eval `python ${PHENOPLIER_CODE_DIR}/libs/conf.py`

# Generate the GTEx ancestry file
python << END
import pandas as pd

import conf

OUTPUT_DIR = conf.EXTERNAL["GTEX_V8_DIR"] / "generated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# read data
data = pd.read_csv(conf.EXTERNAL["GTEX_V8_DIR"] / "analysis_supplement" / "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_899Indiv_QC_TABLE.tsv", sep="\t")
_tmp = data[["SAMPLE_ID", "InferredAncestry"]].dropna()
assert _tmp.shape == (866, 2)

# save
_tmp.to_csv(OUTPUT_DIR / "GTEx_ethnicity.txt", sep="\t", index=False)
END

# Separate samples by ancestry
Rscript ${CODE_DIR}/01_10-samples_by_ethnicity.r

