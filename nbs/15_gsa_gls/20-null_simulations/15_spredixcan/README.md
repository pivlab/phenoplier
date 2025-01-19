# Overview

This folder has the scripts to run the PrediXcan family of methods on GWAS on randomly generated phenotypes (`../10_gwas_harmonization`).


# Load Alpine-specific paths and PhenoPLIER configuration

You need to run `acompile` to start a new interactive session before running the
commands below.

```bash
# load conda environment
module load mambaforge/23.1.0-1 gnu_parallel/20210322
mamba activate phenoplier_light

# load PhenoPLIER config
. /pl/active/pivlab/projects/mpividori/phenoplier/scripts/alpine/env.sh

# load in bash session all PhenoPLIER environmental variables
eval `python ${PHENOPLIER_CODE_DIR}/libs/conf.py`

# make sure they were loaded correctly
# should output something like /project/...
echo $PHENOPLIER_ROOT_DIR
```

# Desktop computer

```bash
# load conda environment
conda activate phenoplier_light

# load PhenoPLIER config
. scripts/env.sh

# load in bash session all PhenoPLIER environmental variables
eval `python ${PHENOPLIER_CODE_DIR}/libs/conf.py`

# Set executor
export PHENOPLIER_JOBS_EXECUTOR="bash"

# make sure they were loaded correctly
echo $PHENOPLIER_ROOT_DIR
```


# Download the necessary data

```bash
python ~/projects/phenoplier/environment/scripts/setup_data.py \
  --actions \
    download_setup_metaxcan \
    download_predixcan_mashr_prediction_models \
    download_mashr_expression_smultixcan_snp_covariance
```


# Run cluster jobs

The `cluster_jobs/` folder has the job scripts to run on Penn's LPC cluster.
To run the jobs in order, you need to execute the command below.
The `_tmp` folder stores logs and needs to be created.

## S-PrediXcan

Here we need to use some templating, because we run across random phenotypes and tissues.

```bash
cd nbs/15_gsa_gls/20-null_simulations/15_spredixcan

run_job() {
  export pheno_id=$1
  export tissue=$2
  
  cat cluster_jobs/01_spredixcan_job-template.sh | envsubst '${pheno_id} ${tissue}' | ${PHENOPLIER_JOBS_EXECUTOR}
}

export -f run_job

# (optional) export function definition so it's included in the Docker container
export PHENOPLIER_BASH_FUNCTIONS_CODE="$(declare -f run_job)"

# Run
parallel -j10 run_job {} ::: {0..299} ::: ${PHENOPLIER_PHENOMEXCAN_PREDICTION_MODELS_MASHR_TISSUES}
```

Checks:

```bash
bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/spredixcan \
  -p "INFO - Sucessfully processed metaxcan association"

bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/spredixcan \
  -p "INFO - 90 % of model's snps"

bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/spredixcan \
  -p "INFO - 99 % of model's snps"
```

If any job failed, check `../10_gwas_harmonization/README.md`, which has python code to get a list of unfinished jobs.


## S-MultiXcan

```bash
cd nbs/15_gsa_gls/20-null_simulations/15_spredixcan

run_job() {
  export pheno_id=$1
  
  cat cluster_jobs/05_smultixcan_job-template.sh | envsubst '${pheno_id}' | ${PHENOPLIER_JOBS_EXECUTOR}
}

export -f run_job

# (optional) export function definition so it's included in the Docker container
export PHENOPLIER_BASH_FUNCTIONS_CODE="$(declare -f run_job)"

# Run
parallel -j10 run_job {} ::: {0..99}
```

Checks:

```bash
bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/smultixcan \
  -p "INFO - Ran multi tissue"

# Check S-PrediXcan files
bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/smultixcan \
  -p "Level 9 - Loading metaxcan " \
  -c 49

# Check tissues loaded
bash check_job.sh \
  -i ${PHENOPLIER_RESULTS_GLS_NULL_SIMS}/twas/smultixcan \
  -p "Level 9 - Processing " \
  -c 49
```


## Monitoring jobs

Check jobs with command `bjobs`.
Or, for a constantly-updated monitoring (refreshing every 2 seconds):
```bash
watch -n 2 bjobs
```

Logs for `random_pheno0` are in `random_pheno1.*` (indexes are different because LPC arrays cannot start with zero).

To kill running jobs:
```bash
bjobs | grep RUN | cut -d ' ' -f1 | xargs -I {} bkill {}
```


# QQ plots

Notebook `15-twas-qqplot.ipynb` checks that the distribution of pvalues is as expected.




REMEMBER TO RUN QQPLOTS NOTEBOOKS WHEN ALL IS DONE
