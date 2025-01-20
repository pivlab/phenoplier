# Overview

This folder has the scripts to run GLS PhenoPLIER (associations between LVs/gene modules and traits on randomly generated phenotypes (`../15_spredixcan`).

Before running these steps, **it is necessary** to generate a correlation matrix for predicted gene expression _specific_ for these random phenotypes (see `nbs/15_gsa_gls/README.md`).


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
    download_phenomexcan_rapid_gwas_pheno_info \
    download_phenomexcan_rapid_gwas_data_dict_file \
    download_uk_biobank_coding_3 \
    download_uk_biobank_coding_6 \
    download_phenomexcan_gtex_gwas_pheno_info \
    download_gene_map_id_to_name \
    download_gene_map_name_to_id \
    download_biomart_genes_hg38 \
    download_multiplier_model_z_pkl
```


# Run cluster jobs

The `cluster_jobs/` folder has the job scripts to run on Penn's LPC cluster.
To run the jobs in order, you need to execute the command below.

## Run LV-trait associations

```bash
cd nbs/15_gsa_gls/20-null_simulations/20_gls_phenoplier

run_job () {
  cluster_job_file="$1"
  export pheno_id="$2"
  
  cat $cluster_job_file | envsubst '${pheno_id}' | ${PHENOPLIER_JOBS_EXECUTOR}
}

export -f run_job

# (optional) export function definition so it's included in the Docker container
export PHENOPLIER_BASH_FUNCTIONS_CODE="$(declare -f run_job)"
```

### With covariates

```bash
# (optional) Run OLS model
#parallel -j10 run_job cluster_jobs/covars/01_gls-use_ols-template.sh {} ::: {0..999}

# GLS:
parallel -j10 run_job cluster_jobs/covars/10_gls_phenoplier-sub_corr-template.sh {} ::: {0..99}
```

```bash
bash ${PHENOPLIER_CODE_DIR}/scripts/check_job.sh \
    -i _tmp/ \
    -f '*.error' \
    -p "INFO: Writing results to"

bash ${PHENOPLIER_CODE_DIR}/scripts/check_job.sh \
    -i _tmp/ \
    -f '*.error' \
    -p "INFO: Using covariates: \['gene_size', 'gene_size_log', 'gene_density', 'gene_density_log'\]"

# for OLS
bash ${PHENOPLIER_CODE_DIR}/scripts/check_job.sh \
    -i _tmp/gls_phenoplier_ols/ \
    -f '*.error' \
    -p "INFO: Using a Ordinary Least Squares (OLS) model"

# for GLS
bash ${PHENOPLIER_CODE_DIR}/scripts/check_job.sh \
    -i _tmp/gls_phenoplier/ \
    -f '*.error' \
    -p "INFO: Correlation matrix is a directory"
```


## Monitoring jobs

Check jobs with command `bjobs`.
Or, for a constantly-updated monitoring (refreshing every 2 seconds):
```bash
watch -n 2 bjobs
```

To kill running jobs:
```bash
bjobs | grep RUN | cut -d ' ' -f1 | xargs -I {} bkill {}
```

# QQ plots

Notebook `05-twas-qqplot.ipynb` checks that the distribution of pvalues is as expected.
