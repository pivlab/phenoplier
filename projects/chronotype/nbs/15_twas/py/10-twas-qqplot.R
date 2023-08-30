# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown] tags=[]
# # Description

# %% [markdown] tags=[]
# It takes TWAS results and verifies that the QQ-plots look fine (without inflation).

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
library(tidyverse)

# %% tags=[]
library(qqman)

# %% [markdown] tags=[]
# # Settings

# %% tags=["parameters"]
PROJECTS_TRAIT_KEY <- "CHRONOTYPE"

# %% tags=["injected-parameters"]
# Parameters
PHENOPLIER_NOTEBOOK_FILEPATH <- "projects/chronotype/nbs/15_twas/10-twas-qqplot.ipynb"


# %% [markdown] tags=[]
# # Paths

# %% tags=[]
BASE_DIR <- Sys.getenv(paste0("PHENOPLIER_PROJECTS_", PROJECTS_TRAIT_KEY, "_RESULTS_DIR"))

# %% tags=[]
BASE_DIR

# %% tags=[]
SPREDIXCAN_DIR <- file.path(BASE_DIR, "twas", "spredixcan")

# %% tags=[]
SPREDIXCAN_DIR

# %% tags=[]
SMULTIXCAN_DIR <- file.path(BASE_DIR, "twas", "smultixcan")

# %% tags=[]
SMULTIXCAN_DIR

# %% [markdown] tags=[]
# # Load trait info file

# %% tags=[]
trait_info <- read_csv(Sys.getenv(paste0("PHENOPLIER_PROJECTS_", PROJECTS_TRAIT_KEY, "_TRAITS_INFO_FILE")), col_names = TRUE)

# %% tags=[]
trait_info

# %% [markdown] tags=[]
# # Chronotype

# %% tags=[]
gwas_title <- "Chronotype (imputed)"

# %% [markdown] tags=[]
# ## S-PrediXcan

# %% [markdown] tags=[]
# ### Load data

# %% tags=[]
twas <- as.data.frame(read_csv(file.path(SPREDIXCAN_DIR, "chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt-gtex_v8-mashr-Whole_Blood.csv")))

# %% tags=[]
dim(twas)

# %% tags=[]
head(twas)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(twas$pvalue, main = paste0("Q-Q plot of S-PrediXcan for ", gwas_title))

# %% [markdown] tags=[]
# ## S-MultiXcan

# %% [markdown] tags=[]
# ### Load data

# %% tags=[]
twas <- as.data.frame(read_table(file.path(SMULTIXCAN_DIR, "chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt-gtex_v8-mashr-smultixcan.txt")))

# %% tags=[]
dim(twas)

# %% tags=[]
head(twas)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(twas$pvalue, main = paste0("Q-Q plot of S-MultiXcan for ", gwas_title))

# %% tags=[]
