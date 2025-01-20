# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown] tags=[]
# # Description

# %% [markdown] tags=[]
# It takes GLS PhenoPLIER results on a random phenotype and verifies that the QQ-plots look fine (without inflation).

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
library(tidyverse)

# %% tags=[]
library(qqman)

# %% [markdown] tags=[]
# # Paths

# %% tags=[]
GLS_NULL_SIMS_DIR <- Sys.getenv("PHENOPLIER_RESULTS_GLS_NULL_SIMS")

# %% tags=[]
GLS_NULL_SIMS_DIR

# %% tags=[]
GLS_PHENOPLIER_DIR <- file.path(GLS_NULL_SIMS_DIR, "phenoplier", "ukbb_eur", "covars", "gls-gtex_v8_mashr-sub_corr")

# %% tags=[]
GLS_PHENOPLIER_DIR

# %% [markdown] tags=[]
# # Random pheno 0

# %%
pheno_code <- "0"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gls_ph <- as.data.frame(read_tsv(file.path(GLS_PHENOPLIER_DIR, paste0("random.pheno", pheno_code,"-gls_phenoplier.tsv.gz"))))

# %% tags=[]
dim(gls_ph)

# %% tags=[]
head(gls_ph)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gls_ph$pvalue_onesided,
    main = paste0("random pheno ", pheno_code, " - QQ plot of GLS PhenoPLIER")
)

# %% [markdown] tags=[]
# # Random pheno 1

# %%
pheno_code <- "1"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gls_ph <- as.data.frame(read_tsv(file.path(GLS_PHENOPLIER_DIR, paste0("random.pheno", pheno_code,"-gls_phenoplier.tsv.gz"))))

# %% tags=[]
dim(gls_ph)

# %% tags=[]
head(gls_ph)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gls_ph$pvalue_onesided,
    main = paste0("random pheno ", pheno_code, " - QQ plot of GLS PhenoPLIER")
)

# %% [markdown] tags=[]
# # Random pheno 2

# %%
pheno_code <- "2"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gls_ph <- as.data.frame(read_tsv(file.path(GLS_PHENOPLIER_DIR, paste0("random.pheno", pheno_code,"-gls_phenoplier.tsv.gz"))))

# %% tags=[]
dim(gls_ph)

# %% tags=[]
head(gls_ph)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gls_ph$pvalue_onesided,
    main = paste0("random pheno ", pheno_code, " - QQ plot of GLS PhenoPLIER")
)

# %% [markdown] tags=[]
# # Random pheno 3

# %%
pheno_code <- "3"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gls_ph <- as.data.frame(read_tsv(file.path(GLS_PHENOPLIER_DIR, paste0("random.pheno", pheno_code,"-gls_phenoplier.tsv.gz"))))

# %% tags=[]
dim(gls_ph)

# %% tags=[]
head(gls_ph)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gls_ph$pvalue_onesided,
    main = paste0("random pheno ", pheno_code, " - QQ plot of GLS PhenoPLIER")
)

# %% [markdown] tags=[]
# # Random pheno 4

# %%
pheno_code <- "4"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gls_ph <- as.data.frame(read_tsv(file.path(GLS_PHENOPLIER_DIR, paste0("random.pheno", pheno_code,"-gls_phenoplier.tsv.gz"))))

# %% tags=[]
dim(gls_ph)

# %% tags=[]
head(gls_ph)

# %% [markdown] tags=[]
# ### QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gls_ph$pvalue_onesided,
    main = paste0("random pheno ", pheno_code, " - QQ plot of GLS PhenoPLIER")
)

# %% tags=[]
