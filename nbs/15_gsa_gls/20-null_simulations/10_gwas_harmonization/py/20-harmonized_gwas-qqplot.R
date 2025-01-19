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
# It takes a GWAS that was imputed and postprocessed (using the PrediXcan scripts here https://github.com/hakyimlab/summary-gwas-imputation) on a random phenotype and verifies that the Manhattan and QQ-plots look fine (without inflation).

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
GWAS_DIR <- file.path(
    GLS_NULL_SIMS_DIR,
    # "_pvalue",
    "harmonized_gwas"
)

# %% tags=[]
GWAS_DIR

# %%
title_prefix = "harmonized"

# %% [markdown] tags=[]
# # Random pheno 0 (beta/se)

# %%
pheno_code <- "0"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv-harmonized.txt"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %% [markdown] tags=[]
# ### Extract chromosome

# %% tags=[]
unique(gwas$chromosome)

# %% tags=[]
gwas$chrom <- gsub("chr([0-9]+)", "\\1", gwas$chromosome)
gwas <- transform(gwas, chrom = as.numeric(chrom))

# %% tags=[]
unique(gwas$chrom)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(pvalue, na.rm = TRUE),
    sd_value   = sd(pvalue, na.rm = TRUE),
    min_value  = min(pvalue, na.rm = TRUE),
    max_value  = max(pvalue, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(pvalue)) %>%
  summarise(
    mean_neg_log10 = mean(LOGP, na.rm = TRUE),
    sd_neg_log10   = sd(LOGP, na.rm = TRUE),
    min_neg_log10  = min(LOGP, na.rm = TRUE),
    max_neg_log10  = max(LOGP, na.rm = TRUE)
  )

# %%
gwas_mlog_stats

# %%
gwas_max_mlog_noninf <- gwas %>%
  mutate(LOGP = -log10(pvalue)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_zs_stats <- gwas %>%
  summarise(
    mean_value = mean(zscore, na.rm = TRUE),
    sd_value   = sd(zscore, na.rm = TRUE),
    min_value  = min(zscore, na.rm = TRUE),
    max_value  = max(zscore, na.rm = TRUE),
  )

# %%
gwas_zs_stats

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(effect_size, na.rm = TRUE),
    sd_value   = sd(effect_size, na.rm = TRUE),
    min_value  = min(effect_size, na.rm = TRUE),
    max_value  = max(effect_size, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(standard_error, na.rm = TRUE),
    sd_value   = sd(standard_error, na.rm = TRUE),
    min_value  = min(standard_error, na.rm = TRUE),
    max_value  = max(standard_error, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown] tags=[]
# ## Manhattan plot

# %% tags=[]
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "chrom",
  bp = "position",
  p = "pvalue",
  snp = "variant_id",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown] tags=[]
# ## QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$pvalue,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% [markdown] tags=[]
# # Random pheno 0 (pvalue)

# %%
pheno_code <- "0"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv-harmonized.txt"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %% [markdown] tags=[]
# ### Extract chromosome

# %% tags=[]
unique(gwas$chromosome)

# %% tags=[]
gwas$chrom <- gsub("chr([0-9]+)", "\\1", gwas$chromosome)
gwas <- transform(gwas, chrom = as.numeric(chrom))

# %% tags=[]
unique(gwas$chrom)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(pvalue, na.rm = TRUE),
    sd_value   = sd(pvalue, na.rm = TRUE),
    min_value  = min(pvalue, na.rm = TRUE),
    max_value  = max(pvalue, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(pvalue)) %>%
  summarise(
    mean_neg_log10 = mean(LOGP, na.rm = TRUE),
    sd_neg_log10   = sd(LOGP, na.rm = TRUE),
    min_neg_log10  = min(LOGP, na.rm = TRUE),
    max_neg_log10  = max(LOGP, na.rm = TRUE)
  )

# %%
gwas_mlog_stats

# %%
gwas_max_mlog_noninf <- gwas %>%
  mutate(LOGP = -log10(pvalue)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_zs_stats <- gwas %>%
  summarise(
    mean_value = mean(zscore, na.rm = TRUE),
    sd_value   = sd(zscore, na.rm = TRUE),
    min_value  = min(zscore, na.rm = TRUE),
    max_value  = max(zscore, na.rm = TRUE),
  )

# %%
gwas_zs_stats

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(effect_size, na.rm = TRUE),
    sd_value   = sd(effect_size, na.rm = TRUE),
    min_value  = min(effect_size, na.rm = TRUE),
    max_value  = max(effect_size, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(standard_error, na.rm = TRUE),
    sd_value   = sd(standard_error, na.rm = TRUE),
    min_value  = min(standard_error, na.rm = TRUE),
    max_value  = max(standard_error, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown] tags=[]
# ## Manhattan plot

# %% tags=[]
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "chrom",
  bp = "position",
  p = "pvalue",
  snp = "variant_id",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown] tags=[]
# ## QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$pvalue,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% tags=[]
