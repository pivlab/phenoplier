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
# It takes a GWAS on a random phenotype and verifies that the Manhattan and QQ-plots look fine (without inflation).

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
library(tidyverse)

# %%
library(qqman)

# %% [markdown] tags=[]
# # Paths

# %% tags=[]
GLS_NULL_SIMS_DIR <- Sys.getenv("PHENOPLIER_RESULTS_GLS_NULL_SIMS")

# %% tags=[]
GLS_NULL_SIMS_DIR

# %% tags=[]
GWAS_DIR <- file.path(GLS_NULL_SIMS_DIR, "gwas")

# %% tags=[]
GWAS_DIR

# %%
title_prefix = "original"

# %% [markdown] tags=[]
# # Random pheno 0

# %%
pheno_code <- "0"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv.gz"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %%
gwas <- gwas %>% filter(P >= 0 & P <= 1)

# %%
dim(gwas)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(P, na.rm = TRUE),
    sd_value   = sd(P, na.rm = TRUE),
    min_value  = min(P, na.rm = TRUE),
    max_value  = max(P, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(P)) %>%
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
  mutate(LOGP = -log10(P)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(BETA, na.rm = TRUE),
    sd_value   = sd(BETA, na.rm = TRUE),
    min_value  = min(BETA, na.rm = TRUE),
    max_value  = max(BETA, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(SE, na.rm = TRUE),
    sd_value   = sd(SE, na.rm = TRUE),
    min_value  = min(SE, na.rm = TRUE),
    max_value  = max(SE, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown]
# ## Manhattan plot

# %%
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "#CHROM",
  bp = "POS",
  p = "P",
  snp = "ID",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown]
# ## QQ-plot

# %%
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$P,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% [markdown] tags=[]
# # Random pheno 1

# %%
pheno_code <- "1"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv.gz"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %%
gwas <- gwas %>% filter(P >= 0 & P <= 1)

# %%
dim(gwas)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(P, na.rm = TRUE),
    sd_value   = sd(P, na.rm = TRUE),
    min_value  = min(P, na.rm = TRUE),
    max_value  = max(P, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(P)) %>%
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
  mutate(LOGP = -log10(P)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(BETA, na.rm = TRUE),
    sd_value   = sd(BETA, na.rm = TRUE),
    min_value  = min(BETA, na.rm = TRUE),
    max_value  = max(BETA, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(SE, na.rm = TRUE),
    sd_value   = sd(SE, na.rm = TRUE),
    min_value  = min(SE, na.rm = TRUE),
    max_value  = max(SE, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown]
# ## Manhattan plot

# %%
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "#CHROM",
  bp = "POS",
  p = "P",
  snp = "ID",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown]
# ## QQ-plot

# %%
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$P,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% [markdown] tags=[]
# # Random pheno 2

# %%
pheno_code <- "2"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv.gz"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %%
gwas <- gwas %>% filter(P >= 0 & P <= 1)

# %%
dim(gwas)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(P, na.rm = TRUE),
    sd_value   = sd(P, na.rm = TRUE),
    min_value  = min(P, na.rm = TRUE),
    max_value  = max(P, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(P)) %>%
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
  mutate(LOGP = -log10(P)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(BETA, na.rm = TRUE),
    sd_value   = sd(BETA, na.rm = TRUE),
    min_value  = min(BETA, na.rm = TRUE),
    max_value  = max(BETA, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(SE, na.rm = TRUE),
    sd_value   = sd(SE, na.rm = TRUE),
    min_value  = min(SE, na.rm = TRUE),
    max_value  = max(SE, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown]
# ## Manhattan plot

# %%
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "#CHROM",
  bp = "POS",
  p = "P",
  snp = "ID",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown]
# ## QQ-plot

# %%
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$P,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% [markdown] tags=[]
# # Random pheno 3

# %%
pheno_code <- "3"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv.gz"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %%
gwas <- gwas %>% filter(P >= 0 & P <= 1)

# %%
dim(gwas)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(P, na.rm = TRUE),
    sd_value   = sd(P, na.rm = TRUE),
    min_value  = min(P, na.rm = TRUE),
    max_value  = max(P, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(P)) %>%
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
  mutate(LOGP = -log10(P)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(BETA, na.rm = TRUE),
    sd_value   = sd(BETA, na.rm = TRUE),
    min_value  = min(BETA, na.rm = TRUE),
    max_value  = max(BETA, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(SE, na.rm = TRUE),
    sd_value   = sd(SE, na.rm = TRUE),
    min_value  = min(SE, na.rm = TRUE),
    max_value  = max(SE, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown]
# ## Manhattan plot

# %%
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "#CHROM",
  bp = "POS",
  p = "P",
  snp = "ID",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown]
# ## QQ-plot

# %%
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$P,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %% [markdown] tags=[]
# # Random pheno 4

# %%
pheno_code <- "4"

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(GWAS_DIR, paste0("random.pheno", pheno_code,".glm.linear.tsv.gz"))))

# %% tags=[]
dim(gwas)

# %% tags=[]
head(gwas)

# %%
gwas <- gwas %>% filter(P >= 0 & P <= 1)

# %%
dim(gwas)

# %% [markdown] tags=[]
# ## Data stats

# %%
gwas_p_stats <- gwas %>%
  summarise(
    mean_value = mean(P, na.rm = TRUE),
    sd_value   = sd(P, na.rm = TRUE),
    min_value  = min(P, na.rm = TRUE),
    max_value  = max(P, na.rm = TRUE),
  )

# %%
gwas_p_stats

# %%
gwas_mlog_stats <- gwas %>%
  mutate(LOGP = -log10(P)) %>%
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
  mutate(LOGP = -log10(P)) %>%
  summarise(
    max_no_inf = max(LOGP[is.finite(LOGP)], na.rm = TRUE)
  ) %>%
  pull(max_no_inf)

# %%
gwas_max_mlog_noninf

# %%
gwas_beta_stats <- gwas %>%
  summarise(
    mean_value = mean(BETA, na.rm = TRUE),
    sd_value   = sd(BETA, na.rm = TRUE),
    min_value  = min(BETA, na.rm = TRUE),
    max_value  = max(BETA, na.rm = TRUE),
  )

# %%
gwas_beta_stats

# %%
gwas_se_stats <- gwas %>%
  summarise(
    mean_value = mean(SE, na.rm = TRUE),
    sd_value   = sd(SE, na.rm = TRUE),
    min_value  = min(SE, na.rm = TRUE),
    max_value  = max(SE, na.rm = TRUE),
  )

# %%
gwas_se_stats

# %% [markdown]
# ## Manhattan plot

# %%
options(repr.plot.width = 20, repr.plot.height = 10)

manhattan(
  gwas,
  chr = "#CHROM",
  bp = "POS",
  p = "P",
  snp = "ID",
  main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - Manhattan plot"),
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
  ylim = c(0, max(gwas_max_mlog_noninf+1, 10)),
)

# %% [markdown]
# ## QQ-plot

# %%
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
    gwas$P,
    main = paste0("random pheno ", pheno_code, " (", title_prefix, ")", " - QQ plot of GWAS p-values")
)

# %%
