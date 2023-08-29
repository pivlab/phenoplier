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
# It takes a GWAS that was imputed and postprocessed (using the PrediXcan scripts here https://github.com/hakyimlab/summary-gwas-imputation) and verifies that the Manhattan and QQ-plots look fine (without inflation).

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
PHENOPLIER_NOTEBOOK_FILEPATH <- "projects/chronotype/nbs/10_gwas_harmonization/15-gwas-qqplot.ipynb"


# %% [markdown] tags=[]
# # Paths

# %% tags=[]
DATA_DIR <- Sys.getenv(paste0("PHENOPLIER_PROJECTS_", PROJECTS_TRAIT_KEY, "_RESULTS_DIR"))

# %% tags=[]
DATA_DIR

# %% tags=[]
INPUT_GWAS_DIR <- file.path(DATA_DIR, "final_imputed_gwas")

# %% tags=[]
INPUT_GWAS_DIR

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
# ## Get GWAS file

# %% tags=[]
row <- which(trait_info$id == "chronotype")
gwas_filename <- trait_info[[row, "gwas_file"]]

# %% tags=[]
gwas_filename

# %% tags=[]
gwas_filename <- tools::file_path_sans_ext(gwas_filename)

# %% [markdown] tags=[]
# ## Load data

# %% tags=[]
gwas <- as.data.frame(read_table(file.path(INPUT_GWAS_DIR, paste0(gwas_filename, "-harmonized-imputed.txt.gz"))))

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
# ## Stats

# %% tags=[]
summary(gwas)

# %% [markdown] tags=[]
# ## Remove NA pvalues

# %% tags=[]
dim(gwas)

# %% tags=[]
gwas <- gwas %>% filter(pvalue >= 0 & pvalue <= 1)

# %% tags=[]
dim(gwas)

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
  main = gwas_title,
  suggestiveline = F,
  genomewideline = -log10(5e-08),
  cex = 0.6,
  cex.axis = 0.9,
)

# %% [markdown] tags=[]
# ## QQ-plot

# %% tags=[]
options(repr.plot.width = 10, repr.plot.height = 10)

qq(
  gwas$pvalue,
  main = paste0("Q-Q plot - ", gwas_title)
)

# %% tags=[]
