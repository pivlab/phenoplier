# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     formats: ipynb,py//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] tags=[]
# # Description

# %% [markdown] tags=[]
# It shows the pathways enriched (from the MultiPLIER models) given an LV name (in Settings below).
# These "pathways enriched" are a set of limited pathways used during training of this PLIER model.
# More pathways might be enriched if using external and more comprehensive databases such as gProfiler or FUMA as shown below.

# %% [markdown] tags=[]
# # Modules loading

# %% tags=[]
# %load_ext autoreload
# %autoreload 2

# %% tags=[]
import re
from pathlib import Path

import pandas as pd

from entity import Trait
import conf

# %% [markdown] tags=[]
# # Settings

# %% tags=["parameters"]
LV_NAME = "LV844"

# %% [markdown] tags=[]
# # Paths

# %%
OUTPUT_FIGURES_DIR = Path(conf.RESULTS_DIR, "demo", f"{LV_NAME.lower()}").resolve()
display(OUTPUT_FIGURES_DIR)
OUTPUT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown] tags=[]
# # Load MultiPLIER summary

# %% tags=[]
multiplier_model_summary = pd.read_pickle(conf.MULTIPLIER["MODEL_SUMMARY_FILE"])

# %% tags=[]
multiplier_model_summary.shape

# %% tags=[]
multiplier_model_summary.head()

# %% [markdown]
# # LV pathways

# %%
lv_pathways = multiplier_model_summary[
    multiplier_model_summary["LV index"].isin((LV_NAME[2:],))
    & (
        (multiplier_model_summary["FDR"] < 0.05)
                | (multiplier_model_summary["AUC"] >= 0.75)
    )
]

# %%
lv_pathways.shape

# %%
lv_pathways = lv_pathways[["pathway", "AUC", "FDR"]].sort_values("FDR")

# %%
lv_pathways = lv_pathways.assign(AUC=lv_pathways["AUC"].apply(lambda x: f"{x:.2f}"))

# %%
lv_pathways = lv_pathways.assign(FDR=lv_pathways["FDR"].apply(lambda x: f"{x:.2e}"))

# %%
lv_pathways = lv_pathways.rename(
    columns={
        "pathway": "Pathway",
    }
)

# %%
lv_pathways.head()

# %% [markdown] tags=[]
# # Load LV data

# %%
from data.recount2 import LVAnalysis

# %%
lv_obj = LVAnalysis(LV_NAME)

# %% [markdown] tags=[]
# Here I show the top 20 genes for our LV. You can see gene symbols, the LV weight (in column `LV603`) and the cytoband.

# %%
lv_obj.lv_genes.head(20)

# %% [markdown]
# # Pathway enrichment using external databases

# %% [markdown]
# ## gProfiler

# %%
print(" ".join(lv_obj.lv_genes.head(70)["gene_name"].tolist()))

# %% [markdown]
# Copy/paste the list of genes above and use gProfiler: https://biit.cs.ut.ee/gprofiler/gost

# %% [markdown]
# Results URL: 
#
# **Notes**: TODO

# %% [markdown]
# ## FUMA

# %%
# print top genes in module
print("\n".join(lv_obj.lv_genes.head(70)["gene_name"].tolist()))

# %%
# save all genes in model to use as background list of genes
lv_obj.lv_genes["gene_name"].to_csv(OUTPUT_FIGURES_DIR / "all_genes.txt", header=None, index=False)

# %%
OUTPUT_FIGURES_DIR / "all_genes.txt"

# %%
# !head /opt/data/results/demo/lv24/all_genes.txt

# %%
# !wc -l /opt/data/results/demo/lv24/all_genes.txt

# %% [markdown]
# Now go to the FUMA GENE2FUNC module here: https://fuma.ctglab.nl/gene2func
#
# 1. Paste the list of top genes above and then upload the `all_genes.txt` file.
# 2. Use a "Title" and click on "Submit"

# %% [markdown]
# **Notes:** significantly expressed in two brain tissues and testis. The brain one makes sense because it's the same cell type.

# %%
