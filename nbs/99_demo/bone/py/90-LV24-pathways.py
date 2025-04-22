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
# Generates manubot tables for pathways enriched (from the MultiPLIER models) given an LV name (in Settings below).

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
LV_NAME = "LV24"

# %%
# assert (
#     conf.MANUSCRIPT["BASE_DIR"] is not None
# ), "The manuscript directory was not configured"

# OUTPUT_FILE_PATH = conf.MANUSCRIPT["CONTENT_DIR"] / "50.00.supplementary_material.md"
# display(OUTPUT_FILE_PATH)
# assert OUTPUT_FILE_PATH.exists()

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
        #         | (multiplier_model_summary["AUC"] >= 0.75)
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

# %% [markdown]
# ## Split names

# %%
lv_pathways["Pathway"] = lv_pathways["Pathway"].apply(lambda x: " ".join(x.split("_")))

# %%
lv_pathways.head()

# %%
