# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     formats: ipynb,py//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] tags=[]
# # Description

# %% [markdown] tags=[]
# It lists the most significant LV-trait associations for asthma, COPD and ACO.

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from IPython.display import HTML

import conf

# %% [markdown] tags=[]
# # Settings

# %% tags=[]
INPUT_DIR = conf.PROJECTS["ASTHMA_COPD"]["RESULTS_DIR"] / "gls_phenoplier"
display(INPUT_DIR)
assert INPUT_DIR.exists()

# %% [markdown] tags=[]
# # Load data

# %% [markdown] tags=[]
# ## MultiPLIER summary

# %% tags=[]
multiplier_model_summary = pd.read_pickle(conf.MULTIPLIER["MODEL_SUMMARY_FILE"])

# %% tags=[]
multiplier_model_summary.shape

# %% tags=[]
multiplier_model_summary.head()

# %% [markdown] tags=[]
# ## GLS associations for Asthma/COPD

# %%
input_filepath = INPUT_DIR / "gls-summary.pkl.gz"
assert input_filepath.exists()

# %%
gls_asthma_copd = pd.read_pickle(input_filepath)

# %%
gls_asthma_copd.shape

# %%
gls_asthma_copd.head()

# %% [markdown] tags=[]
# ## GLS associations for PhenomeXcan

# %%
input_filepath = INPUT_DIR / "gls-summary-phenomexcan.pkl.gz"
assert input_filepath.exists()

# %%
gls_phenomexcan = pd.read_pickle(input_filepath)

# %%
gls_phenomexcan.shape

# %%
gls_phenomexcan.head()

# %% [markdown] tags=[]
# ## GLS associations for eMERGE

# %%
input_filepath = INPUT_DIR / "gls-summary-emerge.pkl.gz"
assert input_filepath.exists()

# %%
gls_emerge = pd.read_pickle(input_filepath)

# %%
gls_emerge.shape

# %%
gls_emerge.head()

# %% [markdown] tags=[]
# # Top hits

# %% tags=[]
with pd.option_context("display.max_columns", None, "display.max_colwidth", None):
    signif_lv_assocs = gls_asthma_copd.sort_values("fdr")  # .drop(columns="phenotype")
    signif_lv_assocs = signif_lv_assocs[signif_lv_assocs["fdr"] < 0.05]

    # convert back "category" data types into str/object
    signif_lv_assocs["lv"] = signif_lv_assocs["lv"].astype(str)
    signif_lv_assocs["phenotype"] = signif_lv_assocs["phenotype"].astype(str)

    display(signif_lv_assocs.shape)
    display(signif_lv_assocs.head(50))

# %% [markdown] tags=[]
# ## Hits per phenotype

# %% tags=[]
signif_lv_assocs["phenotype"].value_counts()

# %% tags=[]
signif_lv_assocs.groupby("phenotype")["fdr"].describe()

# %% tags=[]
for idx, grp in signif_lv_assocs.groupby("phenotype"):
    display(HTML(f"<h3>{idx}</h3>"))
    display(grp)

# %% [markdown] tags=[]
# ## Hits per LV

# %% tags=[]
lv_assocs_by_size = signif_lv_assocs.groupby("lv").size().sort_values(ascending=False)
display(lv_assocs_by_size)

# %%
lv_assocs_by_size.shape

# %%

# %%
lv_assocs_ids = set(map(lambda x: x.split("LV")[1], lv_assocs_by_size.index))

# %%
len(lv_assocs_ids)

# %%
display(list(lv_assocs_ids)[:5])

# %%
lv_pathways = multiplier_model_summary[
    multiplier_model_summary["LV index"].isin(lv_assocs_ids)
    & (
        (multiplier_model_summary["FDR"] < 0.05)
        #         | (multiplier_model_summary["AUC"] >= 0.75)
    )
]

# %%
lv_pathways.shape

# %%
lv_pathways.sort_values("LV index")

# %%
with pd.option_context(
    "display.max_rows", None, "display.max_columns", None, "display.max_colwidth", None
):
    _tmp = lv_pathways.groupby("LV index").apply(lambda x: "; ".join(x.sort_values("p-value")["pathway"].tolist()))
    display(_tmp.to_frame())

# %%
