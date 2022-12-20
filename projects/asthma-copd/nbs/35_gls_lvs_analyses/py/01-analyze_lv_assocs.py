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

# %% [markdown]
# # Take top LVs and show all associations

# %% [markdown]
# Take the top associated LVs (sorted by number of associations with asthma/COPD traits), and then show all the associations with the traits.
#
# This allows to see whether one LV is strongly associated with some traits and not associated at all with the others.
# For example, it could happen that one LVs is strongly associated with asthma and COPD using some threshold (pvalue/FDR < 0.05), but just slightly above that threshold for ACO (pvalue = 0.06). We would not trust much the "no association" with ACO in this case. However, if the p-value for ACO would be 0.76, we'd have more evidence that it is not associated, and thus the LV is more specific to asthma and COPD.

# %%
for lv_code in lv_assocs_by_size.index:
    display(HTML(f"<h2>{lv_code}</h2>"))

    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.max_colwidth",
        None,
    ):
        display(HTML(f"<h3>Asthma/COPD/ACO</h3>"))
        _tmp = gls_asthma_copd[gls_asthma_copd["lv"] == lv_code].sort_values("pvalue")
        display(_tmp)

        display(HTML(f"<h3>PhenomeXcan</h3>"))
        _tmp = gls_phenomexcan[gls_phenomexcan["lv"] == lv_code].sort_values("pvalue")
        _tmp = _tmp[_tmp["fdr"] < 0.05]
        display(_tmp)

        display(HTML(f"<h3>eMERGE</h3>"))
        _tmp = gls_emerge[gls_emerge["lv"] == lv_code].sort_values("pvalue")
        _tmp = _tmp[_tmp["fdr"] < 0.05]
        display(_tmp)

# %%
