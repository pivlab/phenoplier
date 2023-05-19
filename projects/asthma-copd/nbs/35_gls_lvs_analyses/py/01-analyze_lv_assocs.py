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

# %% [markdown]
# # Summary of results

# %% [markdown]
# * LV101:
#   * asthma and ACO
#   * very clear: tonsil innate lymphoid cells (ILC3), NK, ILC2
#   * PhenomeXcan:
#     * hayfever, allergic rhinitis, eczema, asthma, eosinophils
#   * eMERGE:
#     * None
#
# * LV70:
#   * asthma and COPD
#   * Tissues:
#     * It's very hard to find a single tissue for this LV
#     * at the top are: blastocysts, keratinocytes, lung adenocarcinoma, etc
#     * Human skeletal muscle myoblasts (RA)
#       * there is a lot of variation here, but it is not possible to sepearate after reading the paper.
#   * PhenomeXcan:
#     * rheumatoid arthritis (RA), basophil/neutrophil count, neutrophil count, medication: methotrexate (cancer, RA)
#   * eMERGE:
#     * Emphysema (496.1), diseases of white blood cells (288)
#
# * LV948:
#   * asthma and ACO
#   * very similar to LV101
#   * PhenomeXcan:
#     * very similar to LV101 + Cronh's disease + inflammatory bowel disease + ulcerative colitis
#   * eMERGE:
#     * None
#
# * LV705:
#   * asthma and ACO
#   * very clear T cell signature when stimulated and also tonsil innate lymphoid cells (NK) (no ILC as in LV101)
#   * PhenomeXcan:
#     * very similar to LV948
#   * eMERGE:
#     * None
#
# * LV504:
#   * COPD and ACO
#   * lung, but very weak signal
#   * PhenomeXcan:
#     * number of cigarettes currently smoked daily
#   * eMERGE:
#     * ahterosclerosis of the extremeties (440.*)
#
# * LV207:
#   * asthma
#   * very strong breast cancer signal
#   * PhenomeXcan:
#     * neutrophil count and percentage, age of asthma diagnosed, white blood cell count
#   * eMERGE:
#     * None
#
# * LV61:
#   * asthma (ACO has FDR of 0.60 and unadjusted p-value of 0.026).
#   * whole blood
#   * PhenomeXcan:
#     * basophil count, medication: methotrexate, RA, Adult-onset Still disease, asthma
#   * eMERGE:
#     * None
#
# * LV149:
#   * COPD (VERY INTERESTING!)
#   * mostly expressed in embryonic kidney cell lines with single-site and multi-site mutations (something like that)
#   * PhenomeXcan:
#     * blood clot in the leg (DVT, deep venous thrombosis, number of cigarettes, emphysema
#   * eMERGE:
#     * Coagulation defects (286), other venous embolism and thrombosis (452)
#
# * LV844:
#   * asthma (ACO has FDR of 0.27 and unadjusted p-value of 0.006). 
#   * TISSUES
#   * PhenomeXcan:
#     * polymyalgia rheumatica, type 1 diabetes, age diabetes diagnosed, RA, insulin medication, coeliac disease
#   * eMERGE:
#     * None
#
# * LV803:
#   * asthma (COPD has unadjusted p-value of 0.053)
#   * tonsil (similar to LV101) + brain cell types
#   * PhenomeXcan:
#     * RA, HDL cholesterol, Adult-onset Still disease, keratometry measurements
#   * eMERGE:
#     * None
#
# * LV799:
#   * COPD (ACO unadjusted p-value of 0.043)
#   * Blood, peripheral blood, multiple myeloma
#   * PhenomeXcan:
#     * cigarettes smoked, lung cancer (father), Pericarditis
#   * eMERGE:
#     * None
#
# * LV17:
#   * asthma (ACO has unadjusted p-value of 0.048, and COPD 0.43)
#   * tonsil, simillar to LV101, but here NK are ate the top, and then ILC1 and ILC3, ILC2 comes very behind
#   * PhenomeXcan:
#     * wheeze, eosinophil, asthma, RA, hayfever/allergic rhinitis
#   * eMERGE:
#     * None
#
# * LV180:
#   * COPD (asthma and ACO pvalue > 0.30)
#   * Glima + lung fibrobasts cells, hard to see
#   * PhenomeXcan:
#     * smoking, number of cigarettes
#   * eMERGE:
#     * None
#
# * LV696:
#   * COPD, rest has pvalue > 0.40
#   * glioma + breast and lung cells, some tonsil ILC3 ILC2 ILC1, hard to see
#   * PhenomeXcan:
#     * smoking
#   * eMERGE:
#     * Atherosclerosis of the extremities
#
# * LV563:
#   * ACO (asthma has pvalue <0.02)
#   * tonsil, but weird plot, with outliers, some astrocytes, fetal neocortex and other brain cell types
#   * PhenomeXcan:
#     * age asthma diagnosed, height
#   * eMERGE:
#     * Disorders of refraction and accommodation; blindness and low vision
#
# * LV214:
#   * COPD (rest has pvalue > 0.30)
#   * Glioma + breast cancer, tonsil ILC3 NK ILC1 ILC2
#   * PhenomeXcan:
#     * smoking
#   * eMERGE:
#     * None
#
# * LVXXX:
#   * ASSOCS
#   * TISSUES
#   * PhenomeXcan:
#     * ASSOCS
#   * eMERGE:
#     * ASSOC

# %%
