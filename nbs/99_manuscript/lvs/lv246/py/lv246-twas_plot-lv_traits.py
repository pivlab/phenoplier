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

# %% [markdown]
# # Description

# %% [markdown]
# COMPLETE
# Generates a plot with TWAS associations using the top genes and traits for LV246 (related to the hypertension cluster).

# %% [markdown]
# # Modules loading

# %%
# %load_ext autoreload
# %autoreload 2

# %%
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import conf
from data.cache import read_data
from utils import generate_result_set_name
from data.recount2 import LVAnalysis

# %% [markdown]
# # Settings

# %%
LV_NUMBER_SELECTED = 246
LV_NAME_SELECTED = f"LV{LV_NUMBER_SELECTED}"
display(LV_NAME_SELECTED)

# %%
OUTPUT_FIGURES_DIR = Path(
    conf.MANUSCRIPT["FIGURES_DIR"], "lvs_analysis", f"lv{LV_NUMBER_SELECTED}"
).resolve()
display(OUTPUT_FIGURES_DIR)
OUTPUT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Data loading

# %% [markdown]
# ## MultiPLIER summary

# %%
multiplier_model_summary = read_data(conf.MULTIPLIER["MODEL_SUMMARY_FILE"])

# %%
multiplier_model_summary.shape

# %%
multiplier_model_summary.head()

# %% [markdown]
# ## PhenomeXcan results

# %% [markdown]
# ### S-MultiXcan

# %%
smultixcan_results_filename = conf.PHENOMEXCAN[
    "SMULTIXCAN_EFO_PARTIAL_MASHR_PVALUES_FILE"
]
display(smultixcan_results_filename)

# %%
smultixcan_results = pd.read_pickle(smultixcan_results_filename)

# %%
smultixcan_results.shape

# %%
smultixcan_results.head()

# %% [markdown]
# ### fastENLOC

# %%
fastenloc_results_filename = conf.PHENOMEXCAN["FASTENLOC_EFO_PARTIAL_TORUS_RCP_FILE"]
display(fastenloc_results_filename)

# %%
fastenloc_results = pd.read_pickle(fastenloc_results_filename)

# %%
fastenloc_results.shape

# %%
fastenloc_results.head()

# %% [markdown]
# ## S-MultiXcan projection (`z_score_std`)

# %% tags=[]
INPUT_SUBSET = "z_score_std"

# %% tags=[]
INPUT_STEM = "projection-smultixcan-efo_partial-mashr-zscores"

# %% tags=[]
input_filepath = Path(
    conf.RESULTS["DATA_TRANSFORMATIONS_DIR"],
    INPUT_SUBSET,
    f"{INPUT_SUBSET}-{INPUT_STEM}.pkl",
).resolve()
display(input_filepath)

assert input_filepath.exists(), "Input file does not exist"

input_filepath_stem = input_filepath.stem
display(input_filepath_stem)

# %% tags=[]
data = read_data(input_filepath)

# %% tags=[]
data.shape

# %% tags=[]
data.head()

# %% [markdown]
# ## PhenomeXcan LV-trait associations

# %%
input_filepath = Path(conf.RESULTS["GLS"] / "gls-summary-phenomexcan.pkl.gz")
display(input_filepath)

# %%
phenomexcan_lv_trait_assocs = pd.read_pickle(input_filepath)

# %%
phenomexcan_lv_trait_assocs.shape

# %%
phenomexcan_lv_trait_assocs.head()

# %% [markdown]
# ## LV246-trait associations

# %%
from traits import SHORT_TRAIT_NAMES

# %%
lv_assocs = phenomexcan_lv_trait_assocs[
    (phenomexcan_lv_trait_assocs["lv"] == LV_NAME_SELECTED)
    & (phenomexcan_lv_trait_assocs["fdr"] < 0.05)
].sort_values("fdr")

# %%
lv_assocs.shape

# %%
with pd.option_context(
    "display.max_rows", None, "display.max_columns", None, "display.max_colwidth", None
):
    display(lv_assocs)

# %% [markdown]
# # LV analysis

# %%
lv_obj = lv_exp = LVAnalysis(LV_NAME_SELECTED, data)

# %%
lv_gene_sets = multiplier_model_summary[
    multiplier_model_summary["LV index"].isin((str(LV_NUMBER_SELECTED),))
    & (
        (multiplier_model_summary["FDR"] < 0.05)
        | (multiplier_model_summary["AUC"] >= 0.75)
    )
]
display(lv_gene_sets)

# %% [markdown]
# ## Traits

# %%
lv_obj.lv_traits.shape

# %%
lv_obj.lv_traits.head(20)

# %%
lv_traits = lv_obj.lv_traits.copy()

# %% [markdown]
# ## Genes

# %%
lv_obj.lv_genes.shape

# %%
lv_obj.lv_genes.head(10)

# %%
lv_genes = (
    lv_obj.lv_genes[["gene_name", LV_NAME_SELECTED]].set_index("gene_name").squeeze()
)

# %%
lv_genes

# %%
" ".join(lv_genes.index[:20])

# %% [markdown]
# # Prepare data for plot

# %%
from entity import Gene, Trait

# %% [markdown]
# ## Select top genes

# %%
lv_top_genes = lv_genes.head(74).rename(index=Gene.GENE_NAME_TO_ID_MAP)

# %%
lv_top_genes

# %%
# remove genes not present in S-MultiXcan
lv_top_genes = lv_top_genes.loc[
    [g for g in lv_top_genes.index if g in Gene.GENE_ID_TO_NAME_MAP]
]
# assert lv_top_genes.shape[0] == 30, lv_top_genes.shape

# %% [markdown]
# ## Select top traits

# %%
lv_assocs.head()

# %%
lv_assocs["phenotype_full_code"] = lv_assocs["phenotype"].apply(lambda x: Trait.get_trait(code=x).full_code)

# %%
lv_assocs

# %%
smultixcan_results.columns[pd.Series(smultixcan_results.columns).str.lower().str.contains("alz")]

# %%
Trait.get_trait(full_code="F5_DEMENTIA-Dementia").n_cases

# %%
Trait.get_traits_from_efo("alzheimer's disease")[0].n

# %%
Trait.get_traits_from_efo("alzheimer's disease")[0].n_cases

# %%
Trait.get_traits_from_efo("alzheimer's disease")[0].full_code

# %%
smultixcan_results.columns[pd.Series(smultixcan_results.columns).str.lower().str.contains("cholest")]

# %%
Trait.get_traits_from_efo("hypercholesterolemia")[0].n

# %%
Trait.get_traits_from_efo("hypercholesterolemia")[0].n_cases

# %%
Trait.get_traits_from_efo("hypercholesterolemia")[0].full_code

# %%
lv_assocs = lv_assocs.replace({
    "phenotype_full_code": {
        "IGAP_Alzheimer": "alzheimer's disease",
        "20002_1473-Noncancer_illness_code_selfreported_high_cholesterol": "hypercholesterolemia",
    }
})

# %%
lv_top_traits = lv_assocs[["phenotype_full_code", "fdr"]].set_index("phenotype_full_code").squeeze()

# %%
lv_top_traits

# %%
lv_top_traits.shape

# %%
lv_top_traits.head()

# %% [markdown]
# ## Subset S-MultiXcan results

# %%
data_subset = smultixcan_results.loc[lv_top_genes.index, lv_top_traits.index.to_list()].T

# %%
data_subset = data_subset.apply(lambda x: -np.log10(x))

# %%
data_subset.shape

# %%
data_subset = data_subset.rename(columns=Gene.GENE_ID_TO_NAME_MAP)

# %%
data_subset

# %%
# [START] remove duplicated traits / leave the most important/interesting ones only

# %%
# data_subset.index[pd.Series(data_subset.index).str.lower().str.contains("cholest")]

# %%
_trait_code = "F5_DEMENTIA-Dementia"
_trait_obj = Trait.get_traits_from_efo(_trait_code)
if _trait_obj is None:
    _trait_obj = Trait.get_trait(full_code=_trait_code)
else:
    _trait_obj = _trait_obj[0]
assert _trait_obj is not None
display(_trait_obj.n)
display(_trait_obj.n_cases)

# %%
# smultixcan_results.columns[pd.Series(smultixcan_results.columns).str.lower().str.contains("cholest")]

# %%
# data_subset.loc[[
#     "MAGNETIC_IDL.TG",
#     "MAGNETIC_LDL.C",
#     "hypercholesterolemia",
#     "6153_1-Medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_Cholesterol_lowering_medication",
#     "20003_1140861958-Treatmentmedication_code_simvastatin",
#     "",
    
#     "MAGNETIC_CH2.DB.ratio",
#     "6153_1-Medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_Cholesterol_lowering_medication",
#     "",
#     "20003_1141146234-Treatmentmedication_code_atorvastatin",
# ]]

# %%
# I remove some because of space limits
data_subset = data_subset.drop([
    "2654_2-Nonbutter_spread_type_details_Flora_ProActive_or_Benecol",
    "20107_100-Illnesses_of_father_None_of_the_above_group_1",
    "6153_100-Medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_None_of_the_above",
    "6177_1-Medication_for_cholesterol_blood_pressure_or_diabetes_Cholesterol_lowering_medication",
    "1835-Mother_still_alive",
    "22617_2451-Job_SOC_coding_Librarians",
    "AD-Alzheimers_disease",
    "20110_10-Illnesses_of_mother_Alzheimers_diseasedementia",
    "20107_10-Illnesses_of_father_Alzheimers_diseasedementia",
    "20111_10-Illnesses_of_siblings_Alzheimers_diseasedementia",
    "20107_1-Illnesses_of_father_Heart_disease",
    "F5_DEMENTIA-Dementia",
    "20003_1140861958-Treatmentmedication_code_simvastatin",
    "20003_1141146234-Treatmentmedication_code_atorvastatin",
    "20003_1141192410-Treatmentmedication_code_rosuvastatin",
    "20003_1141146138-Treatmentmedication_code_lipitor_10mg_tablet",
])

# %%
# data_subset.shape

# %%
# data_subset.head()

# %%
# [END] remove duplicated traits / leave the most important/interesting ones only

# %%
data_subset["trait_description"] = data_subset.apply(
    lambda x: Trait.get_trait(full_code=x.name).description
    if not Trait.is_efo_label(x.name)
    else x.name,
    axis=1,
)

# %%
# no need to do this (all traits are different)
# calculate mean of traits with the same description
# data_subset = data_subset.groupby("trait_description").mean()

# %%
# data_subset.shape

# %%
data_subset.head()

# %%
plot_data = data_subset.reset_index(drop=True).set_index("trait_description")

# %%
plot_data

# %%
_old_columns = plot_data.index
display(_old_columns.shape)
display(_old_columns)

# %%
# reorder columns
# plot_data = plot_data.loc[
#     [
#         "",
#         "",
#         "",
#         "",
#     ]
# ]

# %%
assert set(plot_data.index) == set(_old_columns)

# %% [markdown]
# ## Subset fastENLOC results

# %%
data_subset2 = (
    fastenloc_results.loc[lv_top_genes.index, lv_top_traits.index].fillna(0.0).T * 100.0
)

# %%
data_subset2.shape

# %%
data_subset2 = data_subset2.rename(columns=Gene.GENE_ID_TO_NAME_MAP)

# %%
data_subset2["trait_description"] = data_subset2.apply(
    lambda x: Trait.get_trait(full_code=x.name).description
    if not Trait.is_efo_label(x.name)
    else x.name,
    axis=1,
)

# %%
# calculate mean of traits with the same description
data_subset2 = data_subset2.groupby("trait_description").mean()
# _tmp = _tmp.assign(cluster_name=_tmp.apply(lambda x: _cluster_column[x.name], axis=1))

# %%
data_subset2.shape

# %%
data_subset2

# %%
data_subset2 = data_subset2.loc[plot_data.index, plot_data.columns]

# %%
data_subset2.shape

# %%
data_subset2.head()

# %%
plot_data2 = data_subset2  # .reset_index().set_index("trait_description")

# %%
plot_data2.shape

# %%
pd.Series(plot_data2.values.flatten()).quantile(
    [0.01, 0.05, 0.10, 0.15, 0.5, 0.55, 0.60, 0.65, 0.75, 0.80, 0.90, 0.95, 0.97, 0.98, 0.99]
)

# %%
pd.Series(plot_data2.values.flatten()).describe()

# %% [markdown]
# ### Filter some genes

# %%
_smultixcan_tmp = plot_data.sum()
_fastenloc_tmp = plot_data2.sum()

# %%
_smultixcan_tmp.shape

# %%
_smultixcan_tmp.describe()

# %%
_fastenloc_tmp.shape

# %%
_fastenloc_tmp.describe()

# %%
_tmp = (_smultixcan_tmp > 5.57) & (_fastenloc_tmp > 2.44)

# %%
_tmp

# %%
_tmp.sum()

# %%
_tmp.head(20)

# %%
_tmp.loc[["DGAT2", "ACACA"]] = True

# %%
selected_genes = _tmp[_tmp].index.to_list()

# %%
selected_genes

# %%
# select genes
plot_data = plot_data[selected_genes]
display(plot_data.shape)

plot_data2 = plot_data2[selected_genes]
display(plot_data2.shape)

# %% [markdown]
# ### Create fastENLOC heatmap data

# %%
fastenloc_heatmap_data = plot_data2.T  # .drop(columns=["cluster_name", "color"]).T

# %%
fastenloc_heatmap_data.shape

# %%
fastenloc_heatmap_data[fastenloc_heatmap_data > 100.0] = 100.0

# %%
fastenloc_heatmap_data.head()

# %%
fastenloc_heatmap_data = (
    fastenloc_heatmap_data.reset_index(drop=True).T.reset_index(drop=True).T
)

# %%
fastenloc_heatmap_data = (
    fastenloc_heatmap_data.rename_axis(index="genes", columns="traits")
    .unstack()
    .rename("rcp")
    .reset_index()
)

# %%
fastenloc_heatmap_data.head()

# %%
_tmp_heatmap_multixcan2 = plot_data.loc[plot_data2.index, plot_data2.columns].T

# %%
fastenloc_heatmap_data = fastenloc_heatmap_data.assign(
    color=fastenloc_heatmap_data.apply(
        lambda x: "white"
        if _tmp_heatmap_multixcan2.iloc[int(x.genes), int(x.traits)] >= 8
        else "black",
        axis=1,
    )
)

# %%
fastenloc_heatmap_data = fastenloc_heatmap_data.assign(
    marker_size=fastenloc_heatmap_data.apply(
        lambda x: x.rcp if x.rcp >= 1.0 else 0.0, axis=1
    )
)

# %%
n_power = 1.30

fastenloc_heatmap_data["marker_size"] = np.power(
    fastenloc_heatmap_data["marker_size"], n_power
)

# %%
fastenloc_heatmap_data["traits"] += 0.5
fastenloc_heatmap_data["genes"] += 0.5

# %%
fastenloc_heatmap_data  # .groupby('color').count()

# %% [markdown]
# # Plot

# %%
_trait_renames = {
    # "CH2DB NMR": "CH2DB (lipids)",
    "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones: Cholesterol lowering medication": "Cholesterol medication",
    # "Treatment/medication code: simvastatin": "Medication: simvastatin",
    # "Treatment/medication code: atorvastatin": "Medication: atorvastatin",
    # "Treatment/medication code: lipitor 10mg tablet": "Medication: lipitor",
    # "Treatment/medication code: rosuvastatin": "Medication: rosuvastatin",
    "hypercholesterolemia": "High-cholesterol",
    "alzheimer's disease": "Alzheimer's disease",
}

# %%
plot_data = _tmp_heatmap_multixcan2.rename(columns=_trait_renames)

# %%
lipids_traits = [
    "LDL Cholesterol NMR",
    "Triglycerides NMR",
    # "CH2DB NMR",
]

diseases_traits = [
    "High-cholesterol",
    "Alzheimer's disease",
    "Any dementia",
]

medications_traits = [
    "Cholesterol medication",
    # "Medication: simvastatin",
    # "Medication: atorvastatin",
    # "Medication: lipitor",
    # "Medication: rosuvastatin",
]

# %%
plot_data = plot_data.loc[:, lipids_traits + diseases_traits + medications_traits]
display(plot_data.shape)

# %%
rc = {
    #     "font.size": 9,
    #     "xtick.labelsize": 10,
    #     "ytick.labelsize": 19,
}

with sns.plotting_context("paper", font_scale=2.50, rc=rc):
    g = sns.clustermap(
        data=plot_data,
        vmin=0.0,
        vmax=10.0,
        row_cluster=False,
        col_cluster=False,
        dendrogram_ratio=0.20,
        xticklabels=True,
        yticklabels=True,
        figsize=(8, 14),
        linewidths=0.25,
        cmap="BuGn",
        cbar_pos=(0.95, 0.64, 0.05, 0.10),
    )

    g.ax_heatmap.set_xlabel(None)
    g.ax_heatmap.set_ylabel(None)

    g.ax_heatmap.get_xaxis().set_ticklabels(
        g.ax_heatmap.get_xaxis().get_ticklabels(),
        rotation=45,
        horizontalalignment="right",
    )

    for l in g.ax_heatmap.get_yaxis().get_ticklabels():
        l.set_style("italic")
        if l.get_text() in ("DGAT2", "ACACA"):
            l.set_fontweight("bold")

    for l in g.ax_heatmap.get_xaxis().get_ticklabels():
        if l.get_text() in lipids_traits:
            l.set_color("black")
        elif l.get_text() in diseases_traits:
            l.set_color("gray")
        elif l.get_text() in medications_traits:
            l.set_color("black")
            

    g.ax_heatmap.scatter(
        fastenloc_heatmap_data["traits"].tolist(),
        fastenloc_heatmap_data["genes"].tolist(),
        marker=".",
        s=fastenloc_heatmap_data["marker_size"].tolist(),
        color=fastenloc_heatmap_data["color"].tolist(),
    )

    for rcp in [10, 25, 50, 75, 100]:
        plt.scatter(
            [],
            [],
            marker=".",
            c="k",
            alpha=0.3,
            s=np.power(rcp, n_power),
            label=str(rcp) + "%",
        )
    leg = plt.legend(
        scatterpoints=1,
        frameon=False,
        labelspacing=1,
        loc="center right",
        bbox_to_anchor=(1.90, -2.0, 0.5, 0.5),
        title="RCP.",
        title_fontsize=18,  # rc['ytick.labelsize'],
        fontsize=14,  # rc['ytick.labelsize'],
    )
    plt.setp(leg.get_title(), multialignment="center")

    g.ax_cbar.set_title("MultiXcan\n-log($p$)", fontdict={"fontsize": 18})
    g.ax_cbar.tick_params(labelsize=14)

    output_filepath = OUTPUT_FIGURES_DIR / f"lv{LV_NUMBER_SELECTED}-lv_traits-twas_plot.svg"
    display(output_filepath)
    plt.savefig(
        output_filepath,
        #         dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )

# %% [markdown]
# # Show p-values for some genes

# %%
_tmp = plot_data.loc["DGAT2"].sort_values()
_tmp = np.power(10, -_tmp)
display(_tmp)

# %%
_tmp = plot_data.loc["ACACA"].sort_values()
_tmp = np.power(10, -_tmp)
display(_tmp)

# %% [markdown]
# # Show RCP for some genes

# %%
_tmp = plot_data2["DGAT2"].sort_values()
_tmp = np.power(10, -_tmp)
display(_tmp)

# %%
_tmp = plot_data2["ACACA"].sort_values()
_tmp = np.power(10, -_tmp)
display(_tmp)

# %%
