# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     formats: ipynb,py//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown] tags=[]
# # Environment variables

# %% tags=[]
# %load_ext autoreload
# %autoreload 2

# %% tags=[]
import conf

# %% tags=[]
N_JOBS = conf.GENERAL["N_JOBS"]
display(N_JOBS)

# %% tags=[]
# %env MKL_NUM_THREADS=$N_JOBS
# %env OPEN_BLAS_NUM_THREADS=$N_JOBS
# %env NUMEXPR_NUM_THREADS=$N_JOBS
# %env OMP_NUM_THREADS=$N_JOBS

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
from pathlib import Path

import statsmodels.api as sm
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from tqdm import tqdm

from gls import GLSPhenoplier

# %% [markdown] tags=[]
# # Settings

# %% tags=[]
OUTPUT_DIR = conf.RESULTS["GLS"]
display(OUTPUT_DIR)

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# %% [markdown] tags=[]
# # Load data

# %% [markdown] tags=[]
# ## PhenomeXcan (S-MultiXcan)

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

# %% tags=[]
data = pd.read_pickle(input_filepath)

# %% tags=[]
data.shape

# %% tags=[]
data.head()

# %% [markdown] tags=[]
# ## Clustering results

# %% tags=[]
CONSENSUS_CLUSTERING_DIR = Path(
    conf.RESULTS["CLUSTERING_DIR"], "consensus_clustering"
).resolve()

display(CONSENSUS_CLUSTERING_DIR)

# %% tags=[]
input_file = Path(CONSENSUS_CLUSTERING_DIR, "best_partitions_by_k.pkl").resolve()
display(input_file)

# %% tags=[]
best_partitions = pd.read_pickle(input_file)

# %% tags=[]
best_partitions.shape

# %% tags=[]
best_partitions.head()

# %% [markdown] tags=[]
# ## MultiPLIER summary

# %% tags=[]
multiplier_model_summary = pd.read_pickle(conf.MULTIPLIER["MODEL_SUMMARY_FILE"])

# %% tags=[]
multiplier_model_summary.shape

# %% tags=[]
multiplier_model_summary.head()

# %% tags=[]
well_aligned_lvs = multiplier_model_summary[
    (multiplier_model_summary["FDR"] < 0.05) | (multiplier_model_summary["AUC"] >= 0.75)
]

display(well_aligned_lvs.shape)
display(well_aligned_lvs.head())

# %% tags=[]
well_aligned_lv_codes = set([f"LV{lvi}" for lvi in well_aligned_lvs["LV index"]])

# %% tags=[]
len(well_aligned_lv_codes)

# %% tags=[]
list(well_aligned_lv_codes)[:5]

# %% [markdown] tags=[]
# # Select partition / cluster pairs

# %% tags=[]
# key: a tuple (partition_k, cluster_id)
# value: a list of tuples
# this specifies
PHENOTYPES_LVS_CONFIG = {
    # cardiovascular
    (29, 14): [(16, 15)],
    (29, 16): [(26, 13), (16, 15)],
    (29, 11): [(26, 13), (16, 15)],
    (29, 21): [(26, 25), (22, 0), (16, 15)],
    (29, 17): [],
    # autoimmune
    (29, 8): [(25, 18), (22, 8)],
    (29, 26): [(25, 12), (25, 18), (22, 8)],
    (29, 13): [(25, 12), (22, 8)],
    # keratometry
    (29, 10): [],
    # heel bone mineral density
    (29, 15): [],
    # hapiness
    (29, 27): [],
    # asthma
    (29, 12): [],
}

# %% tags=[]
CLUSTER_LV_DIR = conf.RESULTS["CLUSTERING_INTERPRETATION"]["BASE_DIR"] / "cluster_lvs"
assert CLUSTER_LV_DIR.exists()

display(CLUSTER_LV_DIR)


# %% tags=[]
def _get_lvs_data(part_k, cluster_idx):
    cluster_lvs = pd.read_pickle(
        CLUSTER_LV_DIR
        / f"part{part_k}"
        / f"cluster_interpreter-part{part_k}_k{cluster_idx}.pkl"
    )

    #     # keep well-aligned lvs only
    #     cluster_lvs = cluster_lvs[
    #         cluster_lvs.index.astype(str).isin(well_aligned_lvs["LV index"])
    #     ]

    cluster_lvs["cluster"] = cluster_idx

    return list(cluster_lvs["name"])


# %% tags=[]
_get_lvs_data(29, 11)[:5]

# %% [markdown] tags=[]
# # GLSPhenoplier

# %% [markdown] tags=[]
# ## Get list of phenotypes/lvs pairs

# %% tags=[]
phenotypes_lvs_pairs = []

for (part_k, cluster_id), extra_for_lvs in PHENOTYPES_LVS_CONFIG.items():
    part = best_partitions.loc[part_k, "partition"]

    # get traits
    cluster_traits = data.index[part == cluster_id]

    # get extra lvs
    lv_list = _get_lvs_data(part_k, cluster_id)

    for extra_part_k, extra_cluster_id in extra_for_lvs:
        extra_lv_list = _get_lvs_data(part_k, cluster_id)
        lv_list.extend(extra_lv_list)

    for phenotype_code in cluster_traits:
        for lv_code in lv_list:
            phenotypes_lvs_pairs.append(
                {
                    "phenotype_part_k": part_k,
                    "phenotype_cluster_id": cluster_id,
                    "phenotype": phenotype_code,
                    "lv": lv_code,
                }
            )

phenotypes_lvs_pairs = pd.DataFrame(phenotypes_lvs_pairs).drop_duplicates()

# %% tags=[]
phenotypes_lvs_pairs = phenotypes_lvs_pairs.sort_values("phenotype").reset_index(
    drop=True
)

# %% tags=[]
phenotypes_lvs_pairs.shape

# %% tags=[]
phenotypes_lvs_pairs.head()

# %% [markdown] tags=[]
# ## Run

# %% tags=[]
output_file = OUTPUT_DIR / "gls_phenotypes.pkl"
display(output_file)

# %% tags=[]
results = []

pbar = tqdm(total=phenotypes_lvs_pairs.shape[0])

for idx, row in phenotypes_lvs_pairs.iterrows():
    phenotype_code = row["phenotype"]
    lv_code = row["lv"]

    pbar.set_description(f"{phenotype_code} - {lv_code}")

    gls_model = GLSPhenoplier(
        smultixcan_result_set_filepath=conf.PHENOMEXCAN[
            "SMULTIXCAN_EFO_PARTIAL_MASHR_ZSCORES_FILE"
        ]
    )
    gls_model.fit_named(lv_code, phenotype_code)
    res = gls_model.results

    results.append(
        {
            "part_k": row["phenotype_part_k"],
            "cluster_id": row["phenotype_cluster_id"],
            "phenotype": phenotype_code,
            "lv": lv_code,
            "lv_with_pathway": lv_code in well_aligned_lv_codes,
            "coef": res.params.loc["lv"],
            "pvalue": res.pvalues.loc["lv"],
            "summary": gls_model.results_summary,
        }
    )

    if (idx % 10) == 0:
        pd.DataFrame(results).to_pickle(output_file)

    pbar.update(1)

pbar.close()

# %% tags=[]
results = pd.DataFrame(results)

# %% tags=[]
results.shape

# %% tags=[]
results.head()

# %% tags=[]
results.sort_values("pvalue").head(10)

# %% [markdown] tags=[]
# ## Save

# %% tags=[]
results.to_pickle(output_file)

# %% tags=[]