# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill
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
# # Description

# %% [markdown] tags=[]
# This notebook evaluates how a spectral clustering method performs as a consensus function. It takes the coassociation matrix, applies several `delta` values to transform it, and computes different clustering quality measures to assess performance. An optimal `delta` value is chosen, and will be used to perform the full analysis later.
#
# This notebook loads the `z_score_std` data version to compute two of the clustering quality measures (Calinski-Harabasz and Davies-Bouldin). The Silhouette score is computed on the ensemble distance matrix, so it is not affected by the data loaded. There are other two notebooks that perform exactly the same steps here but loading the `pca` and `umap` data versions.

# %% [markdown] tags=[]
# # Environment variables

# %% tags=[] trusted=true
from IPython.display import display

import conf

N_JOBS = conf.GENERAL["N_JOBS"]
display(N_JOBS)

# %% tags=[] trusted=true
# %env MKL_NUM_THREADS=$N_JOBS
# %env OPEN_BLAS_NUM_THREADS=$N_JOBS
# %env NUMEXPR_NUM_THREADS=$N_JOBS
# %env OMP_NUM_THREADS=$N_JOBS

# %% [markdown] tags=[]
# # Modules loading

# %% tags=[] trusted=true
# %load_ext autoreload
# %autoreload 2

# %% tags=[] trusted=true
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from utils import generate_result_set_name

# %% [markdown] tags=[]
# # Settings

# %% tags=[] trusted=true
INITIAL_RANDOM_STATE = 100000

# %% trusted=true
# output dir for this notebook
CONSENSUS_CLUSTERING_DIR = Path(
    conf.RESULTS["CLUSTERING_DIR"], "consensus_clustering"
).resolve()

display(CONSENSUS_CLUSTERING_DIR)

# %% [markdown]
# # Load data

# %% tags=[] trusted=true
INPUT_SUBSET = "z_score_std"

# %% tags=[] trusted=true
INPUT_STEM = "projection-smultixcan-efo_partial-mashr-zscores"

# %% tags=[] trusted=true
input_filepath = Path(
    conf.RESULTS["DATA_TRANSFORMATIONS_DIR"],
    INPUT_SUBSET,
    f"{INPUT_SUBSET}-{INPUT_STEM}.pkl",
).resolve()
display(input_filepath)

assert input_filepath.exists(), "Input file does not exist"

input_filepath_stem = input_filepath.stem
display(input_filepath_stem)

# %% trusted=true
data = pd.read_pickle(input_filepath)

# %% trusted=true
data.shape

# %% trusted=true
data.head()

# %% trusted=true
traits = data.index.tolist()

# %% trusted=true
len(traits)

# %% [markdown] tags=[]
# # Ensemble (coassociation matrix)

# %% trusted=true
input_file = Path(CONSENSUS_CLUSTERING_DIR, "ensemble_coassoc_matrix.npy").resolve()
display(input_file)

# %% trusted=true
coassoc_matrix = np.load(input_file)

# %% trusted=true
coassoc_matrix = pd.DataFrame(
    data=coassoc_matrix,
    index=traits,
    columns=traits,
)

# %% trusted=true
coassoc_matrix.shape

# %% trusted=true
coassoc_matrix.head()

# %% trusted=true
dist_matrix = coassoc_matrix

# %% [markdown] tags=[]
# # Clustering

# %% tags=[] trusted=true
from sklearn.cluster import SpectralClustering
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)

from clustering.utils import compute_performance

# %% [markdown] tags=[]
# ## `delta` parameter (gaussian kernel)

# %% [markdown] tags=[]
# ### `delta=1.0`

# %% trusted=true
delta = 1.0

# %% tags=[] trusted=true
with warnings.catch_warnings():
    warnings.filterwarnings("always")

    sim_matrix = np.exp(-(dist_matrix ** 2) / (2.0 * delta ** 2))

    clus = SpectralClustering(
        eigen_solver="arpack",
        #         eigen_tol=1e-3,
        n_clusters=2,
        n_init=10,
        affinity="precomputed",
        random_state=INITIAL_RANDOM_STATE,
    )

    part = clus.fit_predict(sim_matrix)

# %% tags=[] trusted=true
# show number of clusters and their size
pd.Series(part).value_counts()

# %% tags=[] trusted=true
compute_performance(data, part, data_distance_matrix=dist_matrix)

# %% [markdown] tags=[]
# For `delta=1.0`, the algorithm works/converges fine with this data version.

# %% [markdown] tags=[]
# ### `delta>1.0`

# %% trusted=true
delta = 10.0

# %% tags=[] trusted=true
with warnings.catch_warnings():
    warnings.filterwarnings("always")

    sim_matrix = np.exp(-(dist_matrix ** 2) / (2.0 * delta ** 2))

    clus = SpectralClustering(
        eigen_solver="arpack",
        #         eigen_tol=1e-3,
        n_clusters=2,
        n_init=10,
        affinity="precomputed",
        random_state=INITIAL_RANDOM_STATE,
    )

    part = clus.fit_predict(sim_matrix)

# %% tags=[] trusted=true
# show number of clusters and their size
pd.Series(part).value_counts()

# %% tags=[] trusted=true
compute_performance(data, part, data_distance_matrix=dist_matrix)

# %% [markdown] tags=[]
# For `delta` values larger than `1.0`, quality measures go slightly down.

# %% [markdown] tags=[]
# ### `delta<1.0`

# %% trusted=true
delta = 0.20

# %% tags=[] trusted=true
with warnings.catch_warnings():
    warnings.filterwarnings("always")

    sim_matrix = np.exp(-(dist_matrix ** 2) / (2.0 * delta ** 2))

    clus = SpectralClustering(
        eigen_solver="arpack",
        #         eigen_tol=1e-3,
        n_clusters=2,
        n_init=10,
        affinity="precomputed",
        random_state=INITIAL_RANDOM_STATE,
    )

    part = clus.fit_predict(sim_matrix)

# %% tags=[] trusted=true
# show number of clusters and their size
pd.Series(part).value_counts()

# %% tags=[] trusted=true
compute_performance(data, part, data_distance_matrix=dist_matrix)

# %% [markdown] tags=[]
# For `delta` values smaller than `1.0`, quality measures improve.

# %% [markdown] tags=[]
# ### `delta<<<1.0`

# %% trusted=true
delta = 0.10

# %% tags=[] trusted=true
with warnings.catch_warnings():
    warnings.filterwarnings("always")

    sim_matrix = np.exp(-(dist_matrix ** 2) / (2.0 * delta ** 2))

    clus = SpectralClustering(
        eigen_solver="arpack",
        eigen_tol=1e-4,
        n_clusters=2,
        n_init=10,
        affinity="precomputed",
        random_state=INITIAL_RANDOM_STATE,
    )

    part = clus.fit_predict(sim_matrix)

# %% tags=[] trusted=true
# show number of clusters and their size
pd.Series(part).value_counts()

# %% tags=[] trusted=true
compute_performance(data, part, data_distance_matrix=dist_matrix)

# %% [markdown] tags=[]
# For `delta` values around `0.10` the algorithm does not converge, and I need to force convergence with `eigen_tol=1e-4`.
#
# Quality measures in general go down.

# %% [markdown] tags=[]
# ## Extended test

# %% trusted=true
from clustering.methods import DeltaSpectralClustering

# %% tags=[] trusted=true
CLUSTERING_OPTIONS = {}

CLUSTERING_OPTIONS["K_RANGE"] = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40]
CLUSTERING_OPTIONS["N_REPS_PER_K"] = 5
CLUSTERING_OPTIONS["KMEANS_N_INIT"] = 10
CLUSTERING_OPTIONS["DELTAS"] = [
    5.00,
    2.00,
    1.00,
    0.90,
    0.75,
    0.50,
    0.30,
    0.25,
    0.20,
]

display(CLUSTERING_OPTIONS)

# %% tags=[] trusted=true
CLUSTERERS = {}

idx = 0
random_state = INITIAL_RANDOM_STATE

for k in CLUSTERING_OPTIONS["K_RANGE"]:
    for delta_value in CLUSTERING_OPTIONS["DELTAS"]:
        for i in range(CLUSTERING_OPTIONS["N_REPS_PER_K"]):
            clus = DeltaSpectralClustering(
                eigen_solver="arpack",
                n_clusters=k,
                n_init=CLUSTERING_OPTIONS["KMEANS_N_INIT"],
                affinity="precomputed",
                delta=delta_value,
                random_state=random_state,
            )

            method_name = type(clus).__name__
            CLUSTERERS[f"{method_name} #{idx}"] = clus

            random_state = random_state + 1
            idx = idx + 1

# %% tags=[] trusted=true
display(len(CLUSTERERS))

# %% tags=[] trusted=true
_iter = iter(CLUSTERERS.items())
display(next(_iter))
display(next(_iter))

# %% tags=[] trusted=true
clustering_method_name = method_name
display(clustering_method_name)

# %% [markdown] tags=[]
# ## Generate ensemble

# %% tags=[] trusted=true
import tempfile
from clustering.ensembles.utils import generate_ensemble

# %% tags=[] trusted=true
# generate a temporary folder where to store the ensemble and avoid computing it again
ensemble_folder = Path(
    tempfile.gettempdir(),
    f"pre_cluster_analysis",
    clustering_method_name,
).resolve()
ensemble_folder.mkdir(parents=True, exist_ok=True)

# %% tags=[] trusted=true
ensemble_file = Path(
    ensemble_folder,
    generate_result_set_name(CLUSTERING_OPTIONS, prefix=f"ensemble-", suffix=".pkl"),
)
display(ensemble_file)

# %% tags=[] trusted=true
if ensemble_file.exists():
    display(f"Ensemble file exists")
    ensemble = pd.read_pickle(ensemble_file)
else:
    ensemble = generate_ensemble(
        dist_matrix,
        CLUSTERERS,
        attributes=["n_clusters", "delta"],
    )

# %% tags=[] trusted=true
ensemble.shape

# %% tags=[] trusted=true
ensemble.head()

# %% tags=[] trusted=true
ensemble["delta"] = ensemble["delta"].apply(lambda x: f"{x:.2f}")

# %% tags=[] trusted=true
ensemble["n_clusters"].value_counts()

# %% tags=[] trusted=true
_tmp = ensemble["n_clusters"].value_counts().unique()
assert _tmp.shape[0] == 1
assert _tmp[0] == int(
    CLUSTERING_OPTIONS["N_REPS_PER_K"] * len(CLUSTERING_OPTIONS["DELTAS"])
)

# %% tags=[] trusted=true
ensemble_stats = ensemble["n_clusters"].describe()
display(ensemble_stats)

# %% [markdown] tags=[]
# ### Testing

# %% tags=[] trusted=true
assert ensemble_stats["min"] > 1

# %% tags=[] trusted=true
assert not ensemble["n_clusters"].isna().any()

# %% tags=[] trusted=true
assert ensemble.shape[0] == len(CLUSTERERS)

# %% tags=[] trusted=true
# all partitions have the right size
assert np.all(
    [part["partition"].shape[0] == data.shape[0] for idx, part in ensemble.iterrows()]
)

# %% tags=[] trusted=true
# no partition has negative clusters (noisy points)
assert not np.any([(part["partition"] < 0).any() for idx, part in ensemble.iterrows()])

# %% tags=[] trusted=true
# check that the number of clusters in the partitions are the expected ones
_real_k_values = ensemble["partition"].apply(lambda x: np.unique(x).shape[0])
display(_real_k_values)
assert np.all(ensemble["n_clusters"].values == _real_k_values.values)

# %% [markdown] tags=[]
# ### Add clustering quality measures

# %% tags=[] trusted=true
ensemble = ensemble.assign(
    si_score=ensemble["partition"].apply(
        lambda x: silhouette_score(dist_matrix, x, metric="precomputed")
    ),
    ch_score=ensemble["partition"].apply(lambda x: calinski_harabasz_score(data, x)),
    db_score=ensemble["partition"].apply(lambda x: davies_bouldin_score(data, x)),
)

# %% tags=[] trusted=true
ensemble.shape

# %% tags=[] trusted=true
ensemble.head()

# %% [markdown]
# ### Save

# %% trusted=true
ensemble.to_pickle(ensemble_file)

# %% [markdown] tags=[]
# # Cluster quality

# %% tags=[] trusted=true
with pd.option_context("display.max_rows", None, "display.max_columns", None):
    _df = ensemble.groupby(["n_clusters", "delta"]).mean()
    display(_df)

# %% tags=[] trusted=true
with sns.plotting_context("talk", font_scale=0.75), sns.axes_style(
    "whitegrid", {"grid.linestyle": "--"}
):
    fig = plt.figure(figsize=(14, 6))
    ax = sns.pointplot(data=ensemble, x="n_clusters", y="si_score", hue="delta")
    ax.set_ylabel("Silhouette index\n(higher is better)")
    ax.set_xlabel("Number of clusters ($k$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.grid(True)
    plt.tight_layout()

# %% tags=[] trusted=true
with sns.plotting_context("talk", font_scale=0.75), sns.axes_style(
    "whitegrid", {"grid.linestyle": "--"}
):
    fig = plt.figure(figsize=(14, 6))
    ax = sns.pointplot(data=ensemble, x="n_clusters", y="ch_score", hue="delta")
    ax.set_ylabel("Calinski-Harabasz index\n(higher is better)")
    ax.set_xlabel("Number of clusters ($k$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.grid(True)
    plt.tight_layout()

# %% tags=[] trusted=true
with sns.plotting_context("talk", font_scale=0.75), sns.axes_style(
    "whitegrid", {"grid.linestyle": "--"}
):
    fig = plt.figure(figsize=(14, 6))
    ax = sns.pointplot(data=ensemble, x="n_clusters", y="db_score", hue="delta")
    ax.set_ylabel("Davies-Bouldin index\n(lower is better)")
    ax.set_xlabel("Number of clusters ($k$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.grid(True)
    plt.tight_layout()

# %% [markdown] tags=[]
# # Stability

# %% [markdown] tags=[]
# ## Group ensemble by n_clusters

# %% tags=[] trusted=true
parts = ensemble.groupby(["delta", "n_clusters"]).apply(
    lambda x: np.concatenate(x["partition"].apply(lambda x: x.reshape(1, -1)), axis=0)
)

# %% tags=[] trusted=true
parts.shape

# %% tags=[] trusted=true
parts.head()

# %% tags=[] trusted=true
parts.iloc[0].shape

# %% tags=[] trusted=true
assert np.all(
    [
        parts.loc[k].shape == (int(CLUSTERING_OPTIONS["N_REPS_PER_K"]), data.shape[0])
        for k in parts.index
    ]
)

# %% [markdown] tags=[]
# ## Compute stability

# %% tags=[] trusted=true
from sklearn.metrics import adjusted_rand_score as ari
from scipy.spatial.distance import pdist

# %% tags=[] trusted=true
parts_ari = pd.Series(
    {k: pdist(parts.loc[k], metric=ari) for k in parts.index}, name="n_clusters"
)

# %% tags=[] trusted=true
parts_ari_stability = parts_ari.apply(lambda x: x.mean())
display(parts_ari_stability.sort_values(ascending=False).head(15))

# %% tags=[] trusted=true
parts_ari_df = pd.DataFrame.from_records(parts_ari.tolist()).set_index(
    parts_ari.index.copy()
)
parts_ari_df.index.rename(["gamma", "n_clusters"], inplace=True)

# %% tags=[] trusted=true
parts_ari_df.shape

# %% tags=[] trusted=true
_n_total_parts = int(
    CLUSTERING_OPTIONS["N_REPS_PER_K"]
)  # * len(CLUSTERING_OPTIONS["GAMMAS"]))

assert int(_n_total_parts * (_n_total_parts - 1) / 2) == parts_ari_df.shape[1]

# %% tags=[] trusted=true
parts_ari_df.head()

# %% [markdown] tags=[]
# ## Stability plot

# %% tags=[] trusted=true
parts_ari_df_plot = (
    parts_ari_df.stack().reset_index().rename(columns={"level_2": "idx", 0: "ari"})
)

# %% tags=[] trusted=true
parts_ari_df_plot.dtypes

# %% tags=[] trusted=true
parts_ari_df_plot.head()

# %% tags=[] trusted=true
with pd.option_context("display.max_rows", None, "display.max_columns", None):
    _df = parts_ari_df_plot.groupby(["n_clusters", "gamma"]).mean()
    display(_df)

# %% tags=[] trusted=true
with sns.plotting_context("talk", font_scale=0.75), sns.axes_style(
    "whitegrid", {"grid.linestyle": "--"}
):
    fig = plt.figure(figsize=(14, 6))
    ax = sns.pointplot(data=parts_ari_df_plot, x="n_clusters", y="ari", hue="gamma")
    ax.set_ylabel("Averange ARI")
    ax.set_xlabel("Number of clusters ($k$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.grid(True)
    plt.tight_layout()

# %% [markdown] tags=[]
# **CONCLUSION:** the best values for the `delta` parameter seem to be `0.20`, `0.25` and `0.30`. I will also consider `0.50`, since seem to yield potentially good results when `n_clusters` is high (see `umap` results).

# %% tags=[] trusted=true