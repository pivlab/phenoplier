# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] tags=[]
# # Description

# %% [markdown] tags=[]
# Explore GWAS file structure

# %% [markdown] tags=[]
# # Modules

# %% tags=[]
import re
import subprocess
from pathlib import Path
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

import conf
from utils import chunker

# %% [markdown] tags=[]
# # Settings

# %% tags=["parameters"]
PROJECTS_TRAIT_KEY = "CHRONOTYPE"

# %% tags=["injected-parameters"]
# Parameters
PHENOPLIER_NOTEBOOK_FILEPATH = "projects/chronotype/nbs/05_gwas/03-gwas-explore.ipynb"


# %% [markdown] tags=[]
# # Paths

# %% tags=[]
INPUT_GWAS_DIR = conf.PROJECTS[PROJECTS_TRAIT_KEY]["DATA_DIR"] / "gwas"
display(INPUT_GWAS_DIR)
assert INPUT_GWAS_DIR.exists()

# %% tags=[]
PLINK2 = conf.PLINK["EXECUTABLE_VERSION_2"]
display(PLINK2)
assert PLINK2.exists()

# %% [markdown] tags=[]
# # GWAS results files

# %% tags=[]
# check files in directory
gwas_files = sorted(list(INPUT_GWAS_DIR.glob("*.gz")))
display(len(gwas_files))
display(gwas_files[:10])

# %% tags=[]
# get files from traits info file
traits_info = pd.read_csv(conf.PROJECTS[PROJECTS_TRAIT_KEY]["TRAITS_INFO_FILE"])

# %% tags=[]
traits_info.shape

# %% tags=[]
traits_info

# %% tags=[]
gwas_files = [INPUT_GWAS_DIR / t.gwas_file for _, t in traits_info.iterrows()]

# %% tags=[]
len(gwas_files)

# %% tags=[]
gwas_files

# %% [markdown] tags=[]
# # Load GWAS

# %% tags=[]
df = pd.read_csv(gwas_files[0], sep="\t", nrows=10)

# %% tags=[]
df

# %% tags=[]
df = pd.read_csv(gwas_files[0], sep="\t")

# %% tags=[]
df.shape

# %% tags=[]
df.head()

# %% tags=[]
