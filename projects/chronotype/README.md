# PhenoPLIER analyses of Chronotype

## Contents

 * [Overview](#overview)
 * [Quick demo](#quick-demo)
 * [Code and data](#code-and-data)
 * [Setup](#setup)
 * [Running the code](#running-the-code)


## Overview

**TODO**


## Quick demo

**TODO**


## Code and data

**TODO**


## Setup

```bash
docker pull miltondp/phenoplier:chronotype
```


## Running the code

**TODO**


## Instructions for developers

**You very likely do not need to follow these steps**, unless you are a developer working on PhenoPLIER.

### Setup Docker image

**This only needs to be done once.**

Pull the right Docker image for this project and tag it accordingly:

```bash
bash projects/chronotype/scripts/create_docker_image.sh
```

### Load project-specific configuration

```bash
. projects/chronotype/scripts/env.sh
```

### Download data/software

If you have other PhenoPLIER projects, link the `data` and `software` subfolder under
`base` to save space. Otherwise, you need to download the data and setup software with
this command:

```bash
bash scripts/run_docker_dev.sh \
  python environment/scripts/setup_data.py --mode projects
```

### Start JupyterLab server

```bash
bash scripts/run_docker_dev.sh
```


### Run notebook from command-line

```bash
bash scripts/run_docker_dev.sh \
  bash nbs/run_nbs.sh projects/asthma-copd/nbs/05_gwas/05-gwas-inflation_factor.ipynb
```
