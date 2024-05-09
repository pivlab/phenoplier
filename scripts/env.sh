#!/bin/bash

# This file is project-specific and it exports some common environmental
# variables to run the code.  It has to be customized for your need by changing
# the BASE_DIR and N_JOBS below.

#
# Your settings here
#

# BASE_DIR is the parent directory where the code and manuscript repos are
# located.
BASE_DIR=/home/miltondp/projects/phenoplier/clean_orig_base/phenoplier

# Project name
#PROJECT_NAME=phenoplier

# Number of CPUs to use
export PHENOPLIER_N_JOBS=1

#
# Do not edit below
#
export PHENOPLIER_ROOT_DIR=${BASE_DIR}
echo "PHENOPLIER_ROOT_DIR=${PHENOPLIER_ROOT_DIR}"

export PHENOPLIER_MANUSCRIPT_DIR=${BASE_DIR}/manuscript
echo "PHENOPLIER_MANUSCRIPT_DIR=${PHENOPLIER_MANUSCRIPT_DIR}"

export PHENOPLIER_CODE_DIR=/home/miltondp/projects/phenoplier/phenoplier-gls-nulls
echo "PHENOPLIER_CODE_DIR=${PHENOPLIER_CODE_DIR}"

export PYTHONPATH=${PHENOPLIER_CODE_DIR}/libs/:${PYTHONPATH}
echo "PYTHONPATH=${PYTHONPATH}"

# GTEx v8
export PHENOPLIER_GTEX_V8_DIR="/home/miltondp/projects/phenoplier/clean_orig_base/gtex_v8/"
echo "PHENOPLIER_GTEX_V8_DIR=${PHENOPLIER_GTEX_V8_DIR}"

