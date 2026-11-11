#!/bin/bash

# stop the workflow if any step fails to avoid running CellProfiler on stale inputs
set -eo pipefail

# enable conda activate within this script
eval "$(conda shell.bash hook)"
# activate the CellProfiler environment
conda activate mouse_fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run Python script for segmentation and feature extraction
python nbconverted/0.create_loaddata_with_illum_csv.py
python nbconverted/1.extract_features.py
