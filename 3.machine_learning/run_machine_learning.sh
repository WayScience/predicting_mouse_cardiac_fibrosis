#!/bin/bash

# stop the pipeline if any stage fails to avoid using stale outputs downstream
set -eo pipefail

# run commands in the machine learning environment (uv) without activating it
ml_env="uv run --project ../environments/machine_learning_env"

# convert Jupyter notebook(s) to script
$ml_env jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run Python scripts for machine learning
$ml_env python nbconverted/0.split_data.py
$ml_env python nbconverted/1.train_models.py
$ml_env python nbconverted/2.extract_model_performance.py
$ml_env python nbconverted/4.extract_model_coef.py
$ml_env python nbconverted/6.generate_sc_crop_montage.py

# activate the R environment (conda)
eval "$(conda shell.bash hook)"
conda activate mouse_r_fibrosis_env

# run R scripts for visualizing model performance and coefficients
Rscript nbconverted/3.vis_model_performance.r
Rscript nbconverted/5.vis_model_coef.r
Rscript nbconverted/7.vis_data_split_platemap.r

echo "Machine learning complete."
