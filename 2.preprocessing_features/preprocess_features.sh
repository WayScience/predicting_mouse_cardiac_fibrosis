#!/bin/bash

# stop the pipeline if any stage fails to avoid using stale outputs downstream
set -eo pipefail

# activate the preprocessing environment (uv)
source ../environments/preprocessing_env/.venv/bin/activate

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run the notebooks in order
python nbconverted/0.convert_cytotable.py
python nbconverted/1.sc_quality_control.py
python nbconverted/2.single_cell_processing.py
python nbconverted/3.bulk_processing.py

echo "Preprocessing complete."
