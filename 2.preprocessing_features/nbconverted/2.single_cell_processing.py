#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the single-cell profiles to all samples on the plate (MAD robustize).

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
from pycytominer import annotate, feature_select, normalize
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[2]:


# Set the plate to process
plate_id = "CARD-CelIns-CX7_261022150001"

# Directory with QC-labeled profiles
qc_labeled_dir = pathlib.Path("./data/qc_labeled_profiles/").resolve(strict=True)

# Directory to save single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles/")
output_dir.mkdir(parents=True, exist_ok=True)

# Path to the platemap for the pilot plate
platemap_path = pathlib.Path(
    "../metadata/platemaps/pilot_platemap.csv"
).resolve(strict=True)

# Path to the QC-labeled profile for the pilot plate
profile_path = (qc_labeled_dir / f"{plate_id}_qc_labeled.parquet").resolve(
    strict=True
)

# operations to perform for feature selection
# frequency_threshold removes features dominated by one value (in newer pycytominer versions, e.g., 1.7.1,
# variance_threshold only removes features with low raw variance)
feature_select_ops = [
    "blocklist",
    "drop_na_columns",
    "variance_threshold",
    "frequency_threshold",
    "correlation_threshold",
]


# ## Process data with pycytominer

# In[3]:


print("Performing preprocessing on", plate_id)

output_annotated_file = str(output_dir / f"{plate_id}_sc_annotated.parquet")
output_normalized_file = str(output_dir / f"{plate_id}_sc_normalized.parquet")
output_feature_select_file = str(
    output_dir / f"{plate_id}_sc_feature_selected.parquet"
)

profile_df = pd.read_parquet(profile_path)
platemap_df = pd.read_csv(platemap_path)

# Drop all rows in the profiles that failed any Metadata_cqc columns
cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
if cqc_columns:
    profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

print("Performing annotation for", plate_id, "...")
# Step 1: Annotation
annotated_df = annotate(
    profiles=profile_df,
    platemap=platemap_df,
    join_on=["Metadata_well_position", "Image_Metadata_Well"],
)

# Fix metadata column names (add Metadata_ prefix to site and image file/path name columns
# so they are not treated as features downstream)
annotated_df = annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"})
annotated_df = annotated_df.rename(
    columns={
        col: f"Metadata_{col}"
        for col in annotated_df.columns
        if col.startswith(("Image_FileName_", "Image_PathName_"))
    }
)

# Save annotated profiles
output(
    df=annotated_df,
    output_filename=output_annotated_file,
    output_type="parquet",
)

# Step 2: Normalization
normalized_df = normalize(
    profiles=output_annotated_file,
    method="mad_robustize",
    output_file=output_normalized_file,
    output_type="parquet",
    samples="all",
)

# Step 3: Feature selection
print("Performing feature selection for", plate_id, "...")
feature_select(
    profiles=normalized_df,
    operation=feature_select_ops,
    na_cutoff=0,
    output_file=output_feature_select_file,
    output_type="parquet",
    blocklist_file="./blocklist_features.txt",
)

print(f"Annotation, normalization, and feature selection complete for {plate_id}")

