#!/usr/bin/env python
# coding: utf-8

# # Split single-cell data for TGF-β versus vehicle model
# 
# We split the feature selected single-cell profiles into data splits for one binary model that predicts TGF-β treated (fibrotic) versus vehicle (non-failing) treated cardiac fibroblasts.
# 
# 1. **Holdout:** Hold out one whole biological replicate (randomly selected) and one random well per treatment from each of the remaining biological replicates.
# 2. **Training and testing:** Split the remaining data using 70/30 splitting, stratified by treatment.

# In[1]:


import pathlib

import pandas as pd
from sklearn.model_selection import train_test_split


# ## Set paths and variables

# In[2]:


# Set random state for the whole notebook to ensure reproducibility
random_state = 0

# Plate to process
plate = "CARD-CelIns-CX7_261022150001"

# Metadata column used for prediction class and column with biological replicate
label = "Metadata_Treatment"
replicate_col = "Metadata_N_Replicate"

# Path to feature selected single-cell profiles (QC applied)
feature_selected_path = pathlib.Path(
    f"../2.preprocessing_features/data/single_cell_profiles/{plate}_sc_feature_selected.parquet"
).resolve(strict=True)

# Make directory for split data for this model
model_output_dir = pathlib.Path("./data_splits/model_TGFB_vs_Veh")
model_output_dir.mkdir(parents=True, exist_ok=True)

print(f"Found feature selected file: {feature_selected_path.name}")


# ## Load in feature selected data

# In[3]:


# Load the feature selected file as a DataFrame
feature_selected_df = pd.read_parquet(feature_selected_path)

# Drop Location features (position of the cell in the image), which are not biologically meaningful
location_features = [
    col
    for col in feature_selected_df.columns
    if "_Location_" in col and not col.startswith("Metadata_")
]
feature_selected_df = feature_selected_df.drop(columns=location_features)
print(f"Dropped {len(location_features)} Location features: {location_features}")

print(f"Shape: {feature_selected_df.shape}")
print(
    feature_selected_df.groupby([replicate_col, label])
    .agg(n_cells=("Metadata_Well", "size"), n_wells=("Metadata_Well", "nunique"))
    .to_string()
)
feature_selected_df.head()


# ## Perform splits for TGF-β versus vehicle model
# 
# Holdout one whole biological replicate (randomly selected) and one random well per treatment from each of the remaining biological replicates.

# In[4]:


# Randomly select one biological replicate to hold out (make sampling reproducible)
holdout_replicate = (
    pd.Series(sorted(feature_selected_df[replicate_col].dropna().unique()))
    .sample(1, random_state=random_state)
    .item()
)
print(f"Randomly selected biological replicate for holdout: {holdout_replicate}")

# Hold out all rows from the selected biological replicate
holdout_mask = feature_selected_df[replicate_col] == holdout_replicate

# Randomly hold out one well per treatment per other biological replicate (make sampling reproducible)
random_wells = (
    feature_selected_df[~holdout_mask]
    .groupby([replicate_col, label])["Metadata_Well"]
    .apply(lambda x: x.dropna().sample(1, random_state=random_state))
    .explode()
)
print(f"Randomly selected wells for holdout: {random_wells.tolist()}")

# Combine holdout replicate rows with random well rows
well_mask = feature_selected_df["Metadata_Well"].isin(random_wells)
holdout_idx = feature_selected_df.index[holdout_mask | well_mask]

# Create holdout and remaining dataframes
holdout_df = feature_selected_df.loc[holdout_idx].copy()
model_df = feature_selected_df.drop(holdout_idx).copy()

# Record which holdout type each cell belongs to (whole replicate or random well)
holdout_df["Metadata_holdout_type"] = "well"
holdout_df.loc[
    holdout_df[replicate_col] == holdout_replicate, "Metadata_holdout_type"
] = "replicate"

# Save holdout set
holdout_df.to_parquet(model_output_dir / "holdout_split.parquet")

print(f"Holdout set shape: {holdout_df.shape}")
print(
    holdout_df.groupby(["Metadata_holdout_type", replicate_col, label])["Metadata_Well"]
    .agg(n_cells="size", wells=lambda s: sorted(s.unique()))
    .to_string()
)


# In[5]:


print(f"Model data shape (after dropping holdout rows): {model_df.shape}")

# Sanity checks
assert (
    holdout_df.shape[0] + model_df.shape[0] == feature_selected_df.shape[0]
), "Holdout + remaining does not equal original after splitting"
assert (
    holdout_replicate not in model_df[replicate_col].unique()
), "Holdout replicate found in the remaining data"

# Set the ratio of the test data to 30% (training data will be 70%)
test_ratio = 0.30

# Split data into training and test sets
train_df, test_df = train_test_split(
    model_df,
    test_size=test_ratio,
    stratify=model_df[label],
    random_state=random_state,
)

# Save training and test data
train_df.to_parquet(model_output_dir / "training_split.parquet")
test_df.to_parquet(model_output_dir / "testing_split.parquet")

print(f"Training data shape: {train_df.shape}")
print(f"Testing data shape: {test_df.shape}")

# Print out the proportion of samples per treatment in each split
print("\nTreatment distribution in training set:")
print(train_df[label].value_counts(normalize=True))
print("\nTreatment distribution in testing set:")
print(test_df[label].value_counts(normalize=True))

