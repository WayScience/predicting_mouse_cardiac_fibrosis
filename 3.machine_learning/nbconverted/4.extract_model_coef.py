#!/usr/bin/env python
# coding: utf-8

# # Extract logistic regression model coefficients per CellProfiler feature
# 
# The coefficients from the final model will either be positive (x > 0), where if the feature value increases, the more likely the cell is vehicle treated (non-failing), negative (x < 0), where if the feature value increases, the more likely the cell is TGF-β treated (fibrotic), or zero (x = 0), where that feature has no impact on the prediction.

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
from joblib import load


# ## Set paths and variables

# In[2]:


# Name of the model to extract coefficients from
model_name = "model_TGFB_vs_Veh"

# Path to the final model
path_to_final_model = pathlib.Path(
    f"./models/{model_name}_final_downsample.joblib"
).resolve(strict=True)

# Path to the label encoder to confirm which class is positive
path_to_encoder = pathlib.Path(
    "./encoder_results/label_encoder_global.joblib"
).resolve(strict=True)

# Output dir for coefficients data
coeff_dir = pathlib.Path("./coeff_data")
coeff_dir.mkdir(exist_ok=True)

# Channels in this dataset (used to find the channel in each feature name)
channels = ["DNA", "ER", "Golgi", "Mito", "Actin"]


# ## Load in the final model and label encoder

# In[3]:


# Load in final model and label encoder
final_model = load(path_to_final_model)
label_encoder = load(path_to_encoder)

# Confirm the positive class (positive coefficients) is the vehicle class
assert list(label_encoder.classes_) == [
    "TGFB",
    "Veh",
], f"Unexpected classes: {label_encoder.classes_}"
print(
    f"Negative coefficients predict {label_encoder.classes_[0]}, positive coefficients predict {label_encoder.classes_[1]}"
)


# ## Collect coefficients from the model and match with the correct feature in a dataframe
# 
# The feature names are stored in the model, so we use them directly to match each coefficient to its feature.
# We also split out the compartment, feature group, and channel from each feature name for visualization.
# The channel is the first channel name found in the feature name (for correlation features, this is the first of the two channels), and features without a channel (e.g., AreaShape and Neighbors) are labelled as `other`.

# In[4]:


# Get the coefficients and feature names from the model
coefficients = final_model.coef_.flatten()
feature_columns = final_model.feature_names_in_

print(f"Number of coefficients: {len(coefficients)}")
assert len(coefficients) == len(
    feature_columns
), "The number of coefficients does not match the number of features"


def find_channel(feature: str) -> str:
    """Find the first channel in a feature name.

    Args:
        feature (str): CellProfiler feature name

    Returns:
        str: the first channel found in the feature name, or "other" if there is no channel
    """
    for part in feature.split("_")[2:]:
        if part in channels:
            return part
    return "other"


# Create a DataFrame with the coefficients and features
coefficients_df = pd.DataFrame({"Feature": feature_columns, "Coefficient": coefficients})

# Split out components of the feature name for visualization
coefficients_df["Compartment"] = coefficients_df["Feature"].str.split("_").str[0]
coefficients_df["Feature_Group"] = coefficients_df["Feature"].str.split("_").str[1]
coefficients_df["Channel"] = coefficients_df["Feature"].apply(find_channel)

# Save the all coefficient data into a CSV file prior to splitting
coefficients_df.to_csv(coeff_dir / "all_coeffs.csv", index=False)

print(coefficients_df["Channel"].value_counts())
coefficients_df.head()


# ## Split the data frame by negative, positive, and zero coefficients which relate to different class importance

# In[5]:


# Split into DataFrames with positive, negative, and zero coefficients
veh_coeffs_df = coefficients_df[coefficients_df["Coefficient"] > 0].copy()
tgfb_coeffs_df = coefficients_df[coefficients_df["Coefficient"] < 0].copy()
zero_coeffs_df = coefficients_df[coefficients_df["Coefficient"] == 0].copy()

# Make the values in tgfb_coeffs_df positive to prevent issues during plotting
tgfb_coeffs_df["Coefficient"] = tgfb_coeffs_df["Coefficient"].abs()

# Rename the coefficient columns
veh_coeffs_df = veh_coeffs_df.rename(columns={"Coefficient": "Veh_Coeffs"})
tgfb_coeffs_df = tgfb_coeffs_df.rename(columns={"Coefficient": "TGFB_Coeffs"})
zero_coeffs_df = zero_coeffs_df.rename(columns={"Coefficient": "Zero_Coeffs"})

# Save the coefficient data into the "coeff_data" folder
veh_coeffs_df.to_csv(coeff_dir / "veh_coeffs.csv", index=False)
tgfb_coeffs_df.to_csv(coeff_dir / "tgfb_coeffs.csv", index=False)
zero_coeffs_df.to_csv(coeff_dir / "zero_coeffs.csv", index=False)

print("Positive (Veh) coefficients:", veh_coeffs_df.shape[0])
print("Negative (TGFB) coefficients:", tgfb_coeffs_df.shape[0])
print("Zero coefficients:", zero_coeffs_df.shape[0])


# ## Explore the coefficients

# In[6]:


# Show the top ten ranking features for predicting the vehicle (Veh) class
coefficients_df.sort_values(by="Coefficient", ascending=False).head(10)


# In[7]:


# Show the top ten ranking features for predicting the TGF-β (TGFB) class
coefficients_df.sort_values(by="Coefficient", ascending=True).head(10)


# ## Add ranking column with sorted descending values and save the CSV for visualization
# 
# Rank is based on the highest positive coefficient, which will have rank one, and then descending from there.

# In[8]:


# Sort coefficients_df by descending order
ranked_coeffs_df = coefficients_df.sort_values(by="Coefficient", ascending=False).copy()

# Add a new column 'Rank'
ranked_coeffs_df["Rank"] = range(1, len(ranked_coeffs_df) + 1)

# Save the ranked df
ranked_coeffs_df.to_csv(coeff_dir / "ranked_coeffs.csv", index=False)

print(ranked_coeffs_df.shape)
ranked_coeffs_df.head()

