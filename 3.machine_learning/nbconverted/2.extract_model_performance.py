#!/usr/bin/env python
# coding: utf-8

# # Extract model performance metrics
# 
# In this notebook, we extract metrics to evaluate performance such as:
# 
# 1. Precision-recall
# 2. Predicted probabilities

# ## Import libraries

# In[1]:


import pathlib
import sys

import pandas as pd
from joblib import load
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve
from sklearn.preprocessing import LabelEncoder

sys.path.append("../utils")
from training_utils import get_X_y_data


# ## Helper function to collect precision-recall results and predicted probabilities for binary models

# In[2]:


def get_pr_curve_results(
    model: LogisticRegression, df: pd.DataFrame, label: str, label_encoder: LabelEncoder
) -> pd.DataFrame:
    """Collect the precision-recall curve results from a model and dataset.

    Args:
        model (LogisticRegression): loaded in logistic regression model to collect results from
        df (pd.DataFrame): dataframe containing the data to apply model to
        label (str): label with the class being predicted
        label_encoder (LabelEncoder): encoder to transform the labels to integers

    Returns:
        pd.DataFrame: dataframe with the PR curve results for that data and model
    """
    # Get X and y data for the model
    X, y = get_X_y_data(df=df, label=label, shuffle=False)

    assert all(
        col in model.feature_names_in_ for col in X
    ), "Features in the model do not match the columns in the dataset"

    # Transform labels and ensure binary problem
    y_encoded = label_encoder.transform(y)
    assert (
        len(set(y_encoded)) == 2
    ), f"Expected binary labels after encoding, got classes: {sorted(set(y_encoded))}"

    # Use probability for the positive class
    y_scores = model.predict_proba(X)[:, 1]
    assert len(y_scores) == len(
        y_encoded
    ), "Length mismatch between predicted scores and labels"

    precision, recall, _ = precision_recall_curve(y_encoded, y_scores)

    return pd.DataFrame({"precision": precision, "recall": recall})


# In[3]:


def get_predicted_probabilities(
    model: LogisticRegression, df: pd.DataFrame, label: str
) -> pd.DataFrame:
    """Collect predicted probabilities per single-cell from the model and dataset.

    Args:
        model (LogisticRegression): loaded in logistic regression model to collect results from
        df (pd.DataFrame): dataframe containing the data to apply model to
        label (str): label with the class being predicted

    Returns:
        pd.DataFrame: dataframe with the predicted probabilities per single-cell
    """
    # Validate required metadata columns exist for output
    for col in ("Metadata_N_Replicate", "Metadata_Well"):
        assert col in df.columns, f"Required column '{col}' not found in dataframe"

    # Get X and y for the model
    X, y = get_X_y_data(df=df, label=label, shuffle=False)

    # Ensure model features match dataset columns
    assert all(
        col in model.feature_names_in_ for col in X
    ), "Features in the model do not match the columns in the dataset"

    # Use probability for the positive class
    y_scores = model.predict_proba(X)[:, 1]
    assert len(y_scores) == len(
        y
    ), "Length mismatch between predicted probabilities and labels"

    # Holdout type only exists in the holdout split
    holdout_type = (
        df["Metadata_holdout_type"].values
        if "Metadata_holdout_type" in df.columns
        else None
    )

    return pd.DataFrame(
        {
            "actual_label": y,
            "predicted_probability": y_scores,
            "Metadata_N_Replicate": df["Metadata_N_Replicate"].values,
            "Metadata_Well": df["Metadata_Well"].values,
            "Metadata_holdout_type": holdout_type,
        }
    )


# ## Set paths

# In[4]:


# Directory with the training, testing, and holdout datasets
data_dir = pathlib.Path("data_splits")

# Directory with the trained models
model_dir = pathlib.Path("models")

# Directory with encoder
encoder_dir = pathlib.Path("encoder_results")

# Output directory the performance metrics
performance_metrics_dir = pathlib.Path("performance_metrics")
performance_metrics_dir.mkdir(exist_ok=True)

# Label being predicted
label = "Metadata_Treatment"


# ## Create dictionary with all relevant paths per model to extract metrics

# In[5]:


# Extract model names from model filenames
model_names = sorted(
    f.stem.replace("_final_downsample", "")
    for f in model_dir.glob("*_final_downsample.joblib")
)

# Create a nested dictionary with info per model
models_dict = {}
for model in model_names:
    models_dict[model] = {
        "training_data": pathlib.Path(
            data_dir / model / "downsample_training_split.parquet"
        ).resolve(strict=True),
        "testing_data": pathlib.Path(
            data_dir / model / "testing_split.parquet"
        ).resolve(strict=True),
        "holdout_data": pathlib.Path(
            data_dir / model / "holdout_split.parquet"
        ).resolve(strict=True),
        "final_model": pathlib.Path(
            model_dir / f"{model}_final_downsample.joblib"
        ).resolve(strict=True),
        "shuffled_model": pathlib.Path(
            model_dir / f"{model}_shuffled_downsample.joblib"
        ).resolve(strict=True),
        "encoder_result": pathlib.Path(
            encoder_dir / "label_encoder_global.joblib"
        ).resolve(strict=True),
    }

# Print out dictionary keys and paths for verification
for model, paths in models_dict.items():
    lines = [f"Model: {model}"] + [f"  {key}: {path}" for key, path in paths.items()]
    print("\n".join(lines))


# In[6]:


# For each model, print the number of cells per biological replicate and treatment per data split
for model_name, paths in models_dict.items():
    print(f"Model: {model_name}")
    for split_name in ["training_data", "testing_data", "holdout_data"]:
        df = pd.read_parquet(paths[split_name])
        counts = df.groupby(["Metadata_N_Replicate", label]).size().unstack(fill_value=0)
        print(f"  {split_name}:\n{counts.to_string()}\n")


# ## Extract metrics from the data splits applied to their respective models

# In[7]:


# Initialize results list
test_train_pr_results = []
test_train_probability_results = []

# Run through each model and get the PR results
for model_name, paths in models_dict.items():
    # Load the models and data
    final_model = load(paths["final_model"])
    shuffled_model = load(paths["shuffled_model"])
    label_encoder = load(paths["encoder_result"])
    downsample_train_df = pd.read_parquet(paths["training_data"])
    test_df = pd.read_parquet(paths["testing_data"])
    holdout_df = pd.read_parquet(paths["holdout_data"])

    # Confirm the positive class (predicted probability) is the vehicle class
    assert (
        label_encoder.classes_[1] == "Veh"
    ), f"Expected 'Veh' as positive class, got {label_encoder.classes_}"

    print(f"Processing model: {model_name}")

    # Set dictionary with the data splits
    datasets = {"train": downsample_train_df, "test": test_df, "holdout": holdout_df}

    # Loop through all datasets and models
    for dataset_name, dataset in datasets.items():
        for model_type, model in [("final", final_model), ("shuffled", shuffled_model)]:
            # Get per-sample predicted probabilities
            prob_df = get_predicted_probabilities(model=model, df=dataset, label=label)
            prob_df["model_type"] = model_type
            prob_df["dataset"] = dataset_name
            prob_df["model_name"] = model_name
            test_train_probability_results.append(prob_df)

            # Get PR curve results (global)
            pr_df = get_pr_curve_results(
                model=model, df=dataset, label=label, label_encoder=label_encoder
            )
            pr_df["model_type"] = model_type
            pr_df["dataset"] = dataset_name
            pr_df["model_name"] = model_name
            test_train_pr_results.append(pr_df)

            print(f"{model_name} | {model_type} | {dataset_name} → Done")

# Combine all results into one dataframe
all_models_pr_results_df = pd.concat(test_train_pr_results, ignore_index=True)
all_models_probabilities_df = pd.concat(
    test_train_probability_results, ignore_index=True
)

# Save the results
all_models_pr_results_df.to_parquet(
    performance_metrics_dir / "all_models_pr_curve_results.parquet", index=False
)
all_models_probabilities_df.to_parquet(
    performance_metrics_dir / "all_models_predicted_probabilities.parquet", index=False
)

# Check output
print(all_models_pr_results_df.shape)
all_models_pr_results_df.head(2)

