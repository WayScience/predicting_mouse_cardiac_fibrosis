#!/usr/bin/env python
# coding: utf-8

# # Generate single-cell crop montage of top predicted cells per class and biological replicate
# 
# We generate a montage of single-cell crops to visualize how TGF-β treated (fibrotic) and vehicle treated (non-failing) cells differ across biological replicates.
# 
# 1. **Cells:** We only use cells the model was not trained on (testing and holdout data). Biological replicate N3 is the whole held out replicate.
# 2. **Selection:** Per biological replicate and treatment, we select the correctly predicted cells with the highest predicted probability for their ground truth class. To get a good sense of the cells across the replicate, we take at most one cell per image (site) and at most two cells per well.
# 3. **Crops:** All crops are the same size (1 pixel = 1 µm) and are centered on the cell bounding box. If the crop would go outside the image, it is shifted to stay inside the image (the cell still fits completely in the crop). We only use cells that fit completely in the crop and do not touch the image border.
# 4. **Debris:** We skip all cells from images with large out-of-focus debris in the DNA channel, found as a connected DNA-positive area in a crop much larger than a cluster of nuclei. We also skip images with debris found during visual review.
# 5. **Channels:** We show F-actin, Mito, and DNA. F-actin and DNA have the most important features in the model. Of ER and Mito, Mito has the stronger top features.
# 6. **Display:** Images are illumination corrected (same as feature extraction) and the same intensity range and gamma per channel are used across all crops, so intensity differences between classes are kept. The lower limit is the background and the gamma brightens dim signal (e.g., F-actin in vehicle cells) without saturating bright signal.

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from joblib import load
from matplotlib.patches import Rectangle
from scipy import ndimage as ndi

sys.path.append("../utils")
from training_utils import get_X_y_data


# ## Set paths and variables

# In[2]:


# Plate and model names
plate = "CARD-CelIns-CX7_261022150001"
model_name = "model_TGFB_vs_Veh"

# Label being predicted
label = "Metadata_Treatment"

# Paths to model, data splits, QC labeled profiles (bounding boxes and image paths), and illumination functions
model_path = pathlib.Path(f"./models/{model_name}_final_downsample.joblib").resolve(strict=True)
data_split_dir = pathlib.Path(f"./data_splits/{model_name}").resolve(strict=True)
qc_labeled_path = pathlib.Path(
    f"../2.preprocessing_features/data/qc_labeled_profiles/{plate}_qc_labeled.parquet"
).resolve(strict=True)
illum_dir = pathlib.Path(f"../0.illumination_correction/illum_directory/{plate}").resolve(strict=True)

# Output directories
crop_output_dir = pathlib.Path("./sc_crops")
crop_output_dir.mkdir(exist_ok=True)
montage_output_dir = pathlib.Path("./figures/sc_crop_montage")
montage_output_dir.mkdir(parents=True, exist_ok=True)

# Crop and selection parameters
crop_size = 400  # pixels (1 pixel = 1 µm)
pixel_size_um = 1.0
scale_bar_um = 100
grid_rows, grid_cols = 3, 3  # crops per biological replicate and treatment
n_crops = grid_rows * grid_cols
max_cells_per_well = 2

# Debris filter: skip crops where the largest connected area of DNA signal above the intensity
# threshold is larger than the max area (clusters of nuclei are < 6,000 pixels, debris is > 25,000 pixels)
debris_dna_intensity_threshold = 2000
debris_max_area = 10000

# Channels to show with display colors (RGB) for the composite
channel_colors = {
    "Actin": (0.0, 1.0, 0.0),  # green
    "Mito": (1.0, 0.0, 1.0),  # magenta
    "DNA": (0.0, 1.0, 1.0),  # cyan
}

# Display range per channel (same for all crops):
# lower limit = background (median across crops of each crop's lower percentile),
# upper limit = upper percentile of all pixels across all crops
background_percentile = 5
upper_percentile = 99.8

# Gamma applied after scaling to brighten dim signal (same for all crops)
display_gamma = 0.6

# Images to exclude after visual review as (well, site) with reason
exclude_images = {
    ("G04", 11): "faint DNA debris streak that the debris filter misses",
}

# Order and labels for the montage
replicates = ["N1", "N2", "N3", "N4"]
holdout_replicate = "N3"
treatments = ["TGFB", "Veh"]
treatment_labels = {"TGFB": "TGF-β", "Veh": "Vehicle"}
treatment_colors = {"TGFB": "#eb6834", "Veh": "#2a78d6"}
holdout_color = "#1baf7a"  # same as holdout replicate in model performance figures


# ## Get predicted probabilities for cells the model was not trained on

# In[3]:


# Load in final model
final_model = load(model_path)

# Load testing and holdout data and get predicted probability of vehicle (positive class)
split_dfs = []
for split in ["testing", "holdout"]:
    df = pd.read_parquet(data_split_dir / f"{split}_split.parquet")
    X, _ = get_X_y_data(df=df, label=label, shuffle=False)
    df = df.loc[:, df.columns.str.startswith("Metadata_")].copy()
    df["predicted_probability"] = final_model.predict_proba(X)[:, 1]
    df["dataset"] = split
    split_dfs.append(df)

pred_df = pd.concat(split_dfs, ignore_index=True)

# Probability of the ground truth class and whether the prediction is correct
pred_df["true_class_probability"] = np.where(
    pred_df[label] == "Veh",
    pred_df["predicted_probability"],
    1 - pred_df["predicted_probability"],
)
pred_df["correct"] = pred_df["true_class_probability"] >= 0.5

print(pred_df.shape)
print(pred_df.groupby(["Metadata_N_Replicate", label])["correct"].agg(["size", "mean"]))


# ## Add bounding boxes and image paths for each cell

# In[4]:


channel_names = list(channel_colors.keys())
bbox_cols = [
    f"Cells_AreaShape_BoundingBox{m}_{a}" for m in ["Minimum", "Maximum"] for a in ["X", "Y"]
]
image_cols = [f"Image_{t}_Orig{c}" for c in channel_names for t in ["FileName", "PathName"]]
key_cols = ["Metadata_ImageNumber", "Metadata_Cells_Number_Object_Number"]

# Load only the needed columns from the QC labeled profiles
qc_df = pd.read_parquet(qc_labeled_path, columns=key_cols + bbox_cols + image_cols)

# Merge to the predictions (each cell should match exactly one row)
pred_df = pred_df.merge(qc_df, on=key_cols, how="left", validate="one_to_one")
assert pred_df[bbox_cols + image_cols].notna().all().all(), "Some cells are missing bounding boxes or image paths"

# Get the image size from one image
first = pred_df.iloc[0]
image_height, image_width = tifffile.imread(
    pathlib.Path(first["Image_PathName_OrigActin"]) / first["Image_FileName_OrigActin"]
).shape
print(f"Image size: {image_height} x {image_width}")

# Bounding box size and center
pred_df["bbox_width"] = pred_df["Cells_AreaShape_BoundingBoxMaximum_X"] - pred_df["Cells_AreaShape_BoundingBoxMinimum_X"]
pred_df["bbox_height"] = pred_df["Cells_AreaShape_BoundingBoxMaximum_Y"] - pred_df["Cells_AreaShape_BoundingBoxMinimum_Y"]
pred_df["bbox_center_x"] = (pred_df["Cells_AreaShape_BoundingBoxMaximum_X"] + pred_df["Cells_AreaShape_BoundingBoxMinimum_X"]) / 2
pred_df["bbox_center_y"] = (pred_df["Cells_AreaShape_BoundingBoxMaximum_Y"] + pred_df["Cells_AreaShape_BoundingBoxMinimum_Y"]) / 2

# Cells must fit completely in the crop and not touch the image border
pred_df["fits_crop"] = (pred_df["bbox_width"] < crop_size) & (pred_df["bbox_height"] < crop_size)
pred_df["touches_border"] = (
    (pred_df["Cells_AreaShape_BoundingBoxMinimum_X"] <= 0)
    | (pred_df["Cells_AreaShape_BoundingBoxMinimum_Y"] <= 0)
    | (pred_df["Cells_AreaShape_BoundingBoxMaximum_X"] >= image_width - 1)
    | (pred_df["Cells_AreaShape_BoundingBoxMaximum_Y"] >= image_height - 1)
)

print(pred_df.groupby(label)[["fits_crop", "touches_border"]].mean())


# ## Load illumination functions and set functions to load crops and find debris

# In[5]:


# Load illumination functions per channel
illum_functions = {
    c: np.load(illum_dir / f"{plate}_Illum{c}.npy") for c in channel_names
}

# Cache of loaded illumination corrected images
image_cache = {}


def load_corrected_image(row: pd.Series, channel: str) -> np.ndarray:
    """Load an illumination corrected image (divide, same as CellProfiler pipeline)."""
    path = pathlib.Path(row[f"Image_PathName_Orig{channel}"]) / row[f"Image_FileName_Orig{channel}"]
    if path not in image_cache:
        image_cache[path] = tifffile.imread(path).astype(np.float32) / illum_functions[channel]
    return image_cache[path]


def load_crop(row: pd.Series, channel: str) -> np.ndarray:
    """Load an illumination corrected crop centered on the cell bounding box.

    If the crop would go outside the image, it is shifted to stay inside the image.

    Args:
        row (pd.Series): row with the image path and bounding box center for the cell
        channel (str): name of the channel to load

    Returns:
        np.ndarray: crop of size crop_size x crop_size
    """
    image = load_corrected_image(row, channel)
    x0 = int(np.clip(round(row["bbox_center_x"]) - crop_size // 2, 0, image.shape[1] - crop_size))
    y0 = int(np.clip(round(row["bbox_center_y"]) - crop_size // 2, 0, image.shape[0] - crop_size))
    return image[y0 : y0 + crop_size, x0 : x0 + crop_size]


def has_debris(row: pd.Series) -> bool:
    """Check if the crop has a connected area of DNA signal larger than the debris max area."""
    labels, n_objects = ndi.label(load_crop(row, "DNA") > debris_dna_intensity_threshold)
    if n_objects == 0:
        return False
    return np.bincount(labels.ravel())[1:].max() > debris_max_area


# ## Select the top predicted cells per biological replicate and treatment

# In[6]:


excluded_mask = pd.Series(
    list(zip(pred_df["Metadata_Well"], pred_df["Metadata_Site"])), index=pred_df.index
).isin(exclude_images.keys())
print(f"Excluding {excluded_mask.sum()} cells from {len(exclude_images)} images after visual review")

eligible_df = pred_df[pred_df["correct"] & pred_df["fits_crop"] & ~pred_df["touches_border"] & ~excluded_mask]

selected = []
skipped_debris = []
debris_images = set()  # once a crop from an image has debris, skip all cells from that image
for replicate in replicates:
    for treatment in treatments:
        group_df = eligible_df[
            (eligible_df["Metadata_N_Replicate"] == replicate) & (eligible_df[label] == treatment)
        ].sort_values("true_class_probability", ascending=False)

        used_sites, well_counts, picks = set(), {}, []
        for idx, row in group_df.iterrows():
            site_key = (row["Metadata_Well"], row["Metadata_Site"])
            if (
                site_key in used_sites
                or site_key in debris_images
                or well_counts.get(row["Metadata_Well"], 0) >= max_cells_per_well
            ):
                continue
            if has_debris(row):
                debris_images.add(site_key)
                skipped_debris.append((replicate, treatment, row["Metadata_Well"], row["Metadata_Site"]))
                continue
            used_sites.add(site_key)
            well_counts[row["Metadata_Well"]] = well_counts.get(row["Metadata_Well"], 0) + 1
            picks.append(idx)
            if len(picks) == n_crops:
                break

        assert len(picks) == n_crops, f"Only found {len(picks)} cells for {replicate} {treatment}"
        group_selected = eligible_df.loc[picks].copy()
        group_selected["crop_rank"] = range(1, n_crops + 1)
        selected.append(group_selected)

selected_df = pd.concat(selected, ignore_index=True)

# Save selected cells
selected_df[
    ["Metadata_N_Replicate", label, "crop_rank", "dataset", "Metadata_Well", "Metadata_Site",
     "Metadata_ImageNumber", "Metadata_Cells_Number_Object_Number", "bbox_center_x", "bbox_center_y",
     "bbox_width", "bbox_height", "predicted_probability", "true_class_probability"]
].to_csv(crop_output_dir / "selected_cells.csv", index=False)

print(f"Skipped {len(skipped_debris)} crops with debris: {skipped_debris}")
print(selected_df.groupby(["Metadata_N_Replicate", label]).agg(
    n_wells=("Metadata_Well", "nunique"),
    min_true_class_probability=("true_class_probability", "min"),
))


# ## Load illumination corrected crops for the selected cells

# In[7]:


crops = {
    idx: {c: load_crop(row, c) for c in channel_names} for idx, row in selected_df.iterrows()
}
print(f"Loaded {len(crops)} crops")


# ## Create composite crops with the same display range per channel
# 
# The lower limit per channel is the background (median across crops of each crop's lower percentile) and the upper limit is computed from all pixels across all crops, so the same range and gamma are used for every cell.

# In[8]:


display_ranges = {}
for c in channel_names:
    background = np.median([np.percentile(crops[idx][c], background_percentile) for idx in crops])
    upper = np.percentile(np.concatenate([crops[idx][c].ravel() for idx in crops]), upper_percentile)
    display_ranges[c] = (background, upper)
print(display_ranges)


def make_composite(crop_channels: dict) -> np.ndarray:
    """Create an RGB composite from channel crops using the global display ranges."""
    composite = np.zeros((crop_size, crop_size, 3), dtype=np.float32)
    for c, color in channel_colors.items():
        low, high = display_ranges[c]
        scaled = np.clip((crop_channels[c] - low) / (high - low), 0, 1) ** display_gamma
        composite += scaled[..., None] * np.array(color, dtype=np.float32)
    return np.clip(composite, 0, 1)


composites = {idx: make_composite(crops[idx]) for idx in crops}

# Save each composite crop
for idx, row in selected_df.iterrows():
    file_name = (
        f"{row[label]}_{row['Metadata_N_Replicate']}_rank{row['crop_rank']}_"
        f"{row['Metadata_Well']}_site{row['Metadata_Site']}_cell{row['Metadata_Cells_Number_Object_Number']}.png"
    )
    plt.imsave(crop_output_dir / file_name, composites[idx])


# ## Create montage
# 
# Columns are biological replicates and rows are treatments, with a grid of crops per biological replicate and treatment.
# The held out biological replicate is outlined.

# In[9]:


# Spacing (in pixels) between crops within a block and between blocks, and margins for labels
crop_gap = 10
block_gap = 70
margin_left, margin_top, margin_bottom = 260, 140, 140

block_width = grid_cols * crop_size + (grid_cols - 1) * crop_gap
block_height = grid_rows * crop_size + (grid_rows - 1) * crop_gap
canvas_width = margin_left + len(replicates) * block_width + (len(replicates) - 1) * block_gap + block_gap // 2
canvas_height = margin_top + len(treatments) * block_height + (len(treatments) - 1) * block_gap + margin_bottom

canvas = np.zeros((canvas_height, canvas_width, 3), dtype=np.float32)
block_origins = {}
for i, treatment in enumerate(treatments):
    for j, replicate in enumerate(replicates):
        bx = margin_left + j * (block_width + block_gap)
        by = margin_top + i * (block_height + block_gap)
        block_origins[(treatment, replicate)] = (bx, by)
        group = selected_df[
            (selected_df["Metadata_N_Replicate"] == replicate) & (selected_df[label] == treatment)
        ].sort_values("crop_rank")
        for k, idx in enumerate(group.index):
            r, c = divmod(k, grid_cols)
            x = bx + c * (crop_size + crop_gap)
            y = by + r * (crop_size + crop_gap)
            canvas[y : y + crop_size, x : x + crop_size] = composites[idx]

dpi = 300
fig = plt.figure(figsize=(canvas_width / dpi, canvas_height / dpi), dpi=dpi, facecolor="black")
ax = fig.add_axes([0, 0, 1, 1])
ax.imshow(canvas, interpolation="none")
ax.set_axis_off()

# Column labels (biological replicates)
for j, replicate in enumerate(replicates):
    bx, _ = block_origins[(treatments[0], replicate)]
    is_holdout = replicate == holdout_replicate
    ax.text(
        bx + block_width / 2, margin_top - 45,
        f"{replicate} (held out)" if is_holdout else replicate,
        color=holdout_color if is_holdout else "white",
        ha="center", va="bottom", fontsize=16, fontweight="bold",
    )

# Row labels (treatments)
for i, treatment in enumerate(treatments):
    _, by = block_origins[(treatment, replicates[0])]
    ax.text(
        margin_left - 50, by + block_height / 2, treatment_labels[treatment],
        color=treatment_colors[treatment], ha="center", va="center",
        rotation=90, fontsize=18, fontweight="bold",
    )

# Outline the held out biological replicate for each treatment
for treatment in treatments:
    bx, by = block_origins[(treatment, holdout_replicate)]
    pad = 18
    ax.add_patch(
        Rectangle(
            (bx - pad, by - pad), block_width + 2 * pad, block_height + 2 * pad,
            fill=False, edgecolor=holdout_color, linewidth=2.5,
        )
    )

# Channel legend
legend_x = margin_left
legend_y = canvas_height - margin_bottom / 2
for c, color in channel_colors.items():
    ax.text(legend_x, legend_y, "F-actin" if c == "Actin" else c, color=color,
            ha="left", va="center", fontsize=16, fontweight="bold")
    legend_x += 330

# Scale bar on the bottom right crop
bx, by = block_origins[(treatments[-1], replicates[-1])]
bar_length = scale_bar_um / pixel_size_um
bar_x1 = bx + block_width - 30
bar_y = by + block_height - 40
ax.plot([bar_x1 - bar_length, bar_x1], [bar_y, bar_y], color="white", linewidth=4, solid_capstyle="butt")
ax.text(bar_x1 - bar_length / 2, bar_y - 15, f"{scale_bar_um} µm", color="white",
        ha="center", va="bottom", fontsize=11, fontweight="bold")

fig.savefig(montage_output_dir / "sc_crop_montage_top_predicted.png", dpi=dpi, facecolor="black")
fig.savefig(montage_output_dir / "sc_crop_montage_top_predicted.pdf", dpi=dpi, facecolor="black")
plt.close(fig)
print(f"Montage saved to {montage_output_dir}")

