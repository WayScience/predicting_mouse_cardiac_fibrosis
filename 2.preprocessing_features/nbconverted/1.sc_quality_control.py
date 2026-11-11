#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using coSMicQC. We filter the single cells by identifying outliers with z-scores, and use either combinations of features or one feature for each condition. 

# In[1]:


import json
import pathlib

import matplotlib

matplotlib.use("Agg")  # Ensures no GUI windows are opened
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cosmicqc import find_outliers, label_outliers
from cytodataframe import CytoDataFrame


# ## Set paths and variables

# In[2]:


# Set parameters for papermill to use for processing
plate_id = "CARD-CelIns-CX7_261022150001"


# In[3]:


# Directory with data
data_dir = pathlib.Path("./data/converted_profiles/").resolve(strict=True)

# Directory to save qc labeled data
qc_labeled_dir = pathlib.Path("./data/qc_labeled_profiles/")
qc_labeled_dir.mkdir(parents=True, exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path("./qc_figures/")
qc_fig_dir.mkdir(parents=True, exist_ok=True)

# Load only the thresholds for the plate_id in the sc_qc_thresholds.json file
qc_thresholds_path = pathlib.Path("./sc_qc_thresholds.json").resolve(strict=True)
with open(qc_thresholds_path) as f:
    qc_thresholds = json.load(f)
# stop if the plate has no thresholds to avoid an unclear error from cosmicqc downstream
if plate_id not in qc_thresholds:
    raise KeyError(f"No QC thresholds found for plate {plate_id} in {qc_thresholds_path}")
qc_thresholds = qc_thresholds[plate_id]
print(
    f"Loaded QC thresholds for {len(qc_thresholds)} conditions from {qc_thresholds_path}"
)

# Metadata columns to include in output data frame (based on compartment)
nuclei_metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
    "Image_PathName_OrigDNA",
    "Image_FileName_OrigDNA",
    "Image_PathName_OrigActin",
    "Image_FileName_OrigActin",
    "Nuclei_AreaShape_BoundingBoxMaximum_X",
    "Nuclei_AreaShape_BoundingBoxMaximum_Y",
    "Nuclei_AreaShape_BoundingBoxMinimum_X",
    "Nuclei_AreaShape_BoundingBoxMinimum_Y",
]
cells_metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Cells_Location_Center_X",
    "Metadata_Cells_Location_Center_Y",
    "Image_PathName_OrigDNA",
    "Image_FileName_OrigDNA",
    "Image_PathName_OrigActin",
    "Image_FileName_OrigActin",
    "Cells_AreaShape_BoundingBoxMaximum_X",
    "Cells_AreaShape_BoundingBoxMaximum_Y",
    "Cells_AreaShape_BoundingBoxMinimum_X",
    "Cells_AreaShape_BoundingBoxMinimum_Y",
]


# ## Load in plate to perform QC on

# In[4]:


# Set plate path
plate_path = pathlib.Path(f"{data_dir}/{plate_id}_converted.parquet").resolve(
    strict=True
)

# Load the converted plate data
plate_df = pd.read_parquet(plate_path)

# Print the shape of the DataFrame for each plate
print(f"Loaded plate: {plate_path.stem}, Shape: {plate_df.shape}")


# ## Confirm image paths

# In[5]:


# Confirm image paths already point to valid local files (no remapping needed for this plate)
example_path = pathlib.Path(
    plate_df["Image_PathName_OrigDNA"].dropna().iloc[0]
) / plate_df["Image_FileName_OrigDNA"].dropna().iloc[0]
if not example_path.is_file():
    raise FileNotFoundError(f"Example image path does not exist: {example_path}")

# Print example image path
print(plate_df["Image_PathName_OrigDNA"].dropna().iloc[0])


# In[6]:


# Print example image file
print(plate_df["Image_FileName_OrigDNA"].dropna().iloc[0])


# ## Set mapping for outlines

# In[7]:


# Define available compartments
compartments = ["Nuclei", "Cells"]

# Create a dictionary of mappings for each compartment.
# The outline TIFFs (e.g. "NucleiOutlines_PLATE_WELL_SITE.tiff") use the raw,
# unpadded site number. The original images (e.g. "PLATE_WELLf09d4.TIF") use an
# "f"-prefixed, zero-padded site token, a channel digit, and an uppercase ".TIF"
# extension, so the original-image regex has to be built to match that format
# (case-insensitively) rather than reusing the outline's site formatting.
mapping_dict = {
    comp: {
        rf"{comp}Outlines_{record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}_{record['Image_Metadata_Site']}.tiff": (
            rf"(?i){record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}"
            rf"f{int(record['Image_Metadata_Site']):02d}d\d\.tif"
        )
        for record in plate_df[
            ["Image_Metadata_Plate", "Image_Metadata_Well", "Image_Metadata_Site"]
        ].to_dict(orient="records")
    }
    for comp in compartments
}


# ## Detect over-segmented nuclei using shape and intensity

# In[8]:


# Set mapping for Nuclei compartment
outline_to_orig_mapping = mapping_dict["Nuclei"]


# In[9]:


# Find large nuclei outliers for the current plate
oversegmented_nuclei_outliers = find_outliers(
    df=plate_df,
    metadata_columns=nuclei_metadata_columns,
    feature_thresholds=qc_thresholds.get("oversegmented_nuclei", {}),
)

# Convert to CytoDataFrame for outline viewing
oversegmented_nuclei_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(oversegmented_nuclei_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": False,
        "brightness": 5,
        "outline_color": (180, 30, 180),
    },
)[
    [
        "Nuclei_AreaShape_Solidity",
        "Nuclei_Intensity_MassDisplacement_DNA",
        "Image_FileName_OrigDNA",
    ]
]

oversegmented_nuclei_outliers_cdf.sort_values(
    by="Nuclei_Intensity_MassDisplacement_DNA", ascending=False
).head(5)


# ## Detect mis-segmented background as nuclei

# In[10]:


# Find low intensity nuclei outliers for the current plate
low_intensity_outliers = find_outliers(
    df=plate_df,
    metadata_columns=nuclei_metadata_columns,
    feature_thresholds=qc_thresholds.get("low_intensity_nuclei", {}),
)

# Convert to CytoDataFrame for outline viewing
low_intensity_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(low_intensity_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
        # expand crop around the nucleus center to see surrounding context
        "offset_bounding_box": {"x_min": -100, "y_min": -100, "x_max": 100, "y_max": 100},
    },
)[
    [
        "Nuclei_Intensity_MeanIntensity_DNA",
        "Image_FileName_OrigDNA",
    ]
]

low_intensity_outliers_cdf.sort_values(
    by="Nuclei_Intensity_MeanIntensity_DNA", ascending=False
).head(5)


# ## Detect mis-segmented debris (bright) as nuclei

# In[11]:


# Find high intensity nuclei outliers for the current plate
high_intensity_outliers = find_outliers(
    df=plate_df,
    metadata_columns=nuclei_metadata_columns,
    feature_thresholds=qc_thresholds.get("high_intensity_nuclei", {}),
)

# Convert to CytoDataFrame for outline viewing
high_intensity_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(high_intensity_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
        # expand crop around the nucleus center to see surrounding context
        "offset_bounding_box": {"x_min": -100, "y_min": -100, "x_max": 100, "y_max": 100},
    },
)[
    [
        "Nuclei_Intensity_MeanIntensity_DNA",
        "Image_FileName_OrigDNA",
    ]
]

high_intensity_outliers_cdf.sort_values(
    by="Nuclei_Intensity_MeanIntensity_DNA", ascending=True
).head(5)


# ## Detect under-segmented cells

# ### Set up outlines for cells

# In[12]:


# Set mapping for Cells compartment
outline_to_orig_mapping = mapping_dict["Cells"]


# In[13]:


# Find under-segmented cells outliers for the current plate
small_cells_outliers = find_outliers(
    df=plate_df,
    metadata_columns=cells_metadata_columns,
    feature_thresholds=qc_thresholds.get("undersegmented_cells", {}),
)

# Convert to CytoDataFrame for outline viewing
small_cells_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(small_cells_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
        # expand crop around the cell center to see surrounding context
        "offset_bounding_box": {"x_min": -100, "y_min": -100, "x_max": 100, "y_max": 100},
    },
)[
    [
        "Cells_AreaShape_Area",
        "Image_FileName_OrigActin",
        "Image_FileName_OrigDNA",
    ]
]

small_cells_outliers_cdf.sort_values(
    by="Cells_AreaShape_Area", ascending=False
).head(10)


# ## Detect out of focus/blurry or mis-segmented cells
# 
# Using low levels of granularity will detect cells without sharp features and contrast (smoother) that we expect from in-focus cells.
# Using abnormally high levels detect mis-segmented cells of background due to issues with nuclei segmentations.

# In[14]:


# Find blurry cells outliers for the current plate
blurry_cells_outliers = find_outliers(
    df=plate_df,
    metadata_columns=cells_metadata_columns,
    feature_thresholds=qc_thresholds.get("blurry_cells", {}),
)

# Convert to CytoDataFrame for outline viewing
blurry_cells_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(blurry_cells_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
        # expand crop around the cell center to see surrounding context
        "offset_bounding_box": {"x_min": -200, "y_min": -200, "x_max": 200, "y_max": 200},
    },
)[
    [
        "Cells_Granularity_1_Actin",
        "Image_FileName_OrigActin",
        "Image_FileName_OrigDNA",
    ]
]

blurry_cells_outliers_cdf.sort_values(
    by="Cells_Granularity_1_Actin", ascending=True
).head(10)


# In[15]:


# Find misegmented background cells outliers for the current plate
missegmented_background_cells_outliers = find_outliers(
    df=plate_df,
    metadata_columns=cells_metadata_columns,
    feature_thresholds=qc_thresholds.get("missegmented_background_cells", {}),
)

# Convert to CytoDataFrame for outline viewing
missegmented_background_cells_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(missegmented_background_cells_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
        # expand crop around the cell center to see surrounding context
        "offset_bounding_box": {"x_min": -200, "y_min": -200, "x_max": 200, "y_max": 200},
    },
)[
    [
        "Cells_Granularity_1_Actin",
        "Image_FileName_OrigActin",
        "Image_FileName_OrigDNA",
    ]
]

missegmented_background_cells_outliers_cdf.sort_values(
    by="Cells_Granularity_1_Actin", ascending=True
).head(10)


# ## Detect mis-segmented cells based on low intensity and how irregular the boundaries are
# 
# We have seen that due to low staining in some Actin FOVs that segmentation will mis-segment to include more background than it should.

# In[16]:


# Find mis-segmented cells outliers for the current plate
missegmented_cells_outliers = find_outliers(
    df=plate_df,
    metadata_columns=cells_metadata_columns,
    feature_thresholds=qc_thresholds.get("missegmented_cells", {}),
)

# Convert to CytoDataFrame for outline viewing
missegmented_cells_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(missegmented_cells_outliers),
    data_outline_context_dir=f"../1.cellprofiler_processing/cp_output/{plate_id}/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 1,
        "outline_color": (180, 30, 180),
        # expand crop around the cell center to see surrounding context
        "offset_bounding_box": {"x_min": -200, "y_min": -200, "x_max": 200, "y_max": 200},
        
    },
)[
    [
        "Image_FileName_OrigActin",
        "Image_FileName_OrigDNA",
        "Cells_Intensity_MeanIntensity_Actin",
        "Cells_AreaShape_FormFactor",
    ]
]

missegmented_cells_outliers_cdf.sort_values(
    by="Cells_Intensity_MeanIntensity_Actin", ascending=False
).head(10)


# ## Label outliers and save data for the plate

# In[17]:


# Label every single cell with each QC condition for this plate from the
# thresholds file, adding Metadata_cqc_<condition>_is_outlier columns
plate_cleaned_name = plate_df["Image_Metadata_Plate"].iloc[0]
qc_labeled_df = label_outliers(
    df=plate_df,
    feature_thresholds=qc_thresholds,
    export_path=f"{qc_labeled_dir}/{plate_cleaned_name}_qc_labeled.parquet",
)

# Determine cells that failed any QC condition
qc_fail_cols = [
    col
    for col in qc_labeled_df.columns
    if col.startswith("Metadata_cqc_") and col.endswith("_is_outlier")
]
failed_mask = qc_labeled_df[qc_fail_cols].any(axis=1)

# Calculate the total percentage of single cells that failed QC
total_cells = qc_labeled_df.shape[0]
total_failed = failed_mask.sum()
percent_failed = (total_failed / total_cells) * 100 if total_cells > 0 else 0

# Print the number of cells that failed each QC condition
print(qc_labeled_df[qc_fail_cols].sum())

# Verify the result and include the percentage of failed QC
print(
    f"QC labeled data saved for plate {plate_cleaned_name}. "
    f"Shape: {qc_labeled_df.shape}. Total percent single cells failed QC: {percent_failed:.2f}%"
)


# In[18]:


# Set the default value to 'Passed nuclei QC' (only nuclei QC conditions are shown
# in this plot, so cells failing a cell QC condition are still labeled as passed)
plate_df["Outlier_Status"] = "Passed nuclei QC"

# Label each nuclei outlier type
plate_df.loc[
    plate_df.index.isin(oversegmented_nuclei_outliers.index), "Outlier_Status"
] = "Over-segmented nuclei"

plate_df.loc[plate_df.index.isin(low_intensity_outliers.index), "Outlier_Status"] = (
    "Low intensity nuclei"
)

plate_df.loc[plate_df.index.isin(high_intensity_outliers.index), "Outlier_Status"] = (
    "High intensity nuclei"
)

# Create scatter plot
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=plate_df,
    x="Nuclei_AreaShape_Solidity",
    y="Nuclei_Intensity_MassDisplacement_DNA",
    hue="Outlier_Status",
    palette={
        "Passed nuclei QC": "#006400",
        "Over-segmented nuclei": "#990090",
        "Low intensity nuclei": "#D55E00",
        "High intensity nuclei": "#0072B2",
    },
    alpha=0.6,
)

# Add threshold lines
plt.axvline(
    x=oversegmented_nuclei_outliers["Nuclei_AreaShape_Solidity"].max(),
    color="r",
    linestyle="--",
    label="Min. threshold for Nuclei Solidity",
)
plt.axhline(
    y=oversegmented_nuclei_outliers["Nuclei_Intensity_MassDisplacement_DNA"].min(),
    color="b",
    linestyle="--",
    label="Min. threshold for Nuclei Mass Displacement",
)

# Add title and labels
plt.title(f"Nuclei Solidity vs. Mass Displacement ({plate_id})")
plt.xlabel("Nuclei Solidity")
plt.ylabel("Nuclei Mass Displacement (DNA)")
plt.tight_layout()

# Show the legend
plt.legend(loc="upper left", bbox_to_anchor=(0, 1.0), prop={"size": 10})

# Save the figure with plate_id in filename
plt.savefig(qc_fig_dir / f"{plate_id}_nuclei_outliers.png", dpi=500)
plt.close()


# In[19]:


# Create KDE plot of all cells
plt.figure(figsize=(10, 6))
sns.kdeplot(
    x=plate_df["Cells_AreaShape_Area"],
    fill=True,
    color="#4682B4",  # steel blue
    alpha=0.6,
)

# Add threshold line for small cells
plt.axvline(
    x=small_cells_outliers["Cells_AreaShape_Area"].max(),
    color="r",
    linestyle="--",
    label=f'Small cells threshold: < {small_cells_outliers["Cells_AreaShape_Area"].max():.2f}',
)

# Labels and title
plt.ylabel("Density")
plt.xlabel("Cells Area")
plt.title(f"Distribution of Cell Areas ({plate_id})")
plt.legend()
plt.tight_layout()

# Save figure
plt.savefig(qc_fig_dir / f"{plate_id}_small_cell_outliers.png", dpi=500)
plt.close()

