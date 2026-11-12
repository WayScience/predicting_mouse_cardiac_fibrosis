# Preprocess features

In this module, we convert the CellProfiler outputs from the [CellProfiler processing module](../1.cellprofiler_processing/) into single-cell profiles, remove poor-quality segmentations, and generate normalized and feature-selected single-cell and well-level (bulk) profiles.

## Convert SQLite to parquet

The notebook [0.convert_cytotable.ipynb](./0.convert_cytotable.ipynb) uses [CytoTable](https://github.com/cytomining/CytoTable) to merge the per-compartment tables (nuclei, cells, and cytoplasm) from the CellProfiler SQLite output into one row per single cell, saved as a parquet file in `data/converted_profiles`.
We add site, cell count, and image path columns to the default CellProfiler preset, and add the `Metadata_` prefix to the location centroid and cell count columns.

## Single-cell quality control (QC)

The notebook [1.sc_quality_control.ipynb](./1.sc_quality_control.ipynb) uses [coSMicQC](https://github.com/cytomining/coSMicQC) to find single cells that are outliers (by z-score) for one feature or a combination of features.
The thresholds per plate are set in [sc_qc_thresholds.json](./sc_qc_thresholds.json), and we view example outliers with their segmentation outlines using [CytoDataFrame](https://github.com/cytomining/CytoDataFrame) to confirm each threshold is removing poor-quality segmentations.

| QC condition | Feature(s) | Purpose |
| --- | --- | --- |
| Over-segmented nuclei | `Nuclei_AreaShape_Solidity`, `Nuclei_Intensity_MassDisplacement_DNA` | Nuclei segmentations that are irregularly shaped and include more than one nucleus |
| Low intensity nuclei | `Nuclei_Intensity_MeanIntensity_DNA` | Background mis-segmented as nuclei |
| High intensity nuclei | `Nuclei_Intensity_MeanIntensity_DNA` | Bright debris mis-segmented as nuclei |
| Under-segmented cells | `Cells_AreaShape_Area` | Cells segmentations that are too small |
| Blurry cells | `Cells_Granularity_1_Actin` (low) | Out-of-focus cells without sharp features |
| Mis-segmented background cells | `Cells_Granularity_1_Actin` (high) | Background mis-segmented as cells due to issues with nuclei segmentation |
| Mis-segmented cells | `Cells_Intensity_MeanIntensity_Actin`, `Cells_AreaShape_FormFactor` | Cells segmentations with low Actin staining that include more background than they should |

We then use the coSMicQC `label_outliers` function to label each single cell with a `Metadata_cqc_<condition>_is_outlier` column for every QC condition in the thresholds file, and the labeled profiles are saved in `data/qc_labeled_profiles`.
Single cells that fail any QC condition are removed in the single-cell and bulk processing notebooks.
Figures showing the QC thresholds are saved in the [qc_figures](./qc_figures/) folder.

## Single-cell profiles

The notebook [2.single_cell_processing.ipynb](./2.single_cell_processing.ipynb) uses [pycytominer](https://github.com/cytomining/pycytominer) to process the single cells that passed QC:

1. **Annotate:** Add the plate map metadata (treatment and biological replicate) to each single cell, and add the `Metadata_` prefix to the site and image file/path name columns so they are not treated as features.
2. **Normalize:** Normalize the features to all single cells on the plate (`mad_robustize`).
3. **Feature selection:** Remove features that are in the [blocklist](./blocklist_features.txt), have any missing values, have low variance, are dominated by one value (frequency threshold), or are highly correlated.

The single-cell profiles are saved in `data/single_cell_profiles`.

> [!NOTE]
> We use an updated [blocklist](./blocklist_features.txt) instead of the default pycytominer blocklist.
> The default blocklist uses the standard Cell Painting channel names (`AGP`, `DNA`, `ER`, `Mito`, and `RNA`), but our channels are named `Actin`, `DNA`, `ER`, `Golgi`, and `Mito`.
> Our blocklist contains the same nuclei features as the default (`Correlation_Manders`, `Correlation_RWC`, and `Granularity_14` to `Granularity_16`), renamed to match our channels.

> [!IMPORTANT]
> We include the `frequency_threshold` operation in feature selection.
> In older pycytominer versions (e.g., 1.3.1), `variance_threshold` removed features dominated by one value (e.g., Euler number, which is the same for 97% of cells).
> In newer versions (e.g., 1.7.1, used here), `variance_threshold` only removes features with low raw variance, and the old behavior moved to `frequency_threshold`.
> Without it, features that are mostly one value pass feature selection with extreme values (~10^18) after `mad_robustize` (their MAD is 0), which breaks model training.

## Bulk profiles

The notebook [3.bulk_processing.ipynb](./3.bulk_processing.ipynb) uses pycytominer to create well-level profiles from the single cells that passed QC:

1. **Aggregate:** Take the median of each feature per well.
2. **Annotate:** Add the plate map metadata to each well.
3. **Normalize:** Normalize the features to the vehicle-treated wells (negative controls) (`mad_robustize`).
4. **Feature selection:** Remove features that are in the blocklist, have any missing values, have low variance, are dominated by one value, or are highly correlated (> 0.95), then remove features with low variance or dominated by one value within the vehicle-treated wells.
5. **Spherize:** Whiten the profiles using the vehicle-treated wells as the reference (`ZCA-cor`).

The bulk profiles are saved in `data/bulk_profiles`.

## Perform preprocessing on data

To run all of the notebooks in order, run the bash script [preprocess_features.sh](./preprocess_features.sh) using the command below:

```bash
bash preprocess_features.sh
```

This script activates the preprocessing environment (uv), converts the notebooks to Python scripts in the [nbconverted](./nbconverted/) folder, and runs them in order.
If any stage fails, the script stops and exits with an error so later stages do not run on stale outputs.

> [!NOTE]
> The `data` folder is not tracked on GitHub because the profiles are too large.
> Run this module to generate the profiles locally.
