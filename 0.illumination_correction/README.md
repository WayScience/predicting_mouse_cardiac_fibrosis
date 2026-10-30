# Illumination correction & whole image quality control (QC)

In this module, we calculate illumination correction (IC) functions for each channel using a modified version of the CellProfiler illumination correction pipeline from the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.
The IC functions are applied to the images in the [CellProfiler processing module](../1.cellprofiler_processing/) prior to segmentation and feature extraction.

## Create LoadData CSV

The notebook [0.create_loaddata_csv.ipynb](./0.create_loaddata_csv.ipynb) finds the raw images for each plate and creates a CellProfiler `LoadData` CSV (one row per image set, grouped by plate, well, and site) in the [loaddata_csvs](./loaddata_csvs/) folder.
Each image set contains five channels, which are mapped to the following CellProfiler image names:

| Image channel | CellProfiler image name |
| --- | --- |
| d0 | OrigActin |
| d1 | OrigMito |
| d2 | OrigGolgi |
| d3 | OrigER |
| d4 | OrigDNA |

## CellProfiler pipeline

The [IC pipeline](./pipeline/illum.cppipe) performs two tasks:

1. **Whole image QC:** We measure the percentage of pixels at the maximum intensity value (`PercentMaximal`) in each channel. If more than 1% of the pixels in any channel are at the maximum intensity, the whole image set is skipped. We remove these over-saturated images so they do not skew the IC calculation. We do not filter for blur here as blurry images do not impact the IC calculation as much, and we remove them downstream during single-cell QC.
2. **Calculate IC functions:** We apply the same [IC pipeline](https://github.com/broadinstitute/imaging-platform-pipelines/blob/master/JUMP_production/JUMP_illum_LoadData_v1.cppipe) that the Broad Institute used for the JUMP Cell Painting Consortium, where per channel, images are downsampled (0.2x), then an IC function is calculated across all passing image sets in the plate (`Regular` method, median filter smoothing with a filter size of 20). The function is upsampled back to the original image size and saved as a `.npy` file per plate and channel in the `illum_directory` folder.

## Perform IC on data

Before running, update the image paths in [0.create_loaddata_csv.ipynb](./0.create_loaddata_csv.ipynb) and [1.illum_correct.ipynb](./1.illum_correct.ipynb) to where the raw images are located on your machine.

To create the LoadData CSV and calculate the IC functions, run the bash script [perform_ic.sh](./perform_ic.sh) using the command below:

```bash
source perform_ic.sh
```

This script activates the CellProfiler environment (`mouse_fibrosis_cp_env`), converts the notebooks to Python scripts in the [nbconverted](./nbconverted/) folder, and runs them in order.
CellProfiler runs in parallel per plate, and logs are saved to the `logs` folder.
