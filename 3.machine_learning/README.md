# Machine learning

In this module, we train a binary logistic regression model to predict whether single cells are TGF-β treated (fibrotic) or vehicle treated (non-failing) using the feature selected single-cell profiles from the [preprocessing module](../2.preprocessing_features/).

## Split data

The notebook [0.split_data.ipynb](./0.split_data.ipynb) generates the data splits:

1. **Location features:** We drop the two Location features that pass feature selection (`Cytoplasm_Location_MaxIntensity_Y_Actin` and `Cells_Location_MaxIntensity_X_Actin`), since they describe where a cell sits in the image, not its biology.
2. **Holdout:** We hold out one whole biological replicate (randomly selected) and one random well per treatment from each remaining biological replicate.
3. **Training and testing:** We split the remaining data 70/30 into training and testing data, stratified by treatment.

All random selections use `random_state = 0`.
With this seed, we hold out all of N3 plus wells B07 and G04 (N1), D10 and F05 (N2), and B06 and F03 (N4), listed as TGFB then Veh.

### Data split plate map

The notebook [7.vis_data_split_platemap.ipynb](./7.vis_data_split_platemap.ipynb) plots the data splits on the plate map ([figures/data_split_platemap.png](./figures/data_split_platemap.png)).
Held out wells show one color, and the other wells show pie charts of their cells in training data, removed by downsampling, and testing data, since we split single cells (not wells).
A bar plot below shows the cells per biological replicate, treatment, and data split, and [data_split_summary/cells_per_well_per_split.csv](./data_split_summary/cells_per_well_per_split.csv) has the counts per well.

## Train models

The notebook [1.train_models.ipynb](./1.train_models.ipynb) downsamples the training data to the smallest class and trains two binary logistic regression models:

1. **Final:** Trained on the training data.
2. **Shuffled:** Trained on the training data with each feature column shuffled independently, used as a baseline.

Both models use an elastic net penalty (`saga` solver, `max_iter = 1000`, `class_weight = "balanced"`).
Hyperparameters (`C` from 0.01 to 10 and `l1_ratio` from 0 to 1) are optimized using `RandomizedSearchCV` with stratified 5-fold cross validation and weighted F1 score.
The vehicle (`Veh`) class is the positive class (encoded as 1), so predicted probabilities are the probability of a vehicle (non-failing) prediction.

## Extract and visualize model performance

The notebook [2.extract_model_performance.ipynb](./2.extract_model_performance.ipynb) applies the final and shuffled models to the training (downsampled), testing, and holdout datasets and extracts the precision-recall curves and per-cell predicted probabilities to the [performance_metrics](./performance_metrics/) folder.

The notebook [3.vis_model_performance.ipynb](./3.vis_model_performance.ipynb) uses R to plot the precision-recall curves, the distribution of predicted probabilities (overall and per biological replicate), and the accuracy per biological replicate and treatment.
For accuracy, the holdout data is split into the held out biological replicate (`holdout_replicate`) and the held out wells (`holdout_well`).
Figures are saved to the [figures](./figures/) folder.

## Extract and visualize model coefficients

The notebook [4.extract_model_coef.ipynb](./4.extract_model_coef.ipynb) extracts the coefficient for each feature from the final model and saves them to the [coeff_data](./coeff_data/) folder.
Positive coefficients predict vehicle (non-failing) cells and negative coefficients predict TGF-β (fibrotic) cells.
We split out the compartment, feature group, and channel from each feature name, where the channel is the first channel found in the feature name (for correlation features, this is the first of the two channels) and features without a channel (e.g., AreaShape and Neighbors) are labelled as `other`.

The notebook [5.vis_model_coef.ipynb](./5.vis_model_coef.ipynb) uses R to visualize the coefficients, following the figures from [cellpainting_predicts_cardiac_fibrosis](https://github.com/jenna-tomkinson/cellpainting_predicts_cardiac_fibrosis):

1. Heatmaps of the top absolute coefficient per compartment, feature group, and channel, for all features and for the TGF-β and vehicle features separately
2. A ranked plot of all coefficients from most positive (vehicle) to most negative (TGF-β), with the top five features per class labelled
3. A bar plot of the top 15 features per class

Figures are saved to the [figures/coefficients](./figures/coefficients/) folder.

## Single-cell crop montage

The notebook [6.generate_sc_crop_montage.ipynb](./6.generate_sc_crop_montage.ipynb) builds a montage of single-cell crops that shows how TGF-β and vehicle cells differ across biological replicates.
Columns show biological replicates, rows show treatments, and each block holds nine crops.
An outline marks the held out biological replicate (N3).

- **Cells:** We use only cells the model did not train on (testing and holdout data) and pick the correctly predicted cells with the highest probability for their ground truth class. We take at most one cell per image and two per well.
- **Crops:** We crop 400 x 400 pixels (1 pixel = 1 µm, from the image metadata) around the center of the cell bounding box, not the nucleus, so the whole cell fits. We shift crops at the image edge to stay inside the image and skip cells that do not fit.
- **Debris:** We skip whole images when a crop has a connected DNA area (> 2,000 intensity) over 10,000 pixels (nuclei clusters stay under 6,000 pixels). We also skip one image (G04 site 11) with a faint streak the filter misses.
- **Channels:** We show F-actin (green), Mito (magenta), and DNA (cyan). F-actin and DNA carry the top model features, and Mito has stronger top features than ER.
- **Display:** We apply illumination correction and use one display range and gamma (0.6) per channel for every crop, so intensity differences between treatments remain. The lower limit is the background, which keeps it black. Vehicle cells have about half the F-actin intensity of TGF-β cells, so the gamma keeps their F-actin visible.

The notebook saves the selected cells to [sc_crops/selected_cells.csv](./sc_crops/selected_cells.csv), each composite crop to [sc_crops](./sc_crops/), and the montage as PNG and PDF to [figures/sc_crop_montage](./figures/sc_crop_montage/).
We only track the selected cells CSV and the montage PDF on GitHub; the crop PNGs and the montage PNG are too large, so run the notebook to generate them.

> [!NOTE]
> This notebook reads the raw images from the mounted image directory, so mount it before running.

## Perform machine learning on data

To run all of the notebooks in order, run the bash script [run_machine_learning.sh](./run_machine_learning.sh) using the command below:

```bash
bash run_machine_learning.sh
```

The script runs each Python step with `uv run` in the machine learning environment, so it does not depend on which environment the shell has active.
It converts the notebooks to Python and R scripts in the [nbconverted](./nbconverted/) folder, runs the Python scripts, then activates the R environment (conda) to run the R scripts.
If any stage fails, the script stops so later stages do not run on stale outputs.

> [!NOTE]
> The `data_splits` folder is not tracked on GitHub because the data splits are too large.
> Run this module to generate the data splits locally.
