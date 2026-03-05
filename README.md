# Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the code and data for the paper:

**"Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex"**  
Authors: Ao Ma, Zhengjue Ren, Yuqing Li, Shui Yu, Changsong Zhou, Claus C. Hilgetag, Yuhan Chen  
Journal/Conference: Under review

## 📋 Overview

This study investigates how state-dependent brain activity in macaque cortex anchors to three natural structural axes. The repository includes all processed data and analysis code used in the paper.

## 📊 Data

All data files are organized in the `data/` directory.

```
data/
├──  BrainMesh_Monkey_F99.nv
├──  distance_matrix.mat
├──  Empirical_and_Predicted_Rank_Difference.mat
├──  Empirical_Axis3_Afferent_WR.mat
├──  Empirical_Efferent_WR.mat
├──  Empirical_Spectral_Gradient_during_Anesthesia_and_Awake.mat
├──  fmri_Gradient1_and2.mat
├──  Function_Domain.txt
├──  label_467.mat
├──  name91.mat
├──  Power_for_all_frequency_bands.mat
├──  Power_for_all_frequency_bands.txt
├──  Predicted_Axis2_T1_T2.mat
├──  Predicted_Axis3_Afferent_WR.mat
├──  Predicted_Efferent_WR.mat
├──  Predicted_Spectral_Gradient_during_Anesthesia.txt
├──  Predicted_Spectral_Gradient_during_Anesthesia_and_Awake.mat
├──  Predicted_Spectral_Gradient_during_Awake.txt
├──  Receptor.mat
└──  Structural_Connectome.mat       

```

### Data Description

| File                                         | Dimensions | Description                                                                      |
| -------------------------------------------- | ---------- | -------------------------------------------------------------------------------- |
| `BrainMesh_Monkey_F99.nv`                    |            | Template of macaque cortex for BrainNet Viewer                                   |
| `distance_matrix.mat`                        | 91 * 91    | Pairwise edge distances between brain regions                                    |
| `Empirical_Predicted_Rank_Difference.mat`    | 91 * 2     | Column 1 for empirical, column 2 for predicted rank difference                   |
| `Empirical_Axis3_Afferent_WR.mat`            | 91 * 1     | Empirical axis3,  afferent wiring range                                          |
| `Empirical_Spectral_Gradient_Anes_Awake.mat` | 467 * 2    | Column 1 for anesthesia, column 2 for awake empirical spectral gradient          |
| `fmri_Gradient1_and2.mat`                    | 91 * 2     | The first (column 1), the second (column 2) principal component of fMRI gradient |
| `Function_Domain.txt`                        |            | Map for 5 function domains                                                       |
| `label_467.mat`                              | 467 * 1    | Label of brain region for each channel                                           |
| `name91.mat`                                 | 91 * 1     | Name for each brain region                                                       |
| `Predicted_Axis2_T1_T2.mat`                  | 91 * 1     | Empirical axis2,  T1W/T2W                                                        |
| `Predicted_Axis3_Afferent_WR.mat`            | 91 * 1     | Predicted axis3,  afferent wiring range                                          |
| `Predicted_Spectral_Gradient_Anes.txt`       |            | Predicted spectral gradient during anesthesia state                              |
| `Predicted_Spectral_Gradient_Anes_Awake.mat` | 91 * 2     | Predicted spectral gradient for anethesia (column 1), for wakefulness (column2)  |
| `Predicted_Spectral_Gradient_Awake.txt`      |            | Predicted spectral gradient during awake state                                   |
| `Receptor.mat`                               | 91 * 1     | First principal component of 14 receptors                                        |
| `Structural_Connectome.mat`                  | 91 * 91    | Directed weighted structural network                                             |



## 💻 Code

All analysis code is organized in the `code/` directory.

```
code/
├── cal_rank_difference/
│   └── cal_rank_difference.m
│
├── Dual_origin_diffusion_model/
│   ├── corr_lp.m
│   ├── Dual_origin_diffusion_model_code.m
│   └── Geodesic_distance_code.m           
│
├── Structure_dimensionality_reduction/
│   ├── dimensionality_reduction.m
│   └── code_for_UMAP.ipynb
│        
├── SVD/
│	└── SVD_code.m
│
└── generative_model/
 	├── 1_first_step/
    │   ├── pre_weight_1_40_40.m
    │   ├── pre_weight_2_40_40.m
    │   └── trans_matrix_40area.m
    │
    ├── 2_second_step/
    │   ├── Cal_posterior_probability.m
    │   ├── data_array_40area.m
    │   ├── data_array2_40area.m
    │   └── trans_matrix_40area.m   
    │  
    └── 3_validation/
        ├── cal_roc.m
        ├── matrix_sort.m
        └── validation_code.m

```

### Key Scripts and Their Outputs
#### cal_rank_difference
1. cal_rank_difference.m
Calculation of rank difference based on spectral gradient. Quantifies the regional rank changes of the spectral gradient between the anesthesia and awake states. 

Mapped to Paper:
- **Data**: `Empirical_Predicted_Rank_Difference.mat`.
- **Figures**: Fig. 6B (Histogram of rank changes), Fig. 6C (Rank changes vs. Afferent WR), Fig. 6D (Map of empirical rank difference), Fig. 6H (Predicted vs. empirical rank changes).

#### Dual_origin_diffusion_model
1. corr_lp.m
A function that cleans invalid data and calculates correlations.
2. Dual_origin_diffusion_model_code.m
Implements the core dual-origin diffusion model. It generates composite diffusion fields from two source points.
3. Geodesic_distance_code.m
Computes the shortest geodesic distance along the M132 cortical surface mesh using Dijkstra's algorithm.

Mapped to Paper:
- **Data**: `BrainMesh_Monkey_F99.nv`, `Predicted_Axis2_T1_T2.mat`.
- **Figures**: Fig. 3C, Fig. 3D (Neuron density prediction); Fig. 3E, Fig. 3G, Fig. 3H (T1w/T2w prediction); Fig. 3I (Correlation scatter plots); Fig. S5C, Fig. S5E, Fig. S5F, Fig. S5H (Source pair validation distributions).

#### Structure_dimensionality_reduction
1. dimensionality_reduction.m
Dimensionality reduction of structured data using t-SNE, CMDS, or UMAP, followed by k-means or GMM classification.
2. code_for_UMAP.ipynb
Dimensionality reduction using UMAP.

Mapped to Paper:
- **Data**: `Predicted_Axis2_T1_T2.mat`, `Predicted_Axis3_Afferent_WR.mat`.
- **Figures**: Fig. 2B (Correlation matrix), Fig. 2C (t-SNE space); Fig. S4C-I (Validation of categorization across multiple reduction/clustering methods).

#### SVD
1. SVD_code.m
Dimensionality reduction was performed on the data matrix using Singular Value Decomposition (SVD).

Mapped to Paper:
- **Data**: `Receptor.mat`, `Empirical_Spectral_Gradient_Anes_Awake.mat`.
- **Figures**: Fig. 1A, Fig. 1C (Spectral gradient); Fig. S1F (Variance explained for power); Fig. S4A, Fig. S4B (Receptor PC1).

#### generative_model
##### 1_first_step(Symmetric connection prediction)
1. pre_weight_1_40_40.m
Code for predicting log(weight) from CD and SD based on experimental connectivity matrices.
2. pre_weight_2_40_40.m
Prediting 91 * 91 connectome.
3. trans_matrix_40area.m
Transform  91 * 91 matrix into 40 * 40 according to experimental indices.

Mapped to Paper: 
- **Figures**: Fig. 4B (SD and CD exponential decay fitting), Fig. 4C (Model diagram).

##### 2_second_step(Bayesian direction inference)
1. Cal_posterior_probability.m
Calculate posterior probability.
2. data_array_40area.m
Vectorize the matrix into a column vector and remove the diagonal elements.
3. data_array2_40area.m
The matrix was retained as 91 × 40 according to experimental indices, then converted into a column vector with diagonal elements removed.
4. trans_matrix_40area.m   
Transform  91 * 91 matrix into 40 * 40 according to experimental indices.

Mapped to Paper: 
- **Data**: `Structural_Connectome.mat`, `Predicted_Axis3_Afferent_WR.mat`.
- **Figures**: Fig. 4D (Predicted afferent WR), Fig. 4E (Predicted connection matrix), Fig. 4F, Fig. 4G; Fig. S6A (Posterior probability); Fig. S6H,I (Predicted efferent WR).

##### 3_validation(Model evaluation)
1. cal_roc.m
Calculate the ROC curve.
2. matrix_sort.m
Sort the vector along a specific dimension.
3. validation_code.m
Randomly generate train-test splits for validation.

Mapped to Paper: 
- **Figures**: Fig. S6C (ROC curve), Fig. S6D (AUC distribution), Fig. S6E (Sensitivity vs random benchmark), Fig. S6F, Fig. S6G (Sensitivity in train/test splits).

## 🚀 Quick Start

To run the analysis code and reproduce the results, the following software environments and toolboxes are required:

**1. MATLAB Environment**
* Recommended version: `MATLAB R2023b`.
* **Chronux Toolbox**: Required for ECoG multi-taper power spectra analysis. Download at: [http://chronux.org/](http://chronux.org/)
* **BrainNet Viewer**: Required for 3D cortical surface visualization. Download at: [https://www.nitrc.org/projects/bnv/](https://www.nitrc.org/projects/bnv/)

**2. Python Environment**
* Required for running the UMAP dimensionality reduction notebook (`code_for_UMAP.ipynb`).
* Recommended version: `Python 3.8+
* Install required dependencies via pip:
### Installation

```bash

# Install dependencies

pip install pandas 2.2.3
pip install numpy 2.2.5
pip install scipy 1.15.2
pip install umap-learn 0.5.7

```

## 📝 Notes

- **Data size**: Most files are small CSV/text files (<1MB each)
- **Spatial coordinates**: All region-level data uses M132 atlas parcellation (91 regions)

## 📚 Citation

If you use this code or data in your research, please cite:

```
Ma, A., Ren, Z., Li, Y., Yu, S., Zhou, C., Hilgetag, C.C., & Chen, Y. (2026). 
Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex.
```

## 📄 License

This project is licensed under the MIT License.  
Copyright (c) 2026 Yuhan Chen, Ao Ma, Zhengjue Ren, Yuqing Li, Shui Yu, Changsong Zhou, Claus C. Hilgetag
See the [LICENSE](https://github.com/sofia1423/cortical-gradients-anchoring/blob/main/LICENSE) file for the full license text.

## 📧 Contact
Ao Ma: 202531250002@mail.bnu.edu.cn

Zhengjue Ren: zhengjue.ren@gmail.com

Yuqing Li:  202321250018@mail.bnu.edu.cn

Yuhan Chen: yhchen@bnu.edu.cn
