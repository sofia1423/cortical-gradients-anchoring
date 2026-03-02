# Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the code and data for the paper:

**"Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex"**  
Authors: Ao Ma, Zhengjue Ren, Yuqing Li, Shui Yu, Changsong Zhou, Claus C. Hilgetag, Yuhan Chen  
Journal/Conference: 

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

| File                                                          | Dimensions | Description |
| ------------------------------------------------------------- | ---------- | ----------- |
| `BrainMesh_Monkey_F99.nv`                                     |            |             |
| `distance_matrix.mat`                                         |            |             |
| `Empirical_and_Predicted_Rank_Difference.mat`                 |            |             |
| `Empirical_Axis3_Afferent_WR.mat`                             |            |             |
| `Empirical_Efferent_WR.mat`                                   |            |             |
| `Empirical_Spectral_Gradient_during_Anesthesia_and_Awake.mat` |            |             |
| `fmri_Gradient1_and2.mat`                                     |            |             |
| `Function_Domain.txt`                                         |            |             |
| `label_467.mat`                                               |            |             |
| `name91.mat`                                                  |            |             |
| `Power_for_all_frequency_bands.mat`                           |            |             |
| `Power_for_all_frequency_bands.txt`                           |            |             |
| `Predicted_Axis2_T1_T2.mat`                                   |            |             |
| `Predicted_Axis3_Afferent_WR.mat`                             |            |             |
| `Predicted_Efferent_WR.mat`                                   |            |             |
| `Predicted_Spectral_Gradient_during_Anesthesia.txt`           |            |             |
| `Predicted_Spectral_Gradient_during_Anesthesia_and_Awake.mat` |            |             |
| `Predicted_Spectral_Gradient_during_Awake.txt`                |            |             |
| `Receptor.mat`                                                |            |             |
| `Structural_Connectome.mat`                                   |            |             |



## 💻 Code

All analysis code is organized in the `code/` directory.

```
code/
├── cal_rank_difference/
│   └── cal_rank_difference.m        
│
├── Structure_dimensionality_reduction/
│   ├── code_for_UMAP.ipynb
│   └── dimensionality_reduction.m
│        
├── SVD
│	└── SVD_code.m
│
├── generative_model
│   ├──readme.txt
│   │
│	├── 1_first_step/
│   │   ├── pre_weight_1_40_40.m
│   │   ├── pre_weight_2_40_40.m
│   │   └── trans_matrix_40area.m
│   │
│   ├── 2_second_step/
│   │   ├── Cal_posterior_probability.m
│   │   ├── data_array_40area.m
│   │   ├── data_array2_40area.m
│   │   └── trans_matrix_40area.m   
│   │  
│   ├── 3_validation
│   │   ├── cal_roc.m
│   │   ├── matrix_sort.m
│   │   └── validation_code.m

```

### Key Scripts and Their Outputs
#### cal_rank_difference
1. cal_rank_difference.m        

#### Structure_dimensionality_reduction
1. code_for_UMAP.ipynb

2. dimensionality_reduction.m

#### SVD
1. SVD_code.m

#### generative_model
1. readme.txt
##### 1_first_step
1. pre_weight_1_40_40.m

2. pre_weight_2_40_40.m

3. trans_matrix_40area.m

##### 2_second_step
1. Cal_posterior_probability.m

2. data_array_40area.m

3. data_array2_40area.m

4. trans_matrix_40area.m   

##### 3_validation
1. cal_roc.m

2. matrix_sort.m

3. validation_code.m


## 🚀 Quick Start

### Installation

```bash

# Install dependencies
pip install -r requirements.txt
```


## 📝 Notes

- **Data size**: Most files are small CSV/text files (<1MB each)
- **Spatial coordinates**: All region-level data uses M132 atlas parcellation (91 regions)

## 📚 Citation

If you use this code or data in your research, please cite:

```
Ma, A., Ren, Z., Li, Y., Yu, S., Zhou, C., Hilgetag, C.C., & Chen, Y. (2026). 
Anchoring of anesthetized and awake brain states to natural structural axes of the cerebral cortex. 
[Journal Name], [Volume], [Pages].
```

## 📄 License

This project is licensed under the MIT License.  
Copyright (c) 2026 Yuhan Chen, Ao Ma, Zhengjue Ren, Yuqing Li, Shui Yu, Changsong Zhou, Claus C. Hilgetag
See the [LICENSE](https://github.com/sofia1423/cortical-gradients-anchoring/blob/main/LICENSE) file for the full license text.

## 📧 Contact
Ao Ma: 202531250002@mail.bnu.edu.cn
Zhengjue Ren: zhengjue.ren@gmail.com
Yuqing Li: 202321250018@mail.bnu.edu.cn
Yuhan Chen: yhchen@bnu.edu.cn
