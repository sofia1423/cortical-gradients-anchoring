Great — a **group-level directed, weighted 91×91 connectome on the published M132 macaque atlas** is exactly the kind of _derived resource_ that is appropriate to share.

Here’s a **copy-paste package** you can use.

## What to share (minimum set)

1. `M132_group_connectome_directed_weighted.csv` (or `.npy/.mat`)
    
2. `M132_labels.csv` (columns like: `index, region_name, hemisphere, atlas_label_id`)
    
3. `README_connectome.txt` (short, see below)
    

## README (you can paste and edit)

**README_connectome.txt**

- **Atlas / parcellation:** M132 macaque atlas (91 cortical regions)
    
- **File:** `M132_group_connectome_directed_weighted.csv`
    
- **Matrix definition:** a **91 × 91 directed, weighted adjacency matrix**. Entry (W_{ij}) denotes the connection weight **from source region (i) to target region (j)**.
    
- **Weights:** (add one phrase) “Weights represent [streamline count / connection probability / normalized tract strength / log-transformed weight], group-averaged across subjects.”
    
- **Diagonal:** “Diagonal entries are [0 / NaN / retained]”.
    
- **Preprocessing:** “Matrix was [symmetry not enforced / normalized by … / thresholded at …] before gradient estimation.”
    
- **Region order:** Row/column order corresponds to `M132_labels.csv` (index 1–91).
    
- **Intended use:** “These matrices were used to construct the connectome-derived gradients and to parameterize the large-scale dynamical model.”
    

(That’s enough for reviewers; you can expand later for the public release.)

## What to select on the submission “Standardized datasets” question

This connectome file is a **derived dataset**, not typically a Cell Press “standardized dataset type.” So you would usually answer:

- **Standardized datasets:** **No**  
    and then you share the connectome under **Data and code availability** (and/or as a reviewer-access file link).
    

## Text for the submission portal (the 2000-character box if needed)

If you need to provide reviewer access instructions somewhere (often in that box or in “Data and code availability”), you can paste:

> **Derived connectome data access (for peer review):** We provide a group-level directed, weighted white-matter connectome matrix (91×91) in M132 atlas space, along with region label files. These derived data were used to construct the gradients and to parameterize the computational model. Anonymous reviewer access is available via [PRIVATE REVIEW LINK / TOKEN], and the files will be deposited in a public repository (e.g., Zenodo/OSF) and made publicly available as of the date of publication.

## One detail to be careful about

Because it’s **directed**, reviewers will ask what “direction” means for white-matter tractography (which is often undirected unless you used tracer-based directionality or imposed direction via another rule). So in the README and methods, add **one clear line** explaining _how directionality was defined_.

If you tell me **how you defined direction (tractography? tracer? asymmetry rule based on neuron density? etc.)** and what the **weight** is (streamlines, log, normalization), I’ll finalize the exact README text + a very polished 2–3 line “Data availability” statement that matches Neuron STAR Methods tone.