# sci-Comm-Plex

Analysis workflows and notebooks for the paper:
"Mapping kinase-dependent tumor immune adaptation with multiplexed single-cell CRISPR screens."

## Preprint

The preprint describing this work is available on bioRxiv:
https://www.biorxiv.org/content/10.64898/2026.01.08.698516v1



## Preprocessing pipeline

Raw data processing (hash demultiplexing and CROP-seq alignment) is based on the following pipelines:

- **Hash + CROP-seq processing**: adapted from [mcfaline-figueroa-lab/sci-Plex-EGFRi](https://github.com/mcfaline-figueroa-lab/sci-Plex-EGFRi/tree/main/process_from_raw)
- **CROP-seq processing**: adapted from [cole-trapnell-lab/sci-Plex-GxE](https://github.com/cole-trapnell-lab/sci-Plex-GxE/tree/main/process_from_raw)

Our version of the preprocessing pipeline will be added to this repository shortly.

## Repository structure

- `Kinome_sceen_T_cells_coculture/` – Kinome screen T cell co-culture analyses (R Markdown, Python notebooks, Decipher and MrVI sub-analyses).
- `PDN_T_cells_coculture/` – PDN T cell co-culture analyses and quality control.
- `Published_data_analysis/` – Analyses of published datasets, including spatial data.
- `T_cells_alone_analysis/` – T cell–only single-cell and clonotype analyses.
- `Tcell_killing_image_analysis/` – Image-based analyses of T cell killing and cell confluency.

## Data availability

Data will be deposited to GEO and linked here when available.