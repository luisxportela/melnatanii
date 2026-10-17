# melnatanii

Supplementary data and analysis code for:

> Portela LHX, Brunoro I, Nepomuceno A, De Toni KLG, Gonçalves-Oliveira RC, Bovini MG.
> *Overlooked in the herbarium, lost in the field: the discovery of* Mitracarpus elnatanii *sp. nov. (Rubiaceae, Spermacoceae) highlights the vulnerability of herbaceous endemics in the Brazilian restinga.*
> Submitted to **Anais da Academia Brasileira de Ciências**.

*Mitracarpus elnatanii* is described from the restinga of Setiba, Espírito Santo, Brazil, on
the basis of a single 1987 gathering long misidentified as *M. lhotzkyanus*. Beyond the
morphological description, the paper places the new species in the climatic space of four
related species — *M. eichleri*, *M. frigidus*, *M. lhotzkyanus* and *M. salzmannianus* — as an
independent line of evidence for the species hypothesis.

This repository holds the occurrence data, the outputs of those analyses, and the script that
produces them.

---

## Contents

| File | What it is |
|---|---|
| `R_script_climatic_niche.R` | The complete climatic-niche analysis, from the raw rasters and occurrence records to Figure 5 and Supplementary Materials 3–5. |
| `Suppl_material_1_occurrence_records.xlsx` | All 383 compiled occurrence records, with a flag marking the 318 used in the analyses and, for each of the 65 excluded, the reason (no climate data in the cell, coordinates with fewer than two decimal places, or a duplicate within a 30-arc-second cell). Input to the script. |
| `Suppl_material_2_VIF_selection.xlsx` | Variance inflation factors for the 19 bioclimatic variables plus elevation, and the 12 predictors retained. |
| `Suppl_material_3_PCA_all_axes.xlsx` | Loadings and variance explained for all 12 PCA-env axes. |
| `Suppl_material_4_LDA_MANOVA.xlsx` | MANOVA (Pillai's trace) and Linear Discriminant Analysis results. |
| `Suppl_material_5_LDA_confusion_matrix.xlsx` | Leave-one-out cross-validation confusion matrix, with accuracy and Cohen's κ. |

Supplementary Materials 3–5 are written by the script. They are included here so the published
tables can be checked without re-running anything.

---

## The climate layers are not in this repository

The script needs 20 raster layers — the 19 bioclimatic variables plus elevation — and those
are **not deposited here**. They are third-party data with their own terms, and they are large.

Download them from **WorldClim v. 2.1 at 30-arc-second resolution** (Fick & Hijmans 2017,
https://worldclim.org/data/worldclim21.html). Put all 20 `.tif` files in one folder. The script
identifies them by the strings `bio` and `elev` in the filename, resamples them to a common
grid, and crops them to the extent of the occurrence records plus a 1° buffer.

Nothing else is missing: with those rasters and `Suppl_material_1_occurrence_records.xlsx`,
the script reproduces every number and figure in the paper.

---

## Running the analysis

```r
install.packages(c("terra", "readxl", "dplyr", "ggplot2", "ggrepel",
                   "patchwork", "ecospat", "usdm", "MASS", "openxlsx"))

source("R_script_climatic_niche.R")
```

The script installs any missing package itself, then opens **interactive file dialogs**: first
for the folder of climate rasters, then for the occurrence spreadsheet. It is meant to be run
in an interactive session, not in batch.

Tested on **R 4.6.0**.

### What it does, in order

1. Reads and checks the rasters, renames them, crops them to the study area
2. Reads the occurrence records and derives species names
3. **VIF selection** — `vifstep`, threshold 80, on Spearman correlations; 12 of 20 predictors retained
4. Extracts climate values and removes records with missing data, imprecise coordinates, or duplicate species–cell combinations
5. **PCA-env** calibrated on 10,000 background points drawn from the study area, with the records projected onto the first two axes
6. **MANOVA** (Pillai's trace) and **LDA** with leave-one-out cross-validation, reporting accuracy, Cohen's κ and the LD1 coefficients
7. **Niche overlap** — Schoener's *D* and niche equivalency tests via `ecospat`, with Holm-adjusted *P* values
8. Tests whether the *M. elnatanii* record falls inside each species' 95% concentration ellipse, and lists its 10 nearest records in climatic space
9. Builds Figure 5 — the PCA biplot with ellipses, marginal densities and variable-contribution bars

### Outputs

Everything is written to a `resultados_nicho` folder:

```
VIF_results.csv
records_used.csv
Table2_D_P.csv
Table2.csv
closest_records_elnatanii.csv
Figure5.png
Figure5.pdf
Suppl_material_3_PCA_all_axes.xlsx
Suppl_material_4_LDA_MANOVA.xlsx
Suppl_material_5_LDA_confusion_matrix.xlsx
```

The main results are also printed to the console.

### One caveat on the single record

*M. elnatanii* is known from one gathering, so it has no occurrence density and cannot enter
the overlap tests. It is **projected** onto the climatic space defined by the four related
species and tested for membership in their ellipses. The overlap and equivalency statistics
concern the four related species only.

---

## License

- **Code** (`R_script_climatic_niche.R`) — [MIT](LICENSE-MIT)
- **Data** (`.xlsx` files) — [CC BY 4.0](LICENSE-CC-BY-4.0)

Occurrence records were compiled from Herbário Virtual Reflora and speciesLink; the underlying
specimen data belong to the holding herbaria and are subject to their own terms.

---

## Contact

Luís Henrique Ximenes Portela — Escola Nacional de Botânica Tropical, Instituto de Pesquisas
Jardim Botânico do Rio de Janeiro (JBRJ).

Please open an issue for anything that does not reproduce.
