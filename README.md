# Interpretable Machine Learning on Soybean Multi-Omics Data Reveals Drought-Driven Shifts of Plant–Microbe Interactions

This repository contains all scripts and data used for the analysis in **“Interpretable Machine Learning on Soybean Multi-Omics Data Reveals Drought-Driven Shifts of Plant–Microbe Interactions.”**


📄 Preprint available at: [https://doi.org/10.1101/2025.08.13.670005](https://doi.org/10.1101/2025.08.13.670005)


---

<img src="abstract.png" width="1000">

---

## 📂 script_ShapPlantMicro
All scripts 

| File | Description |
|------|--------------|
| `c0_RF-markers.Rmd` | Random Forest marker evaluation |
| `c1_RF_prediction.Rmd` | Random Forest prediction |
| `c2_BLUP_prediction.Rmd` | BLUP (Best Linear Unbiased Prediction) |
| `c3_GWAS.Rmd` | GWAS (Genome-Wide Association Study) analysis |
| `c4_SNP_Mapping.Rmd` | SNP selection by RF and GWAS |
| `c5_SHAP_single.Rmd` | SHAP analysis for single feature selection |
| `c6_SHAP_interaction.Rmd` | SHAP for interaction feature selection |
| `c7_SHAP_difference.Rmd` | SHAP difference matrix visualization |
| `c8_hyperparameter_xgboost.Rmd` | Hyperparameter Tuning for XGBoost |
| `c9_hyperparameter_ranger.Rmd` | Hyperparameter Tuning for ranger |

---

## 📂 data 
Multi-omics data (Genome, Microbiome, Metabolome, Phenotype)

| File | Description |
|------|--------------|
| `SoyData_Drought2.RDS` | Multi omics data - Drought |
| `SoyData_Control2.RDS` | Multi omics data - Control |

---

## 📂 data/genome
Genome marker data (with different LD thresholds)

| File | Description |
|------|--------------|
| `genoMarker_LD0.001_SNP3078.RDS`  | SNP markers filtered at LD = 0.001 |
| `genoMarker_LD0.01_SNP10143.RDS`  | SNP markers filtered at LD = 0.01  |
| `genoMarker_LD0.1_SNP16419.RDS`   | SNP markers filtered at LD = 0.1   |
| `genoMarker_LD0.3_SNP34632.RDS`   | SNP markers filtered at LD = 0.3   |
---


## 📂 data/omics_csv
Genome marker data (with different LD thresholds)

| File | Description |
|------|--------------|
| `genoMarker_LD0.001_SNP3078.csv`  | SNP markers filtered at LD = 0.001 |
| `genoMarker_LD0.01_SNP10143.csv`  | SNP markers filtered at LD = 0.01  |
| `genoMarker_LD0.1_SNP16419.csv`   | SNP markers filtered at LD = 0.1   |
| `genoMarker_LD0.3_SNP34632.csv`   | SNP markers filtered at LD = 0.3   |
| `metabolome_namelist.csv` | Metabolome name list (265 features) |
| `microbiome_namelist.csv` | Microbiome name list (16457 features)|
| `met_Control.csv`   | Metabolome (Control) |
| `met_Drought.csv`   | Metabolome (Drought) |
| `micro_Control.csv` | Microbiome (Control) |
| `micro_Drought.csv` | Microbiome (Drought) |
| `pheno_Control.csv` | Phenotypic traits (Control) |
| `pheno_Drought.csv` | Phenotypic traits (Drought) |
---


for snp full data, also refer to Kanegae et al. (2021). "Whole-genome sequence diversity and associa-
tion analysis of 198 soybean accessions in mini-core collections."


