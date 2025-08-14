# Interpretable Machine Learning on Soybean Multi-Omics Data Reveals Drought-Driven Shifts of Plant–Microbe Interactions

This repository contains all scripts and data used for the analysis in **“Interpretable Machine Learning on Soybean Multi-Omics Data Reveals Drought-Driven Shifts of Plant–Microbe Interactions.”**



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

for snp full data, also refer to Kanegae et al. (2021). "Whole-genome sequence diversity and associa-
tion analysis of 198 soybean accessions in mini-core collections."


