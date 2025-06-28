# immune-cell-eQTL-T2D-complicaitons

This repository contains code for reproducing the transcriptome-wide Mendelian randomization examined in ' Integrating immune-cell transcriptomic data with Mendelian randomization reveals novel causal genes for type 2 diabetes and its complications'. In this study, MR and colocalization analyses revealed the expression of 425 and 123 unique genes associated with T2DM and its complications respectively. We further quantified the impacts of cell-type-related pleiotropy, demonstrating the percentage of pleiotropic eQTLs increased from 46% (classic pleiotropy) to 91.4%, which implied widespread pleiotropy across almost all eQTLs. Applying the MV-IVW-PCA method could substantially attenuate the cell-type-related pleiotropy for the top findings. Our study supports a key role for immune mechanisms in diabetic complications and highlights promising therapeutic targets by distinguishing and minimizing the influence of cell-type-related pleiotropy.

To start using the code, you need to install TwoSampleMR package:

install.packages("remotes")

remotes::install_github("MRCIEU/TwoSampleMR")

devtools::install_github("mrcieu/ieugwasr")

The following code were used in this study:

1. MR. R: To identify the causal effects for T2DM and its complications.
2. Steiger. R: To test the directionality of the eQTL-outcome associations.
3. Colocalization. R: To obtain more reliable evidence for immune-related gene associations with T2DM and its complications (FDR < 0.05).
4. LD check. R: To assess approximate colocalization evidence for causal genes with robust MR evidence (FDR < 0.05).
5-7. WeakIV (MV_IVW_PCA, MVMR, TWMR): To address cell-type-related pleiotropy with limited cis-eQTLs.
