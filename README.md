# immune-cell-eQTL-T2D-complicaitons

This repository contains code for reproducing the transcriptome-wide Mendelian randomization examined in ' Integrating immune-cell transcriptomic data with Mendelian randomization reveals novel causal genes for type 2 diabetes and its complications'. In this study, MR and colocalization analyses revealed the expression of 425 and 123 unique genes associated with T2DM and its complications respectively. We further quantified the impacts of cell-type-related pleiotropy, demonstrating the percentage of pleiotropic eQTLs increased from 14.4% (classic pleiotropy) to 56.6%. Applying the six MVMR methods could substantially attenuate the cell-type-related pleiotropy for the top findings. Our study supports a key role for immune mechanisms in diabetic complications and highlights promising therapeutic targets by distinguishing and minimizing the influence of cell-type-related pleiotropy.

To start using the code, you need to install TwoSampleMR package:

install.packages("remotes")

remotes::install_github("MRCIEU/TwoSampleMR")

devtools::install_github("mrcieu/ieugwasr")

The following code were used in this study:

1. Immune_complication_UKB.R:To explore associations between immunity and tendency to different diabetic complications from the UKB.

2. MR. R: To estimated the putative causal effects of 18,611 gene expression traits on T2DM and its complications.

3. Steiger. R: To test the directionality of the eQTL-outcome associations.

4. Colocalization. R: To obtain more reliable evidence for immune-related gene associations with T2DM and its complications (FDR < 0.05).

5. LD check. R: To assess approximate colocalization evidence for causal genes with robust MR evidence (FDR < 0.05).

6-11. WeakIV: To address cell-type-related pleiotropy with limited cis-eQTLs.
