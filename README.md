This repository contains the raw data and code to perform the meta-analyses in Senior et al. Effects of diet-mixing on life-history traits: meta-analytic tests of predictions from the geometric framework for nutrition.

The code, data and output are organised into four directories.
Running through the scripts in order 1.tree.R with the raw dataset full_data_phil_trans.csv will generate all other output.

The contents of the directories are as follows:

code
1. tree.R - an R script to generate the phylogenetic covariance matrices for the analyses.
2. data_preperation.R - an R script to calculate pairwise effect sizes and the covariance matrices for sampling errors, as well as other data processing steps.
3. models.R - an R script to run the final models and generate tables and figures. 

data
full_data_phil_trans.csv - a csv file with the raw extracted data.
DevTime_tree.tre - a tre file containing the tree for species in the developmental time analyses (created by the script 1. tree.R)
Longevity_tree.tre - a tre file containing the tree for species in the lifespan analyses (created by the script 1. tree.R)
Repro_tree.tre - a tre file containing the tree for species in the reproductive function analyses (created by the script 1. tree.R)
Size_tree.tre - a tre file containing the tree for species in the body size analyses (created by the script 1. tree.R)
PhyloMatrix.Rds - an rds file containing a list of phylogenetic covariance matrices in the analyses (created by the script 1. tree.R)
EffectSizes.Rds - an rds file containing all the processed and estimated effect sizes (created by the script 2. data_preperation.R)
lnCVR_vcv.Rds - an rds file containing a list of the estimated variance-covariance matrices for sampling errors on lnCVR (created by the script 2. data_preperation.R)
lnRR_vcv.Rds - an rds file containing a list of the estimated variance-covariance matrices for sampling errors on lnRR (created by the script 2. data_preperation.R)

plots
figure_2.pdf - a pdf file for figure_2 (created by the script 3. models.R).
figure_3.pdf - a pdf file for figure_3 (created by the script 3. models.R).
figure_4.pdf - a pdf file for figure_4 (created by the script 3. models.R).
figure_S2.pdf - a pdf file for figure_S2 (created by the script 2. data_preperation.R).
figure_S3.pdf - a pdf file for figure_S3 (created by the script 1. tree.R).

tables
MA_table.csv - a csv file containing meta-analytic estimates (created by the script 3. models.R).
MR_table.csv - a csv file containing meta-regression estimates (created by the script 3. models.R).
pb_table.csv - a csv file containing estimates from analyses of publication bias (created by the script 3. models.R).
taxonomic_data.csv - a csv file giving the taxonomy of each species (created by the script 1. tree.R)
