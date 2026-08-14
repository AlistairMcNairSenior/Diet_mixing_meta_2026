
# Script written by AMS at the University of Sydney.
# Formats data, by calculating effect sizes for meta-analysis

# Clean up
rm(list=ls())

# Load packages
library(plyr)
library(ggplot2)
library(Cairo)
library(gridExtra)
library(metafor)
library(ape)
library(corpcor)

# Load data
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Phil_Trans_Meta"
setwd(wd)
data<-read.csv("data/full_data_phil_trans.csv", stringsAsFactors=T)
str(data)
head(data)
summary(data)

# Add a group ID
data$group.ID<-paste0("g", seq(1, nrow(data), 1))

# Get all the different experimental units
comparisons<-unique(data$Comparison.ID)

# Now for every experiment create all pairwise comparisons between single and mixed food groups
for(i in 1:length(comparisons)){
  
  # Pull out the ith experiment and separate the mix and single food groups
  dat_i<-data[which(data$Comparison.ID == comparisons[i]),]
  i_mixed<-dat_i[which(dat_i$treat == "mix"),]
  i_single<-dat_i[which(dat_i$treat == "single"),]
  
  # Add in whether the single food is the single food with max fitness - note for DevTime its minimal level
  if(i_single$Trait_type[1] == "DevTime"){
    i_single$Max.Single<-i_single$mean == min(i_single$mean)
  }else{
    i_single$Max.Single<-i_single$mean == max(i_single$mean)
  }
  
  # Now for each instance of a mixed-food group create a cope of the single-food data set and insert that instance of mixed diet data	
  for(j in 1:nrow(i_mixed)){
    j_comb<-i_single
    j_comb$n_mix<-i_mixed$n[j]
    j_comb$sd_mix<-i_mixed$sd[j]
    j_comb$mean_mix<-i_mixed$mean[j]
    j_comb$n_foods_mix<-i_mixed$n_foods[j]
    j_comb$group.ID_mix<-i_mixed$group.ID[j]
    if(j == 1){
      dat_i_lnRR<-j_comb
    }else{
      dat_i_lnRR<-rbind(dat_i_lnRR, j_comb)
    }
  }
  
  # Save the ith experiment - on the first pass save the first instance
  # Else build on to it.
  if(i == 1){
    dat_contrast<-dat_i_lnRR
  }else{
    dat_contrast<-rbind(dat_contrast, dat_i_lnRR)
  }
}

# Add in an effect size level ID
dat_contrast$Effect.ID<-paste0("e", seq(1, nrow(dat_contrast), 1))

# Calculate the Effect Sizes
lnRR_data<-escalc(measure="ROM", n1i=dat_contrast$n, n2i=dat_contrast$n_mix, m1i=dat_contrast$mean, m2i=dat_contrast$mean_mix, sd1i=dat_contrast$sd, sd2i=dat_contrast$sd_mix)
lnCVR_data<-escalc(measure="CVR", n1i=dat_contrast$n, n2i=dat_contrast$n_mix, m1i=dat_contrast$mean, m2i=dat_contrast$mean_mix, sd1i=dat_contrast$sd, sd2i=dat_contrast$sd_mix)
names(lnRR_data)<-c("lnRR", "v_lnRR")
names(lnCVR_data)<-c("lnCVR", "v_lnCVR")
dat_contrast<-cbind(dat_contrast, lnRR_data, lnCVR_data)
head(dat_contrast)

# Use the log of n mixed foods for analysis
dat_contrast$ln_n_foods_mix<-log(dat_contrast$n_foods_mix)

###################################################
# Create the Variance-co-variance matrix for each #
###################################################

# Work out the covariance for effect sizes that feature a given groupID
single_correction<-ddply(dat_contrast, .(group.ID), summarise, n_effects=length(group.ID), n=n[1], sd=sd[1], mean=mean[1])
mix_correction<-ddply(dat_contrast, .(group.ID_mix), summarise, n_effects=length(group.ID_mix), n=n_mix[1], sd=sd_mix[1], mean=mean_mix[1])
single_correction$cov_lnRR<-(single_correction$sd/single_correction$mean)^2 / single_correction$n + (single_correction$sd/single_correction$mean)^4 / (2 * single_correction$n^2)
single_correction$cov_lnCVR<-single_correction$cov_lnRR + 1/2 * (single_correction$n / (single_correction$n - 1)^2)
mix_correction$cov_lnRR<-(mix_correction$sd/mix_correction$mean)^2 / mix_correction$n + (mix_correction$sd/mix_correction$mean)^4 / (2 * mix_correction$n^2)
mix_correction$cov_lnCVR<-mix_correction$cov_lnRR + 1/2 * (mix_correction$n / (mix_correction$n - 1)^2)

# Start with size
vcv_lnRR<-matrix(0, nrow=nrow(dat_contrast), ncol=nrow(dat_contrast))
rownames(vcv_lnRR)<-dat_contrast$Effect.ID
colnames(vcv_lnRR)<-dat_contrast$Effect.ID
vcv_lnCVR<-vcv_lnRR

# Loop to fill the off diags
for(i in 1:nrow(dat_contrast)){
  
  # For single group matches - get the group ID and right covariances
  group_i<-dat_contrast$group.ID[i]
  cov_i_lnRR<-single_correction$cov_lnRR[match(group_i, single_correction$group.ID)]
  cov_i_lnCVR<-single_correction$cov_lnCVR[match(group_i, single_correction$group.ID)]
  
  # Get all instances with that group ID and remove the focal instance
  matches<-which(dat_contrast$group.ID == group_i)
  matches<-matches[-which(matches == i)]
  
  # If there are matches then add the covariances
  if(length(matches) > 0){
    vcv_lnRR[i,matches]<-cov_i_lnRR
    vcv_lnCVR[i,matches]<-cov_i_lnCVR
  }
  
  # Clean up
  rm(matches)
  rm(group_i)
  rm(cov_i_lnRR)
  rm(cov_i_lnCVR)
  
  # Repeat for mix group ID
  group_i<-dat_contrast$group.ID_mix[i]
  cov_i_lnRR<-mix_correction$cov_lnRR[match(group_i, mix_correction$group.ID_mix)]
  cov_i_lnCVR<-mix_correction$cov_lnCVR[match(group_i, mix_correction$group.ID_mix)]
  matches<-which(dat_contrast$group.ID_mix == group_i)
  matches<-matches[-which(matches == i)]
  if(length(matches) > 0){
    vcv_lnRR[i,matches]<-cov_i_lnRR
    vcv_lnCVR[i,matches]<-cov_i_lnCVR
  }
  rm(matches)
  rm(group_i)
  rm(cov_i_lnRR)
  rm(cov_i_lnCVR)
  
}

# Add the sampling errors in the diagonals
diag(vcv_lnRR)<-dat_contrast$v_lnRR
diag(vcv_lnCVR)<-dat_contrast$v_lnCVR

# Subset by trait, check for positive definite
size_lnRR_vcv<-vcv_lnRR[which(dat_contrast$Trait_type == "Size"),which(dat_contrast$Trait_type == "Size")]
is.positive.definite(size_lnRR_vcv)
size_lnCVR_vcv<-vcv_lnCVR[which(dat_contrast$Trait_type == "Size"),which(dat_contrast$Trait_type == "Size")]
is.positive.definite(size_lnCVR_vcv)

dev_lnRR_vcv<-vcv_lnRR[which(dat_contrast$Trait_type == "DevTime"),which(dat_contrast$Trait_type == "DevTime")]
is.positive.definite(dev_lnRR_vcv)
dev_lnCVR_vcv<-vcv_lnCVR[which(dat_contrast$Trait_type == "DevTime"),which(dat_contrast$Trait_type == "DevTime")]
is.positive.definite(size_lnCVR_vcv)

repro_lnRR_vcv<-vcv_lnRR[which(dat_contrast$Trait_type == "Repro"),which(dat_contrast$Trait_type == "Repro")]
is.positive.definite(repro_lnRR_vcv)
repro_lnCVR_vcv<-vcv_lnCVR[which(dat_contrast$Trait_type == "Repro"),which(dat_contrast$Trait_type == "Repro")]
is.positive.definite(repro_lnCVR_vcv)

long_lnRR_vcv<-vcv_lnRR[which(dat_contrast$Trait_type == "Longevity"),which(dat_contrast$Trait_type == "Longevity")]
is.positive.definite(repro_lnRR_vcv)
long_lnCVR_vcv<-vcv_lnCVR[which(dat_contrast$Trait_type == "Longevity"),which(dat_contrast$Trait_type == "Longevity")]
is.positive.definite(repro_lnCVR_vcv)

# Does forcing PD retain sampling errors? - check correlations
size_lnRR_vcv_PD<-make.positive.definite(size_lnRR_vcv)
size_lnCVR_vcv_PD<-make.positive.definite(size_lnCVR_vcv)
cor(diag(size_lnRR_vcv_PD), diag(size_lnRR_vcv))
cor(diag(size_lnCVR_vcv_PD), diag(size_lnCVR_vcv))

dev_lnRR_vcv_PD<-make.positive.definite(dev_lnRR_vcv)
dev_lnCVR_vcv_PD<-make.positive.definite(dev_lnCVR_vcv)
cor(diag(dev_lnRR_vcv_PD), diag(dev_lnRR_vcv))
cor(diag(dev_lnCVR_vcv_PD), diag(dev_lnCVR_vcv))

repro_lnRR_vcv_PD<-make.positive.definite(repro_lnRR_vcv)
repro_lnCVR_vcv_PD<-make.positive.definite(repro_lnCVR_vcv)
cor(diag(repro_lnRR_vcv_PD), diag(repro_lnRR_vcv))
cor(diag(repro_lnCVR_vcv_PD), diag(repro_lnCVR_vcv))

long_lnRR_vcv_PD<-make.positive.definite(long_lnRR_vcv)
long_lnCVR_vcv_PD<-make.positive.definite(long_lnCVR_vcv)
cor(diag(long_lnRR_vcv_PD), diag(long_lnRR_vcv))
cor(diag(long_lnCVR_vcv_PD), diag(long_lnCVR_vcv))

# In all cases the correlation is > 0.99, suggesting we will not over/underweight effect sizes in downstream analyses, so we will use the PD

###################################################
############ SUBSET AND SAVE THE DATASETS #########
###################################################

# Subset the effect sizes and meta-data and save in a list
effect_list<-list(dat_contrast[which(dat_contrast$Trait_type == "Size"),], 
              dat_contrast[which(dat_contrast$Trait_type == "DevTime"),],
              dat_contrast[which(dat_contrast$Trait_type == "Repro"),],
              dat_contrast[which(dat_contrast$Trait_type == "Longevity"),])
names(effect_list)<-c("Size", "DevTime", "Repro", "Longevity")

# Save as an RDS file
save(effect_list, file="data/EffectSizes.Rds")

# save the PD VCVs
lnRR_vcv_list<-list(size_lnRR_vcv_PD,
                    dev_lnRR_vcv_PD,
                    repro_lnRR_vcv_PD,
                    long_lnRR_vcv_PD)
names(lnRR_vcv_list)<-c("Size", "DevTime", "Repro", "Longevity")

lnCVR_vcv_list<-list(size_lnCVR_vcv_PD,
                    dev_lnCVR_vcv_PD,
                    repro_lnCVR_vcv_PD,
                    long_lnCVR_vcv_PD)

names(lnCVR_vcv_list)<-c("Size", "DevTime", "Repro", "Longevity")

# Save as RDS files
save(lnRR_vcv_list, file="data/lnRR_vcv.Rds")
save(lnCVR_vcv_list, file="data/lnCVR_vcv.Rds")

############################################
########## DESCRIPTIVE STUFF ###############
############################################

head(data)

# Number of articles and experiments
dim(data)
length(unique(data$Consumer.Sp))

# Observations and experiments by trait
length(which(data$Trait_type == "Size"))
length(which(data$Trait_type == "Repro"))
length(which(data$Trait_type == "DevTime"))
length(which(data$Trait_type == "Longevity"))

# Number of experiments per trait type
ddply(data, .(Trait_type), summarise, n_exp=length(unique(Comparison.ID)), n_obs=length(Comparison.ID))

# Number of effect sizes in each trait
dim(size_lnRR_vcv_PD)
dim(repro_lnRR_vcv_PD)
dim(dev_lnRR_vcv_PD)
dim(long_lnRR_vcv_PD)

##################################################### 
############ MEAN-VARIANCE CORRELATION ##############
#####################################################

# Group centre the means and SDs
data$lnX<-log(data$mean)
data$lnSD<-log(data$sd)
data$z_lnX<-NA
data$z_lnSD<-NA
comparisons<-unique(data$Comparison.ID)
for(i in 1:length(comparisons)){
	tag_i<-which(data$Comparison.ID == comparisons[i])
	data$z_lnX[tag_i]<-data$lnX[tag_i] - mean(data$lnX[tag_i])
	data$z_lnSD[tag_i]<-data$lnSD[tag_i] - mean(data$lnSD[tag_i])
}


# Check centered size data
a<-ggplot(data[which(data$Trait_type == "Size"),], aes(x = z_lnX, y = z_lnSD)) +
  geom_point(size=0.3) + theme_bw() +
  xlab("Centered log Mean") + ylab("Centered log SD") +
  ggtitle("A. Body Size")

# Check correlation
cor.test(data[which(data$Trait_type == "Size"),]$z_lnX, data[which(data$Trait_type == "Size"),]$z_lnSD, use='pairwise.complete.obs')

# Check centered reproductive data
b<-ggplot(data[which(data$Trait_type == "Repro"),], aes(x = z_lnX, y = z_lnSD)) +
		geom_point(size=0.3) + theme_bw() +
		xlab("Centered log Mean") + ylab("Centered log SD") +
		ggtitle("B. Reproductive Function")

# Check correlation
cor.test(data[which(data$Trait_type == "Repro"),]$z_lnX, data[which(data$Trait_type == "Repro"),]$z_lnSD, use='pairwise.complete.obs')

# Check centered longevity data
c<-ggplot(data[which(data$Trait_type == "DevTime"),], aes(x = z_lnX, y = z_lnSD)) +
		geom_point(size=0.3) + theme_bw() +
		xlab("Centered log Mean") + ylab("Centered log SD") +
		ggtitle("C. Development Time")

# Check correlation
cor.test(data[which(data$Trait_type == "DevTime"),]$z_lnX, data[which(data$Trait_type == "DevTime"),]$z_lnSD, use='pairwise.complete.obs')

# Check centered longevity data
d<-ggplot(data[which(data$Trait_type == "Longevity"),], aes(x = z_lnX, y = z_lnSD)) +
  geom_point(size=0.3) + theme_bw() +
  xlab("Centered log Mean") + ylab("Centered log SD") +
  ggtitle("D. Lifespan")

# Check correlation
cor.test(data[which(data$Trait_type == "Longevity"),]$z_lnX, data[which(data$Trait_type == "Longevity"),]$z_lnSD, use='pairwise.complete.obs')

# Figure 1
CairoPDF("plots/figure_S2.pdf", height=4*2, width=4*2)

grid.arrange(a,b,c,d, nrow=2, ncol=2)

dev.off()
