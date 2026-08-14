
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
library(ggbeeswarm)

# Load data
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Phil_Trans_Meta"
setwd(wd)
load("data/EffectSizes.Rds")
# Check it out
names(effect_list)
head(effect_list[[1]])

# Load the VCVs
load("data/lnRR_vcv.Rds")
names(lnRR_vcv_list)
load("data/lnCVR_vcv.Rds")
names(lnCVR_vcv_list)

# Load the phylogenetic covariance matrices
phyloM<-readRDS("data/PhyloMatrix.Rds")
names(phyloM)

#########################################
############# META-ANALYSES #############
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models
repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_mean)
repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_var)

# Fit the longevity models
long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_mean)
long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_var)

# Create the table of MA results
MA_table<-data.frame(Trait=c("Size", "", "Developmental Time", "", "Reproductive Function", "", "Lifespan", ""),
                     Effect=rep(c("lnRR", "lnCVR"), 4),
                     Estimate=c(size_MA_mean$b, size_MA_var$b, dev_MA_mean$b, dev_MA_var$b, repro_MA_mean$b, repro_MA_var$b, long_MA_mean$b, long_MA_var$b),
                     CI_lower=c(size_MA_mean$ci.lb, size_MA_var$ci.lb, dev_MA_mean$ci.lb, dev_MA_var$ci.lb, repro_MA_mean$ci.lb, repro_MA_var$ci.lb, long_MA_mean$ci.lb, long_MA_var$ci.lb),
                     CI_upper=c(size_MA_mean$ci.ub, size_MA_var$ci.ub, dev_MA_mean$ci.ub, dev_MA_var$ci.ub, repro_MA_mean$ci.ub, repro_MA_var$ci.ub, long_MA_mean$ci.ub, long_MA_var$ci.ub),
                     M2_total=NA,
                     M2_exp=NA,
                     M2_phylo=NA,
                     M2_within=NA)

# Calculate the M2 values and add in
Tau2_total_mean<-c(sum(size_MA_mean$sigma2), sum(size_MA_var$sigma2), sum(dev_MA_mean$sigma2), sum(dev_MA_var$sigma2), sum(repro_MA_mean$sigma2), sum(repro_MA_var$sigma2), sum(long_MA_mean$sigma2), sum(long_MA_var$sigma2))
MA_table$M2_total<-Tau2_total_mean / (Tau2_total_mean + MA_table$Estimate^2)
MA_table$M2_exp<-c(size_MA_mean$sigma2[1], size_MA_var$sigma2[1], dev_MA_mean$sigma2[1], dev_MA_var$sigma2[1], repro_MA_mean$sigma2[1], repro_MA_var$sigma2[1], long_MA_mean$sigma2[1], long_MA_var$sigma2[1]) / (Tau2_total_mean + MA_table$Estimate^2)
MA_table$M2_phylo<-c(size_MA_mean$sigma2[2], size_MA_var$sigma2[2], dev_MA_mean$sigma2[2], dev_MA_var$sigma2[2], repro_MA_mean$sigma2[2], repro_MA_var$sigma2[2], long_MA_mean$sigma2[2], long_MA_var$sigma2[2]) / (Tau2_total_mean + MA_table$Estimate^2)
MA_table$M2_within<-c(size_MA_mean$sigma2[3], size_MA_var$sigma2[3], dev_MA_mean$sigma2[3], dev_MA_var$sigma2[3], repro_MA_mean$sigma2[3], repro_MA_var$sigma2[3], long_MA_mean$sigma2[3], long_MA_var$sigma2[3]) / (Tau2_total_mean + MA_table$Estimate^2)

# Save it
write.table(MA_table, file="tables/MA_table.csv", sep=",", col.names=names(MA_table), row.names=F)

# Checking the older studies, for effect magnitude for comparison with Senior et al. 2015
rma(yi=lnRR, vi=v_lnRR, data=effect_list$Size[which(effect_list$Size$Year < 2010),])
rma(yi=lnCVR, vi=v_lnCVR, data=effect_list$Size[which(effect_list$Size$Year < 2010),])
rma(yi=lnRR, vi=v_lnRR, data=effect_list$Repro[which(effect_list$Repro$Year < 2010),])
rma(yi=lnCVR, vi=v_lnCVR, data=effect_list$Repro[which(effect_list$Repro$Year < 2010),])
rma(yi=lnRR, vi=v_lnRR, data=effect_list$Longevity[which(effect_list$Longevity$Year < 2010),])
rma(yi=lnCVR, vi=v_lnCVR, data=effect_list$Longevity[which(effect_list$Longevity$Year < 2010),])

rma(yi=lnRR, vi=v_lnRR, mods=~Year, data=effect_list$Longevity)
rma(yi=lnCVR, vi=v_lnCVR, mods=~Year, data=effect_list$Longevity)

#########################################
########### ORCHARD PLOTS ###############
#########################################

# Calculate Prediction Intervals
PI_lower<-MA_table$Estimate - 1.96 * sqrt(Tau2_total_mean)
PI_upper<-MA_table$Estimate + 1.96 * sqrt(Tau2_total_mean)

# Colours for points
cols<-c("gold", "cornflowerblue")
size_pes<-0.025
axis.text_size<-12

# Size Plot
rows<-c(1,2) # Rows in the MA table containing the mean and var estimates
plot_plot<-data.frame(yi=c(effect_list$Size$lnRR, effect_list$Size$lnCVR),
                      Precision=1/sqrt(c(effect_list$Size$v_lnRR, effect_list$Size$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Size)), rep("lnCVR", nrow(effect_list$Size))))

A<-ggplot(plot_plot, aes(x=yi, y=effect, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="A. Body Size") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=PI_lower[rows[1]], xend=PI_upper[rows[1]], y=2, yend=2, size=0.5) +
  annotate("segment", x=PI_lower[rows[2]], xend=PI_upper[rows[2]], y=1, yend=1, size=0.5) +
  annotate("segment", x=MA_table$Estimate[rows[1]], xend=MA_table$Estimate[rows[1]], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$Estimate[rows[2]], xend=MA_table$Estimate[rows[2]], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=MA_table$CI_lower[rows[1]], xend=MA_table$CI_upper[rows[1]], y=2, yend=2, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$CI_lower[rows[2]], xend=MA_table$CI_upper[rows[2]], y=1, yend=1, size=1.4, col="firebrick") 

# Devtime Plot
rows<-c(3,4) # Rows in the MA table containing the mean and var estimates
plot_plot<-data.frame(yi=c(effect_list$DevTime$lnRR, effect_list$DevTime$lnCVR),
                      Precision=1/sqrt(c(effect_list$DevTime$v_lnRR, effect_list$DevTime$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$DevTime)), rep("lnCVR", nrow(effect_list$DevTime))))

B<-ggplot(plot_plot, aes(x=yi, y=effect, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="B. Development Time") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=PI_lower[rows[1]], xend=PI_upper[rows[1]], y=2, yend=2, size=0.5) +
  annotate("segment", x=PI_lower[rows[2]], xend=PI_upper[rows[2]], y=1, yend=1, size=0.5) +
  annotate("segment", x=MA_table$Estimate[rows[1]], xend=MA_table$Estimate[rows[1]], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$Estimate[rows[2]], xend=MA_table$Estimate[rows[2]], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=MA_table$CI_lower[rows[1]], xend=MA_table$CI_upper[rows[1]], y=2, yend=2, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$CI_lower[rows[2]], xend=MA_table$CI_upper[rows[2]], y=1, yend=1, size=1.4, col="firebrick") 


# Reproduction Plot
rows<-c(5,6) # Rows in the MA table containing the mean and var estimates
plot_plot<-data.frame(yi=c(effect_list$Repro$lnRR, effect_list$Repro$lnCVR),
                      Precision=1/sqrt(c(effect_list$Repro$v_lnRR, effect_list$Repro$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Repro)), rep("lnCVR", nrow(effect_list$Repro))))

C<-ggplot(plot_plot, aes(x=yi, y=effect, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="C. Reproductive Function") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.15,0.2), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=PI_lower[rows[1]], xend=PI_upper[rows[1]], y=2, yend=2, size=0.5) +
  annotate("segment", x=PI_lower[rows[2]], xend=PI_upper[rows[2]], y=1, yend=1, size=0.5) +
  annotate("segment", x=MA_table$Estimate[rows[1]], xend=MA_table$Estimate[rows[1]], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$Estimate[rows[2]], xend=MA_table$Estimate[rows[2]], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=MA_table$CI_lower[rows[1]], xend=MA_table$CI_upper[rows[1]], y=2, yend=2, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$CI_lower[rows[2]], xend=MA_table$CI_upper[rows[2]], y=1, yend=1, size=1.4, col="firebrick") 

# Reproduction Plot
rows<-c(7,8) # Rows in the MA table containing the mean and var estimates
plot_plot<-data.frame(yi=c(effect_list$Longevity$lnRR, effect_list$Longevity$lnCVR),
                      Precision=1/sqrt(c(effect_list$Longevity$v_lnRR, effect_list$Longevity$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Longevity)), rep("lnCVR", nrow(effect_list$Longevity))))

D<-ggplot(plot_plot, aes(x=yi, y=effect, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="D. Lifespan") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=PI_lower[rows[1]], xend=PI_upper[rows[1]], y=2, yend=2, size=0.5) +
  annotate("segment", x=PI_lower[rows[2]], xend=PI_upper[rows[2]], y=1, yend=1, size=0.5) +
  annotate("segment", x=MA_table$Estimate[rows[1]], xend=MA_table$Estimate[rows[1]], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$Estimate[rows[2]], xend=MA_table$Estimate[rows[2]], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=MA_table$CI_lower[rows[1]], xend=MA_table$CI_upper[rows[1]], y=2, yend=2, size=1.4, col="firebrick") +
  annotate("segment", x=MA_table$CI_lower[rows[2]], xend=MA_table$CI_upper[rows[2]], y=1, yend=1, size=1.4, col="firebrick") 

pdf("plots/figure_2.pdf", height=10, width=10)

grid.arrange(A, B, C, D, ncol=2, nrow=2)

dev.off()

#########################################
###### META-REGRESSIONS TABLE ###########
#########################################

MR_table<-data.frame(Trait=c("Size", "", "", "", "", "Development Time", "", "", "", "", "Reproductive Function", "", "", "", "", "Lifespan", "", "", "", ""),
                     Moderator=rep(c("Max Single Food", "Habitat", "Trophic Level", "Defence", "log N Foods"), 4),
                     lnRR_Estimate=NA,
                     lnRR_CI_lower=NA,
                     lnRR_CI_upper=NA,
                     lnCVR_Estimate=NA,
                     lnCVR_CI_lower=NA,
                     lnCVR_CI_upper=NA)

# Now fit the models to fill this in

#########################################
##### MAX SINGLE FOOD FITNESS ###########
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models
repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_mean) # Significant
repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_var)

# Fit the longevity models
long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_mean)
long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~Max.Single, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_var)

# Fill the in MR table
MR_table$lnRR_Estimate[c(1,6,11,16)]<-c(size_MA_mean$b[2], dev_MA_mean$b[2], repro_MA_mean$b[2], long_MA_mean$b[2])
MR_table$lnRR_CI_lower[c(1,6,11,16)]<-c(size_MA_mean$ci.lb[2], dev_MA_mean$ci.lb[2], repro_MA_mean$ci.lb[2], long_MA_mean$ci.lb[2])
MR_table$lnRR_CI_upper[c(1,6,11,16)]<-c(size_MA_mean$ci.ub[2], dev_MA_mean$ci.ub[2], repro_MA_mean$ci.ub[2], long_MA_mean$ci.ub[2])
MR_table$lnCVR_Estimate[c(1,6,11,16)]<-c(size_MA_var$b[2], dev_MA_var$b[2], repro_MA_var$b[2], long_MA_var$b[2])
MR_table$lnCVR_CI_lower[c(1,6,11,16)]<-c(size_MA_var$ci.lb[2], dev_MA_var$ci.lb[2], repro_MA_var$ci.lb[2], long_MA_var$ci.lb[2])
MR_table$lnCVR_CI_upper[c(1,6,11,16)]<-c(size_MA_var$ci.ub[2], dev_MA_var$ci.ub[2], repro_MA_var$ci.ub[2], long_MA_var$ci.ub[2])

#########################################
##### HABITAT META-REGRESSION ###########
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models
repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_mean) # Significant
repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_var)

# Fit the longevity models
long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_mean)
long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~Habitat, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_var)

# Fill the in MR table
MR_table$lnRR_Estimate[c(1,6,11,16)+1]<-c(size_MA_mean$b[2], dev_MA_mean$b[2], repro_MA_mean$b[2], long_MA_mean$b[2])
MR_table$lnRR_CI_lower[c(1,6,11,16)+1]<-c(size_MA_mean$ci.lb[2], dev_MA_mean$ci.lb[2], repro_MA_mean$ci.lb[2], long_MA_mean$ci.lb[2])
MR_table$lnRR_CI_upper[c(1,6,11,16)+1]<-c(size_MA_mean$ci.ub[2], dev_MA_mean$ci.ub[2], repro_MA_mean$ci.ub[2], long_MA_mean$ci.ub[2])
MR_table$lnCVR_Estimate[c(1,6,11,16)+1]<-c(size_MA_var$b[2], dev_MA_var$b[2], repro_MA_var$b[2], long_MA_var$b[2])
MR_table$lnCVR_CI_lower[c(1,6,11,16)+1]<-c(size_MA_var$ci.lb[2], dev_MA_var$ci.lb[2], repro_MA_var$ci.lb[2], long_MA_var$ci.lb[2])
MR_table$lnCVR_CI_upper[c(1,6,11,16)+1]<-c(size_MA_var$ci.ub[2], dev_MA_var$ci.ub[2], repro_MA_var$ci.ub[2], long_MA_var$ci.ub[2])

#########################################
##### TROPHIC META-REGRESSION ###########
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models
repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_mean)
repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_var)

# Fit the longevity models
long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_mean)
long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~Trophic.Level, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_var)

# Fill the in MR table
MR_table$lnRR_Estimate[c(1,6,11,16)+2]<-c(size_MA_mean$b[2], dev_MA_mean$b[2], repro_MA_mean$b[2], long_MA_mean$b[2])
MR_table$lnRR_CI_lower[c(1,6,11,16)+2]<-c(size_MA_mean$ci.lb[2], dev_MA_mean$ci.lb[2], repro_MA_mean$ci.lb[2], long_MA_mean$ci.lb[2])
MR_table$lnRR_CI_upper[c(1,6,11,16)+2]<-c(size_MA_mean$ci.ub[2], dev_MA_mean$ci.ub[2], repro_MA_mean$ci.ub[2], long_MA_mean$ci.ub[2])
MR_table$lnCVR_Estimate[c(1,6,11,16)+2]<-c(size_MA_var$b[2], dev_MA_var$b[2], repro_MA_var$b[2], long_MA_var$b[2])
MR_table$lnCVR_CI_lower[c(1,6,11,16)+2]<-c(size_MA_var$ci.lb[2], dev_MA_var$ci.lb[2], repro_MA_var$ci.lb[2], long_MA_var$ci.lb[2])
MR_table$lnCVR_CI_upper[c(1,6,11,16)+2]<-c(size_MA_var$ci.ub[2], dev_MA_var$ci.ub[2], repro_MA_var$ci.ub[2], long_MA_var$ci.ub[2])

#########################################
##### DEFENCE META-REGRESSION ###########
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models - NA, Data Deficient
# repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
# summary(repro_MA_mean) 
# repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
# summary(repro_MA_var)

# Fit the longevity models  - NA, Data Deficient
# long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
# summary(long_MA_mean)
# long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~Defence, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
# summary(long_MA_var)

# Fill the in MR table
MR_table$lnRR_Estimate[c(1,6,11,16)+3]<-c(size_MA_mean$b[2], dev_MA_mean$b[2], NA, NA)
MR_table$lnRR_CI_lower[c(1,6,11,16)+3]<-c(size_MA_mean$ci.lb[2], dev_MA_mean$ci.lb[2], NA, NA)
MR_table$lnRR_CI_upper[c(1,6,11,16)+3]<-c(size_MA_mean$ci.ub[2], dev_MA_mean$ci.ub[2], NA, NA)
MR_table$lnCVR_Estimate[c(1,6,11,16)+3]<-c(size_MA_var$b[2], dev_MA_var$b[2], NA, NA)
MR_table$lnCVR_CI_lower[c(1,6,11,16)+3]<-c(size_MA_var$ci.lb[2], dev_MA_var$ci.lb[2], NA, NA)
MR_table$lnCVR_CI_upper[c(1,6,11,16)+3]<-c(size_MA_var$ci.ub[2], dev_MA_var$ci.ub[2], NA, NA)

#########################################
##### N FOODS META-REGRESSION ###########
#########################################

# Fit the size mean and variation model
size_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_mean)
size_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
summary(size_MA_var)

# Fit the development time models
dev_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_mean)
dev_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
summary(dev_MA_var)

# Fit the repro models
repro_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_mean) # Significant
repro_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
summary(repro_MA_var)

# Fit the longevity models
long_MA_mean<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_mean)
long_MA_var<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
summary(long_MA_var)

# Fill the in MR table
MR_table$lnRR_Estimate[c(1,6,11,16)+4]<-c(size_MA_mean$b[2], dev_MA_mean$b[2], repro_MA_mean$b[2], long_MA_mean$b[2])
MR_table$lnRR_CI_lower[c(1,6,11,16)+4]<-c(size_MA_mean$ci.lb[2], dev_MA_mean$ci.lb[2], repro_MA_mean$ci.lb[2], long_MA_mean$ci.lb[2])
MR_table$lnRR_CI_upper[c(1,6,11,16)+4]<-c(size_MA_mean$ci.ub[2], dev_MA_mean$ci.ub[2], repro_MA_mean$ci.ub[2], long_MA_mean$ci.ub[2])
MR_table$lnCVR_Estimate[c(1,6,11,16)+4]<-c(size_MA_var$b[2], dev_MA_var$b[2], repro_MA_var$b[2], long_MA_var$b[2])
MR_table$lnCVR_CI_lower[c(1,6,11,16)+4]<-c(size_MA_var$ci.lb[2], dev_MA_var$ci.lb[2], repro_MA_var$ci.lb[2], long_MA_var$ci.lb[2])
MR_table$lnCVR_CI_upper[c(1,6,11,16)+4]<-c(size_MA_var$ci.ub[2], dev_MA_var$ci.ub[2], repro_MA_var$ci.ub[2], long_MA_var$ci.ub[2])

# Save it
write.table(MR_table, file="tables/MR_table.csv", sep=",", col.names=names(MR_table), row.names=F)

#########################################
######## FIGURE1 FOR META_REG ###########
#########################################

# Colours for points
cols<-c("gold", "cornflowerblue")
size_pes<-0.05
axis.text_size<-12

# First Size by max food 
plot_model_lnRR<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
plot_model_lnCVR<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Size, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
plot_data<-data.frame(yi=c(effect_list$Size$lnRR, effect_list$Size$lnCVR),
                      Precision=1/sqrt(c(effect_list$Size$v_lnRR, effect_list$Size$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Size)), rep("lnCVR", nrow(effect_list$Size))),
                      max_food=c(effect_list$Size$Max.Single, effect_list$Size$Max.Single))
plot_data$cat<-paste0(plot_data$effect, "
                      ", "Max Single ", plot_data$max_food)

A<-ggplot(plot_data, aes(x=yi, y=cat, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="A. Body Size") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model_lnCVR$b[1], xend=plot_model_lnCVR$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$b[2], xend=plot_model_lnCVR$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnCVR$ci.lb[1], xend=plot_model_lnCVR$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$ci.lb[2], xend=plot_model_lnCVR$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") +  
  annotate("segment", x=plot_model_lnRR$b[1], xend=plot_model_lnRR$b[1], y=3-size_pes, yend=3+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$b[2], xend=plot_model_lnRR$b[2], y=4-size_pes, yend=4+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnRR$ci.lb[1], xend=plot_model_lnRR$ci.ub[1], y=3, yend=3, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$ci.lb[2], xend=plot_model_lnRR$ci.ub[2], y=4, yend=4, size=1.4, col="firebrick")

# First Development Time by max food 
plot_model_lnRR<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
plot_model_lnCVR<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$DevTime, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
plot_data<-data.frame(yi=c(effect_list$DevTime$lnRR, effect_list$DevTime$lnCVR),
                      Precision=1/sqrt(c(effect_list$DevTime$v_lnRR, effect_list$DevTime$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$DevTime)), rep("lnCVR", nrow(effect_list$DevTime))),
                      max_food=c(effect_list$DevTime$Max.Single, effect_list$DevTime$Max.Single))
plot_data$cat<-paste0(plot_data$effect, "
                      ", "Max Single ", plot_data$max_food)

B<-ggplot(plot_data, aes(x=yi, y=cat, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="B. Development Time") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model_lnCVR$b[1], xend=plot_model_lnCVR$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$b[2], xend=plot_model_lnCVR$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnCVR$ci.lb[1], xend=plot_model_lnCVR$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$ci.lb[2], xend=plot_model_lnCVR$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") +  
  annotate("segment", x=plot_model_lnRR$b[1], xend=plot_model_lnRR$b[1], y=3-size_pes, yend=3+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$b[2], xend=plot_model_lnRR$b[2], y=4-size_pes, yend=4+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnRR$ci.lb[1], xend=plot_model_lnRR$ci.ub[1], y=3, yend=3, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$ci.lb[2], xend=plot_model_lnRR$ci.ub[2], y=4, yend=4, size=1.4, col="firebrick")

# First Repro by max food 
plot_model_lnRR<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
plot_model_lnCVR<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
plot_data<-data.frame(yi=c(effect_list$Repro$lnRR, effect_list$Repro$lnCVR),
                      Precision=1/sqrt(c(effect_list$Repro$v_lnRR, effect_list$Repro$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Repro)), rep("lnCVR", nrow(effect_list$Repro))),
                      max_food=c(effect_list$Repro$Max.Single, effect_list$Repro$Max.Single))
plot_data$cat<-paste0(plot_data$effect, "
                      ", "Max Single ", plot_data$max_food)

C<-ggplot(plot_data, aes(x=yi, y=cat, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="C. Reroductive Function") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.1,0.2), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model_lnCVR$b[1], xend=plot_model_lnCVR$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$b[2], xend=plot_model_lnCVR$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnCVR$ci.lb[1], xend=plot_model_lnCVR$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$ci.lb[2], xend=plot_model_lnCVR$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") +  
  annotate("segment", x=plot_model_lnRR$b[1], xend=plot_model_lnRR$b[1], y=3-size_pes, yend=3+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$b[2], xend=plot_model_lnRR$b[2], y=4-size_pes, yend=4+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnRR$ci.lb[1], xend=plot_model_lnRR$ci.ub[1], y=3, yend=3, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$ci.lb[2], xend=plot_model_lnRR$ci.ub[2], y=4, yend=4, size=1.4, col="firebrick")

# First Longevity by max food 
plot_model_lnRR<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Longevity, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
plot_model_lnCVR<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Longevity, mods=~Max.Single-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Longevity), data=effect_list$Longevity)
plot_data<-data.frame(yi=c(effect_list$Longevity$lnRR, effect_list$Longevity$lnCVR),
                      Precision=1/sqrt(c(effect_list$Longevity$v_lnRR, effect_list$Longevity$v_lnCVR)),
                      effect=c(rep("lnRR", nrow(effect_list$Longevity)), rep("lnCVR", nrow(effect_list$Longevity))),
                      max_food=c(effect_list$Longevity$Max.Single, effect_list$Longevity$Max.Single))
plot_data$cat<-paste0(plot_data$effect, "
                      ", "Max Single ", plot_data$max_food)

D<-ggplot(plot_data, aes(x=yi, y=cat, color=effect, fill=effect, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="Effect Magnitude", y="", subtitle="D. Lifespan") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("lnRR" = cols[1], "lnCVR" = cols[2])) +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model_lnCVR$b[1], xend=plot_model_lnCVR$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$b[2], xend=plot_model_lnCVR$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnCVR$ci.lb[1], xend=plot_model_lnCVR$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnCVR$ci.lb[2], xend=plot_model_lnCVR$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") +  
  annotate("segment", x=plot_model_lnRR$b[1], xend=plot_model_lnRR$b[1], y=3-size_pes, yend=3+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$b[2], xend=plot_model_lnRR$b[2], y=4-size_pes, yend=4+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model_lnRR$ci.lb[1], xend=plot_model_lnRR$ci.ub[1], y=3, yend=3, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model_lnRR$ci.lb[2], xend=plot_model_lnRR$ci.ub[2], y=4, yend=4, size=1.4, col="firebrick")


pdf("plots/figure_3.pdf", height=10, width=15)

grid.arrange(A, B, C, D, layout_matrix=rbind(c(1,1,2,2),
                                                c(3,3,4,4)))

dev.off()


#########################################
######$## FIGURE FOR META_REG ###########
#########################################

# Colours for points
cols<-c("gold", "cornflowerblue")
size_pes<-0.025
axis.text_size<-12

# First Size lnRR by n_foods
plot_model<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Size, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Size), data=effect_list$Size)
plot_data<-effect_list$Size
plot_data$Precision<-1/sqrt(plot_data$v_lnRR)

A<-ggplot(plot_data, aes(x=ln_n_foods_mix, y=lnRR, size=Precision)) +
  geom_point(alpha=0.3) +
  geom_hline(yintercept=0, lty=2) +
  geom_abline(intercept=plot_model$b[1], slope=plot_model$b[2]) +
  theme_bw() + 
  labs(x="log N Foods", y="lnRR", subtitle="A. Body Size") +
  theme(legend.position = c(0.85,0.8), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size))

# Repro lnRR by n_foods
plot_model<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
plot_data<-effect_list$Repro
plot_data$Precision<-1/sqrt(plot_data$v_lnRR)

B<-ggplot(plot_data, aes(x=ln_n_foods_mix, y=lnRR, size=Precision)) +
  geom_point(alpha=0.3) +
  geom_hline(yintercept=0, lty=2) +
  geom_abline(intercept=plot_model$b[1], slope=plot_model$b[2]) +
  theme_bw() + 
  labs(x="log N Foods", y="lnRR", subtitle="B. Reproductive Function") +
  theme(legend.position = c(0.85,0.15), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size))

# Repro lnCVR by n_foods
plot_model<-rma.mv(yi=lnCVR, V=lnCVR_vcv_list$Repro, mods=~ln_n_foods_mix, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
plot_data<-effect_list$Repro
plot_data$Precision<-1/sqrt(plot_data$v_lnCVR)

C<-ggplot(plot_data, aes(x=ln_n_foods_mix, y=lnCVR, size=Precision)) +
  geom_point(alpha=0.3) +
  geom_hline(yintercept=0, lty=2) +
  geom_abline(intercept=plot_model$b[1], slope=plot_model$b[2]) +
  theme_bw() + 
  labs(x="log N Foods", y="lnCVR", subtitle="C. Reproductive Function") +
  theme(legend.position = c(0.85,0.15), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size))

# DevTime by Trophic Level
plot_model<-rma.mv(yi=lnRR, V=lnRR_vcv_list$DevTime, mods=~Trophic.Level-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$DevTime), data=effect_list$DevTime)
plot_data<-effect_list$DevTime
plot_data$Precision<-1/sqrt(plot_data$v_lnCVR)

D<-ggplot(plot_data, aes(x=lnRR, y=Trophic.Level, color=Trophic.Level, fill=Trophic.Level, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="lnRR", y="", subtitle="D. Development Time") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("Primary.Consumer" = cols[1], "Secondary.Consumer" = cols[2])) +
  theme(legend.position = c(0.85,0.2), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model$b[1], xend=plot_model$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model$b[2], xend=plot_model$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model$ci.lb[1], xend=plot_model$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model$ci.lb[2], xend=plot_model$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") 

# Repro by Habitat
plot_model<-rma.mv(yi=lnRR, V=lnRR_vcv_list$Repro, mods=~Habitat-1, random=list(~1|Comparison.ID, ~1|Consumer.Sp, ~1|Effect.ID), R=list(Consumer.Sp = phyloM$Repro), data=effect_list$Repro)
plot_data<-effect_list$Repro
plot_data$Precision<-1/sqrt(plot_data$v_lnCVR)

E<-ggplot(plot_data, aes(x=lnRR, y=Habitat, color=Habitat, fill=Habitat, size=Precision)) +
  geom_beeswarm(alpha=0.3, corral="random") +
  geom_vline(xintercept=0, lty=2) +
  theme_bw() + 
  labs(x="lnRR", y="", subtitle="E. Reproductive Function") +
  guides(color="none", fill="none") +
  scale_color_manual(values=c("Terrestrial" = cols[1], "Marine" = cols[2])) +
  theme(legend.position = c(0.9,0.15), legend.background = element_blank(), axis.text = element_text(size = axis.text_size), plot.subtitle = element_text(size = axis.text_size)) +
  annotate("segment", x=plot_model$b[1], xend=plot_model$b[1], y=1-size_pes, yend=1+size_pes, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model$b[2], xend=plot_model$b[2], y=2-size_pes, yend=2+size_pes, size=1.4, col="firebrick") + 
  annotate("segment", x=plot_model$ci.lb[1], xend=plot_model$ci.ub[1], y=1, yend=1, size=1.4, col="firebrick") +
  annotate("segment", x=plot_model$ci.lb[2], xend=plot_model$ci.ub[2], y=2, yend=2, size=1.4, col="firebrick") 


pdf("plots/figure_4.pdf", height=10, width=15)

grid.arrange(A, B, C, D, E, layout_matrix=rbind(c(1,1,2,2,3,3),
                                          c(4,4,4,5,5,5)))

dev.off()

#########################################
######### PUBLICATION BIAS ##############
#########################################

# Table to hold results
pb_table<-data.frame(Trait=c("Body Size", "Development Time", "Reproductive Function", "Lifespan"),
                     Reg_test=NA,
                     Reg_test.p=NA,
                     Trim_fill=NA,
                     Trim_fill_adj=NA)

# Size
pb<-rma(yi=effect_list$Size$lnRR, vi=diag(lnRR_vcv_list$Size))
regtest(pb) # Regtest indicates some evidence
pb_table$Reg_test[1]<-regtest(pb)$est
pb_table$Reg_test.p[1]<-regtest(pb)$pval
tf<-trimfill(pb) # Trimfill says no missing studies 
pb_table$Trim_fill[1]<-paste0(tf$k0, " ", tf$side)
pb_table$Trim_fill_adj[1]<-tf$b[1] - pb$b[1] # Trimfill estimates no adjustment

# Development Time
pb<-rma(yi=effect_list$DevTime$lnRR, vi=diag(lnRR_vcv_list$DevTime))
regtest(pb) # Regtest indicates some evidence
pb_table$Reg_test[2]<-regtest(pb)$est
pb_table$Reg_test.p[2]<-regtest(pb)$pval
tf<-trimfill(pb) # Trimfill says no missing studies 
pb_table$Trim_fill[2]<-paste0(tf$k0, " ", tf$side)
pb_table$Trim_fill_adj[2]<-tf$b[1] - pb$b[1] # Trimfill estimates no adjustment

# Repro
pb<-rma(yi=effect_list$Repro$lnRR, vi=diag(lnRR_vcv_list$Repro))
regtest(pb) # Regtest indicates some evidence
pb_table$Reg_test[3]<-regtest(pb)$est
pb_table$Reg_test.p[3]<-regtest(pb)$pval
tf<-trimfill(pb) # Trimfill says no missing studies 
pb_table$Trim_fill[3]<-paste0(tf$k0, " ", tf$side)
pb_table$Trim_fill_adj[3]<-tf$b[1] - pb$b[1] # Trimfill estimates no adjustment

# Lifespan
pb<-rma(yi=effect_list$Longevity$lnRR, vi=diag(lnRR_vcv_list$Longevity))
regtest(pb) # Regtest indicates some evidence
pb_table$Reg_test[4]<-regtest(pb)$est
pb_table$Reg_test.p[4]<-regtest(pb)$pval
tf<-trimfill(pb) # Trimfill says no missing studies 
pb_table$Trim_fill[4]<-paste0(tf$k0, " ", tf$side)
pb_table$Trim_fill_adj[4]<-tf$b[1] - pb$b[1] # Trimfill estimates no adjustment

# Save the table
write.table(pb_table, file="tables/pb_table.csv", row.names=F, col.names=names(pb_table), sep=",")

#########################################
#### CHANGES IN PUBLICATION TRENDS ######
#########################################

for(i in 1:length(effect_list)){
	effect_list[[i]]$new_data<-as.factor(effect_list[[i]]$Year > 2011)
}

size_MA_mean<-rma(yi=lnRR, vi=v_lnRR, mods=~new_data, scale=~new_data, data=effect_list$Size)
summary(size_MA_mean)

dev_MA_mean<-rma(yi=lnRR, vi=v_lnRR, mods=~new_data, scale=~new_data, data=effect_list$DevTime)
summary(dev_MA_mean)

repro_MA_mean<-rma(yi=lnRR, vi=v_lnRR, mods=~new_data, scale=~new_data, data=effect_list$Repro)
summary(repro_MA_mean)

long_MA_mean<-rma(yi=lnRR, vi=v_lnRR, mods=~new_data, scale=~new_data, data=effect_list$Longevity)
summary(long_MA_mean)
