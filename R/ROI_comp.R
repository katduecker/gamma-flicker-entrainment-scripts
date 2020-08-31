### Phd project gamma entrainment visual flicker

# dependent sample t-test on peak mni coordinates

# [c] PGR: K. Duecker,  kxd888@student.bham.ac.uk
# supervisor: O. Jensen, Centre for Human Brain Health, Birmingham UK

# script: last updated 31.08.2020

# Comparison peak coordinates

library(compute.es)
library(ez)
library(ggplot2)
library(multcomp)
library(nlme)
library(pastecs)
library(reshape)
library(WRS2)
source("http://www-rcf.usc.edu/~rwilcox/Rallfun-v14")
library(Hmisc)
library(DSUR)
library(dplyr)
library(svglite)
install.packages("pwr");
library(pwr)

BEAMdir <- "Z:\\results\\beamformer MNI\\LCMV"
setwd(BEAMdir)

roi <- read.table("ROI_comp4mm.csv",
                     header = TRUE,
                     sep = ";")

# t-test
roi_x <- t.test(roi$x_IGF,roi$x_FL, alternative = "less", paired = TRUE)
roi_y <- t.test(roi$y_IGF,roi$y_FL, alternative = "less", paired = TRUE)
roi_z <- t.test(roi$z_IGF,roi$z_FL, alternative = "less", paired = TRUE)

# corrected p values
pvalscor <- p.adjust(c(roi_x$p.value,roi_y$p.value,roi_z$p.value), method = "bonferroni", n = 3)

# power analysis: effect size
t <- roi_z$statistic[[1]]
df <- roi_z$parameter[[1]]
r <- t^2/(t^2+df)

