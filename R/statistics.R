### Phd project gamma entrainment visual flicker

# ANOVA & dependent sample t-test: power change at IGF during flicker

# [c] PGR: K. Duecker,  kxd888@student.bham.ac.uk
# supervisor: O. Jensen, Centre for Human Brain Health, Birmingham UK

# script: last updated 31.08.2020


install.packages("compute.es")
install.packages("ez")
install.packages("ggplot2")
install.packages("multcomp")
install.packages("nlme")
install.packages("pastecs")
install.packages("reshape")
install.packages("WRS2")
install.packages("Hmisc")
install.packages("Rcmdr")
install.packages("DSUR")
install.packages("Rtools")
install.packages("dplyr")
install.packages("svglite")
#Initiate packages
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
#library(Rcmdr)

## ANOVA ####################################

ANOVAdir <- "Z:\\results\\statistics\\ANOVA\\Batch_3\\"
setwd(ANOVAdir)

# load in power values (averaged TFRs)
TFR1_3 <- read.table("TFRT1vsT3_int.csv",
                     header = TRUE,
                     sep = ",",
                     dec =".")

# number of subjects
numSUB = nrow(TFR1_3);
subID <- factor(seq(1,numSUB,1))

# data frame
TFR <- data.frame( "participant" = subID,
                   "belT1" = TFR1_3$belT1,
                   "belT3" = TFR1_3$belT3,
                   "aboT1" = TFR1_3$aboT1,
                   "aboT3" = TFR1_3$aboT3)

rm(TFR1_3)

# melt to long data frame
longTFR <- melt(TFR, id = "participant", measured = c("belT1", "belT3", "aboT1", "aboT3"))
names(longTFR) <- c("participant", "timexfreq", "powchan")
longTFR$freq <- gl(2,numSUB*2, labels = c("below", "above"))
longTFR$time <- gl(2,numSUB,numSUB*4, labels =c("T1","T2"))

by(longTFR$powchan, list(longTFR$freq, longTFR$time), stat.desc)

# plot & identify outliers
TFRBoxplot <- ggplot(longTFR, aes(freq, powchan))
TFRBoxplot + geom_boxplot() + facet_wrap(~time, nrow = 1) + labs(x = "Frequency relative to IGF", y = "relative power change [a.u.]")
imageFile <- paste(ANOVAdir,"boxplot_belvsaboT1T3.emf",sep="\\")
ggsave(file = imageFile, device = "emf", dpi = 1000)


# contrasts
BelvsAbo <- c(-1,1)    # below vs above IGF
contrasts(longTFR$freq) <- BelvsAbo

T1vsT3 <- c(-1,1)      # T1 vs T3
contrasts(longTFR$time) <- T1vsT3


# Anova
TFRModel <- ezANOVA(data = longTFR, dv = .(powchan), wid = .(participant), within = .(freq,time), type = 3, detailed = TRUE)
TFRModel


##### T-Test relative power change T1 vs T2 below vs above ########################################################


Tdir <- "Z:\\results\\statistics\\T-Test\\Batch_3"
setwd(Tdir)

# load data
TFR1_3 <- read.table("TFRChangebelvsabovsout.csv",
                     header = TRUE,
                     sep = ",",
                     dec =".")

# number of participants
numSUB = nrow(TFR1_3);
subID <- factor(seq(1,numSUB,1))

# data frame
TFR <- data.frame( "participant" = subID,
                   "below IGF" = TFR1_3$belRC,
                   "above IGF" = TFR1_3$aboRC)

rm(TFR1_3)

# melt into long format
longTFR <- melt(TFR, id = "participant", measured = c("belRC", "aboRC"))
names(longTFR) <- c("participant", "freq", "powchan")

# plot and check
TFRBoxplot <- ggplot(longTFR, aes(freq, powchan))
TFRBoxplot + geom_boxplot() + theme(axis.text.x= element_text(size = 16)) + theme(axis.text.y= element_text(size = 16)) +labs(x = "Frequency", y = "relative power change T1 vs T2")
imageFile <- paste(Tdir,"boxplot_RCbelabo.svg",sep="\\")
ggsave(file = imageFile, device = "svg", dpi = 600)
imageFile <- paste(Tdir,"boxplot_RCbelabo.emf",sep="\\")

# n<30 -> SHapiro-Wilk test for normality of differences
d <- with(longTFR, powchan[freq == "below.IGF"] - powchan[freq == "above.IGF"])
shapiro.test(d)   # normally distributed

dep.t.test <- t.test(powchan ~ freq, data = longTFR, paired = TRUE,alternative = "less",conf.level=0.95)

dep.t.test
stat.desc(TFR, basic = FALSE)

# effect size
t <- dep.t.test$statistic[[1]]
df <- dep.t.test$parameter[[1]]
r <- t^2/(t^2+df)

# Yuen's test based on trimmed means
yuend(TFR$below.IGF, TFR$above.IGF, tr = .2 )     # robust
Dqcomhd(TFR$below.IGF, TFR$above.IGF,nboot = 200, q = c(0.25, 0.5, 0.75))
