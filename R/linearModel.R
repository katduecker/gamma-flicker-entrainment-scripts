### Phd project gamma entrainment visual flicker

# fit regression model to power (and plv as function of frequency

# [c] PGR: K. Duecker,  kxd888@student.bham.ac.uk
# supervisor: O. Jensen, Centre for Human Brain Health, Birmingham UK

# script: last updated 31.08.2020

install.packages("car");
install.packages("QuantPsyc")
install.packages("reshape2")
install.packages("GeneNet")
install.packages("sandwich")    # solve the problem of heteroscedasticity in the data
install.packages("Hmisc")
library(car)
library(QuantPsyc)
library(reshape2)
library(GeneNet)
library(sandwich)
library(BayesFactor)
library(Hmisc)
setwd("z:\\results\\statistics\\Linear regression\\Batch_3")

leftLim = -6          # minimum IGF distance
rightLim = 12         # maximum IGF distance

##### power as a function of frequency ##########################################

cond = "entrainment"    # condition

# load data: distance frequency to IGF
DISTIGF <-  read.csv(paste("pow_fun_IGF_",cond,"_-0.75_RFT.csv",sep=""),
                     sep = ",",
                     header = TRUE)

numSub  = nrow(DISTIGF);
names(DISTIGF) <- seq(leftLim, rightLim, by =2)
DISTIGF$subj <- seq(1,numSub,1)

# melt into long format and change to numeric
DISTIGFlong <- melt(DISTIGF,id="subj", measured ="distIGF")
names(DISTIGFlong) <- c("subj","freq","powchan")
DISTIGFlong$powchan <- as.numeric(sub(",", ".", DISTIGFlong$powchan))
DISTIGFlong$freq <- as.numeric(as.character(DISTIGFlong$freq))

# power as a function of distance to the IGF
IGFmodel <- lm(powchan~freq, data = DISTIGFlong)
summary(IGFmodel)

## test: Homoscedasticity
#leveneTest(FREQdatlong$powchan,FREQdatlong$freq, center = mean) # violated!
leveneTest(DISTIGFlong$powchan,DISTIGFlong$freq,center = mean) # violated!

# violated! Use ln(power) instead! (not reported in manus)
DISTIGFlog <-  read.csv(paste("pow_fun_IGF_",cond,"_RFT_log.csv",sep=""),
                     sep = ",",
                     header = TRUE)

numSub  = nrow(DISTIGFlog);
names(DISTIGFlog) <- seq(leftLim, rightLim, by =2)
DISTIGFlog$subj <- seq(1,numSub,1)

DISTIGFlog <- melt(DISTIGFlog,id="subj", measured ="distIGF")
names(DISTIGFlog) <- c("subj","freq","log_pow")
DISTIGFlog$freq <- as.numeric(as.character(DISTIGFlog$freq))

IGFmodel <- lm(log_pow~freq, data = DISTIGFlog)
summary(IGFmodel)

### plv as function of frequency #################################################################################
cond = "resonance"    # condition

# load data: distance frequency to IGF
DISTIGF <-  read.csv(paste("plv_fun_IGF_",cond,"_RFT_sliwin0.5.csv",sep=""),
                     sep = ",",
                     header = TRUE)

# melt into long format
names(DISTIGF) <- seq(leftLim,rightLim,by=2)
numSub  = nrow(DISTIGF);
DISTIGF$subj <- seq(1,numSub,1)
DISTIGFlong <- melt(DISTIGF,id="subj", measured ="distIGF")
names(DISTIGFlong) <- c("subj","freq","plv")

# change to numeric
DISTIGFlong$plv <- as.numeric(sub(",", ".", DISTIGFlong$plv))
DISTIGFlong$freq <- as.numeric(as.character(DISTIGFlong$freq))

# power as a function of distance to the IGF
IGFmodel <- lm(plv~freq, data = DISTIGFlong)
summary(IGFmodel)


## test: Homoscedasticity
#leveneTest(PLVdatlong$plv,PLVdatlong$freq, center = mean) # not violated!

leveneTest(DISTIGFlong$plv,DISTIGFlong$freq)


##### Binomial Test: IGF > resonance frequency
setwd("Z:\\results\\power")
resfreq <-  read.csv("res_props_sli_0.5.csv",
                     sep = ",",
                     header = TRUE)

# resonance frequency flicker&gratings
igf_fligam <- resfreq$ï..IGF > resfreq$res_gamma

btest = binom.test(sum(igf_fligam),length(igf_fligam),p=0.5)
btest 

btest_bf = proportionBF(sum(igf_fligam),length(igf_fligam),0.5)
btest_bf

# resonance frequency flicker only
igf_fli <- resfreq$ï..IGF > resfreq$res_rft
igf_fli

btest = binom.test(sum(igf_fli),length(igf_fli),p=0.5)
btest 

btest_bf = proportionBF(sum(igf_fli),length(igf_fli),0.5)
btest_bf

# correlation between resoance frequencies
fli_fligam=rcorr(resfreq$res_gamma,resfreq$res_rft,type="pearson")
fli_fligam

fligam_fli <- resfreq$res_gamma > resfreq$res_rft
fligam_fli
# -uncorrelated

# is flicker freq higher when gratings are on?
# -> not significant
# -> bayes goes in that direction, but inconclusive
btest = binom.test(sum(fligam_fli),length(fligam_fli),p=0.5)
btest 
btest_bf = proportionBF(sum(fligam_fli),length(fligam_fli),0.5)
btest_bf