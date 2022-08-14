## Perform PCA on 3D MNI coordinates
install.packages("effsize")
library(BayesFactor)
library(effsize)
lcmv_dir <- "Z:\\results\\beamformer MNI\\LCMV\\"
setwd(lcmv_dir)

# load in power value
coord <- read.table("ROI_comp4mm_all.csv",
                     header = TRUE,
                     sep = ",",
                     dec =".")


# data frames

#IGF
igfdf <- data.frame("x" = coord$x_IGF,
                   "y" = coord$y_IGF,
                   "z" = coord$z_IGF)

# flicker
flidf <- data.frame("x" = coord$x_FL,
                    "y" = coord$y_FL,
                    "z" = coord$z_FL)

#flicker+grating
fgdf <- data.frame("x" = coord$x_FG,
                    "y" = coord$y_FG,
                    "z" = coord$z_FG)

# bind
x = rbind(igfdf,flidf,fgdf)

#pca
pc_mni = prcomp(x,scale=TRUE)
biplot(pc_mni)

pc_mni$var = pc_mni$sdev^2
pc_mni$var

# proportion of variance explained
pve = pc_mni$var/sum(pc_mni$var)
pve

# t-test along first PC
l_igf <- pc_mni$x[1:22,1]
l_fli <- pc_mni$x[23:44,1]
l_fg  <- pc_mni$x[45:66,1]
t_igf_fli <- t.test(l_igf,l_fli,paired=TRUE)
t_igf_fli 
cohd_igf_fli = cohen.d(l_igf,l_fli,paired=TRUE)
t_igf_g <- t.test(l_igf,l_fg,paired=TRUE)
t_igf_g
cohd_igf_fg = cohen.d(l_igf,l_fg,paired=TRUE)
t_fli_fg <- t.test(l_fli,l_fg,paired=TRUE)
t_fli_fg

# correct for multiple comparisons
pvalscor <- p.adjust(c(t_igf_fli$p.value,t_igf_g$p.value,t_fli_fg$p.value), method = "BH", n = 3)
pvalscor


# Bayes
# igf fli
diffsc = l_igf - l_fli
t_igf_fli_bayes = ttestBF(x=diffsc)
t_igf_fli_bayes

# igf fligrat
diffsc = l_igf - l_fg
t_igf_fg_bayes = ttestBF(x=diffsc)
t_igf_fg_bayes

# fli fligrat
diffsc = l_fli - l_fg
t_fli_fg_bayes = ttestBF(x=diffsc)
t_fli_fg_bayes

############# OTHER PCs

##t-test along 2nd PC
l_igf <- pc_mni$x[1:22,2]
l_fli <- pc_mni$x[23:44,2]
l_fg  <- pc_mni$x[45:66,2]
t_igf_fli <- t.test(l_igf,l_fli,paired=TRUE)
t_igf_fli 
cohd_igf_fli = cohen.d(l_igf,l_fli,paired=TRUE)

t_igf_g <- t.test(l_igf,l_fg,paired=TRUE)
t_igf_g
t_fli_fg <- t.test(l_fli,l_fg,paired=TRUE)
t_fli_fg
pvalscor <- p.adjust(c(t_igf_fli$p.value,t_igf_g$p.value,t_fli_fg$p.value), method = "BH", n = 3)
pvalscor


# t-test along 3rd PC
l_igf <- pc_mni$x[1:22,3]
l_fli <- pc_mni$x[23:44,3]
l_fg  <- pc_mni$x[45:66,3]
t_igf_fli <- t.test(l_igf,l_fli,paired=TRUE)
t_igf_fli 
t_igf_g <- t.test(l_igf,l_fg,paired=TRUE)
t_igf_g
t_fli_fg <- t.test(l_fli,l_fg,paired=TRUE)
t_fli_fg
pvalscor <- p.adjust(c(t_igf_fli$p.value,t_igf_g$p.value,t_fli_fg$p.value), method = "BH", n = 3)
pvalscor

# save pc's
write.csv(pc_mni$x,'pc_mni.csv')
# means
write.csv(pc_mni$center,'pc_mni_m.csv')
# std
write.csv(pc_mni$scale,'pc_mni_std.csv')






###########################

# PC's per comparison?

# igf vs flicker
igf_fli = rbind(igfdf,flidf)

#pca
pc_igf_fli = prcomp(igf_fli,scale=TRUE)
biplot(pc_igf_fli)

pc_igf_fli$var = pc_igf_fli$sdev^2
pc_igf_fli$var

# proportion of variance explained
pve = pc_igf_fli$var/sum(pc_igf_fli$var)
pve

# t-test along first PC
l_igf <- pc_igf_fli$x[1:22]
l_fli <- pc_igf_fli$x[23:44]
t_igf_fli <- t.test(l_igf,l_fli)


# IGF vs flicker&gratings
igf_fg = rbind(igfdf,fgdf)

#pca
pc_igf_fg = prcomp(igf_fg,scale=TRUE)
biplot(pc_igf_fli)

pc_igf_fg$var = pc_igf_fg$sdev^2
pc_igf_fg$var

# proportion of variance explained
pve = pc_igf_fg$var/sum(pc_igf_fg$var)
pve

# t-test along first PC
l_igf <- pc_igf_fg$x[1:22]
l_fg <- pc_igf_fg$x[23:44]
t_igf_fg <- t.test(l_igf,l_fli)
t_igf_fg

# flicker vs flicker&gratings
fli_fg = rbind(flidf,fgdf)

#pca
pc_fli_fg = prcomp(fli_fg,scale=TRUE)
biplot(pc_fli_fg)

pc_fli_fg$var = pc_fli_fg$sdev^2
pc_fli_fg$var

# proportion of variance explained
pve = pc_fli_fg$var/sum(pc_fli_fg$var)
pve

# t-test along first PC
l_fli <- pc_fli_fg$x[1:22]
l_fg <- pc_fli_fg$x[23:44]
t_fli_fg <- t.test(l_fli,l_fg)
t_fli_fg

######

# PCA igf
pc_igf = prcomp(igfdf)
biplot(pc_igf)

pc_igf$var = pc_igf$sdev^2
pc_igf$var

# proportion of variance explained
pve = pc_igf$var/sum(pc_igf$var)
pve

# PCA flicker
pc_fli = prcomp(igfdf)
biplot(pc_igf)

pc_igf$var = pc_igf$sdev^2
pc_igf$var

# PCA flicker+gratings


## old data 

coord <- read.table("ROI_comp4mm.csv",
                    header = TRUE,
                    sep = ";",
                    dec =",")

igfdf <- data.frame("x" = coord$x_IGF,
                    "y" = coord$y_IGF,
                    "z" = coord$z_IGF)

# flicker
flidf <- data.frame("x" = coord$x_FL,
                    "y" = coord$y_FL,
                    "z" = coord$z_FL)

x = rbind(igfdf,flidf)

# PCA
pc_igf_fli = prcomp(x, scale=TRUE)
biplot(pc_igf_fli)

pc_igf_fli$var = pc_igf_fli$sdev^2
pc_igf_fli$var

# proportion of variance explained
pve = pc_igf_fli$var/sum(pc_igf_fli$var)
pve

# PCA on differences

diff_igf_fli$x = 