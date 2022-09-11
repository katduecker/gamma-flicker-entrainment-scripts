# Entrainment of neuronal gamma oscillations using Rapid Frequency Tagging

This repository contains all scripts developed for the analysis of Magnetoencephalography (MEG) data for [Duecker et al., 2021, J Neurosci](https://www.jneurosci.org/content/41/31/6684.abstract).

## MEG data analysis in MATLAB

### Preprocessing 

1. (a1) semi-automatic artefact rejection
2. (a2) Independent Component analysis & component suppression
3. (a3) estimating source models and leadfield matrix for Linearly Constrained Minimum Variance (LCMV) beamformer [(van Veen et al., 1997)](https://pubmed.ncbi.nlm.nih.gov/9282479/)
4. (c) separate cleaned data into conditions (tagging frequencies) based on signal in photodiode

### Gamma oscillations & flicker responses - power
2. (e) Identify individual gamma frequency: Time-Frequency analysis of power for the moving grating (without flicker) interval
3. (f) Gamma oscillations and flicker response: Time-Frequency analysis of power for all trials in flicker&grating (in code "entrainment") condtion
4. (g) Response to invisible flicker: TF analysis of power for all trials in flicker (in code "resonance") condition

### Gamma oscillations & flicker responses - coherence
5. (h) time frequency analysis of coherence between MEG sensors of interest and photodiode signal

6. (i) plot average power and coherence as a function of frequency for both conditions

### LCMV Beamformer
7. (k1) estimate spatial filter based on the covariance matrix over all data (including all conditions)
8. (k2) project flicker response and gamma oscillations into source space (40,000 virtual channels are broken up into small subsets)
9. (k3) concatenate subsets of virtual channels & plot results
10. (k4) project flicker&grating interval into source space and contrast to grating only interval to isolate flicker response

### Phase analysis
11. (n) calculate phase difference between MEG sensors and photodiode
12. (n2) algorithm to find plateaus

## Statistical analysis in R

1. statistics: ANOVA on power at IGF before and during flicker, for flicker frequencies above and below IGF
2. Linear Model fit to power and coherence as a function of frequency
3. PCA_MNI: principal component analysis of peak MNI coordinates of gamma oscillations and flicker response in conditions "flicker" and "flicker&grating";
comparison of coordinates along first principal component
