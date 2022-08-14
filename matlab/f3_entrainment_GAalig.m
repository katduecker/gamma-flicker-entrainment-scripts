%% Entrainment gamma rFT
% PhD project 1
%
% Grandaverage TFR in flicker+gamma condition

% Function
% - loads in the results of the TFR in the gammatron+flicker condition
% (f_entrainment.m) per subject
% - computes relative power change using baseline interval identified in
% TOIBSL
% - aligns frequency vector such that the IGF is = 0
% - averages over all subjects

% INPUT
% BATCH: sample
% TOI: time window of interest (flicker+gamma on)
% TOIBSL: time window baseline
% oneof: 1/f corrected data (1 or 0)
% sliwin: length of sliding window in seconds
%
% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

%clear all; close all; clc;
%BATCH = 3;

function f3_entrainment_GAalig(BATCH, TOIBSL,oneof,sliwin)

flickfreq = 52:2:90;                % flicker frequencies

%% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';

if oneof
    addstr = 'oneof';
else
    addstr = '';
end

addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;

% paths
PATHCOND = fullfile(MAINPATH, 'results','preprocessing','conditions',['Batch_',num2str(BATCH)],addstr);
PATHENT = fullfile(MAINPATH,'results','power','entrainment',['Batch_',num2str(BATCH)],addstr);
PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);

if sliwin > .5
    addstr = 'long slidewin';
else
    addstr = '';
end

PATHCOND = fullfile(PATHCOND,addstr);
PATHENT = fullfile(PATHENT,addstr);

% read in subjects
folds = dir(PATHENT);
for f = 1:length(folds)
    SUBJ = [SUBJ;folds(f).name];
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

% store SOIs and IGFs
for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end
% keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
exclSUBJ = SUBJ(find(IGF < 56));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;
%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

% flicker frequencies to be considered
leftLim = minIGF - 52;
rightLim = 90-maxIGF;

% width of frequency range that will be shown in TFR (y-axis)
speclim = [30-minIGF 100-maxIGF];

%% Grandaverage

% preparation: cell array to store TFRs for averaging
% rows: frequencies
% columns: subjects
TFRGAprep = cell(length(-leftLim:2:rightLim),length(SUBJ));

for s = 1:length(SUBJ)
    
    SUBJPATH = fullfile(PATHENT,SUBJ{s});
    
    % load in TFR result
    if oneof || sliwin > .5
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma_',num2str(sliwin),'.mat']))
    else
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma.mat']))
    end
    
    % cut out IGF - leftLim, +rightLim
    freqvec = [IGF(s)-leftLim:2:IGF(s)+rightLim];
   
    % cut out the respective frequency conditions
    TFRal = TFR(find(ismember(flickfreq,freqvec)));
    clear TFR
    
    
    freqIGF = TFRal{1}.freq - IGF(s);               % align frequency vector to IGF
    for t = 1:length(TFRal)
        TFRal{t}.freq = freqIGF;
        % combine planar
        cfg = [];
        cfg.channel = 'MEGGRAD';
        TFRal{t} = ft_selectdata(cfg, TFRal{t});
        cfg = [];
        cfg.method = 'sum';
        TFRal{t} = ft_combineplanar(cfg, TFRal{t});
    end
    
    % relative power change compared to baseline
    cfg = [];
    cfg.baseline = TOIBSL;
    cfg.baselinetype = 'relchange';
    cfg.parameter = 'powspctrm';
    for t = 1:length(TFRal)
        TFR{t} = ft_freqbaseline(cfg, TFRal{t});
    end
    clear TFRal
    
    % select frequencies of interest (width of spectrum, y-axis) and average over channel
    cfg = [];
    cfg.frequency = [speclim(1) speclim(2)];
    cfg.channel = SOI_all{s}';
    cfg.avgoverchan = 'yes';
    for f = 1:length(freqvec)
        TFRGAprep{f,s} = ft_selectdata(cfg,TFR{f});
        TFRGAprep{f,s}.label ={'MEG2032+2033'};   % replace sensor label to be able to average (workaround)
    end
    clear TFRal
end

% Grandaverage
cfg = [];
cfg.paramter = 'powspctrm';
for f = 1:size(TFRGAprep,1)
    TFRGA{f} = ft_freqgrandaverage(cfg, TFRGAprep{f,:});
end

plTFR = TFRGAprep{1};
clear TFRGAprep

% store in smaller structure (to reduce size of file)
for f = 1:length(TFRGA)
    TFRGAent{f} = plTFR;
    TFRGAent{f}.powspctrm = TFRGA{f}.powspctrm;
end
clear TFRGA


save(fullfile(PATHENT,['GA_align_bsl',num2str(TOIBSL(1)),'_sliwi_',num2str(sliwin),'.mat']),'TFRGAent','-v7.3')

end
