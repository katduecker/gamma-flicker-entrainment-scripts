%% Entrainment gamma rFT
% Phd project 1
% 
%
% Coherence and phase-locking between diode and MEG signal
% (only plv shown in manus)
% Time-Frequency-Analysis

% INPUTS: 
% - s: subject ID
% - BATCH: sample
% - unit: coh or plv
% - condit: 1: flicker&gratings, 2: flicker
% - oneof: 1/f corrected data? (0,1)
% - sliWin: length of sliding window in seconds


% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function h_coherence(s, BATCH, unit, condit,oneof,sliWin)
%clear all; close all; clc;

%% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH, 'matlab','kd fun'))

ft_defaults;

RESULTPATH = fullfile(MAINPATH, 'results');
PATHCOND = fullfile(RESULTPATH, 'preprocessing','conditions',['Batch_',num2str(BATCH)]);
PATHCOH = fullfile(RESULTPATH, 'coherence',['Batch_',num2str(BATCH)]);
PATHPLV = fullfile(RESULTPATH, 'phase_locking',['Batch_',num2str(BATCH)]);

% coherence or phase locking?
if strcmp(unit,'coh')
    coh = 1;
    plv = 0;
elseif strcmp(unit,'plv')
    coh = 0;
    plv = 1;
end

% entrainment or resonance? = flicker&gratings or flicker?
if condit == 1
    ent = 1;
    res = 0;
elseif condit == 2
    ent = 0;
    res = 1;
end


folds = dir(PATHCOND);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(1:2) = [];


if coh && ent
    PATHOUT = fullfile(PATHCOH, 'entrainment');
elseif coh && res
    PATHOUT = fullfile(PATHCOH, 'resonance');
elseif plv && ent
    PATHOUT = fullfile(PATHPLV, 'entrainment');
elseif plv && res
    PATHOUT = fullfile(PATHPLV, 'resonance');
end

if oneof
   PATHOUT = fullfile(PATHOUT,'oneof');
   PATHCOND = fullfile(PATHCOND,'oneof');

end

PATHIN = fullfile(PATHCOND, SUBJ{s});
PATHOUT = fullfile(PATHOUT,SUBJ{s});
mkdir(PATHOUT)
if ent
    load(fullfile(PATHIN, 'entrainment.mat'))
elseif res
    load(fullfile(PATHIN, 'resonance.mat'))
end



%% separate into freq conditons & adjust sample info

flickfreq = 52:2:90;
if ent
    [COND, ~] = kd_freqconditions(entrainment, 'MISC004', flickfreq);
elseif res
    [COND, ~] = kd_freqconditions(resonance, 'MISC004', flickfreq);
end

clear entrainment resonance

%% Coherence/Phase-locking

% cross-spectral density
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'powandcsd';           % csd: cross-spectral density
cfg.taper = 'hanning';
cfg.foi          = 4:2:120;                                                % frequencies of interest
cfg.numfoi       = length(cfg.foi);
cfg.t_ftimwin    = ones(length(cfg.foi),1).* sliWin;                      % length of time window in seconds

% select time windows: entire trial
if ent
    cfg.toi = [-1+0.5:0.05:6-0.5];
elseif res
    cfg.toi = [-1+0.5:0.05:4-0.5];
end
cfg.keeptrials = 'yes';
cfg.channel = {'MEG', 'MISC004'};
cfg.channelcmb =  {'MEG', 'MISC004'};

for e = 1:length(COND)
    CSD{e} = ft_freqanalysis(cfg, COND{e});
end

clear COND

% coherence
if coh
    cfg = [];
    cfg.method = 'coh';
    cfg.channelcmb ={'MEG', 'MISC004'};
    
    for e = 1:length(CSD)
        COH{e} = ft_connectivityanalysis(cfg, CSD{e});
    end
    
    save(fullfile(PATHOUT, ['coherence_',num2str(sliWin),'.mat']), 'COH')
    clear CSD COH
    
% phase-locking
elseif plv
    cfg = [];
    cfg.method = 'plv';
    cfg.channel = {'MEG','MISC004'};
    for e = 1:length(CSD)
        PLV{e} = ft_connectivityanalysis(cfg, CSD{e});
        
    end
    save(fullfile(PATHOUT, ['phase_lock_',num2str(sliWin),'.mat']), 'PLV')
    clear CSD PLV
end
        
end
