%% Entrainment gamma rFT
% PhD project 1

% flicker+gamma condition: Time-Frequency Analysis, using Hanning Taper

% INPUT
% - s: subject index
% - BATCH: sample
% - oneof: 1/f corrected data?
% - sliWin: length of sliding window in seconds

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

function f_entrainment(s, BATCH,oneof,sliWin)

%% settings

% paths
ft_defaults;
MAINPATH = '/rds/projects/2018/jenseno-entrainment';

if oneof
    addstr = 'oneof';
else
    addstr = '';
end
addpath(fullfile(MAINPATH,'matlab','kd fun'))
PATHCOND = fullfile(MAINPATH, 'results','preprocessing','conditions',['Batch_',num2str(BATCH)],addstr);


% list subjects (linux friendly)
folds = dir(PATHCOND);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));

% subject path in
PATHIN = fullfile(PATHCOND, SUBJ{s});

% subject path out
if sliWin > .5
    PATHOUT = fullfile(MAINPATH, 'results', 'power','entrainment',['Batch_',num2str(BATCH)],addstr,'long slidewin',SUBJ{s});
else
    PATHOUT = fullfile(MAINPATH,'results','power','entrainment',['Batch_',num2str(BATCH)],addstr,SUBJ{s});
end
mkdir(PATHOUT)

dfreq            = 1/sliWin;

% fliceker frequencies
flickfreq = 52:2:90;


% load in clean subject flicker+gamma data
load(fullfile(PATHIN, 'entrainment.mat'))

%% Frequency conditions

% separate into freq conditons & adjust sample info
[ENTR, ~] = kd_freqconditions(entrainment, 'MISC004', flickfreq);
clear entrainment

%% Time-Frequency Analysis

TFR = cell(1,length(ENTR));

cfg = [];
cfg.output       = 'pow';
cfg.channel      = {'MISC004', 'MEG'};
cfg.taper        = 'hanning';
cfg.method       = 'mtmconvol';                                            
cfg.foi          = 4:dfreq:100;                                            % frequencies of interest
cfg.numfoi       = length(cfg.foi);
cfg.pad          = 8;                                                      % length in seconds to which data is padded out
cfg.padtype      = 'zero';

cfg.t_ftimwin    = ones(length(cfg.foi),1).* sliWin;                      % length of time window in seconds
cfg.toi          = [-0.75:0.05:6-0.25];                                   % sliding over trial in steps of 0.05 
cfg.keeptrials = 'no';

% TFR separately for each flicker frequencies
for d = 1:length(ENTR)
    TFR{d} = ft_freqanalysis(cfg, ENTR{d});
    
    % number of trials
    numTrl{d} = length(ENTR{d}.trial);
end
clear ENTR


save(fullfile(PATHOUT, ['TFR_ent_alpha_gamma_',num2str(sliWin),'.mat']), 'TFR', 'numTrl')
clear TFR numTrlc

clearvars -except *PATH PATH* SUBJ SID sliWin dfreq BATCH flickfreq
end

%run('f3_entrainment_aligFreq.m')
