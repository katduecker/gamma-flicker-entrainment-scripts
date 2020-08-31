%% Entrainment gamma rFT
% PhD project 1
% 
% Function runs Time-Frequency Analysis (Hanning taper) on first 3 seconds
% of trials in flicker+gamma condition (Grating stimulus without flicker
% on)
% INPUT: s: subject id
% - BATCH: sample
% - oneof: use 1/f corrected data?
% -sliWin: length of sliding window in seconds


% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

function e_gammatron(s,BATCH,oneof,sliWin)

% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
ft_defaults;

% load in 1/f corrected data?
if oneof
    addstr = 'oneof';
else
    addstr = '';
end

PATHCOND = fullfile(MAINPATH,'results','preprocessing','conditions',['Batch_',num2str(BATCH)],addstr);
PATHGAM = fullfile(MAINPATH, 'results', 'power','gammatron',['Batch_',num2str(BATCH)],addstr);
mkdir(PATHGAM)
folds = dir(PATHCOND);                                       % read in subjects
SUBJ = {};

% list subject files
for f = 1:length(folds)
    SUBJ = [SUBJ; folds(f).name];
end
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));
dfreq            = 1/sliWin;                                % frequency steps: depend on sliding window

if sliWin > .5
    PATHGAM = fullfile(MAINPATH, 'results', 'power','gammatron',['Batch_',num2str(BATCH)],addstr,'long slidewin');
end
PATHIN = fullfile(PATHCOND, SUBJ{s});

PATHOUT = fullfile(PATHGAM, SUBJ{s});
mkdir(PATHOUT)

%% Time-frequency analysis
load(fullfile(PATHIN, 'entrainment.mat'))


%% TFR using hanning taper

cfg = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.taper        = 'hanning';
cfg.method       = 'mtmconvol';
cfg.foi          = 30:dfreq:100;
cfg.t_ftimwin    = ones(length(cfg.foi),1).* sliWin;                      % length of time window in seconds
cfg.toi          = [-0.75:0.05:2-0.25];     %sliding over trial in steps of 0.05 % times on which the analysis window should be centered: percentage: overlap
cfg.keeptrials = 'no';
cfg.polyremoval = -1;
cfg.pad          = 8;                                                      % length in seconds to which data can be padded out
cfg.padtype      = 'zero';
TFRgam = ft_freqanalysis(cfg, entrainment);


%% IGF peak

% relative baseline
cfg = [];
cfg.baseline = [-0.75 -0.25];
cfg.baselinetype = 'relchange';
cfg.parameter = 'powspctrm';

rlChanGam = ft_freqbaseline(cfg, TFRgam);

% average over TOI (gratings on, after evoked broadband gamma response)
cfg = [];
cfg.latency = [0.25 1.75];
cfg.avgovertime = 'yes';
rlChanGam = ft_selectdata(cfg,rlChanGam);

save(fullfile(PATHOUT, ['TFR_gammatron_slide_',num2str(sliWin),'.mat']), 'TFRgam')
save(fullfile(PATHOUT, ['IGF_peak_slide_',num2str(sliWin),'.mat']), 'rlChanGam')

clearvars -except *PATH PATH* SUBJ dfreq sliWin cl
end