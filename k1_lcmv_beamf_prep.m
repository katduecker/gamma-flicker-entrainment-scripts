%% Entrainment gamma rFT
% PhD project 1
% 
%
% LCMV beamformer preparation
% Band-pass filter data (50-92 Hz)
% estimate LCMV spatial filter based on leadfield, headmodel


% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%clear all; close all; clc;

function k3_lcmv_beamf_prep(s)

% set up paths
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH, 'matlab','narrowtag'))
addpath(fullfile(MAINPATH, 'matlab','kd fun'))
addpath(fullfile(MAINPATH,'fieldtrip'));

PATHBEAM = fullfile(MAINPATH, 'results', 'beamformer MNI','LCMV');
ft_defaults;

RAWPATH = fullfile(MAINPATH, 'subjects','Batch_3');
PATHVOL = fullfile(MAINPATH, 'results', 'anatomy - vol, leadf, sourcemodel','Batch_3');
PATHCOND = fullfile(MAINPATH, 'results', 'preprocessing', 'conditions','Batch_3');
PATHGAM = fullfile(MAINPATH, 'results','power','gammatron','Batch_3');
PATHENT = fullfile(MAINPATH, 'results','power','entrainment','Batch_3');
PATHRES = fullfile(MAINPATH, 'results','power','resonance','Batch_3');
PATHPLOT = fullfile(MAINPATH, 'results','plots','beamformer');

%% read in subjects - linux friendly
folds = dir(PATHBEAM);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

for sc = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{sc}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(sc) = gamFreq;
    clear SOI gamFreq
end

% is subject one to be excluded?
exclSUBJ = SUBJ(find(IGF < 58));
if ismember(SUBJ{s},exclSUBJ)
    error('Subject excluded due to IGF <58');
end

SUBJPATH = fullfile(PATHBEAM, SUBJ{s});
mkdir(SUBJPATH)
% Time windows of interest: set in advance
fsample = 1000;
toibsl = [-.75 -.25-1/fsample];
toistim = [.75 1.25-1/fsample];
% Load Template
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
%  load 8 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d4mm'));
%% Source localization IGF
load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))            % load IGF
clear SOI
load(fullfile(PATHVOL, SUBJ{s}, 'vol_lf_source_meggrad4mm.mat'))       % load sourcemodel etc



%% Prepare data

% broadband filter
% append
load(fullfile(PATHCOND, SUBJ{s}, 'entrainment.mat'))
load(fullfile(PATHCOND, SUBJ{s}, 'resonance.mat'))

data = ft_appenddata([],entrainment,resonance);
clear entrainment resonance

disp('Filtering data..')
% bp filter
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 50;
cfg.lpfilter = 'yes';
cfg.lpfreq = 90;
data = ft_preprocessing(cfg,data);
trl = cell2mat(cellfun(@length,data.trial,'UniformOutput',false));

cfg = [];
cfg.trials = trl == 7000;
ent = ft_selectdata(cfg,data);
cfg.trials = trl == 5000;
res = ft_selectdata(cfg,data);
clear data
% load in flicker+gamma condition
cfg = [];
cfg.latency = toibsl;
cfg.channel = 'meggrad';
entbsl = ft_selectdata(cfg,ent);
cfg.latency = toistim;
entstim = ft_selectdata(cfg,ent);
clear ent

% load flicker cond
cfg = [];
cfg.latency = toibsl;
cfg.channel = {'meggrad'};
resbsl = ft_selectdata(cfg,res);
cfg.latency = toistim;
resstim = ft_selectdata(cfg,res);
% load photodiode (important for frequency conditions later
cfg.channel = 'MISC004';
phdiode = ft_selectdata(cfg,res);
clear res

data = ft_appenddata([],entbsl,entstim,resbsl,resstim);

%% Covariance matrix

% replace time vector to not confuse ft_timelock
[data.time{:}] = deal(entbsl.time{1});
% joint covariance matrix
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'no';
cov_data = ft_timelockanalysis(cfg,data);

%% Estimate spatial filter
disp('Estimating spatial filter')

cfg = [];
cfg.method            = 'lcmv';
cfg.grid              = leadf;
cfg.headmodel         = vol;
cfg.keeptrials        = 'no';
cfg.lcmv.projectnoise = 'yes';  % estimate noise distribution
cfg.lcmv.lambda       = '5%';   % data rank deficient ?(insufficient information to estimate model) lambda: regulation
cfg.lcmv.keepfilter   = 'yes';
cfg.lcmv.fixedori     = 'yes';     % fixed dipole orientation
cfg.lcmv.projectmom   = 'yes'; % project dipole time series in direction of maximal power 
cfg.grad              = mGrad;
cfg.senstype          = 'meg';
sourcetime_data       = ft_sourceanalysis(cfg, cov_data);

% channels inside brain
chanin = find(sourcetime_data.inside);

%% separate data again into gamma oscillation and flicker response

cfg = [];
cfg.trials = [1:length(entbsl.trial)];
gambsl = ft_selectdata(cfg,data);
cfg.trials = [length(entbsl.trial)+1:length(entbsl.trial)+length(entstim.trial)];
gam = ft_selectdata(cfg,data);
cfg.trials = [length(entbsl.trial)+length(entstim.trial)+1:...
    length(entbsl.trial)+length(entstim.trial)+length(resbsl.trial)];

flibsl = ft_selectdata(cfg,data);
cfg.trials = [length(entbsl.trial)+length(entstim.trial)+length(resbsl.trial)+1:...
    length(entbsl.trial)+length(entstim.trial)+length(resbsl.trial)+length(resstim.trial)];
fli = ft_selectdata(cfg,data);
clear ent* res*

save(fullfile(SUBJPATH,'filtdat.mat'),'gambsl','gam','flibsl','fli','phdiode')
save(fullfile(SUBJPATH,'spatfilt.mat'),'sourcetime_data','chanin')

% % chanin as csv file for python submission
% if s == 1
%     csvwrite(fullfile(PATHBEAM,'chanidx.csv'),chanin)
% end
