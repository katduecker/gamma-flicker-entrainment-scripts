%% Entrainment gamma rFT
% PhD project 1
% 
% Condition: & gratings
% project gamma oscillations from sensor into source space

% Inputs
% s: subject id
% varargin: virtual channels (submitted automatically with python script)

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function k3_lcmv_gamma(s,varargin)

%% set up paths
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

% read in subj
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

% % is subject one to be excluded?
% exclSUBJ = SUBJ(find(IGF < 58));
% if ismember(SUBJ{s},exclSUBJ)
%     return
% end
SUBJPATH = fullfile(PATHBEAM, SUBJ{s});

%% Load filtered data and spatial filters
load(fullfile(SUBJPATH,'filtdat.mat'),'gambsl','gam')
load(fullfile(SUBJPATH,'spatfilt.mat'))


% load in data and spatial filter
disp('Localizing gamma')

%% LCMV beamformer
% set dummy variables
bslsrc = gambsl;
bslsrc.label = {'vc'};
bslsrc.trial = [];

gamsrc = gam;
gamsrc.label = {'vc'};
gamsrc.trial = [];
bslpow = [];
stimpow = [];
% freq analysis
cfg = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';
cfg.output      = 'pow';
cfg.foilim      = [IGF(s) IGF(s)];
cfg.keeptrials  = 'no';
cfg.tapsmofrq   = 1.5;

% channel indices
vc = cell2mat(varargin);
% multiply with filter
for v = 1:length(vc)
    disp(['Channel ', num2str(vc(v))])
    for tr = 1:length(gambsl.trial)
        bslsrc.trial{tr} = sourcetime_data.avg.filter{vc(v)}*gambsl.trial{tr};
        gamsrc.trial{tr} = sourcetime_data.avg.filter{vc(v)}*gam.trial{tr};
    end
    
    % power at IGF
    freq = ft_freqanalysis(cfg,bslsrc);
    bslpow = [bslpow;freq.powspctrm];
    clear freq
    freq = ft_freqanalysis(cfg,gamsrc);
    stimpow = [stimpow;freq.powspctrm];
    clear freq

end
if exist(fullfile(SUBJPATH,['bslgam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']))
    delete(fullfile(SUBJPATH,['bslgam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']))
end
save(fullfile(SUBJPATH,['bslgam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']),'bslpow','vc')
if exist(fullfile(SUBJPATH,['gam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']))
    delete(fullfile(SUBJPATH,['gam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']))
end
save(fullfile(SUBJPATH,['gam [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']),'stimpow','vc')