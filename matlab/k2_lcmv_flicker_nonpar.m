%% Entrainment gamma rFT
% PhD project 1: entrainment
% 
% Condition: flicker
% project flicker response from sensor into source space

% Inputs
% s: subject id
% varargin: virtual channels (submitted automatically with python script)

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function k3_lcmv_flicker_nonpar(s,varargin)


%% settings
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

SUBJPATH = fullfile(PATHBEAM, SUBJ{s});

% store SOIs and IGFs
for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    clear SOI gamFreq
end

%% load filtered data and spatial filter
load(fullfile(SUBJPATH,'filtdat.mat'),'flibsl','fli','phdiode')
load(fullfile(SUBJPATH,'spatfilt.mat'))

% virtual channel indices (submitted with python script)
vc = cell2mat(varargin);

%% Beamforming
% Push flicker response through spatial filter
flickfreq = 52:2:90;

[bslcond, flickcond] = freqcond(phdiode,flibsl,fli,flickfreq);

% only consider frequencies up to 78 Hz
flickfreq = 52:2:78;

% empty variables to fill in
bslpow = zeros(length(vc),length(flickfreq));                % baseline
stimpow = zeros(length(vc),length(flickfreq));               % stimulation

% set up empty cells for parallel computing

% freqanalysis cfg
cfg = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';
cfg.output      = 'pow';
cfg.keeptrials  = 'no';
cfg.tapsmofrq   = 1.5;

% empty cells for spatial filter
bslsrc      = cell(1,length(flickfreq));
flisrc      = bslsrc;

for f = 1:length(flickfreq)
    % place holder cells
    bslsrc{f} = bslcond{f};
    bslsrc{f}.label = {'vc'};
    bslsrc{f}.trial = [];
    flisrc{f} = bslsrc{f};
end

% for each channel, loop over all frequencies (up to 78 Hz)
for v = 1:length(vc)
    for f = 1:length(flickfreq)
        % filter
        bslsrc{f}.trial = [];
        flisrc{f}.trial = [];
        for tr = 1:length(bslcond{f}.trial)
            bslsrc{f}.trial{tr} = sourcetime_data.avg.filter{vc(v)}*bslcond{f}.trial{tr};
            flisrc{f}.trial{tr} = sourcetime_data.avg.filter{vc(v)}*flickcond{f}.trial{tr};
        end
        
        cfg.foilim      = [flickfreq(f) flickfreq(f)];

        % fft
        freq = ft_freqanalysis(cfg,bslsrc{f});
        bslpow(v,f) = freq.powspctrm;
        clear freq
        freq = ft_freqanalysis(cfg,flisrc{f});
        stimpow(v,f) = freq.powspctrm;
        clear freq
    end
end

save(fullfile(SUBJPATH,['bslfli [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']),'bslpow','vc')
save(fullfile(SUBJPATH,['fli [',num2str(vc(1)),' ', num2str(vc(end)),'].mat']),'stimpow','vc')

%% Function: separate into frequency conditions
function [bslcond, flickcond] = freqcond(phdiode,flibsl,fli,flickfreq)
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.taper = 'hanning';
cfg.foi = 30:2:100;
cfg.tapsmofrq = 2;

% photodiode
MISC = ft_freqanalysis(cfg, phdiode);

% plot power spectra and save frequencies in frequency vector
freqvec = [];
for p = 1:size(MISC.powspctrm,1)
    [pow pos] = max(squeeze(MISC.powspctrm(p,:,:)));
    freq = round(MISC.freq(pos));
    freqvec = [freqvec freq];
end

bslcond = cell(1, length(flickfreq));             % condition cells
flickcond = bslcond;
cfg = [];
cfg.avgoverrpt = 'no';                             % average over trial = avg over freq
for f = 1:length(flickfreq)
    findFRQ = logical(freqvec == flickfreq(f));        % find trials in which this frequency was used
    cfg.trials = findFRQ;
    bslcond{f} = ft_selectdata(cfg, flibsl);
    flickcond{f} = ft_selectdata(cfg, fli);
end