%% Entrainment gamma rFT
% Phd project 1
%
%
% Flicker condition
% relative power change in response to flicker compared to baseline


% INPUT
% - s: subject index
% - BATCH: sample
% - which_resp: 'induced' (used in published data) or 'evoked'
% - oneof: 1/f corrected data?
% - sliWin: length of sliding window in seconds

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function g_resonance(s, BATCH,which_resp,oneof,sliWin)

% settings
ft_defaults;
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH,'matlab','kd fun'))
PATHCOND = fullfile(MAINPATH, 'results','preprocessing','conditions',['Batch_',num2str(BATCH)]);

folds = dir(PATHCOND);                                       % read in subjects
SUBJ = {};
for f = 1:length(folds)
    SUBJ = [SUBJ; folds(f).name];
end
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));

PATHIN = fullfile(PATHCOND, SUBJ{s});

if strcmp(which_resp,'induced')
    induced = 1;
elseif strcmp(which_resp,'evoked')
    evoked = 0;
else
    error('Define whether you want to look at the induced or evoked response.')
end

if oneof
    addstr = 'oneof';
else
    addstr = '';
end
if sliWin > .5
    PATHOUT = fullfile(MAINPATH, 'results', 'power','resonance',['Batch_',num2str(BATCH)],addstr,'long slidewin',SUBJ{s});
else
    PATHOUT = fullfile(MAINPATH,'results','power','resonance',['Batch_',num2str(BATCH)],addstr,SUBJ{s});
end
mkdir(PATHOUT)

dfreq            = 1/sliWin;
    
load(fullfile(PATHIN, 'resonance.mat'));
    
    
% separate into frequency conditons & adjust sample info
flickfreq = 52:2:90;
[RESON, ~] = kd_freqconditions(resonance, 'MISC004', flickfreq);

clear resonance

%% Time frequency analysis

TFR = cell(1,length(RESON));

cfg = [];
cfg.output       = 'pow';
cfg.channel      = {'MISC004','MEG'};
cfg.taper        = 'hanning';
cfg.method       = 'mtmconvol';                                                   % multitaper: creates several orthogonal tapering windows
cfg.foi          = 4:dfreq:100;                                                % frequencies of interest
cfg.numfoi       = length(cfg.foi);
cfg.t_ftimwin    = ones(length(cfg.foi),1).* sliWin;                      % length of time window in seconds
cfg.toi          = [-0.75:0.05:4-0.25];     %sliding over trial in steps of 0.02 % times on which the analysis window should be centered: percentage: overlap
cfg.keeptrials = 'no';

for r = 1:length(RESON)
    TFR{r} = ft_freqanalysis(cfg, RESON{r});
    numTrl{r} = length(RESON{r}.trial);
end

clear RESON
% plot & check if data is noisy & find the right sensors
%     for t = 1:length(TFR)
%
%         TFR{t}.time = TFR{t}.time - 1;
%         %     % combine planar
%         %     cfg = [];
%         %     cfg.method = 'sum';
%         %     TFRplan{t} = ft_combineplanar(cfg, TFR{t});
%         %
%         %     figure;
%         %     plot(TFR{t}.freq, squeeze(mean(TFR{t}.powspctrm(find(strcmp(TFR{t}.label, 'MISC004')),:,:),3)))
%         %     xticks(30:4:100)
%         %     title(['Frequency: ', num2str(flickfreq(t))])
%         %
%         %     cfg = [];
%         %     cfg.baseline = [-0.5 0];
%         %     cfg.baselinetype = 'relchange';
%         %     %cfg.zlim = [0 10];
%         %     cfg.layout = 'neuromag306cmb.lay';
%         %     figure; ft_multiplotTFR(cfg,TFRplan{t})
%         %     disp(['current frequency: ', num2str(flickfreq(t))]);
%         %     pause
%         %     yticks(30:2:100)
%         %     close all
%     end
%     % plotTFR = {cnoise, SOI, zval};
%     if induced
%         save(fullfile(PATHOUT, 'TFRresonance_cycfreq.mat'), 'TFR', 'numTrl')
%     elseif evoked
%         save(fullfile(PATHOUT, 'TFR_evoked.mat'), 'TFR', 'plotTFR')
%     end
save(fullfile(PATHOUT, ['TFR_res_alpha_RFT_',num2str(sliWin),'.mat']), 'TFR', 'numTrl')

clearvars -except *PATH PATH* SUBJ SID BATCH dfreq sliWin
end

