%% Entrainment gamma rFT
% PhD project 1
% 
% Grandaverage phase-locking value between photodiode and MEG signal as a function of frequency

% INPUTS: 
% - BATCH: sample (manus: 3)
% - condit: 1: flicker&gratings, 2: flicker
% - TOIBSL: baseline interval (recommended [-0.75 -0.25])
% - resptype: 'RFT': flicker response, 'IGF': power at IGF (shown in manus:
% RFT)
% - oneof: 1/f corrected data? (0,1)
% - sliWin: length of sliding window in seconds

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function i_freq_fun_plv(BATCH,condit,resptype,oneof,sliwin)

%% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;

addstr = {'entrainment','resonance'};               % conditions


flickfreq = [52:2:90];                      % stimulation frequency

TOI = {[2.25 3.75];[0.25 1.75]};                    % stimulation time windows

POWPATH = fullfile(MAINPATH, 'results', 'power');

if condit == 1
    ent = 1;
    res = 0;
elseif condit == 2
    res = 1;
    ent = 0;
end
SUBJ = {};
IGF = [];
SOI_all = {};
numSens = [];

GAMPATH = fullfile(POWPATH,'gammatron',['Batch_',num2str(BATCH)]);
PLVPATH = fullfile(MAINPATH,'results','phase_locking',['Batch_',num2str(BATCH)],addstr{condit});
PATHSTAT = fullfile(MAINPATH,'results', 'statistics','Linear regression',['Batch_',num2str(BATCH)]);

% not in manus 1/f correction and longer sliding windows
if oneof
    PATHSTAT = fullfile(PATHSTAT,'oneof');
    PLVPATH = fullfile(PLVPATH,'oneof');
end

if sliwin > .5
    PLVPATH = fullfile(PLVPATH,'long slidewin');
end

% read in subjects
folds = dir(PLVPATH);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

for s = 1:length(SUBJ)
    load(fullfile(GAMPATH, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end
 
%keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));
minIGF = min(IGF);
maxIGF = max(IGF);

leftLim = minIGF - 52;
rightLim = 20-leftLim-2;
SUBJ = keepSUBJ;
clear keepSUBJ maxIGF minIGF numSens numSens_B

PATHPLOT = fullfile(MAINPATH, 'results','plots', 'arnold tongue',['Batch_',num2str(BATCH)]);
PATHOUT = fullfile(MAINPATH,'arnold tongue',['Batch_',num2str(BATCH)]);
mkdir(PATHOUT)
mkdir(PATHSTAT)

mkdir(fullfile(PATHPLOT,'phase locking'));


% can't combine planars for coherence.. pull SOI apart
for h = 1:length(SOI_all)
    w = 1;              % help index
    for l = 1:length(SOI_all{h})
        SOInc{h}{w} = SOI_all{h}{l}(1:strfind(SOI_all{h}{l},'+')-1);         % first sensor
        w = w + 1;
        SOInc{h}{w} = ['MEG',SOI_all{h}{l}(strfind(SOI_all{h}{l},'+')+1:end)];
        w = w + 1;
    end
end


%% Load TFRs 
arnoldIGF = cell(1,length(SUBJ));                    % store TFRs with IGF aligned frequency vector
arnoldFREQ = cell(1,length(SUBJ));                   % store TFRs with normal frequ vector

for s = 1:length(SUBJ)
    PLV_FREQ = {};

    SUBJPATH = fullfile(PLVPATH,SUBJ{s});
    if oneof || sliwin > .5
        load(fullfile(SUBJPATH,['phase_lock_',num2str(sliwin),'.mat']))
    else
    load(fullfile(SUBJPATH,'phase_lock.mat'))
    end
    % all frequ conditions: relative power change
    for c = 1:length(flickfreq)

        % average over time: flicker on (0-2 seconds after stimulus onset)
        cfg = [];
        cfg.latency = TOI{condit};
        cfg.avgovertime = 'yes';
        PLV_FREQ{c} = ft_selectdata(cfg,PLV{c});

    end
    arnoldFREQ{s} = PLV_FREQ;
    
    % aligned frequency vector
    freqvec = IGF(s)-leftLim:(1/sliwin):IGF(s)+rightLim;
    freq_cond = find(ismember(flickfreq, freqvec));
    PLV_IGF = PLV_FREQ(freq_cond);                          % only select FOI
    for i = 1:length(PLV_IGF)       
        PLV_IGF{i}.plvspctrm = PLV_IGF{i}.plvspctrm(:,find(ismember(floor(PLV_IGF{i}.freq),freqvec)));
        PLV_IGF{i}.freq = freqvec;
    end
    arnoldIGF{s} = PLV_IGF;
    clear PLV PLV_FREQ PLV_IGF freqvec freq_cond
end

% average: IGF aligned
for s = 1:size(arnoldIGF,2)
    for g = 1:size(arnoldIGF{s},2) % frequencies

        % work around: add channel dimension
        arnoldIGF{s}{g}.label = arnoldIGF{s}{g}.labelcmb(:,1);
        arnoldIGF{s}{g}.dimord = 'chan_freq_time';
       % average over SOI
        cfg = [];
        cfg.channel = SOInc{s}';
        cfg.avgoverchan = 'yes';
        arnoldIGFav{s}{g} = ft_selectdata(cfg, arnoldIGF{s}{g});
    end
end
clear arnoldIGF

% average: flicker frequencies
arnoldFREQav = cell(1,length(SUBJ));
for s = 1:size(arnoldFREQ,2)
    for g = 1:size(arnoldFREQ{s},2) % frequencies

        arnoldFREQ{s}{g}.label = arnoldFREQ{s}{g}.labelcmb(:,1);
        arnoldFREQ{s}{g}.dimord = 'chan_freq_time';
        % average over SOI
        cfg = [];
        cfg.channel = SOInc{s}';
        cfg.avgoverchan = 'yes';
        cfg.frequency = [flickfreq(1) flickfreq(end)];
        arnoldFREQav{s}{g} = ft_selectdata(cfg, arnoldFREQ{s}{g});
        arnoldFREQav{s}{g}.freq = floor(arnoldFREQav{s}{g}.freq);
    end
end

clear arnoldFREQ 

%% Prepare for Linear Model in R

% store as table
for f = 1:size(arnoldFREQav{1},2)
    for s = 1:length(arnoldFREQav)
        FREQ{f}{s} = arnoldFREQav{s}{f};
        FREQ{f}{s}.label = {'MEG1922'};
        if strcmp(resptype,'RFT')
            M(s,f) = arnoldFREQav{s}{f}.plvspctrm(f);
        elseif strcmp(resptype,'IGF')
            M(s,f) = arnoldFREQav{s}{f}.plvspctrm(find(arnoldFREQav{s}{f}.freq == IGF(s)));
        end
    end   
    
    % Check whether Grandaverage and average of M are the same
%     cfg = [];
%     cfg.parameter = 'plvspctrm';
%     GAFREQ{f} = ft_freqgrandaverage(cfg,FREQ{f}{:});        
end
clear FREQ arnoldFREQav
T = array2table(M);
% D = array2table(DISTIGF);

if exist(fullfile(PATHSTAT,['plv_fun_',addstr{condit},'_',resptype,'.xlsx']))
    delete(fullfile(PATHSTAT,['plv_fun_',addstr{condit},'_',resptype,'.xlsx']))
end
writetable(T,fullfile(PATHSTAT,['plv_fun_',addstr{condit},'_',resptype,'_sliwin',num2str(sliwin),'.xlsx']),'Sheet',1);          % this will have to be converted to 
% Linear model IGF
for f = 1:size(arnoldIGFav{1},2)
    for s = 1:length(arnoldIGFav)
        %         FREQ{f}{s} = arnoldIGFav{s}{f};
        %         FREQ{f}{s}.label = {'MEG1922'};
        if strcmp(resptype,'RFT')
            M_IGF(s,f) = arnoldIGFav{s}{f}.plvspctrm(f);
        elseif strcmp(resptype,'IGF')
            M_IGF(s,f) = arnoldIGFav{s}{f}.plvspctrm(arnoldIGFav{s}{f}.freq == IGF(s));
        end
    end
    
end

T_IGF = array2table(M_IGF);
if exist(fullfile(PATHSTAT,['plv_fun_IGF_',addstr{condit},'_',resptype,'.xlsx']))
    delete(fullfile(PATHSTAT,['plv_fun_IGF_',addstr{condit},'_',resptype,'.xlsx']))
end
writetable(T_IGF,fullfile(PATHSTAT,['plv_fun_IGF_',addstr{condit},'_',resptype,'_sliwin',num2str(sliwin),'.xlsx']),'Sheet',1);
%clearvars -except PATH* *PATH IGF SUBJ M M_IGF flickfreq SOInc leftLim rightLim BATCH

%% save files

% average plv
plvspct = mean(M,1);
stdplvspct = std(M,1);

% repeat for distance to IGF
freqvec = -leftLim:2:rightLim;
plvspct_IGF = mean(M_IGF,1);
stdplvspct_IGF = std(M_IGF,1);

if exist(fullfile(PLVPATH,['Batch_',num2str(BATCH)','_freq_fun_',addstr{condit},'_',resptype,'.mat']))
    delete(fullfile(PLVPATH,['Batch_',num2str(BATCH)','_freq_fun_',addstr{condit},'_',resptype,'.mat']))
end
filename = ['Batch_',num2str(BATCH)','_',resptype];
if oneof
    filename = [filename,'oneof'];
end
if sliwin > .5
    filename = [filename,'long slidewin'];
end    
save(fullfile(PLVPATH,[filename,'.mat']), 'plvspct', 'stdplvspct','plvspct_IGF', 'stdplvspct_IGF','freqvec', 'M', 'M_IGF')
end