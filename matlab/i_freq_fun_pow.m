%% Entrainment gamma rFT
% PhD project
% 
% Grandaverage power of response to RFT as a function of frequency

% INPUTS: 
% - BATCH: sample
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

% clear all; close all; clc;
function i_freq_fun_pow(BATCH,condit,TOIBSL,resptype,oneof,sliwin)

% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;


POWPATH = fullfile(MAINPATH, 'results', 'power');
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'arnold tongue',['Batch_',num2str(BATCH)]);

addstr = {'entrainment','resonance'};               % condititions
PATHOUT = fullfile(POWPATH,addstr{condit}); mkdir(PATHOUT)
addstrcondi = addstr{condit};                                       % seems redundant but important!!
PATHSTAT = fullfile(MAINPATH,'results', 'statistics','Linear regression',['Batch_',num2str(BATCH)]);
mkdir(PATHSTAT)

mkdir(fullfile(PATHPLOT,'power'));

PATHOUT = fullfile(POWPATH,addstrcondi);
mkdir(PATHOUT)

flickfreq = [52:2:90];                      % stimulation frequencies

allTOI = {[2.25 3.75],[0.25 1.75]};         % stimulation time windows
TOI = allTOI{condit};                       % time window of interest

SUBJ = {};
IGF = [];
SOI_all = {};
numSens = [];

GAMPATH = fullfile(POWPATH,'gammatron',['Batch_',num2str(BATCH)]);
PATHIN = fullfile(POWPATH,addstr{condit},['Batch_',num2str(BATCH)]);

% not shown in manus: 1/f correction and longer sliding windows
if oneof
    PATHIN = fullfile(PATHIN,'oneof');
end

if sliwin > .5
    PATHIN = fullfile(PATHIN,'long slidewin');
end

% read in subjects
folds = dir(PATHIN);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];
    

% read in IGF
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

%% Adjust frequency vector: IGF +- 

% loop: average over time and select frequencies
arnoldIGF = cell(1,length(SUBJ));                    % store TFRs with IGF aligned frequency vector
arnoldFREQav = cell(1,length(SUBJ));                   % store TFRs with normal frequ vector
h = 0;

for s = 1:length(SUBJ)
    
    SUBJPATH = fullfile(PATHIN,SUBJ{s});
    
    TFR_FREQ = {};
    % load TFRs
    if condit == 1
        if oneof || sliwin > .5
            load(fullfile(SUBJPATH,['TFR_',addstr{condit}(1:3),'_alpha_gamma_',num2str(sliwin),'.mat']))
        else
            load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma.mat']))
        end
    elseif condit == 2
        if oneof || sliwin > .5
            load(fullfile(SUBJPATH,['TFR_',addstr{condit}(1:3),'_alpha_RFT_',num2str(sliwin),'.mat']))
        else
            load(fullfile(PATHIN,SUBJ{s},['TFR_',addstr{condit}(1:3),'_alpha_RFT.mat']))
        end
    end
    
    % relative power change per flicker frequency
    for c = 1:length(flickfreq)
        % relative baseline
        cfg = [];
        cfg.baseline = TOIBSL;%[-0.7 -0.2];
        cfg.baselinetype = 'relchange';
        TFR_FREQ{c} = ft_freqbaseline(cfg,TFR{c});      
   
         % average over TOI & SOI
        cfg = [];  
        cfg.latency = TOI;
        disp(TOI)
        cfg.avgovertime = 'yes';
        cfg.channel = SOInc{s}';
        cfg.avgoverchan = 'yes';
        TFR_FREQav{c} = ft_selectdata(cfg, TFR_FREQ{c});
        
        clear arniGRI        
        
        % RFT = flicker
        % select power at current flicker frequency
        if strcmp(resptype,'RFT')
            TFR_FREQav{c}.powspctrm = TFR_FREQav{c}.powspctrm(:,find(ismember(floor(TFR_FREQav{c}.freq),flickfreq(c))));
            TFR_FREQav{c}.freq = flickfreq(c);
        elseif strcmp(resptype,'IGF')
            % select power at IGF during flicker
            TFR_FREQav{c}.powspctrm = TFR_FREQav{c}.powspctrm(:,find(ismember(floor(TFR_FREQav{c}.freq),IGF(s))));
            TFR_FREQav{c}.freq = 0;
        end

    end
    
    % aligned frequency vector (to IGF)
    freqvec = IGF(s)-leftLim:(1/sliwin):IGF(s)+rightLim;
    freq_condit = find(ismember(flickfreq, freqvec));
    TFR_IGF = TFR_FREQav(freq_condit);                          % only select FOI
% 
    for c = 1:length(TFR_IGF)
        TFR_IGF{c}.freq = freqvec(c)-IGF(s);
    end
    arnoldIGF{s} = TFR_IGF;
    arnoldFREQav{s} = TFR_FREQav;
    clearvars TFR* freqvec freq_condit
end


%% Linear Model in R

% store power values in table and save as excel
for f = 1:size(arnoldFREQav{1},2)
    for s = 1:length(arnoldFREQav)
        M(s,f) = arnoldFREQav{s}{f}.powspctrm;
    end         

end
clear arnoldFREGav 
T = array2table(M);
if exist(fullfile(PATHSTAT,['pow_fun_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']))
    delete(fullfile(PATHSTAT,['pow_fun_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']));
end
writetable(T,fullfile(PATHSTAT,['pow_fun_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']),'Sheet',1);
 

% Linear Model IGF-leftLim+rightLim only
for f = 1:size(arnoldIGF{1},2)
    for s = 1:length(arnoldIGF)
        M_IGF(s,f) = arnoldIGF{s}{f}.powspctrm;
    end   
       
end
clear arnoldIGF

T = array2table(M_IGF);
if exist(fullfile(PATHSTAT,['pow_fun_IGF_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']))
    delete(fullfile(PATHSTAT,['pow_fun_IGF_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']));
end
writetable(T,fullfile(PATHSTAT,['pow_fun_IGF_',addstrcondi ,'_',num2str(TOIBSL(1)),'_',resptype,'.xlsx']),'Sheet',1);

% average plv at current frequency
powspct = mean(M,1);                % grandaverage per frequency
stdpowspct = std(M,1);

%% IGF +- leftLim rightLIm

powspct_IGF = mean(M_IGF,1);        % grandaverage per IGF +-
stdpowspct_IGF = std(M_IGF,1);
freqvec = -leftLim:2:rightLim;


%% save files
if exist(fullfile(POWPATH,['Batch_',num2str(BATCH)','_freq_fun_',addstrcondi ,'_',resptype,'.mat']))
    delete(fullfile(POWPATH,['Batch_',num2str(BATCH)','_freq_fun_',addstrcondi ,'_',resptype,'.mat']))
end

filename = ['Batch_',num2str(BATCH)','_freq_fun_basel_',num2str(TOIBSL(1)),'_',resptype];
if oneof
    filename = [filename,'oneof'];
end
if sliwin > .5
    filename = [filename,'long slidewin'];
end    
save(fullfile(PATHOUT,[filename,'.mat']), 'powspct', 'stdpowspct','powspct_IGF', 'stdpowspct_IGF','freqvec','M','M_IGF')
end
