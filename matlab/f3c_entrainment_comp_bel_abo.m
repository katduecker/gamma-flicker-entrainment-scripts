%% Entrainment gamma rFT
% PhD project 1
% 
% Script prepares TFR data for statistically testing power at IGF during
% flicker<IGF vs flicker>IGF (statistics done in R)

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

clear all; close all; clc;

%% settings
BATCH = 3;                          % sample (3 = actual data) 2 & 1: PILOTS
relCh = 1;                          % relative change? (1 or 0 for power)
oneof = 0;                          % 1/f corrected data?
sliwin = .5;                        % length of sliding window in seconds

prep_anova = 1;                     % prepare data for ANOVA (compare relative power change at each time point, for flicker>IGF and flicker<IGF)
prep_t_test = 0;                    % prepare data for t-test > compare relative change T1 to T2 in both conditions
highf_bl = 0;                       % use high frequencies as baseline for t-test

TOI = [2.5 3.5];                    % time window of interest
flickfreq = 52:2:90;                % flicker frequencies

% paths
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
MATPATH = 'C:\Users\katha\Documents\MATLAB';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MATPATH,'fieldtrip'));
ft_defaults;

PATHCOND = fullfile(MAINPATH, 'results','preprocessing','conditions',['Batch_',num2str(BATCH)]);
PATHENT = fullfile(MAINPATH,'results','power','entrainment',['Batch_',num2str(BATCH)]);

PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'power','entrainment',['Batch_',num2str(BATCH)]);
mkdir(PATHPLOT)

if prep_anova
    PATHOUT = fullfile(MAINPATH,'results','statistics','ANOVA',['Batch_',num2str(BATCH)]);
elseif prep_t_test || highf_bl
    PATHOUT = fullfile(MAINPATH,'results','statistics','T-Test',['Batch_',num2str(BATCH)]);
end

if oneof
    PATHENT = fullfile(PATHENT,'oneof');
    PATHPLOT = fullfile(PATHPLOT,'oneof');
    PATHOUT = fullfile(PATHOUT,'oneof');
end
if sliwin > .5
    PATHENT = fullfile(PATHENT,'long slidewin');
    PATHPLOT = fullfile(PATHPLOT,'long slidewin');
    PATHOUT = fullfile(PATHOUT,'long slidewin');
end

mkdir(PATHPLOT)
mkdir(PATHOUT)
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'power','entrainment',['Batch_',num2str(BATCH)]);
mkdir(PATHPLOT)

% read in subjects & keep those with IGF>56
folds = dir(PATHENT);
for f = 1:length(folds)
    SUBJ = [SUBJ;folds(f).name];
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

% find min and max IGF
minIGF = min(IGF(IGF > 0));
maxIGF = max(IGF);

%test: keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
exclSUBJ = SUBJ(find(IGF > 56));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;

%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

% flicker frequencies to be considered
leftLim = minIGF - 52;
rightLim = 90-maxIGF;

%% Prepare for statistics in R

for s = 1:length(SUBJ)
    
    SUBJPATH = fullfile(PATHENT,SUBJ{s});
    
    % load in TFR results
    if oneof || sliwin > .5
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma_',num2str(sliwin),'.mat']))
    else
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma.mat']))
    end

    % combine planar
    cfg = [];
    cfg.channel = 'MEGGRAD';
    for f = 1:length(TFR)
        TFR{f} = ft_selectdata(cfg, TFR{f});
    end
    
    cfg = [];
    cfg.method = 'sum';
    for f = 1:length(TFR)
        TFR{f} = ft_combineplanar(cfg, TFR{f});
    end
      
    if relCh && ~highf_bl
        cfg = [];
        if prep_t_test                  % t-test: T1 is baseline
            cfg.baseline = [0.75 1.75];
        elseif prep_anova               % ANOVA: pre-stimulus baseline
            cfg.baseline = [-0.75 -0.25];
        end
        cfg.baselinetype = 'relchange';
        cfg.parameter = 'powspctrm';
        for t = 1:length(TFR)
            TFR{t} = ft_freqbaseline(cfg, TFR{t});
        end
    else
        error('Set relCh = 0 when using high frequencies as baseline.')
    end
    
    % extract frequency conditons for below (-6 and -4), above (4 and 6)
    % and high (highest 2 flicker frequencies)
    TFRbel = {TFR{find(flickfreq-IGF(s) == -6)},TFR{find(flickfreq-IGF(s) == -4)}}; % below: -6, -4
    TFRabo = {TFR{find(flickfreq-IGF(s) == 6)},TFR{find(flickfreq-IGF(s) == 4)}};   % above: 6, 4
     TFRout = {TFR{find(flickfreq-IGF(s) == rightLim)},TFR{find(flickfreq-IGF(s) == rightLim-2)}}; % highest frequencies

    clear TFR
    
    % average over channel & time
    cfg = [];
    cfg.channel = SOI_all{s}';
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
        
    if prep_anova
        cfg.latency = [0.5 1.5];
        % average over 1-sec interval before stim
        for f = 1:length(TFRbel)
            belSOIT1{f}  = ft_selectdata(cfg, TFRbel{f});
            aboSOIT1{f}  = ft_selectdata(cfg, TFRabo{f});
            outSOIT1{f} = ft_selectdata(cfg, TFRout{f});

        end
    end
    
    % average over 1-second interval during flicker
    cfg.latency = TOI;
    for f = 1:length(TFRbel)
        if prep_anova
            belSOIT3{f} = ft_selectdata(cfg, TFRbel{f});
            aboSOIT3{f} = ft_selectdata(cfg, TFRabo{f});

            outSOIT3{f} = ft_selectdata(cfg, TFRout{f});
            
        elseif prep_t_test || highf_bl
            belSOIRC{f} = ft_selectdata(cfg, TFRbel{f});
            aboSOIRC{f} = ft_selectdata(cfg, TFRabo{f});

            outSOIRC{f} = ft_selectdata(cfg, TFRout{f});

        end
    end
    clear TFR*
    
    for f = 1:length(belSOIRC)
        belSOIRC{f}.powspctrm = (belSOIRC{f}.powspctrm./outSOIRC{f}.powspctrm)-1;
        aboSOIRC{f}.powspctrm = (aboSOIRC{f}.powspctrm./outSOIRC{f}.powspctrm)-1;
    end
    
    
    
    % average over frequency conditions
    cfg = [];
    cfg.keepindividual = 'no';
    cfg.foilim = [IGF(s) IGF(s)];
    if prep_anova
        belT1 = ft_freqgrandaverage(cfg, belSOIT1{:});
        aboT1 = ft_freqgrandaverage(cfg, aboSOIT1{:});
        outT1 = ft_freqgrandaverage(cfg, outSOIT1{:});
        belT3 = ft_freqgrandaverage(cfg, belSOIT3{:});
        aboT3 = ft_freqgrandaverage(cfg, aboSOIT3{:});
        outT3 = ft_freqgrandaverage(cfg, outSOIT3{:});
        % make matrix: below: T1, T3; above: T1, T3; highest frequencies
            M(s,:) = [belT1.powspctrm belT3.powspctrm aboT1.powspctrm aboT3.powspctrm outT1.powspctrm outT3.powspctrm];


    elseif prep_t_test || highf_bl
        belRC = ft_freqgrandaverage(cfg, belSOIRC{:});
        aboRC = ft_freqgrandaverage(cfg, aboSOIRC{:});
        outRC = ft_freqgrandaverage(cfg, outSOIRC{:}); 
        M(s,:) = [belRC.powspctrm aboRC.powspctrm outRC.powspctrm];

    end
  

    clearvars bel* abo* TFR*
end


if ~relCh
    if prep_anova
        if exist(fullfile(PATHOUT,'TFRT1vsT3_int.xlsx'))
            delete(fullfile(PATHOUT,'TFRT1vsT3_int.xlsx'))
        end
        T = table(M(:,1), M(:,2), M(:,3),M(:,4),M(:,5),M(:,6),'VariableNames',{'belT1','belT3','aboT1','aboT3','outT1','outT3'});
        writetable(T,fullfile(PATHOUT,'TFRT1vsT3_int_raw.xlsx'))
        
    elseif prep_t_test
        if exist(fullfile(PATHOUT,'TFRChangebelvsabovsout.xlsx'))
            delete(fullfile(PATHOUT,'TFRChangebelvsabovsout.xlsx'))
        end
        T = table(M(:,1), M(:,2), M(:,3),'VariableNames',{'belRC','aboRC','outRC'});
        writetable(T,fullfile(PATHOUT,'TFRChangebelvsabovsout_raw.xlsx'))
    end
else
    % Make table
    if prep_anova
        if exist(fullfile(PATHOUT,'TFRT1vsT3_int.xlsx'))
            delete(fullfile(PATHOUT,'TFRT1vsT3_int.xlsx'))
        end
        T = table(M(:,1), M(:,2), M(:,3),M(:,4),M(:,5),M(:,6),'VariableNames',{'belT1','belT3','aboT1','aboT3','outT1','outT3'});
        writetable(T,fullfile(PATHOUT,'TFRT1vsT3_int.xlsx'))
        
    elseif prep_t_test
        if exist(fullfile(PATHOUT,'TFRChangebelvsabovsout.xlsx'))
            delete(fullfile(PATHOUT,'TFRChangebelvsabovsout.xlsx'))
        end
        T = table(M(:,1), M(:,2), M(:,3),'VariableNames',{'belRC','aboRC','outRC'});
        writetable(T,fullfile(PATHOUT,'TFRChangebelvsabovsout.xlsx'))
    elseif highf_bl
        if exist(fullfile(PATHOUT,'TFRrc_highf_bsl.xlsx'))
            delete(fullfile(PATHOUT,'TFRrc_highf_bsl.xlsx'))
        end
        T = table(M(:,1), M(:,2), M(:,3),'VariableNames',{'belRC','aboRC','outRC'});
        writetable(T,fullfile(PATHOUT,'TFRrc_highf_bsl.xlsx'))
    end
end
%save(fullfile(MAINPATH,'powerchanpointT1T3.mat'),'M')