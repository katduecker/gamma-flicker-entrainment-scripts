%% Entrainment gamma rFT
% Phd project 1
% 
% flicker condition plots (single subject): power as a function of flicker frequency

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%% settings
clear all; close all; clc;
BATCH = 3;
ticksize = 16;
labelsize = 20;
MAINPATH = 'X:\';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;
PATHIN = fullfile(MAINPATH,'results','power','resonance',['Batch_',num2str(BATCH)]);
PATHGAM = fullfile(MAINPATH, 'results','power', 'gammatron',['Batch_',num2str(BATCH)]);
PATHPLOTS = fullfile(MAINPATH, 'results','plots', 'power','resonance', ['Batch_',num2str(BATCH)]);mkdir(PATHPLOTS)

folds = dir(PATHIN);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;                                       % list subject folders
end
SUBJ(1:2) = [];
SUBJ(~strncmp(SUBJ,'201',3))=[];

for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

% exclude subjects with IGF < 56
SUBJ = SUBJ(IGF>56);
IGF = IGF(IGF>56);
flickfreq = 52:2:90;

PATHPLOTS = fullfile(MAINPATH, 'results','plots', 'power','resonance', ['Batch_',num2str(BATCH)],'relative change');mkdir(PATHPLOTS)

%% single subject plots
ncol = 2;
nrow = 2;
h = -1;
fig = figure;
set(gcf, 'Position', [0, 0, 1920/1.75, 1080/1.75])

for s = [1,3]
    h = h+2;
    SUBJPATH = fullfile(PATHIN, SUBJ{s});
    load(fullfile(SUBJPATH, 'TFR_res_alpha_RFT.mat'));
    %   load(fullfile(SUBJPATH, 'flickerResp.mat'));
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'));
    
    for t = 1:length(TFR)
        cfg = [];
        cfg.channel = 'MEGGRAD';
        
        TFRgrad{t} = ft_selectdata(cfg,TFR{t});
        % combine planar
        cfg = [];
        cfg.method = 'sum';
        TFRcmb{t} = ft_combineplanar(cfg, TFRgrad{t});
        %TFRrescmb{t}.dimord = 'chan_freq';
    end
    
    
    for r = 1:length(TFRcmb)
        
        cfg = [];
        cfg.baseline = [-0.75 -0.25];
        cfg.baselinetype = 'relchange';
        cfg.paramter = 'powspctrm';
        TFRcmb{r} = ft_freqbaseline(cfg, TFRcmb{r});
        
        cfg = [];
        cfg.channel = SOI';
        cfg.avgoverchan = 'yes';
        cfg.latency = [0.25 1.75];
        cfg.avgovertime = 'yes';
        TFRSOI{r} = ft_selectdata(cfg,TFRcmb{r});
    end
    
    subplot(nrow,ncol,h)
    for r = 1:length(TFRcmb)
        plot(TFRSOI{r}.freq,TFRSOI{r}.powspctrm);
        hold on
        xlim([50 100])
        xticks([10:10:100])
        xlabel('Frequency [Hz]')
        ylabel({'relative power', 'change'})
        a = gca;
        a.FontSize = ticksize;
        a.FontName ='Arial';
        a.XLabel.FontSize = labelsize;
        a.YLabel.FontSize = labelsize;
        
        
    end
    
    hold on
    if s == 1
        amp = [-1 2];
    elseif s == 3
        amp = [-2 6];
    end
    line([IGF(s) IGF(s)],amp,'Color', 'black', 'LineStyle', ':','LineWidth',1.5)
    hold on
    text([IGF(s)],TFRSOI{find(flickfreq==IGF(s))}.powspctrm(find(TFRSOI{r}.freq==IGF(s))),num2str(IGF(s)), 'FontSize', ticksize,'FontName', 'Arial');

    clearvars TFR* SOI*
end
print(fullfile('X:\Manuscript','spectra_alpha_RFT'),'-dsvg','-r600')

for s = 1%:length(SUBJ)
    SUBJPATH = fullfile(PATHIN, SUBJ{s});
    
    load(fullfile(SUBJPATH, 'TFR_res_alpha_RFT.mat'));
    %   load(fullfile(SUBJPATH, 'flickerResp.mat'));
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat')); % sensors of interest
    
    
    % Plot T-F-Spectra for Sensors of Interest
    for t = 1:length(TFR)
        % combine planar
        cfg = [];
        cfg.method = 'sum';
        TFRcmb{t} = ft_combineplanar(cfg, TFR{t});
    end
    
    f = figure;
    set(gcf, 'Position', [0, 0, 1920, 1080])
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.channel = SOI';
    cfg.colormap = 'jet';
    
    % zlims defined manually
    if s == 5
        cfg.zlim =[-0.8 0.8];
    elseif s == 7 || s == 10
        cfg.zlim = [-0.5 0.5];
    elseif s == 9
        cfg.zlim = [-2 2];
    elseif s== 1 || s == 12 || s == 16
        cfg.zlim = [-1.5 1.5];
    else
        cfg.zlim = [-1 1];
    end
    for t = 1:length(TFR)
        subplot(5,4,t)
        cfg.colorbar = 'no';
        ft_singleplotTFR(cfg,TFRcmb{t})
        if flickfreq(t) == gamFreq
            title([num2str(flickfreq(t)), ' Hz = IGF rFT, ', num2str(numTrl{t}), ' trials'],'FontSize', 10, 'FontWeight', 'normal', 'FontName', 'Arial');
        else
            title([num2str(flickfreq(t)), ' Hz rFT, ', num2str(numTrl{t}), ' trials'],'FontSize', 10, 'FontWeight', 'normal', 'FontName', 'Arial');
        end
        hold on
        % mark flicker freq
        line(TFRcmb{t}.time, repmat(flickfreq(t),1,length(TFRcmb{t}.time)), 'Color', 'black', 'LineStyle', '-.')
        % mark subharmonic
        line(TFRcmb{t}.time, repmat(flickfreq(t)/2,1,length(TFRcmb{t}.time)), 'Color', 'blue', 'LineStyle', ':')

        if t > 16
            xlabel('time [s]')
            hLabel = get(gca,'XLabel');
            set(hLabel, 'FontSize', 10, 'FontName', 'Arial');
        end

        if find(t == [1, 5, 9, 13, 17])
            ylabel('Frequency [Hz]')
            jLabel = get(gca,'YLabel');
            set(jLabel, 'FontSize',  10, 'FontName', 'Arial');
        end
        if find(t == [4, 8, 12, 16, 20])
            cfg.colorbar = 'yes';
            cb = colorbar;
            cb.Label.String = 'relative power change';
            cb.Label.String = 'relative power change';
            c.Label.FontSize = 9;
            c.Label.FontName = 'Arial';
        else
            cfg.colorbar = 'no';
        end
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'FontName','Arial','FontSize',9)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Arial','FontSize',9)
        text(-0.25,flickfreq(t)+5,num2str(flickfreq(t)), 'FontSize', 7)
        text(-0.25,flickfreq(t)/2+5,num2str(flickfreq(t)/2), 'FontSize', 7)
        %      text(-0.25,10+5,num2str(10))
    end
    saveas(f, fullfile(PATHPLOTS,[SUBJ{s},'reson_alpha_RFT.emf']))
    saveas(f, fullfile(PATHPLOTS,[SUBJ{s},'reson_alpha_RFT.png']))

    close all
    clear relChcmb relChRES TFR TFRcmb
end
