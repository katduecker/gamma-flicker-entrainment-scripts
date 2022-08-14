%% Entrainment gamma rFT
% Phd project 1
%
% Flicker condition: 
% compute power grandaverage separately for each flicker frequency
% plot result

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%% settings
clear all; close all; clc
BATCH = 3;
TOI = [0.5 1.5];                    % time window of interest
ticksize = 16;
labelsize = 20;
MAINPATH = 'X:\';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;

SUBJ = {};
IGF = [];
SOI_all = {};
numSens = [];

PATHIN = fullfile(MAINPATH, 'results','power','resonance',['Batch_',num2str(BATCH)]);

PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);

% read in subjects
folds = dir(PATHIN);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];


%% IGFs & SOIs

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
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;
%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

leftLim = minIGF - 52;
rightLim = 90 - maxIGF;

PATHPLOTS = fullfile(MAINPATH, 'results','plots', 'power','resonance',['Batch_',num2str(BATCH)]);
mkdir(PATHPLOTS)


%% Load in all subjects, baseline etc

for s = 1:length(SUBJ)
    
    SUBJPATH = fullfile(PATHIN,SUBJ{s});
    
    load(fullfile(SUBJPATH, 'TFR_res_alpha_RFT.mat'))
    for t = 1:length(TFR)
        
        % relative baseline
        cfg = [];
        cfg.baseline = [-0.75 -0.25];
        cfg.baselinetype = 'relchange';
        cfg.parameter = 'powspctrm';
        TFRRC = ft_freqbaseline(cfg, TFR{t});
        
        % select grad an combine
        cfg = [];
        cfg.channel = 'MEGGRAD';
        TFRg = ft_selectdata(cfg, TFRRC);
        cfg = [];
        cfg.method = 'sum';
        TFRcmb = ft_combineplanar(cfg, TFRg);
        clear TFRg numTrl
        
        % average over SOI for subject
        cfg = [];
        cfg.channel = SOI_all{s}';
        cfg.avgoverchan = 'yes';
        TFRSOI = ft_selectdata(cfg, TFRcmb);
        
        cfg.latency = [0.25 1.75];
        cfg.avgovertime = 'yes';
        SPEC = ft_selectdata(cfg,TFRcmb);
        clear TFRcmb
        TFRGAprep{t,s} = TFRSOI;                % TFR
        SPECGAprep{t,s} = SPEC;                 % averaged TFR -> power spectrum
        clear TFRSOI SPEC
    end
end

%% Grandaverage
cfg = [];
cfg.parameter = 'powspctrm';
for f = 1:size(SPECGAprep,1)
    % replace with one sensor (for averaging)
    for s = 1:length(SUBJ)
        TFRGAprep{f,s}.label = {'MEG1922+1923'};
        SPECGAprep{f,s}.label = {'MEG1922+1923'};
    end
    TFRGA{f} = ft_freqgrandaverage(cfg, TFRGAprep{f,:});
    SPECGA{f} = ft_freqgrandaverage(cfg, SPECGAprep{f,:});
end
clear SPECGAprep TFR TFRRC

%% PLOTS

% plot power spectra
ncol = 2;
nrow = 2;
fig = figure;
set(gcf, 'Position', [0, 0, 1920/1.75, 1080/1.75])
subplot(nrow,ncol,1)
for f = 1:length(SPECGA)
    plot(SPECGA{f}.freq,SPECGA{f}.powspctrm);
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
print(fullfile('X:\Manuscript','spectra_RFT_GA'),'-dsvg','-r600')

%% Plot TFR (not shown in manus)
% uncomment to plot with colorbar
flickfreq = 52:2:90;
ncol = 5;
nrow = 4;
h = figure;
set(gcf, 'Position', [0, 0, 1920, 1080])
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colormap = 'jet';
cfg.channel = TFRGA{1}.label;
cfg.zlim = [-1 1];
cfg.xlim = [-0.75 3.5];
cfg.ylim = [4 100];
cfg.title = [''];

for t = 1:length(TFRGA)
    subplot(nrow,ncol,t)
    cfg.colorbar = 'no';

    ft_singleplotTFR(cfg,TFRGA{t})
    %      if find(t == [ncol:ncol:length(TFRGA)])
    %
    %         cb = colorbar;        %     cb.Ticks = [floor(cb.Limits(1)):0.5:floor(cb.Limits(2))];
    %         cb.Label.String = {'relative power'; 'change [a.u.]'};
    %         cb.Label.FontSize = 12;
    % %
    %      end
    title(' ')
    hold on
    % mark flicker freq
    line(TFRGA{t}.time(find(TFRGA{t}.time==cfg.xlim(1)):find(TFRGA{t}.time==cfg.xlim(2))), repmat(flickfreq(t),1,length(TFRGA{t}.time(find(TFRGA{t}.time==cfg.xlim(1)):find(TFRGA{t}.time==cfg.xlim(2))))), 'Color', 'white', 'LineStyle', '-.', 'LineWidth', 1.5)
    % mark subharmonic
    line(TFRGA{t}.time(find(TFRGA{t}.time==cfg.xlim(1)):find(TFRGA{t}.time==cfg.xlim(2))), repmat(flickfreq(t)/2,1,length(TFRGA{t}.time(find(TFRGA{t}.time==cfg.xlim(1)):find(TFRGA{t}.time==cfg.xlim(2))))), 'Color', 'white', 'LineStyle', ':','LineWidth', 1.5)
    text(-0.75,flickfreq(t)+8,num2str(flickfreq(t)), 'FontSize', 12, 'FontName','Arial')
    text(-0.75,flickfreq(t)/2+8,num2str(flickfreq(t)/2), 'FontSize', 12, 'FontName','Arial')
    % mark 10 Hz
    %      line(TFRcmb{t}.time, repmat(10,1,length(TFRcmb{t}.time)), 'Color', 'blue', 'LineStyle', '--')
    % set x axis for lower row
    yticks([])
    xticks([])
    % cb.Ticks = [floor(cb.Limits(1)):0.5:floor(cb.Limits(2))];
    %     cb.Label.String = 'power change';
    %     cb.Label.FontName = 'Arial';
    %     cb.Label.FontSize = 12;
    %     cb.TickLabelsMode = 'manual';
    
    
    if t > length(TFRGA)-ncol
        xlabel('time [s]');
        xticks([-0.5:1:3.5])
        b = gca;
        b.FontSize = 12;
        b.XLabel.FontSize = 14;
        b.XAxis.FontName = 'Arial';
        
        
        %         cb = colorbar;
        %       %  cb.Ticks = [floor(cb.Limits(1)):0.5:floor(cb.Limits(2))];
        %         cb.Label.String = 'power change';
        %         cb.Label.FontName = 'Arial';
        %         cb.Label.FontSize = 12;
        %         cb.TickLabelsMode = 'manual';
        %
        %         cb.Label.Position = [2.0306 7.9076e-07 0];
    end
    %     a = get(gca,'XTickLabel');
    %     set(gca,'XTickLabel',a,'FontName','Arial','FontSize',8)
    if find(t == [1:ncol:length(TFRGA)])
        ylabel('Frequency [Hz]');
        yticks([10:15:100])
        a = gca;
        a.FontSize = 12;
        a.FontName = 'Arial';
        a.YLabel.FontSize = 14;
        a.XLabel.FontSize = 14;
        
    end
    
    %     if find(t == [ncol:ncol:length(TFRGA)])
    %         cb = colorbar;
    %         cb.FontSize = 10;
    %         cb.FontName = 'Arial';
    %         %  cb.Ticks = [floor(cb.Limits(1)):0.5:floor(cb.Limits(2))];
    %
    %         cb.Label.String = {'relative power change'};
    %
    %         cb.Label.FontSize = 12;
    %         cb.TickLabelsMode = 'manual';
    %     end
    
    
    %      text(-0.25,10+5,num2str(10))
end

print(h, fullfile(PATHPLOTS,['TFRGA_res_Batch_',num2str(BATCH)]),'-dsvg','-r1000');
print(h, fullfile(PATHPLOTS,['TFRGA_res_Batch_',num2str(BATCH)']),'-dpng','-r1000');
clearvars TFR*
