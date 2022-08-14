%% Entrainment gamma rFT
% PhD project 1

% Script 
% - loads in TFRs of flicker+gamma for all subjects kept for analysis
% - aligns to IGF & computes grand average

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

clear all; close all; clc;

%% settings

% axes label font size
ticksize = 16;
labelsize = 20;

% paths
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
MATPATH = 'C:\Users\katha\Documents\MATLAB';
addpath(fullfile(MATPATH,'fieldtrip'));
ft_defaults;
PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron','Batch_3');
PATHRES = fullfile(MAINPATH,'results', 'power','resonance','Batch_3');
PATHPLOTS = fullfile(MAINPATH,'results','plots', 'power','IGF', 'Batch_3');

% list subjects (not linux friendly)
cd(PATHGAM)
SUBJ = cellstr(ls());
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));


% exclude subjects whose IGF<58
for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

SUBJ = SUBJ(IGF>56);
IGF = IGF(IGF>56);

minIGF = min(IGF);
maxIGF = max(IGF);

leftLim = minIGF - 52;
rightLim = 20-leftLim-2;

freqvec = -leftLim:2:rightLim;
speclim = [30-minIGF 100-maxIGF];                   % frequency band relative to IGFs


%% TFR Grandaverage
TFRGAprep = cell(size(SUBJ));

for s = 1:length(SUBJ)
    % load in single subject data
    load(fullfile(PATHGAM,SUBJ{s},['TFR_gammatron.mat']))

    % combine planar
    cfg = [];
    cfg.channel = 'MEGGRAD';
    TFRgam = ft_selectdata(cfg, TFRgam);
    cfg = [];
    cfg.method = 'sum';
    TFRgamcmb = ft_combineplanar(cfg, TFRgam);

    clear TFRgam
    
    % compute relative change
    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.parameter = 'powspctrm';
    TFRgam = ft_freqbaseline(cfg,TFRgamcmb);
    
    % select frequencies of interest, SOI & store in cell
    subf = IGF(s)+speclim;
    cfg.frequency   = subf;
    cfg.channel     =  SOI_all{s};
    cfg.avgoverchan = 'yes';
    TFRGAprep{s}    = ft_selectdata(cfg,TFRgam);
    TFRGAprep{s}.freq = speclim(1):2:speclim(2);
    TFRGAprep{s}.label = {'MEG2032+2033'};
    clear TFRgam subf
end

% Grandaverage
cfg = [];
TFRGA = ft_freqgrandaverage(cfg,TFRGAprep{:});

%% plot (make sure same size as single subj gamma plots)

ncol = 2;
nrow = 2;
fig = figure;
set(gcf, 'Position', [0, 0, 1920/1.75, 1080/1.75])

% TFR
subplot(nrow,ncol,1)
imagesc(TFRGA.time,TFRGA.freq,squeeze(TFRGA.powspctrm))
axis xy
colormap jet
cb = colorbar;
cb.Limits = [-round(cb.Limits(2)) round(cb.Limits(2))];
caxis([-cb.Limits(2)+.5 cb.Limits(2)-.5])
cb.FontSize = ticksize;
cb.FontName = 'Arial';

cb.Label.String = {'relative power', 'change'};
cb.Label.FontName = 'Arial';
cb.Label.FontSize = labelsize;
cb.Label.Position = [2.0306 7.9076e-07 0];
xlabel ('time [s]')
ylabel('Frequency IGF+ [Hz]')
a = gca;
a.FontSize = ticksize;
a.FontName ='Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
yticks([TFRGA.freq(3):8:TFRGA.freq(end)])

% spectrum
subplot(nrow,ncol,2)
plot(TFRGA.freq,mean(squeeze(TFRGA.powspctrm),2),'LineWidth',1,'Color',[0 0.4470 0.7410])
xlabel('Frequency IGF+ [Hz]')
xlim([TFRGA.freq(1) TFRGA.freq(end)])
xticks([TFRGA.freq(3):8:TFRGA.freq(end)])
a = gca;
a.FontSize = ticksize;
a.FontName ='Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
ylim([0 a.YLim(2)])
yticks([0:0.2:0.8])

print(fig, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','GAMGA'),'-dsvg','-r600')
saveas(fig, fullfile('X:\Manuscript','GAMGA.svg'))
