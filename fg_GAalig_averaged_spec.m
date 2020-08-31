%% Entrainment gamma rFT
% PhD project 1
%

% This script:
% - loads in the results of the single subject TFR (conditons specified in settings)
% - averages over the flicker time interval and SOI per participant
% - computes the Grandaverage
% - plots the result using imagesc()

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

clear all; close all; clc

%% settings
oneof = 0;                      % 1/f corrected data?
sliwin = 0;                     % length of sliding window in seconds
condit = 1;                     % condition: 1 = flicker+gamma; 2: flicker
BATCH = 3;                      % sample

TOIBSL = [-0.75 -0.25];             % baseline interval
TOI = {[2.25 3.75],[0.25 1.75]};    % flicker interval
flickfreq = 52:2:90;                % flicker frequencies

% font size plots
ticksize = 16;
labelsize = 20;

% paths
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
MATPATH = 'C:\Users\katha\Documents\MATLAB';
addstr = {'entrainment','resonance'};
addpath(fullfile(MATPATH,'fieldtrip'))
addpath(fullfile(MAINPATH,'matlab','template function'))
ft_defaults



plotfile = ['av_spectra_',addstr{condit}];              % filename of plot


if oneof
    stroneof = 'oneof';
    plotfile = [plotfile,'_oneof'];
else
    stroneof = '';
end 
PATHIN = fullfile(MAINPATH,'results','power',addstr{condit},['Batch_',num2str(BATCH)],stroneof);
PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'power',addstr{condit},['Batch_',num2str(BATCH)]);
mkdir(PATHPLOT)

if sliwin > .5
    strlwin = 'long slidewin';
    plotfile = [plotfile,'_long slidewin'];

else
     strlwin = '';
end

PATHIN = fullfile(PATHIN, strlwin);
mkdir(PATHPLOT)

SUBJ = {};
% read in subjects & keep IGF>56
folds = dir(PATHIN);
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
%test: keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
exclSUBJ = SUBJ(find(IGF < 58));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;

% find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

% frequencies to be explored with respect to IGF
leftLim = minIGF - 52;
rightLim = 90-maxIGF;
freqvec = -leftLim:1/sliwin:rightLim;                   % frequency vector aligned to IGF


%% Grandaverage

TFRgam_ga = {};                      % help variables to store data for grandaverage
TFR_gam_ga_IGF = {};

for s = 1:length(SUBJ)
    SUBJPATH = fullfile(PATHIN, SUBJ{s});

    if condit == 1
        if oneof || sliwin >.5
            load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma_',num2str(sliwin),'.mat']))
        else
            load(fullfile(SUBJPATH, 'TFR_ent_alpha_gamma.mat'));
        end
    elseif condit == 2
        if oneof || sliwin >.5
            load(fullfile(SUBJPATH, ['TFR_res_alpha_RFT_',num2str(sliwin),'.mat']))
        else
            load(fullfile(SUBJPATH, ['TFR_res_alpha_RFT.mat']))
        end
    end

 
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'));   
    
    for r = 1:length(TFR) 
        cfg = [];
        cfg.channel = 'MEGGRAD';
        TFRgrad{r} = ft_selectdata(cfg,TFR{r});
        
        % combine planar
        cfg = [];
        cfg.method = 'sum';
        TFRgrad{r} = ft_combineplanar(cfg, TFRgrad{r});   
        
        % calculate relative power change
        cfg = [];
        cfg.baseline = [-0.75 -0.25];
        cfg.baselinetype = 'relchange';
        cfg.paramter = 'powspctrm';
        TFRcmb{r} = ft_freqbaseline(cfg, TFRgrad{r});
        
        % average over SOI and flicker interval
        cfg = [];
        cfg.channel = SOI';
        cfg.avgoverchan = 'yes';
        cfg.latency = TOI{condit};
        cfg.avgovertime = 'yes';
        TFRgamma{r} = ft_selectdata(cfg,TFRcmb{r});                 % all frequency conditions

    end
    clear TFRcmb TFRgrad
    
    %% frequencies of interest based on IGF
    freqvec = [IGF(s)-leftLim:1/sliwin:IGF(s)+rightLim];
    % select those frequency conditions
    TFR_IGF_gam = TFRgamma(find(ismember(flickfreq,freqvec)));
    
    % select those frequencies
    cfg = [];
    cfg.frequency = [freqvec(1) freqvec(end)];
    for r = 1:length(TFR_IGF_gam)
        TFR_IGF_gam{r} = ft_selectdata(cfg,TFR_IGF_gam{r});
        TFR_IGF_gam{r}.freq = -leftLim:1/sliwin:rightLim;
    end
    
    % all frequency conditions
    TFRgam_ga = [TFRgam_ga;TFRgamma];
    
    % only frequencies of interest (IGF-6 +16)
    TFR_gam_ga_IGF = [TFR_gam_ga_IGF;TFR_IGF_gam];
    close all;
    clearvars TFRgamma TFRalpha TFRcmb TFRgrad TFR_IGF_gam TFR_IGF_alpha TFR
end


cfg = [];
cfg.parameter = 'powspctrm';
for f = 1:size(TFRgam_ga,2)
    % replace with one sensor (to be able to compute GA)
    for s = 1:length(SUBJ)
        TFRgam_ga{s,f}.label = {'MEG1922+1923'};
     %   TFRalph_ga{s,f}.label = {'MEG1922+1923'};

    end
    TFRGAM{f} = ft_freqgrandaverage(cfg, TFRgam_ga{:,f});
   % TFRALPHA{f} = ft_freqgrandaverage(cfg, TFRalph_ga{:,f});
  
end
clear TFRgam_ga

cfg = [];
cfg.parameter = 'powspctrm';
for f = 1:size(TFR_gam_ga_IGF,2)
    % replace with one sensor (to be able to compute GA)
    for s = 1:length(SUBJ)
        TFR_gam_ga_IGF{s,f}.label = {'MEG1922+1923'};
     %   TFRalph_ga{s,f}.label = {'MEG1922+1923'};
         TFR_gam_ga_IGF{s,f}.freq = -leftLim:1/sliwin:rightLim;
    end
    TFRGAM_IGF{f} = ft_freqgrandaverage(cfg, TFR_gam_ga_IGF{:,f});
   % TFRALPHA{f} = ft_freqgrandaverage(cfg, TFRalph_ga{:,f});
  
end
clear TFR_gam_ga_IGF


% prepare to plot
for f = 1:length(flickfreq)
    gampow(:,f) = TFRGAM{f}.powspctrm';
   % alphpow(f,:) = TFRALPHA{f}.powspctrm;  
end

for f = 1:length(TFRGAM_IGF)
    igfpow(:,f) = TFRGAM_IGF{f}.powspctrm';
   % alphpow(f,:) = TFRALPHA{f}.powspctrm;  
end
freqvec = -leftLim:1/sliwin:rightLim;
clear TFRGAM TFRGAM_IGF
save(fullfile(PATHIN,'imagesc_pow_stimfreq.mat'),'gampow','igfpow','freqvec')



%% PLOT 
h = figure;
set(gcf, 'Position', [0, 0, 1920,465])
subplot(1,2,1)
imagesc(flickfreq,[40:1/sliwin:100],gampow)
axis xy
colormap jet
cb = colorbar;
cb.FontName = 'Arial';
cb.FontSize = ticksize;
cb.Label.String = {'relative power change'};
cb.Label.FontSize = labelsize;
caxis([-ceil(cb.Limits(2)) ceil(cb.Limits(2))])
if mod(ceil(cb.Limits(2)),2)
    cb.Ticks = [-ceil(cb.Limits(2)):ceil(cb.Limits(2))];
end
xlabel('RFT stimulation frequency [Hz]')
xticks([52:6:90])
ylabel({'Response frequency [Hz]'})
b = gca;
b.FontSize = ticksize;
b.FontName = 'Arial';
b.XLabel.FontSize = labelsize;
b.YLabel.FontSize = labelsize;
subplot(1,2,2)
imagesc(freqvec(~mod(freqvec,2)),freqvec,igfpow)
axis xy
colormap jet
cb2 = colorbar;
cb2.FontName = 'Arial';
cb2.FontSize = ticksize;
cb2.Label.String = {'relative power change'};
cb2.Label.FontSize = labelsize;
xlabel('RFT frequency IGF + [Hz]')
xticks([freqvec(2):4:freqvec(end)])
yticks([freqvec(1):4:freqvec(end)])
ylabel({'Response RFT IGF + [Hz]'})
b = gca;
b.FontSize = ticksize;
b.FontName = 'Arial';
b.XLabel.FontSize = labelsize;
b.YLabel.FontSize = labelsize;
caxis([-ceil(cb.Limits(2)) ceil(cb.Limits(2))])
% if mod(ceil(cb.Limits(2)),2)
% caxis([-ceil(cb.Limits(2)) ceil(cb.Limits(2))])
% end
% print(h, fullfile(PATHPLOT,plotfile),'-dpng','-r300');
print('-painters',h, fullfile(PATHPLOT,plotfile),'-dpng','-r600')


% cfg = [];
% cfg.parameter = 'powspctrm';
% for f = 1:size(TFR_gam_ga_IGF,2)
%     for s = 1:length(SUBJ)        
%         cfg = [];
%         cfg.frequency = [IGF(s)-leftLim, IGF(s)+rightLim];
%         TFR_gam_ga_IGF{s,f} = ft_selectdata(cfg,TFR_gam_ga_IGF{s,f});           % select frequencies of interest
%         TFR_gam_ga_IGF{s,f}.label = {'MEG1922+1923'};
%         TFR_gam_ga_IGF{s,f}.freq = TFR_gam_ga_IGF{s,f}.freq-IGF(s);             % subtract IGF (same frequency vector for all)
%         TFR_alpha_ga_IGF{s,f}.label = {'MEG1922+1923'};
%     end
%     TFRGAM_IGF{f} = ft_freqgrandaverage(cfg, TFR_gam_ga_IGF{:,f});
%     TFRALPHA_IGF{f} = ft_freqgrandaverage(cfg, TFR_alpha_ga_IGF{:,f});
% end
% clear TFR_gam_ga_IGF TFR_alpha_ga_IGF
% 
% for f = 1:length( TFRGAM_IGF)
%     gampow_IGF(f,:) =  TFRGAM_IGF{f}.powspctrm;
%     alphpow_IGF(f,:) = TFRALPHA_IGF{f}.powspctrm;  
% end
% freqvec_gam_IGF = TFRGAM_IGF{1}.freq;
% freqvec_alpha_IGF = TFRALPHA_IGF{1}.freq;
% clear TFRGAM_IGF TFRALPHA_IGF
% 
% freqvec = [-leftLim:2:rightLim];
% h = figure;
% set(gcf, 'Position', [0, 0, 1920, 1080])
% subplot(1,2,1)
% imagesc(freqvec_alpha_IGF,freqvec,alphpow_IGF)
% cb = colorbar;
% cb.FontName = 'Arial';
% cb.FontSize = 10;
% cb.Label.String = {'relative power change [a.u.]'};
% cb.Label.FontSize = 12;
% colormap jet
% %ylim([52 90])
% ylabel('RFT frequency condition in IGF + [Hz]')
% xlabel({'Alpha-band activity [Hz]'})
% caxis([-max(max(alphpow_IGF))-0.25 max(max(alphpow_IGF))+0.25])
% subplot(1,2,2)
% imagesc(freqvec_gam_IGF,freqvec,gampow_IGF)
% cb = colorbar;
% cb.FontName = 'Arial';
% cb.FontSize = 10;
% cb.Label.String = {'relative power change [a.u.]'};
% cb.Label.FontSize = 12;
% colormap jet
% ylabel('RFT frequency condition in IGF + [Hz]')
% xlabel({'Response RFT + Gamma in IGF + [Hz]'})
% caxis([-max(max(gampow_IGF))-0.25,max(max(gampow_IGF))+0.25])
% print(h, fullfile(PATHPLOT,['spectra_ent_gamma_alpha_IGF']),'-dpng','-r300');

