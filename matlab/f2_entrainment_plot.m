%% Entrainment gamma rFT
% PhD project 1
%
%
% flicker+gammma condition: plot TFRs for each participant

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

clear all; close all; clc;

%% settings
BATCH = 3;                  % sample

oneof = 1;                  % 1/f corrected data?
sliwin = 1;                 % sliding window length
relCh = 1;                  % relative baseline or not?


% font sizes for plots
ticksize = 12;
labelsize = 15;

do_cb = 0;                  % colorbar?

% paths
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH, 'matlab','template function'))
addpath(fullfile('C:\Users\katha\Documents\MATLAB','fieldtrip'));

ft_defaults;
PATHIN = fullfile(MAINPATH,'results','power','entrainment',['Batch_',num2str(BATCH)]);
GAMPATH = fullfile(MAINPATH, 'results','power', 'gammatron',['Batch_',num2str(BATCH)]);
PATHPLOTS = fullfile(MAINPATH, 'results','plots','power', 'entrainment', ['Batch_',num2str(BATCH)]);

if relCh
    PATHPLOTS = fullfile(PATHPLOTS,'relative change');
else
    PATHPLOTS = fullfile(PATHPLOTS,'spectra');
end
if oneof
    PATHIN = fullfile(PATHIN,'oneof');
    PATHPLOTS = fullfile(PATHPLOTS,'oneof');

end
if sliwin > .5
    PATHIN = fullfile(PATHIN,'long slidewin');
    PATHPLOTS = fullfile( PATHPLOTS,'long slidewin');

end
mkdir(PATHPLOTS)

% read in subjects
folds = dir(PATHIN);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];
flickfreq = 52:2:90;
for s = 1:length(SUBJ)
    load(fullfile(GAMPATH, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

% find min and max IGF
minIGF = min(IGF(IGF > 0));
maxIGF = max(IGF);

leftLim = minIGF - 52;
rightLim = 20-leftLim-2;

% keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));

SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;

%% Exploratory: plot spectra
% average over TFRs: flicker+gamma interval

for s = 1:7
    SUBJPATH = fullfile(PATHIN, SUBJ{s});

    
    % load results of f_entrainment.m
    if oneof || sliwin > .5
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma_',num2str(sliwin),'.mat']))
    else
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma.mat']))
    end
 
    load(fullfile(GAMPATH, SUBJ{s}, 'SOI_freq.mat'));                       % load identified SOI
    
    for t = 1:length(TFR)
        % select gradiometers
        cfg = [];
        cfg.channel = 'MEGGRAD';
        TFRgrad{t} = ft_selectdata(cfg,TFR{t});
        
        % combine planar
        cfg = [];
        cfg.method = 'sum';
        TFRcmb{t} = ft_combineplanar(cfg, TFRgrad{t});
    end
    
    
    for r = 1:length(TFRcmb)
        
        % relative power change
        if relCh
            cfg = [];
            cfg.baseline = [-0.75 -0.25];
            cfg.baselinetype = 'relchange';
            cfg.paramter = 'powspctrm';
            TFRcmb{r} = ft_freqbaseline(cfg, TFRcmb{r});
        end
        
        % average over SOI and flicker+gamma interval
        cfg = [];
        cfg.channel = SOI';
        cfg.avgoverchan = 'yes';
        cfg.latency = [2.5 3.5];
        cfg.avgovertime = 'yes';
        TFRSOI{r} = ft_selectdata(cfg,TFRcmb{r});
    end
    
    % plot
   f =  figure;
    for r = 1:length(TFRcmb)
        plot(TFRSOI{r}.freq,TFRSOI{r}.powspctrm);
        hold on
        xlim([4 100])
        xlabel('Frequency [Hz]')
        ylabel('Relative power change [a.u.]')
       % title(['IGF ',num2str(gamFreq)]);
    end
   
    hold on 
    text([gamFreq],TFRSOI{find(flickfreq==gamFreq)}.powspctrm(find(round(TFRSOI{r}.freq)==gamFreq)),['IGF  = ',num2str(gamFreq)], 'FontSize', 8,'FontName', 'Arial');
    if relCh
        print(fullfile(PATHPLOTS,[SUBJ{s},'_RC_entrain_alpha_gamma_',num2str(sliwin)]),'-dpng','-r300')
    else
        print(fullfile(PATHPLOTS,[SUBJ{s},'_spectra_entrain_alpha_gamma_',num2str(sliwin)]),'-dpng','-r300')
    end
    close all

    close all;
    clearvars TFR* SOI*
end

%% PLOT TFR

% layout of subplots
nrow = 4;
ncol = 5;
for s = 1
   SUBJPATH = fullfile(PATHIN, SUBJ{s});

  % load results of f_entrainment.m

 if oneof || sliwin > .5
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma_',num2str(sliwin),'.mat']))
    else
        load(fullfile(SUBJPATH, ['TFR_ent_alpha_gamma.mat']))
 end
 load(fullfile(GAMPATH,SUBJ{s}, 'SOI_freq.mat'))                % load in SOIs
    % flicker frequencies
    flickfreq = 52:2:90;
    
    % combine planar
    for t = 1:length(TFR)
        % select gradiometers
        cfg = [];
        cfg.channel = 'MEGGRAD';
        TFR{t} = ft_selectdata(cfg, TFR{t});
        
        % combine
        cfg = [];
        cfg.method = 'sum';
        TFRcmb{t} = ft_combineplanar(cfg, TFR{t});
    end
    clear TFR
    
    
    % subplot TFR
    cfg = [];
    cfg.parameter = 'powspctrm';
    % relative power change
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.channel = SOI';
    cfg.colormap = 'jet';
    cfg.ylim = [30 100];
    
    % zlims identified by hand
    if s == 1 || s==26
         cfg.zlim =[-3 3];
     elseif s == 12 
         cfg.zlim = [-6 6];

     else
       cfg.zlim = 'maxmin';
    end

    % subplot
    f = figure;
    set(gcf, 'Position', [0, 0, 1920, 1080])
    % subplot layout identified by hand
    [ha pos] = tight_subplot(nrow,ncol,[.05 .025],[.125 .075],[.1 .08]);
    for t = 1:length(TFRcmb)

        cfg.colorbar = 'no';
        axes(ha(t))
        ft_singleplotTFR(cfg,TFRcmb{t})
        yticks([30:10:100])
        
        % title: number of trials
        if flickfreq(t) == gamFreq
            title(['IGF ', num2str(numTrl{t}), ' trials'],'FontSize', ticksize, 'FontWeight','normal','FontName', 'Arial');
        else
            title([num2str(numTrl{t}), ' trials'],'FontSize', ticksize, 'FontWeight','normal','FontName', 'Arial');
        end
        hold on
        
        % mark flicker freq
        line(TFRcmb{t}.time(TFRcmb{t}.time>=2 & TFRcmb{t}.time<=4), repmat(flickfreq(t),1,length(TFRcmb{t}.time(TFRcmb{t}.time>=2 & TFRcmb{t}.time<=4))), 'Color', 'black', 'LineStyle', '-.','LineWidth', 1)
        % mark gamma freq
        line(TFRcmb{t}.time(TFRcmb{t}.time>=0), repmat(gamFreq,1,length(TFRcmb{t}.time(TFRcmb{t}.time>=0))), 'Color', 'black', 'LineStyle', '--','LineWidth', 1)

        if flickfreq(t) > gamFreq
            text(1.5,flickfreq(t)+4,num2str(flickfreq(t)), 'FontSize', ticksize)
        elseif flickfreq(t) <= gamFreq
            text(1.5,flickfreq(t)-4,num2str(flickfreq(t)), 'FontSize', ticksize)
        end
        
        yticks([])
        xticks([])
        
        % x-label & xticks only in bottom row
        if t > length(TFRcmb) - ncol
            xlabel('time [s]')
            xticks([-0:1:5])
            set(ha(t),'XTickLabel',string([0:1:5]))
        end
        
        % y-label & yticks only in first column
        if find(t == [1:ncol:length(TFRcmb)])
            ylabel('Frequency [Hz]')
            set(ha(t),'YTickLabel',string([30:20:100]))
            yticks([30:20:100])
            
        end
        a = gca;
        a.FontName = 'Arial';
        a.FontSize = ticksize;
        a.XLabel.FontSize = labelsize;
        a.YLabel.FontSize = labelsize;
        a.Title.FontSize = ticksize;
        
        % colorbar (only if specified above)
        if do_cb
            if find(t == [ncol:ncol:length(TFRcmb)])
                cb = colorbar;
                if cb.Limits(2) < 1
                    cb.Limits = [-round(cb.Limits(2),1) round(cb.Limits(2),1)];
                else
                    cb.Limits = [-ceil(cb.Limits(2)) ceil(cb.Limits(2))];
                end
                cb.Ticks = [cb.Limits(1),0,cb.Limits(2)];
                cb.FontSize = ticksize;
                cb.FontName = 'Arial';
                cb.FontName = 'Arial';
                cb.FontSize = ticksize;
                cb.Label.String = {'relative power', 'change'};
                cb.Label.FontSize = labelsize;
                cb.Label.FontName = 'Arial';
            end
        end
        
    end
    
    % save
    if do_cb
        print(f, fullfile(PATHPLOTS,[SUBJ{s},'entrain_alpha_gamma_colorbar_oneof_longeslide']),'-dpng','-r300')
        %print(f, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript',[SUBJ{s},'entrain_alpha_gamma_cb']),'-dsvg','-r600')
    else
        print(f, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript',[SUBJ{s},'entrain_alpha_gamma']),'-depsc','-r600')
        
        
    end
    
    close all;
    
    
    clear TFRcmb numTrl
end



