%% Entrainment gamma rFT
% PhD project 1
%
%
% entrainment condition: Plot Grandaverage TFR

% This script:
% loads in the results of the Grandavergae TFR in the flicker+gamma condition
% aligned to the participant's IGF 
% plots the result

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen


%% settings
clear all; close all; clc
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
MATPATH = 'C:\Users\katha\Documents\MATLAB';

oneof = 0;                  % 1/f corrected data?
sliwin = .5;                % sliding window length?

addpath(fullfile(MATPATH,'fieldtrip'))
addpath(fullfile(MAINPATH,'matlab','template function'))
ft_defaults
BATCH = 3;                  % which group of participants?

% font size plots
ticksize = 12;
labelsize = 16;
do_cb = 1;                  % colorbar?

TOIBSL = [-0.75 -0.25];     % baseline interval
TOI = [2.25 3.75];          % flicker+gamma interval

% load in subjects, discard IGF<58
SUBJ = {};
IGF = [];
SOI_all = {};
numSens = [];

if oneof
    addstr = 'oneof';
else
    addstr = '';
end 

PATHENT = fullfile(MAINPATH,'results','power','entrainment',['Batch_',num2str(BATCH)],addstr);
PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'power','entrainment',['Batch_',num2str(BATCH)],addstr);
mkdir(PATHPLOT)

if sliwin > .5
    addstr = 'long slidewin';
else
    addstr = '';
end

PATHENT = fullfile(PATHENT,addstr);
PATHPLOT = fullfile(PATHPLOT,addstr);
mkdir(PATHPLOT)
% read in subjects
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


%test: keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
exclSUBJ = SUBJ(find(IGF < 56));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));

SUBJ = keepSUBJ;
%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

leftLim = minIGF - 52;
rightLim = 90-maxIGF;
freqvec = -leftLim:1/sliwin:rightLim;


%% Plot
load(fullfile(PATHENT,['GA_align_bsl',num2str(TOIBSL(1)),'_sliwi_',num2str(sliwin),'.mat']))

% set up subplot
ncol = 4;
nrow = 3;
h = figure;
set(gcf, 'Position', [0, 0, 1920, 1080/2])
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colormap = 'jet';
cfg.zlim = [-2 2];
%cfg.xlim = [TFRGA{1}.time(1) TFRGA{1}.time(end)-0.25];
cfg.baseline = 'no';
cfg.title = '';
% set subplot layout (using tight_subplot)
[ha pos] = tight_subplot(nrow,ncol,[.05 .025],[.125 .075],[.1 .08]);

for t = 1:length(TFRGAent)

    if find(t == [ncol:ncol:length(TFRGAent)])
       cfg.colorbar = 'no';
    else
     cfg.colorbar = 'no';
    end
    axes(ha(t))
    ft_singleplotTFR(cfg,TFRGAent{t})           % plot TFR
    xlim([-0.5 5.5])
    xticks([])
    
    % add line to mark IGF, RFT 
    if freqvec(t) == 0
        title([''])
        hold on  
        line(TFRGAent{t}.time(find(TFRGAent{t}.time>=0.1)), repmat(freqvec(t),1,length(TFRGAent{t}.time(find(TFRGAent{t}.time>=0.1)))), 'Color', 'black', 'LineStyle', '-.', 'LineWidth',1)
        text(-0.55,0-1,'IGF', 'FontSize', ticksize, 'FontName', 'Arial')

    elseif freqvec(t) ~= 0
    title([''])
        %   title(['Entrainment at IGF + ', num2str(TFRGAcmb{t}.freq(t)), ', Grandaverage'],'FontSize', 10', 'FontWeight', 'normal', 'FontName', 'Arial');
    hold on
        % mark flicker freq
    line(TFRGAent{t}.time(find(TFRGAent{t}.time>=2 & TFRGAent{t}.time<=4)), repmat(freqvec(t),1,length(TFRGAent{t}.time(find(TFRGAent{t}.time>=2& TFRGAent{t}.time<=4)))), 'Color', 'black', 'LineStyle', '-.', 'LineWidth',1)
    % mark gamma freq
    line(TFRGAent{t}.time(find(TFRGAent{t}.time>=0.1)), repmat(0,1,length(TFRGAent{t}.time(find(TFRGAent{t}.time>=0.1)))), 'Color', 'black', 'LineStyle', '--', 'LineWidth',1)
    
    %      line(TFRcmb{t}.time, repmat(10,1,length(TFRcmb{t}.time)), 'Color', 'blue', 'LineStyle', '--')
    % set x axis for lower row
    text(-0.5,0-1,'IGF', 'FontSize', ticksize,'FontName', 'Arial')
    
    if freqvec(t) > 0 && freqvec(t) < 10
        text(2-.3,freqvec(t)+1,num2str(freqvec(t)), 'FontSize', ticksize,'FontName', 'Arial')
    elseif freqvec(t) > 0 && freqvec(t) >= 10
        text(2-.5,freqvec(t)+1,num2str(freqvec(t)), 'FontSize', ticksize,'FontName', 'Arial')
    elseif freqvec(t) < 0
        text(2-.4,freqvec(t)-1,num2str(freqvec(t)), 'FontSize', ticksize,'FontName', 'Arial')   
    end
    end

    % set axes
    % x axis label
    if t > length(TFRGAent)-ncol
        xlabel('time [s]')
        xticks([-0:1:5])
        set(ha(t),'XTickLabel',string([0:1:5]))

    end

    if find(t == [1:ncol:length(TFRGAent)])
        ylabel({'Frequency';'IGF +[Hz]'})
        yticks([-24:8:24]) 
       set(ha(t),'YTickLabel',string([-24:8:24]))
    end
    b = gca;
    b.FontSize = ticksize;
    b.FontName = 'Arial';
    b.XLabel.FontSize = labelsize;
    b.YLabel.FontSize = labelsize;
    
    % colorbar?
    if do_cb
        if find(t == [ncol:ncol:length(TFRGAent)])
            cb = colorbar;
            cb.FontSize = ticksize;
            cb.FontName = 'Arial';
            cb.Label.String = {'relative power'; 'change'};
            cb.Label.FontSize = labelsize;
        end
    end
    
end


if do_cb
    print(h, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','TFRGA_ent_IGF_cb'),'-dsvg','-r600');
else
    print(h, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','TFRGA_ent_IGF'),'-dsvg','-r600');
  %  print(h, fullfile(PATHPLOT,['TFRGA_ent_IGF']),'-depsc','-r1000');
end

%savefig(h, fullfile(PATHPLOT,['TFR_GA_ent.fig']))


clear TFRGAcmb TFRGAgrad TFR TFRal sens_map freqIGF 

close all

