%% Entrainment gamma rFT
% PhD project 1
% 
% Identification of Individual Gamma Frequencies (IGFs)
% - Sensors of Interest (SOI) identified by hand
% - IGF identified by averaging over 0.25 - 1.75 after stimulus onset

% [c] Katharina Duecker
% PGR Centre for Human Brain Health, University of Birmingham
% k.duecker@bham.ac.uk

% supervisor: Ole Jensen

clear all; close all; clc;

%% settings

% plotting
ticksize = 16;
labelsize = 20;
bluecode = [0 0.4470 0.7410];

% data
BATCH = 3;              % sample
oneof = 0;              % 1/f corrected data?
sli = 0.5;              % sliding window length

% set folders
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;
PATHIN = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOTS = fullfile(MAINPATH,'results','plots', 'power','IGF', ['Batch_',num2str(BATCH)]);


mkdir(PATHPLOTS)
cd(PATHIN)
SUBJ = cellstr(ls());
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));

%% Identify IGFs
for s = 1:length(SUBJ)
    
    if ~exist(fullfile(PATHIN,SUBJ{s},'SOI_freq.mat'))
        close all
        load(fullfile(PATHIN,SUBJ{s},['TFR_gammatron_slide_',num2str(sli),'.mat']))
       
        % combine planar
        cfg = [];
        cfg.method = 'sum';
        TFRgamcmb = ft_combineplanar(cfg, TFRgam);
        
        % plot planar & identify SOI by hand
        cfg = [];
        cfg.baseline       = [-0.75 -0.25];
        cfg.baselinetype   = 'relchange';
        cfg.xlim           = [-0.75 1.75];
        cfg.zlim           = 'maxabs';
        cfg.layout         = 'neuromag306cmb.lay';
        cfg.colormap       = 'jet';
        figure; ft_multiplotTFR(cfg, TFRgamcmb);
        
        SOI = input('Sensors of interest? {}');
        
        load(fullfile(PATHIN,SUBJ{s},['IGF_peak_slide_',num2str(sli),'.mat']))
        
        cfg = [];
        cfg.method = 'sum';
        powChanGam = ft_combineplanar(cfg, rlChanGam);
        
        
        % extract powspctrm SOIs
        spctr = powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:,:);
        
        % find maximum
        [amp1 p] = max(mean(spctr,1));
        
        % find second maximum (replace first maximum with zero)
        mockspec = spctr;
        mockspec(:,p) = zeros(length(SOI),1);
        [amp2 p2] = max(mean(mockspec,1));
        
        % if second maximum follows first maximum and is not the very next
        % point, gamFreq is the mean between the two frequencies
        if p2 > p && p2 ~= p +1 && powChanGam.freq(p2) < 80
            gamFreq = mean([powChanGam.freq(p),powChanGam.freq(p2)]);
            if mod(gamFreq,2)
                gamFreq = gamFreq +1;
            end
        else
            % else gamma is the frequency at maximum
            gamFreq = powChanGam.freq(p);
        end
        
        save(fullfile(PATHIN,SUBJ{s},'SOI_freq.mat'), 'SOI', 'gamFreq')
        
        clear TFRgam TFRgamcmb relCh powChanGam
    end
end

%% Keep subjects whose IGF > 56 Hz
for s = 1:length(SUBJ)
    load(fullfile(PATHIN, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

keepSUBJ = SUBJ(IGF>56);
keepIGF = IGF(IGF>56);
exclSUBJ = setdiff(SUBJ,keepSUBJ);


%% Example participants

% setting for subplot
ncol = 2;
nrow = 2;
h = 0;
fig = figure;
set(gcf, 'Position', [0, 0, 1920/1.75, 1080/1.75])

for s = [1,3]
    h = h+1;
    
    load(fullfile(PATHIN,keepSUBJ{s},['TFR_gammatron.mat']))
    load(fullfile(PATHIN,keepSUBJ{s},'SOI_freq.mat'))
    
    %% plot TFR
    subplot(ncol,nrow,h)
    cfg = [];
    cfg.channel = 'MEGGRAD';
    TFRgam = ft_selectdata(cfg, TFRgam);
    cfg = [];
    cfg.method = 'sum';
    TFRgamcmb = ft_combineplanar(cfg, TFRgam);
    
    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.xlim = [-0.75 1.75];
    cfg.zlim = 'maxabs';
    cfg.channel = SOI';
    cfg.colormap = 'jet';
    cfg.ylim = [30 90];
    ft_singleplotTFR(cfg, TFRgamcmb);
    title('');
    yticks([])
    xticks([])
    
    % mark flicker freq
    line(TFRgamcmb.time(find(TFRgamcmb.time == -0.75):end), repmat(gamFreq,1,length(TFRgamcmb.time)-find(TFRgamcmb.time == -0.75)+1), 'Color', 'black', 'LineStyle', '-.','LineWidth',1)
    text(-0.75,gamFreq+4,num2str(gamFreq),'FontSize', ticksize,'FontName', 'Arial')
    
    c = colorbar;
    if c.Limits(2) < 1
        c.Limits = [-round(c.Limits(2),1) round(c.Limits(2),1)];
    else
        c.Limits = [-ceil(c.Limits(2)) ceil(c.Limits(2))];
    end
    c.Ticks = [c.Limits(1),0,c.Limits(2)];
    c.FontSize = ticksize;
    c.FontName = 'Arial';
    
    c.Label.String = {'relative power', 'change'};
    c.Label.FontName = 'Arial';
    c.Label.FontSize = labelsize;
    c.Label.Position = [2.0306 7.9076e-07 0];
    
    ylabel('Frequency [Hz]')
    yticks([30:30:90])
    if h <3
        xlabel(' ')
        xticks([ ])
    else
        xlabel('time [s]')
        xticks([-0.5:0.5:1.5])
    end
          
    a = gca;
    a.FontSize = ticksize;
    a.FontName ='Arial';
    a.XLabel.FontSize = labelsize;
    a.YLabel.FontSize = labelsize;
    
    %% plot spectrum ( = averaged TFR)
    
    h = h+1;
    subplot(ncol,nrow,h)

    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.parameter = 'powspctrm';
    powChanGam = ft_freqbaseline(cfg, TFRgamcmb);
    
    cfg = [];
    cfg.latency = [0.25 1.75];
    cfg.avgovertime = 'yes';
    powChanGam = ft_selectdata(cfg,powChanGam);
    
    clear TFRgamcmb
    % find gamma frequency and relative power change
    [amp p] = max(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1));
    plot(powChanGam.freq, mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1),'LineWidth',1,'Color',bluecode)
    ylim([0 amp+0.4]);
    xlim([30 90])
    set(gca, 'YTick',linspace(0,round(amp+0.3,1),3))
    ytickformat('%.1f')
    title('');
    text(gamFreq, std(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1))/3, num2str(gamFreq), 'FontSize', ticksize,'FontName', 'Arial');
    
    line([gamFreq gamFreq],[0 amp+0.4],'Color', 'black', 'LineStyle', ':','LineWidth',1.5)
    %  ylabel({'relative power', 'change [a.u.]'})
    
    if h <3
        xlabel(' ')
        xticks([ ])
    else
        xlabel('Frequency [Hz]')
        xticks([30:10:100])
    end
    
    a = gca;
    a.FontSize = ticksize;
    a.FontName ='Arial';
    a.XLabel.FontSize = labelsize;
    a.YLabel.FontSize = labelsize;
    clear powChanGam   

end
print(fig, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','IGF_ex_SUBJ'),'-dsvg','-r1000')


%% Sanity checks: plot TFRs of all subjects that were kept in analysis
close all

f = figure;
for s = 1:length(SUBJ)
    load(fullfile(PATHIN,SUBJ{s},'TFR_gammatron.mat'))
    % combine planar
    cfg = [];
    cfg.method = 'sum';
    TFRgamcmb = ft_combineplanar(cfg, TFRgam);
    load(fullfile(PATHIN,SUBJ{s},'SOI_freq.mat'))
    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.xlim = [0 2];
    cfg.zlim = 'maxabs';
    cfg.marker = 'on';
    cfg.layout = 'neuromag306cmb.lay';
    cfg.comment = 'no';
    cfg.colormap = 'jet';
    cfg.colorbar = 'yes';
    cfg.ylim = [gamFreq gamFreq];
    subplot(5,3,s)
    ft_topoplotTFR(cfg, TFRgamcmb);
    c = colorbar;
    c.Label.String = 'relative power change';
    c.Label.FontSize = 9;
    c.Label.FontName = 'Arial';
    title(['IGF subject ', num2str(s)],'FontSize', 10, 'FontWeight', 'normal', 'FontName', 'Arial');
end
savefig(fullfile('X:\results\PLOTS_GAM', 'topo_all'))

% plot spectra of subjects kept
ncol = 5;
nrow = 5;
f = figure;
set(gcf, 'Position', [0, 0, 1920, 1080])
for s = 1:length(keepSUBJ)
    if sli > .5
        
    load(fullfile(PATHIN,keepSUBJ{s},['TFR_gammatron_slide_',num2str(sli),'.mat'])) 
    else
        load(fullfile(PATHIN,keepSUBJ{s},['TFR_gammatron.mat']))
    end
    load(fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)],keepSUBJ{s},'SOI_freq.mat'))
      
    % combine planar
    cfg = [];
    cfg.channel = 'MEGGRAD';
    TFRgam = ft_selectdata(cfg, TFRgam);
    cfg = [];
    cfg.method = 'sum';
    TFRgamcmb = ft_combineplanar(cfg, TFRgam);
    
    % plot planar
    subplot(nrow,ncol,s)
    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.xlim = [-0.75 1.75];
    cfg.zlim = 'maxabs';
    cfg.channel = SOI';
    cfg.colormap = 'jet';
    cfg.ylim = [30 90];

    ft_singleplotTFR(cfg, TFRgamcmb);
    title('');
    yticks([])
    xticks([])
    
    % mark flicker freq
    line(TFRgamcmb.time(find(TFRgamcmb.time == -0.75):end), repmat(gamFreq,1,length(TFRgamcmb.time)-find(TFRgamcmb.time == -0.75)+1), 'Color', 'black', 'LineStyle', '-.','LineWidth',1)
    text(-0.75,gamFreq+4,num2str(gamFreq),'FontSize', ticksize,'FontName', 'Arial')

    c = colorbar;
     if c.Limits(2) < 1
            c.Limits = [-round(c.Limits(2),1) round(c.Limits(2),1)];
        else
            c.Limits = [-ceil(c.Limits(2)) ceil(c.Limits(2))];
        end
        c.Ticks = [c.Limits(1),0,c.Limits(2)];
        c.FontSize = ticksize;
        c.FontName = 'Arial';
    if find(s >= length(keepSUBJ)-ncol+1)
        xlabel('time [s]')  
        xticks([-0.5:0.5:1.5])
    end

    if find(s == [1:ncol:length(SUBJ)-ncol])
        ylabel('Frequency [Hz]')
        yticks([30:30:90])
    elseif find(s == [ncol:ncol:length(keepSUBJ)])
        c.Label.String = {'relative power', 'change'};
        c.Label.FontName = 'Arial';
        c.Label.FontSize = labelsize;
        c.Label.Position = [2.0306 7.9076e-07 0];
    end
    if s == length(keepSUBJ)
        c.Label.String = {'relative power', 'change'};
        c.Label.FontSize = labelsize;
        c.Label.Position = [2.0306 7.9076e-07 0];              
    end 
    a = gca;
    a.FontSize = ticksize;
    a.FontName ='Arial';
    a.XLabel.FontSize = labelsize;
    a.YLabel.FontSize = labelsize;
end
print(f, fullfile(PATHPLOTS,['subplot_IGF_Batch_3_slide_',num2str(sli)]),'-dpng','-r300')
print(f, fullfile('X:\Manuscript','subplot_IGF_slide'),'-dsvg','-r1000')


%% plot excluded subj
ncol = 4;
nrow = 2;
f = figure;
set(gcf, 'Position', [0, 0, 1920/1.2, 1080/2])
for s = 1:length(exclSUBJ)
    if sli > .5
        
    load(fullfile(PATHIN,exclSUBJ{s},['TFR_gammatron_slide_',num2str(sli),'.mat'])) 
    else
        load(fullfile(PATHIN,exclSUBJ{s},['TFR_gammatron.mat']))
    end
    load(fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)],exclSUBJ{s},'SOI_freq.mat'))

    cfg = [];
    cfg.channel = 'MEGGRAD';
    TFRgam = ft_selectdata(cfg, TFRgam);
    cfg = [];
    cfg.method = 'sum';
    TFRgamcmb = ft_combineplanar(cfg, TFRgam);
    
    %plot planar
    subplot(nrow,ncol,s)
    cfg = [];
    cfg.baseline = [-0.75 -0.25];
    cfg.baselinetype = 'relchange';
    cfg.xlim = [-0.75 1.75];
    cfg.zlim = 'maxabs';
    cfg.channel = SOI';
    cfg.colormap = 'jet';
    cfg.ylim = [30 90];

    ft_singleplotTFR(cfg, TFRgamcmb);
   title('');
   yticks([])
   xticks([])
   
    % mark flicker freq
    line(TFRgamcmb.time(find(TFRgamcmb.time == -0.75):end), repmat(gamFreq,1,length(TFRgamcmb.time)-find(TFRgamcmb.time == -0.75)+1), 'Color', 'black', 'LineStyle', '-.','LineWidth',1)
    text(-0.75,gamFreq+4,num2str(gamFreq),'FontSize', ticksize,'FontName', 'Arial')

    c = colorbar;
     if c.Limits(2) < 1
            c.Limits = [-round(c.Limits(2),1) round(c.Limits(2),1)];
        else
            c.Limits = [-ceil(c.Limits(2)) ceil(c.Limits(2))];
        end
        c.Ticks = [c.Limits(1),0,c.Limits(2)];
        c.FontSize = ticksize;
        c.FontName = 'Arial';
    if find(s >= length(exclSUBJ)-ncol+1)
        xlabel('time [s]')  
        xticks([-0.5:0.5:1.5])
    end

    if find(s == [1:ncol:length(exclSUBJ)-ncol+1])
        ylabel('Frequency [Hz]')
        yticks([30:30:90])
    elseif find(s == [ncol:ncol:length(exclSUBJ)])
        c.Label.String = {'relative power', 'change'};
        c.Label.FontName = 'Arial';
        c.Label.FontSize = labelsize;
        c.Label.Position = [2.0306 7.9076e-07 0];
    end
    if s == length(keepSUBJ)
        c.Label.String = {'relative power', 'change'};
        c.Label.FontSize = labelsize;
        c.Label.Position = [2.0306 7.9076e-07 0];
       
                
    end
    
    a = gca;
    a.FontSize = ticksize;
    a.FontName ='Arial';
    a.XLabel.FontSize = labelsize;
    a.YLabel.FontSize = labelsize;


end
print(f, fullfile('X:\Manuscript','subplot_IGF_excl'),'-dsvg','-r1000')


close all
