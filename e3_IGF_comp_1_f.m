%% Entrainment gamma rFT
% PhD project 1
%
%
% compare IGF of 500 ms sliding window approach with 1/f corrected & 1
% second sliding window data

% [c] student: K. Duecker
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

clear all; close all; clc;

BATCH = 3;

% main folder
MAINPATH = 'X:\';
addpath(fullfile(MAINPATH,'fieldtrip'));
%rmpath('rds/bear-apps/handbuilt/software/matlab/ThirdPartyToolboxes');
ft_defaults;
PATHIN = fullfile(MAINPATH,'results', 'power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOTS = fullfile(MAINPATH,'results','plots', 'power','IGF', ['Batch_',num2str(BATCH)]);

mkdir(PATHPLOTS)
cd(PATHIN)
SUBJ = cellstr(ls());
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));
PATHOOF = fullfile(PATHIN,'oneof');
PATHLSW = fullfile(PATHIN,'long slidewin');
PATHOL = fullfile(PATHOOF,'long slidewin');

PATHS = {PATHOOF,PATHLSW,PATHOL};
sli = [.5, 1, 1];
IGFcomp = zeros(length(SUBJ),4);
for s = 1:length(SUBJ)
    load(fullfile(PATHIN,SUBJ{s},'SOI_freq.mat'))
    IGFcomp(s,1) = gamFreq;
    % one over f corrected
    
    for pths = 1:length(PATHS)
        load(fullfile(PATHS{pths},SUBJ{s},['IGF_peak_slide_',num2str(sli(pths)),'.mat']))
        
        cfg = [];
        cfg.method = 'sum';
        powChanGam = ft_combineplanar(cfg, rlChanGam);
        
        
        % copy spectrum
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
            %else gamma is the frequency at maximum
            gamFreq = powChanGam.freq(p);
        end
        
        IGFcomp(s,pths+1) = gamFreq;
        
        clear powChanGam mocjspec spctr amp* p* gamFreq
    end
    
end

%% plot subjects for which a different IGF was found
ncol = 6;
nrow = 5;
h = figure;
set(gcf, 'Position', [0, 0, 1920, 1080])
for s = 1:length(SUBJ)
    if IGFcomp(s,1) > IGFcomp(s,3)
    textpos1 = 2;
    textpos2 = -6;
    elseif IGFcomp(s,1) < IGFcomp(s,3)
        textpos1 = -6;
    textpos2 = 2;
    end
   load(fullfile(PATHIN,SUBJ{s},'SOI_freq.mat'))
    load(fullfile(PATHIN,SUBJ{s},['IGF_peak.mat']))
    % combine planar
    cfg = [];
    cfg.method = 'sum';
    powChanGam = ft_combineplanar(cfg, rlChanGam);
    clear rlChanGam
    subplot(nrow,ncol,s)
    plot(powChanGam.freq, mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
       [amp p] = max(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1));

    ylim([0 amp+0.4]);
    xlim([30 90])
    set(gca, 'YTick',linspace(0,round(amp+0.3,1),3))
    ytickformat('%.1f')
    title('');
  %  yticks([])
    xticks([])
    text(IGFcomp(s,1)+textpos1, std(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1))/3, num2str(IGFcomp(s,1)), 'FontSize', 8,'FontName', 'Arial','Color',[0 0.4470 0.7410]);

    line([IGFcomp(s,1) IGFcomp(s,1)],[0 amp+0.4],'Color', 'black', 'LineStyle', ':')
    
    clear powChanGam amp
    
    %% long sliding window
    hold on
    load(fullfile(PATHLSW,SUBJ{s},['IGF_peak_slide_',num2str(1),'.mat']))
    % combine planar
    cfg = [];
    cfg.method = 'sum';
    powChanGam = ft_combineplanar(cfg, rlChanGam);
    clear rlChanGam
    p2 = plot(powChanGam.freq, mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
        [amp p] = max(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1));
    p2.Color(4) = 0.75;
    ylim([0 amp+0.4]);
    xlim([30 90])
    set(gca, 'YTick',linspace(0,round(amp+0.3,1),3))
    ytickformat('%.1f')
    title('');
    %  yticks([])
    xticks([])
    text(IGFcomp(s,3)+textpos2, std(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1))/3, num2str(IGFcomp(s,3)), 'FontSize', 8,'FontName', 'Arial','Color',[0.8500 0.3250 0.0980]);
    [amp p] = max(mean(powChanGam.powspctrm(find(ismember(powChanGam.label, SOI)),:),1));
    
    line([IGFcomp(s,3) IGFcomp(s,3)],[0 amp+0.4],'Color', 'black', 'LineStyle', ':')
    
    if find(s == [1:ncol:length(SUBJ)])
        ylabel('relative power change [a.u.]')
        jLabel = get(gca,'YLabel');
        
        set(jLabel, 'FontSize', 10, 'FontName', 'Arial');
    else
       
    end
    if find(s == [nrow*ncol-ncol+1:length(SUBJ)])
        xlabel('frequency [Hz]')  
        hLabel = get(gca,'XLabel');
        set(hLabel, 'FontSize', 10, 'FontName', 'Arial');
        set(hLabel, 'Units', 'pixels');
        xticks([30:10:100])
    end
    jLabel = get(gca,'YLabel');
    set(jLabel, 'FontName', 'Arial', 'FontSize', 9);
    hLabel = get(gca,'XLabel');
    set(hLabel, 'FontName', 'Arial', 'FontSize', 9);
    clear powChanGam   
    if s == 30
    legend('500 ms','line','1 sec','line')
    end
    clear powChanGam amp
end
print(h, fullfile(PATHPLOTS,'IGF_peak_500ms_vs_1s'),'-dpng','-r300')
print(h, fullfile(PATHPLOTS,'IGF_peak_500ms_vs_1s'),'-dsvg','-r600')
