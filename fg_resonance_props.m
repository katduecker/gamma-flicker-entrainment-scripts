%% Entrainment gamma rFT
% PhD project 3
% 
% Resonance properties: 
% identify frequency that induces strongest flicker response

% [c] PGR: K. Duecker
%          
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

clear all; close all; clc

% settings
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')
ft_defaults;
ticksize = 12;
labelsize = 16;
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
PATHGAM = fullfile(MAINPATH, 'results','power','gammatron','Batch_3');
PATHENT = fullfile(MAINPATH, 'results','power','entrainment','Batch_3');
PATHRES = fullfile(MAINPATH, 'results','power','resonance','Batch_3');
PATHPLOT = fullfile(MAINPATH,'results','plots','power');
flickfreq = 52:2:90;
folds = dir(PATHGAM);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];


for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end


% keep subjects whose IGF is > 56
SUBJ = SUBJ(IGF > 56);
SOI_all = SOI_all(IGF > 56);
IGF = IGF(IGF > 56);
IGFpow_gamma = zeros(size(SUBJ));
IGFpow_rft = zeros(size(SUBJ));

RES_props = zeros(length(SUBJ),5);

% store IGF and max flicker frequency
for s = 1:length(SUBJ)
    RES_props(s,1) = IGF(s);            % store IGF
    
    %% resonance properties flicker&gratings condition
     load(fullfile(PATHENT,SUBJ{s}, 'TFR_ent_alpha_gamma.mat'))
     
     % get spectra
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
    clear TFRgrad TFR
    
    for r = 1:length(TFRcmb)

        cfg = [];
        cfg.baseline = [-0.75 -0.25];
        cfg.baselinetype = 'relchange';
        cfg.paramter = 'powspctrm';
        TFRcmb{r} = ft_freqbaseline(cfg, TFRcmb{r});

        % latency: flicker&gratings on
        cfg = [];
        cfg.channel = SOI_all{s}';
        cfg.avgoverchan = 'yes';
        cfg.latency = [2.25 3.75];
        cfg.avgovertime = 'yes';
        TFRSOI{r} = ft_selectdata(cfg,TFRcmb{r});
    end
       
    
    TFRpow = [];
    for r = 1:length(TFRSOI)
        TFRpow = [TFRpow,TFRSOI{r}.powspctrm(TFRSOI{r}.freq == flickfreq(r))];
    end
    [maxpow,idx] = max(TFRpow);
    RES_props(s,2) = flickfreq(idx);
    RES_props(s,3) = maxpow;
    
    % power at IGF during RFT
     for r = find(flickfreq == RES_props(s,2))
        cfg.frequency = IGF(s);
        TFRSOI = ft_selectdata(cfg,TFRcmb{r});
     end
      clear TFRcmb
    IGFpow_gamma(s) = TFRSOI.powspctrm;
    clear idx TFRpow TFRSOI
    
    
    %% resonance properties flicker without grating
     load(fullfile(PATHRES,SUBJ{s}, 'TFR_res_alpha_RFT.mat'));

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
        
        % latency: flicker only
        cfg = [];
        cfg.channel = SOI_all{s}';
        cfg.avgoverchan = 'yes';
        cfg.latency = [0.25 1.75];
        cfg.avgovertime = 'yes';
        TFRSOI{r} = ft_selectdata(cfg,TFRcmb{r});
    end
    
    
    TFRpow = [];
    for r = 1:length(TFRSOI)
        TFRpow = [TFRpow,TFRSOI{r}.powspctrm(TFRSOI{r}.freq == flickfreq(r))];
    end
    [maxpow,idx] = max(TFRpow);
    RES_props(s,4) = flickfreq(idx);
    RES_props(s,5) = maxpow;
    % power at IGF during RFT
     for r = find(flickfreq == RES_props(s,4))
        cfg.frequency = IGF(s);
        TFRSOI = ft_selectdata(cfg,TFRcmb{r});
     end
        IGFpow_rft(s) = TFRSOI.powspctrm;
    clear idx TFRpow TFRSOI TFRcmb

end

% save as excel
save(fullfile(MAINPATH, 'results','power','res_props.mat'),'RES_props','IGFpow_gamma','IGFpow_rft')
T = table(RES_props(:,1), RES_props(:,2), RES_props(:,3),RES_props(:,4),RES_props(:,5),'VariableNames',{'IGF','res_gamma','pow_res_gamma','res_rft','pow_res_rft'});
    writetable(T,fullfile(MAINPATH, 'results','power','res_props.xlsx'))

 %% Pearson correlation
 % IGF & preferred freq RFTonly
 
 [Rrft,prft] = corrcoef(IGF,RES_props(:,4))
 [Rgam,pgam] = corrcoef(IGF,RES_props(:,2))

    %% Resonance relative to IGF
fig = figure;
set(gcf, 'Position', [0, 0, 480,465])
s = scatter(RES_props(:,1),RES_props(:,2),'filled','MarkerFaceColor',[0 0.4470 0.7410]);
% s.MarkerFaceAlpha = .5;
% hold on
% s = scatter(RES_props(:,1),RES_props(:,4),'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
% s.MarkerFaceAlpha = .5;
hold on
plot(52:2:90,52:2:90,'LineWidth',1.5,'LineStyle','-.','Color',[0 0 0])
xlabel('IGF [Hz]')
ylabel('Resonance Frequency [Hz]')
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
xlim([52 90])
xticks([52:6:90])
ylim([52 90])
yticks([52:6:90])
% lgd = legend('res gammaRFT', 'res RFTonly');

print(fig,fullfile(PATHPLOT,'IGF_resfreq_IGF'),'-dsvg','-r600')

    
    %% scatterplot
fig = figure;
set(gcf, 'Position', [0, 0, 1920/2, 1080/2])

scatter(1:length(SUBJ),RES_props(:,1)','filled','MarkerFaceColor',[0 0 0])
hold on
s = scatter(1:length(SUBJ),RES_props(:,2)','filled','MarkerFaceColor',[0 0.4470 0.7410]);
s.MarkerFaceAlpha = .75;
hold on
s = scatter(1:length(SUBJ),RES_props(:,4)','filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
s.MarkerFaceAlpha = .75;
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,2)),length(SUBJ),1),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,2))-(std(RES_props(:,2))./sqrt(length(SUBJ))),length(SUBJ),1),'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle', '-.')
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,2))+(std(RES_props(:,2))./sqrt(length(SUBJ))),length(SUBJ),1),'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle', '-.')
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,4)),length(SUBJ),1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,4))-(std(RES_props(:,4))./sqrt(length(SUBJ))),length(SUBJ),1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'LineStyle', '-.')
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,4))+(std(RES_props(:,4))./sqrt(length(SUBJ))),length(SUBJ),1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'LineStyle', '-.')
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,1)),length(SUBJ),1),'Color','k','LineWidth',1.5)
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,1))-(std(RES_props(:,1))./sqrt(length(SUBJ))),length(SUBJ),1),'Color','k','LineWidth',1,'LineStyle', '-.')
hold on
line(1:length(SUBJ),repmat(mean(RES_props(:,1))+(std(RES_props(:,1))./sqrt(length(SUBJ))),length(SUBJ),1),'Color','k','LineWidth',1,'LineStyle', '-.')
lgd = legend('IGF','res gammaRFT', 'res RFTonly');
xlim([1 22])
ylim([52 90])
xlabel('subject')
ylabel('Frequency Hz')
yticks(52:4:90)
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

print(fig,fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','res_prop_subj'),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,'res_prop_subj'),'-dpng','-r600')

close all;


%% Power change at Resonance frequency & power at IGF

% gamma during RFT in both conditions - senseless probably
fig = figure;
set(gcf, 'Position', [0, 0, 1920/6, 1080/2])

scatter(repmat(1,1,length(SUBJ)),IGFpow_gamma,'filled','MarkerFaceColor',[0 0.4470 0.7410])
hold on
scatter(repmat(2,1,length(SUBJ)),IGFpow_rft,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980])
xticks([1,2])
xlim([0.6 2.4])
xlabel('condition')
ylabel('relative power change at IGF')
xticklabels({'gammaRFT', 'RFTonly'})
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

close all;

%% power at resonance
fig = figure;
set(gcf, 'Position', [0, 0, 1920/2, 1080/2])
RES_props(:,3) = RES_props(:,3)./2;
RES_props(:,5) = RES_props(:,5)./2;

s = scatter(RES_props(:,2),RES_props(:,3),'filled','MarkerFaceColor',[0 0.4470 0.7410]);
hold on
s = scatter(RES_props(:,4),RES_props(:,5),'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
hold on
freqrange = ceil(mean(RES_props(:,2)./2)-std(RES_props(:,2)./2)):2:ceil(mean(RES_props(:,2))+std(RES_props(:,2)));
line(freqrange,repmat(mean(RES_props(:,3)),1,length(freqrange)),'Color',[0 0.4470 0.7410],'LineWidth',1.5)

hold on
powrange = [mean(RES_props(:,3))-(std(RES_props(:,3)./sqrt(length(SUBJ)))),...
    mean(RES_props(:,3))+(std(RES_props(:,3)./sqrt(length(SUBJ))))];
line(freqrange,repmat(powrange(1),1,length(freqrange)),'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle', '-.')
hold on
line(freqrange,repmat(powrange(2),1,length(freqrange)),'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle', '-.')

hold on
freqrange = ceil(mean(RES_props(:,4))-std(RES_props(:,4))):2:ceil(mean(RES_props(:,4))+std(RES_props(:,4)));
line(freqrange,repmat(mean(RES_props(:,5)),1,length(freqrange)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)

hold on
powrange = [mean(RES_props(:,5))-(std(RES_props(:,5)./sqrt(length(SUBJ)))),...
    mean(RES_props(:,5))+(std(RES_props(:,5)./sqrt(length(SUBJ))))];
line(freqrange,repmat(powrange(1),1,length(freqrange)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'LineStyle', '-.')
hold on
line(freqrange,repmat(powrange(2),1,length(freqrange)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1,'LineStyle', '-.')

ylabel('Relative power change')
xlabel('Resonance frequency [Hz]')
lgd = legend('res gammaRFT', 'res RFTonly');
xlim([52 90])
xticks([52:4:90])

a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
print(fig,fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript','res_prop_subj_power'),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,'res_prop_subj_power'),'-dpng','-r600')




[RHOgam Pgam] = corr(RES_props(:,1),RES_props(:,2),'type','Spearman')
[RHOrft Prft] = corr(RES_props(:,1),RES_props(:,4),'type','Spearman')
