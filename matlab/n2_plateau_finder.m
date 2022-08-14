%% Entrainment gamma RFT
% Phd project 1

% script identifies plateaus in phase angle time series
%
% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health


%% settings
clear all; close all;clc;
BATCH = 3;                  
condit = 2;
% condit = 1: flicker&gamma condition
% condit = 2: flicker condition

method = 'hann';        

MAINPATH = 'Z:\';
flickfreq = 52:2:90;
% font sizes 
ticksize = 18;
labelsize = 20;

% define & add paths
addpath(genpath(fullfile(MAINPATH,'matlab')));
addstr = {'entrainment', 'resonance'};

PATHIN = fullfile(MAINPATH,'results','phase slips'); 
PATHOUT = fullfile(PATHIN,method,addstr{condit});
PATHPLOT = fullfile(MAINPATH,'results','plots','phase slips');
mkdir(PATHPLOT)
PATHGAM = fullfile(MAINPATH,'results','power','gammatron',['Batch_',num2str(BATCH)]);

% list/load in subjects, discard subjects with IGF < 58
folds = dir(PATHGAM);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];            % files not starting with '201' are not subjects

for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

SUBJ = SUBJ(IGF > 56);
SOI_all = SOI_all(IGF > 56);
IGF = IGF(IGF > 56);
%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);
leftLim = minIGF - 52;


%% Plateaus: get gradient of noisy surrogate data 
% ... to find appropriate limit & to check phase angle script

Fs = 1000;                      % sampling rate
t =  (1:2000)/Fs;               % length of time window
foi = 58;                       % frequency of interest 
width = 3;                      % width of taper in cycles 
Nfoi = floor(width*Fs/foi);     % length of taper in samples

% surrogate signals
s1 = rand(1,2000) + sin(2*pi*foi*t);
s2 = sin(2*pi*foi*t+pi) + .2*rand(1,2000);

% complex hanning taper
tap = hanning(Nfoi)'.*exp(i*2*pi*foi.*(1:Nfoi)/1000);

% convolve
s1f = conv(s1,tap);
s2f = conv(s2,tap);
% positive frequencies
s1f = s1f(ceil(Nfoi/2):length(s1f)-floor(Nfoi/2));
s2f = s2f(ceil(Nfoi/2):length(s2f)-floor(Nfoi/2));

% subtract unwrapped phase angles
s_diff = unwrap(angle(s1f))-unwrap(angle(s2f));

% plot & check
subplot(121)
plot(s_diff)
title('58 Hz sine + noise')
ylabel('phase angle')
subplot(122)
s2 = sin(2*pi*(foi+1)*t+pi) + .2*rand(1,2000);
s2f = conv(s2,tap);
s2f = s2f(ceil(Nfoi/2):length(s2f)-floor(Nfoi/2));
s_diff = unwrap(angle(s1f))-unwrap(angle(s2f));
plot(s_diff)
title('58 Hz sine + noise, 50 Hz sine + noise')
ylabel('phase angle')

% get mean slope and standard deviation as plateau criterion
slope = mean(abs(gradient(s_diff)));
stdslope = std(gradient(s_diff));

%% Find plateaus in data

soi = 1;                % sensor of interest index

% helper variables
fplat = {};
fanout = [];
ntrials = [];

for condit = 1:2                % loop over conditions
    PATH = fullfile(PATHIN,method,addstr{condit});
    for s = 1:length(SUBJ)
        load(fullfile(PATH,[SUBJ{s},'mod_s_diff.mat']))
        % only check IGF
%         if strcmp(method,'hann')
%             f = find(flickfreq == IGF(s));
%         elseif strcmp(method,'hilb')
%             f = 4;
%         end
        freqal = IGF(s)-6:2:IGF(s)+6;
        for fi = 1:length(freqal)
            f = find(flickfreq == freqal(fi));
            platlen = ceil(1000/freqal(fi));                % minimum length of plateau
            ntrials(condit,fi,s) = length(s_diff{soi}{f});
            % fanninf out: std at end of trials
            alltr = vertcat(s_diff{soi}{f}{:});
            
            % fanning out in cycles
            fanout(condit,fi,s) = std(alltr(:,end))/(2*pi);
            % fanning out in ms
            fanout_ms(condit,fi,s) = platlen*fanout(condit,fi,s);
            clear alltr
            
            % plateau finder algorithm
            for tr = 1:length(s_diff{soi}{f})
                
                slilen = 1;                             % length of signal sliding window has passed (in samples)
                p = 0;                                  % plateau counter
                
                
                while slilen <= length(s_diff{soi}{f}{tr})-slilen
                    % compute mean absolute gradient of current interval
                    % with length of one flicker frequency cycle
                    curgrad = mean(abs(gradient(s_diff{soi}{f}{tr}(slilen:slilen+platlen-1))));
                    
                    % if gradient smaller than criterion -> save as plateau
                    if curgrad <= slope + stdslope
                        p = p+1;
                        fplat{condit,fi,s,tr,p} = [slilen slilen+platlen];
                        slilen = slilen+platlen-1;
                    else
                        slilen = slilen+1;
                    end
                end
%                 
%                 % check visually
%                 for pl=1:p
%                     plot(s_diff{soi}{f}{tr})
%                     xlim([fplat{condit,s,tr,pl}(1)-100 fplat{condit,s,tr,pl}(2)+100])
%                     hold on
%                     text(fplat{condit,s,tr,pl}(1),s_diff{soi}{f}{tr}(fplat{condit,s,tr,pl}(1)),num2str(fplat{condit,s,tr,pl}(1)))
%                     text(fplat{condit,s,tr,pl}(2),s_diff{soi}{f}{tr}(fplat{condit,s,tr,pl}(2)),num2str(fplat{condit,s,tr,pl}(2)))
%                     
%                     pause
%                     close all
%                 end
            end
        end
            clear slilen platlen s_diff
    end
end
   

%% Fanning out per condition

m_fanout = mean(fanout_ms,3);
ste_fanout = std(fanout_ms,[],3)/sqrt(length(SUBJ));

%% Number of plateaus

for condit = 1:2
    for fi = 1:size(fanout,2)
        for s = 1:length(SUBJ)
            for tr = 1:ntrials(condit,s)
                c = 0;
                for h = 1:length(fplat(condit,fi,s,tr,:))
                    c = c + size(fplat{condit,fi,s,tr,h},1);
                end
                numplat(condit,fi,s,tr) = c;
                
            end
            nplattr_rat(condit,fi,s) = sum(numplat(condit,fi,s,:))/ntrials(condit,fi,s);
        end
    end
end

m_plattr = mean(nplattr_rat,3);
ste_plattr = std(nplattr_rat,[],3)/sqrt(length(SUBJ));

%mkdir(fullfile(MAINPATH,'results','phase plats'))
save(fullfile(MAINPATH,'results','phase plats',['plats_grad_',num2str(slope),'.mat']),'slope','stdslope','m_fanout','m_plattr','nplattr_rat','fanout')

% diff in number of plateaus
diffplat = mean(diff(nplattr_rat,1));
stediffplat = std(diff(nplattr_rat,1))/sqrt(length(SUBJ));
M = array2table(horzcat([1:length(SUBJ)]',nplattr_rat'),'VariableNames',{'subj','fligam','fli'});
writetable(M,fullfile(MAINPATH,'results','phase plats',['plats_grad_',num2str(slope),'.xlsx']))

%% Plot shaded error bar
bluecode = [0 0.4470 0.7410]; % blue
orcode = [0.8500 0.3250 0.0980]; %orange

fig = figure;
set(gcf, 'Position', [0, 0, 1920/1.5, 1080],'renderer','Painters')

subplot(221)
h = shadedErrorBar(-6:2:6,m_fanout(1,:),ste_fanout(1,:),'lineprops',{'w-','markerfacecolor',bluecode});
h.mainLine.Color = bluecode; % blue
h.patch.FaceColor = bluecode;
hold on
h = shadedErrorBar(-6:2:6,m_fanout(2,:),ste_fanout(2,:),'lineprops',{'w-','markerfacecolor',orcode});
h.mainLine.Color = orcode; % blue
h.patch.FaceColor = orcode;
ylabel('std(phase) at t = 2s [ms]')
xlabel('frequency IGF + [Hz]')
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = ticksize;
a.YLabel.FontSize = ticksize;

subplot(223)
h = shadedErrorBar(-6:2:6,m_plattr(1,:),ste_plattr(1,:),'lineprops',{'w-','markerfacecolor',bluecode});
h.mainLine.Color = bluecode; % blue
h.patch.FaceColor = bluecode;
hold on
h = shadedErrorBar(-6:2:6,m_plattr(2,:),ste_plattr(2,:),'lineprops',{'w-','markerfacecolor',orcode});
h.mainLine.Color = orcode; % blue
h.patch.FaceColor = orcode;
ylabel('plateaus/trial')
xlabel('frequency IGF + [Hz]')
yticks(3:1:5)
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = ticksize;
a.YLabel.FontSize = ticksize;

print(fig,fullfile(PATHPLOT,['slope_',num2str(slope),'_plats_IGFpm.svg']),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,['slope_',num2str(slope),'_plats_IGFpm.png']),'-dpng','-r600')








%% Plot results
bluecode = [0 0.4470 0.7410]; % blue
orcode = [0.8500 0.3250 0.0980]; %orange

fig = figure;
set(gcf, 'Position', [0, 0, 1920/4, 1080])
subplot(211)
scatter(ones(size(SUBJ)),m_fanout(2,:),'filled','MarkerFaceColor',orcode)
hold on 
scatter(ones(size(SUBJ))*2,fanout_GA(1,:),'filled','MarkerFaceColor',bluecode)
hold on
% mark mean flicker
line([0.75:0.25:1.25],repmat(mean(fanout_GA(2,:)),1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5)
hold on
line([0.75:0.25:1.25],repmat(mean(fanout_GA(2,:))-(std(fanout_GA(2,:),[],2)/sqrt(length(fanout_GA(2,:)))),...
    1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5,'LineStyle', '-.')
hold on
line([0.75:0.25:1.25],repmat(mean(fanout_GA(2,:))+(std(fanout_GA(2,:),[],2)/sqrt(length(fanout_GA(2,:)))),...
    1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5,'LineStyle', '-.')
hold on
% mark mean gamma+flicker
line([1.75:0.25:2.25],repmat(mean(fanout_GA(1,:)),1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5)
hold on
line([1.75:0.25:2.25],repmat(mean(fanout_GA(1,:))-(std(fanout_GA(1,:),[],2)/sqrt(length(fanout_GA(1,:)))),...
    1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5,'LineStyle', '-.')
hold on
line([1.75:0.25:2.25],repmat(mean(fanout_GA(1,:))+(std(fanout_GA(1,:),[],2)/sqrt(length(fanout_GA(1,:)))),...
    1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5,'LineStyle', '-.')
xticks([1 2])
ylabel('phase angle')
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = ticksize;
a.YLabel.FontSize = ticksize;
a.XTickLabel = {'flicker','flicker+gamma'};

% number of plateas
subplot(212)
scatter(ones(size(SUBJ)),nplattr_rat(2,:),'filled','MarkerFaceColor',orcode)
hold on 
scatter(ones(size(SUBJ))*2,nplattr_rat(1,:),'filled','MarkerFaceColor',bluecode)
hold on
% mark mean flicker
line([0.75:0.25:1.25],repmat(mean(nplattr_rat(2,:)),1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5)
hold on
line([0.75:0.25:1.25],repmat(mean(nplattr_rat(2,:))-(std(nplattr_rat(2,:),[],2)/sqrt(length(nplattr_rat(2,:)))),...
    1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5,'LineStyle', '-.')
hold on
line([0.75:0.25:1.25],repmat(mean(nplattr_rat(2,:))+(std(nplattr_rat(2,:),[],2)/sqrt(length(nplattr_rat(2,:)))),...
    1,length([0.75:0.25:1.25])),'Color',orcode,'LineWidth',1.5,'LineStyle', '-.')
hold on
% mark mean gamma+flicker
line([1.75:0.25:2.25],repmat(mean(nplattr_rat(1,:)),1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5)
hold on
line([1.75:0.25:2.25],repmat(mean(nplattr_rat(1,:))-(std(nplattr_rat(1,:),[],2)/sqrt(length(nplattr_rat(1,:)))),...
    1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5,'LineStyle', '-.')
hold on
line([1.75:0.25:2.25],repmat(mean(nplattr_rat(1,:))+(std(nplattr_rat(1,:),[],2)/sqrt(length(nplattr_rat(1,:)))),...
    1,length([0.75:0.25:1.25])),'Color',bluecode,'LineWidth',1.5,'LineStyle', '-.')
xticks([1 2])
ylabel({'plateaus/trial'})
a = gca;
a.FontSize = ticksize;
a.FontName = 'Arial';
a.XLabel.FontSize = ticksize;
a.YLabel.FontSize = ticksize;
a.XTickLabel = {'flicker','flicker+gamma'};
print(fig,fullfile(PATHPLOT,['slope_',num2str(slope),'_plats.svg']),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,['slope_',num2str(slope),'_plats.png']),'-dpng','-r600')