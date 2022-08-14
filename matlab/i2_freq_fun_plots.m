%% Entrainment gamma rFT
% PhD project 1
% 
% power and plv as a function of frequency

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%% Settings
clear all; close all; clc;
BATCH = 3;
MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
ticksize = 18;
labelsize = 20;
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH, 'matlab','template function'))
addpath(fullfile(MAINPATH, 'matlab','kd fun'))
addpath(fullfile('C:\Users\katha\Documents\MATLAB','fieldtrip'));
ft_defaults;

oneof = 0;
sliwin = .5;
TOIBSL = [-.75 -.25];
flickfreq = [52:2:90];                      % stimulation frequency
addstr = {'entrainment','resonance'}; 

POWPATH = fullfile(MAINPATH, 'results', 'power');

% colours
bluecode = [0 0.4470 0.7410]; % blue
orcode = [0.8500 0.3250 0.0980]; %orange
yellcode = [0.9290  0.4980 0.1250]; % yellow
cycode = [0.3010 0.7450 0.9330]; %cyan

GAMPATH = fullfile(POWPATH,'gammatron',['Batch_',num2str(BATCH)]);
PATHENT = fullfile(POWPATH,'entrainment',['Batch_',num2str(BATCH)]);
RESPATH = fullfile(POWPATH,'resonance',['Batch_',num2str(BATCH)]);

COHPATH = fullfile(MAINPATH, 'results', 'coherence',['Batch_',num2str(BATCH)]);
PLVPATH = fullfile(MAINPATH,'results','phase_locking',['Batch_',num2str(BATCH)]);
    
PATHPLOT = fullfile(MAINPATH, 'results','plots', 'arnold tongue',['Batch_',num2str(BATCH)]);
PATHOUT = fullfile(MAINPATH,'freq fun',['Batch_',num2str(BATCH)]);
PATHSTAT = fullfile(MAINPATH,'results', 'statistics','Linear regression',['Batch_',num2str(BATCH)]);

% filenames
savepow = ['pow_fun_batch_',num2str(BATCH)];
loadpow = ['Batch_',num2str(BATCH)','_freq_fun_basel_',num2str(TOIBSL(1)),'_RFT'];
loadplvcoh = ['Batch_',num2str(BATCH)','_RFT'];
saveplv = ['plv_fun_batch_',num2str(BATCH)];
savecoh = ['coh_fun_batch_',num2str(BATCH)];

% not in manus: 1/f corrected, long sliding window
if oneof
    loadpow = [loadpow,'oneof'];
    loadplvcoh = [loadplvcoh,'oneof'];
    savepow = [savepow, 'oneof'];
    saveplv = [saveplv, 'oneof'];
    savecoh = [savecoh, 'oneof'];
end

if sliwin > .5
    COHPATH = fullfile(COHPATH,'long slidewin');
    PLVPATH = fullfile(PLVPATH,'long slidewin');
    POWPATH = fullfile(POWPATH,'long slidewin');

end

mkdir(PATHOUT)
folds = dir(GAMPATH);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = []; 

for s = 1:length(SUBJ)
    load(fullfile(GAMPATH, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end
 
%keep subjects whose IGF is > 56
keepSUBJ = SUBJ(find(IGF > 56));
SOI_all = SOI_all(find(IGF > 56));
IGF = IGF(find(IGF > 56));
mIGF = mean(IGF);

for h = 1:length(SOI_all)
    w = 1;              % help index
    for l = 1:length(SOI_all{h})
        SOInc{h}{w} = SOI_all{h}{l}(1:strfind(SOI_all{h}{l},'+')-1);         % first sensor
        w = w + 1;
        SOInc{h}{w} = ['MEG',SOI_all{h}{l}(strfind(SOI_all{h}{l},'+')+1:end)];
        w = w + 1;
    end
end
%% Relative power change as a function of freq

% compute standard error
for a = 1:length(addstr)-1
    condit{a} = load(fullfile(POWPATH,addstr{a},[loadpow,'.mat']));
    condit{a}.sepowspct = condit{a}.stdpowspct./sqrt(size(condit{a}.M,2));
    condit{a}.sepowspct_IGF = condit{a}.stdpowspct_IGF./sqrt(size(condit{a}.powspct_IGF,2));
end    

% IGF models (parameters copied from R)
% res: flicker; ent: flicker&gratings
igfresmodel_pow = .92785 -.07281 .* condit{1}.freqvec;
igfresmodel_pow_Rsqu = .2712;
igfentmodel_pow = 2.3837 -.1638 .* condit{1}.freqvec;
igfentmodel_pow_Rsqu = .2688;

% plot
f = figure;
set(gcf, 'Position', [0, 0, 1920, 1080/3.5])
% plot1 = pow ~ frequency
subplot(121)
h = shadedErrorBar(flickfreq,condit{1}.powspct,condit{1}.stdpowspct,'lineprops',{'w-','markerfacecolor',cycode});
h.mainLine.Color = cycode; % cyan
h.patch.FaceColor = cycode;
hold on

h = shadedErrorBar(flickfreq,condit{2}.powspct,condit{2}.stdpowspct,'lineprops',{'w-','markerfacecolor',yellcode});
h.mainLine.Color = yellcode; % yellow
h.patch.FaceColor = yellcode;
xlim([flickfreq(1) flickfreq(end)])
xticks([flickfreq(1):4:flickfreq(end)])
ylim([0 5])
a = gca;
a.FontName = 'Arial';
a.XAxis.FontSize = ticksize;
a.YAxis.FontSize = ticksize;
ylabel('relative power change')
xlabel('Frequency [Hz]')
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

% lgd = legend('SE1','SE2','SE3','gammaRFT','SE1','SE2','SE3','RFTonly');
% lgd.FontName = 'Arial';
% lgd.FontSize = ticksize;
clear h

% plot1 = pow ~ distance to IGF
subplot(122)
h = shadedErrorBar(condit{1}.freqvec,condit{1}.powspct_IGF,condit{1}.stdpowspct_IGF,'lineprops',{'w-','markerfacecolor',cycode});
h.mainLine.Color = cycode; 
h.patch.FaceColor = cycode;

hold on
plot(condit{1}.freqvec,igfentmodel_pow,'LineStyle','-.','Color',bluecode,'LineWidth',1.5)
hold on
h = shadedErrorBar(condit{1}.freqvec,condit{2}.powspct_IGF,condit{2}.stdpowspct_IGF,'lineprops',{'w-','markerfacecolor',yellcode});
h.mainLine.Color = yellcode; % yellow
h.patch.FaceColor = yellcode;
hold on
plot(condit{1}.freqvec,igfresmodel_pow,'LineStyle','-.','Color',orcode,'LineWidth',1.5)

ylim([0 5])
xlim([condit{1}.freqvec(1) condit{1}.freqvec(end)])
a = gca;
a.FontName = 'Arial';
a.XAxis.FontSize = ticksize;
a.YAxis.FontSize = ticksize;
ylabel('relative power change')
xlabel('Frequency [Hz]')
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

% lgd = legend('SE1','SE2','SE3','gammaRFT','lin model gammaRFT','SE1','SE2','SE3','RFTonly','lin model RFTonly');
% lgd.FontName = 'Arial';
% lgd.FontSize = ticksize;
clear h
print(f,  fullfile(PATHPLOT,'power',[savepow,'_sd']),'-dpng','-r300')
%print(f, fullfile(PATHPLOT,'power',savepow),'-dsvg','-r600')


%% Phase locking as a function of freq

% compute standard error
for a = 1:2
    condit{a} = load(fullfile(PLVPATH,addstr{a},[loadplvcoh,'.mat']));
    condit{a}.seplvspct = condit{a}.stdplvspct./sqrt(size(condit{a}.M,2));
    condit{a}.seplvspct_IGF = condit{a}.stdplvspct_IGF./sqrt(size(condit{a}.plvspct_IGF,2));
end    

% Linear regression models (estimated in R)
igfresmodel_plv = .472761 - .011408.*condit{1}.freqvec;
igfresmodel_plv_Rsqu = .2299;
igfentmodel_plv = .423507 - .009877.*condit{1}.freqvec;
igfentmodel_plv_Rsqu = .1939;


% plot
f = figure;
set(gcf, 'Position', [0, 0, 1920, 1080/2.5])
% plot1: plv ~ freq
subplot(121)
h = shadedErrorBar(flickfreq,condit{1}.plvspct,condit{1}.seplvspct,'lineprops',{'w-','markerfacecolor',bluecode});
h.mainLine.Color = cycode; % cyan
h.patch.FaceColor = cycode;
hold on
h = shadedErrorBar(flickfreq,condit{2}.plvspct,condit{2}.seplvspct,'lineprops',{'w-','markerfacecolor',orcode});
h.mainLine.Color = yellcode; % yellow
h.patch.FaceColor = yellcode;
xlim([flickfreq(1) flickfreq(end)])
xticks([flickfreq(1):4:flickfreq(end)])
ylim([0 .6])
a = gca;
a.FontName = 'Arial';
a.XAxis.FontSize = ticksize;
a.YAxis.FontSize = ticksize;
ylabel('phase locking value [a.u.]')
xlabel('Frequency [Hz]')
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

% lgd = legend('SE1','SE2','SE3','gammaRFT','SE1','SE2','SE3','RFTonly');
% lgd.FontName = 'Arial';
% lgd.FontSize = ticksize;

% plot1: plv ~ freq
subplot(122)
h= shadedErrorBar(condit{1}.freqvec,condit{1}.plvspct_IGF,condit{1}.seplvspct_IGF,'lineprops',{'w-','markerfacecolor',bluecode});
h.mainLine.Color = cycode; % cyan
h.patch.FaceColor = cycode;


hold on
plot(condit{1}.freqvec,igfentmodel_plv,'LineStyle','-.','Color',bluecode,'LineWidth',1.5)
hold on
h=shadedErrorBar(condit{1}.freqvec,condit{2}.plvspct_IGF,condit{2}.seplvspct_IGF,'lineprops',{'w-','markerfacecolor',orcode});
h.mainLine.Color = yellcode; % yellow
h.patch.FaceColor = yellcode;
xlim([condit{1}.freqvec(1) condit{1}.freqvec(end)])
hold on
plot(condit{1}.freqvec,igfresmodel_plv,'LineStyle','-.','Color',orcode,'LineWidth',1.5)
ylim([0 .6])
a = gca;
a.FontName = 'Arial';
a.XAxis.FontSize = ticksize;
a.YAxis.FontSize = ticksize;
ylabel('phase locking value [a.u.]')
xlabel('Frequency [Hz]')
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;

% lgd = legend('SE1','SE2','SE3','gammaRFT','lin model gammaRFT','SE1','SE2','SE3','RFTonly','lin model RFTonly');
% lgd.FontName = 'Arial';
% lgd.FontSize = ticksize;
% 
% print(f, fullfile(PATHPLOT,'phase locking',saveplv),'-dpng','-r300')
% print(f, fullfile(PATHPLOT,'phase locking',saveplv),'-dsvg','-r300')
print(f,  fullfile('X:\Manuscript',saveplv),'-dsvg','-r1000')
