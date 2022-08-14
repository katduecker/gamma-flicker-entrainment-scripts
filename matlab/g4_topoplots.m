%% Entrainment gamma rFT
% PhD project 1

% topoplots of gamma oscillations and flicker response

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%% settings
clear all; close all; clc;

MAINPATH = 'E:\University of Birmingham\Proj1_Entrainment';
ticksize = 16;
labelsize = 20;
addpath(fullfile('C:\Users\katha\Documents\MATLAB','fieldtrip'));
ft_defaults;
PATHGAM = fullfile(MAINPATH,'results', 'power','gammatron','Batch_3');
PATHRES = fullfile(MAINPATH,'results', 'power','resonance','Batch_3');
PATHPLOTS = fullfile(MAINPATH,'results','plots', 'power','IGF', 'Batch_3');

cd(PATHGAM)
SUBJ = cellstr(ls());
SUBJ = SUBJ(logical(cell2mat(cellfun(@(x) strncmp(x,'201',3),SUBJ,'UniformOutput',false))));

for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end

SUBJ = SUBJ(IGF>56);
IGF = IGF(IGF>56);

s = 3;                  % representative subject ID (used 1,3)

%% topo gamma oscillations

load(fullfile(PATHGAM,SUBJ{s},['TFR_gammatron.mat']))
load(fullfile(PATHGAM,SUBJ{s},'SOI_freq.mat'))

% select grads
cfg = [];
cfg.channel = 'MEGGRAD';
TFRgam = ft_selectdata(cfg, TFRgam);

% combine planar
cfg = [];
cfg.method = 'sum';
TFRgamcmb = ft_combineplanar(cfg, TFRgam);

% relative power change settings
cfg = [];
cfg.xlim = [0.25 1.75];
cfg.ylim = [IGF(s) IGF(s)];
cfg.baseline = [-.75 -0.25];
cfg.baselinetype = 'relchange';
cfg.layout = 'neuromag306cmb.lay';
cfg.comment = ' ';
%cfg.marker = 'labels';

f = figure; ft_topoplotTFR(cfg,TFRgamcmb)
c = colorbar;
c.Limits = [-c.Limits(2) c.Limits(2)];
c.FontSize = ticksize;
c.FontName = 'Arial';
c.Label.FontSize = labelsize;
c.Label.String = 'absolute power change';
colormap jet
caxis([-c.Limits(2) c.Limits(2)])

print(f, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript',[SUBJ{s},'_topo_IGF']),'-depsc','-r600')
print(f, fullfile(PATHPLOTS,[SUBJ{s},'_topo_IGF']),'-dpng','-r600')

close all; clear TFR*

%% topo flicker without gratings (evoked response)

% load TFR data (separated in 20 flicker frequencies)
load(fullfile(PATHRES,SUBJ{s},['TFR_res_alpha_RFT.mat']))

flickfreq = 52:2:78;

for f = 1:length(flickfreq)
    cfg = [];
    cfg.frequency = flickfreq(f);
    TFR_s{f} = ft_selectdata(cfg,TFR{f});
    TFR_s{f}.freq = 0;
end

% average over all frequencies to get common sensor
TFR_rft = ft_freqgrandaverage([],TFR_s{:});

% select grads
cfg = [];
cfg.channel = 'MEGGRAD';
TFR_rft = ft_selectdata(cfg,  TFR_rft);
TFR_rft.grad = TFR{1}.grad;

% combine planar
cfg = [];
cfg.method = 'sum';
TFR_rft = ft_combineplanar(cfg, TFR_rft);

% relative power change and topoplot
cfg = [];
cfg.xlim = [0.25 1.75];
cfg.ylim = [0 0];
cfg.baseline = [-.75 -0.25];
cfg.baselinetype = 'relchange';
cfg.layout = 'neuromag306cmb.lay';
cfg.comment = ' ';
%cfg.marker = 'labels';
f = figure; ft_topoplotTFR(cfg,TFR_rft)
c = colorbar;
c.Limits = [-round(c.Limits(2)) round(c.Limits(2))];
c.FontSize = ticksize;
c.FontName = 'Arial';
c.Label.FontSize = labelsize;
c.Label.String = 'relative power change';
colormap jet
caxis([-round(c.Limits(2)) round(c.Limits(2))])

print(f, fullfile('C:\Users\katha\Desktop\PhD\1_entrainment\Manuscript',[SUBJ{s},'_topo_RFT']),'-dsvg','-r600')
print(f, fullfile(PATHPLOTS,[SUBJ{s},'_topo_RFT']),'-dpng','-r600')
