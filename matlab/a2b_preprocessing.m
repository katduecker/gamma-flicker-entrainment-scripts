
%% Entrainment gamma rFT
% Masterthesis/PhD project
%
% preprocessing II b
% reject ICA components

% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

clear all; close all; clc;


MAINPATH = 'X:\';
addpath(fullfile(MAINPATH,'fieldtrip'))
ft_defaults;
BATCH = 3;
PATHIN = fullfile(MAINPATH,'subjects',['Batch_',num2str(BATCH)]);
PATHPRE = fullfile(MAINPATH,'results', 'preprocessing','filtered',['Batch_',num2str(BATCH)]);
PATHOUT = fullfile(MAINPATH,'results', 'preprocessing','ICA',['Batch_',num2str(BATCH)]);

% list subjects (Linux friendly)
folds = dir(PATHIN);
SUBJ = {};
for f = 1:length(folds)
    SUBJ = [SUBJ; folds(f).name];
end

% delete files that are not SUBJ data
SUBJ(1:2) = [];
%% Reject components
s = 28;
while s <=length(SUBJ)
    if exist(fullfile(PATHOUT, [SUBJ{s},'_badcomp.mat']))
        s = s + 1;
        continue
    else
        disp(['Loading subject',num2str(s)])
        load(fullfile(PATHOUT, [SUBJ{s},'_ICAcomp.mat']))
        cfg = [];
        cfg.channel = [1:10];
        cfg.continuous='no';
        cfg.viewmode = 'component';
        cfg.ylim = 'maxmin';
        cfg.layout = 'neuromag306mag.lay';
        cfg.compscale = 'local';
        ft_databrowser(cfg, dataICA);
        pause
        badComp = input('Bad components?:[]');
        
        save(fullfile(PATHOUT, [SUBJ{s},'_badcomp.mat']), 'badComp')
        close all; clear dataICA badComp
    end
end