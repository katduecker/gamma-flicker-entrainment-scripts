%% Entrainment gamma rFT
% PhD project 1
%
% preprocessing II
% - ICA & store artefactual components

% note: resonance: flicker
% entrainment: flicker&gratings

% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

function a2_preprocessing(s, BATCH)
%INPUT: s: subject ID, BATCH: sample (3: experiment, 1-2: PILOT)

ft_defaults;
MAINPATH = '/rds/projects/2018/jenseno-entrainment/';
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
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

SUBJPATH = fullfile(PATHIN, SUBJ{s});
% all MEG files per subject
files = dir(fullfile(SUBJPATH, 'rawdata'));
for f = 1:length(files)
    testfiles{f} =files(f).name;
end
testfiles(1:2) = [];
clear files

load(fullfile(PATHPRE, [SUBJ{s},'_trl_chan.mat']))


%% read in data
% read in resonance condition
for t = 1:length(testfiles)
    cfg = [];
    cfg.dataset = fullfile(SUBJPATH, 'rawdata', testfiles{t});
    cfg.trialdef.eventtype = 'Trigger';
    cfg.trialdef.eventvalue = 8;
    cfg.trialdef.prestim = 0;
    cfg.trialdef.poststim = 5;
    cfg = ft_definetrial(cfg);
    % remove trials with wrong trigger values
    if rejTRLRes{1}
        cfg.trl(rejTRLRes{t},:) = [];
    end
    %         cfg.demean = 'yes';
    %         cfg.baselinewindow = [0.25 0.75];
    cfg.detrend = 'yes';
    resonance{t} = ft_preprocessing(cfg);
end

for t = 1:length(testfiles)
    % read in entrainment
    cfg = [];
    cfg.dataset = fullfile(SUBJPATH, 'rawdata', testfiles{t});
    cfg.trialdef.eventtype = 'Trigger';
    cfg.trialdef.eventvalue = 16;
    cfg.trialdef.prestim = 0;
    cfg.trialdef.poststim = 7;
    cfg = ft_definetrial(cfg);
    if rejTRLEnt{1}
        cfg.trl(rejTRLEnt{t},:) = [];
    end
    %         cfg.demean = 'yes';
    %         cfg.baselinewindow = [0.25 0.75];
    cfg.detrend = 'yes';
    entrainment{t} = ft_preprocessing(cfg);
end


% append to one data set
data = ft_appenddata([],resonance{1:end}, entrainment{1:end});
clear entrainment resonance


% load grad structure
grad = {};
for t = 1:length(testfiles)
    grad = [grad; ft_read_sens(fullfile(SUBJPATH, 'rawdata', testfiles{t}))];
end


% average grad positions
mGrad = grad{1};
for g = 2:length(grad)
    mGrad.chanpos = mGrad.chanpos + grad{g}.chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
data.grad = mGrad;

clear grad
%% Reject identified artifactual trials

keepTrl = true([1, length(data.time)]);
keepTrl(rejTrl) = 0;

cfg = [];
cfg.trials = keepTrl;
dataCL = ft_selectdata(cfg, data);
clear data
%% Fix sample info

sampleinfo = [1 size(dataCL.trial{1},2)];
for c = 2:size(dataCL.trial,2)
    sampleinfo = [sampleinfo; [sampleinfo(c-1,1), sampleinfo(c-1,2)]+[size(dataCL.trial{c-1},2),size(dataCL.trial{c},2)]];
end

dataCL.sampleinfo = sampleinfo;

% hack: replace sensors with neighbouring sensors
if ~isempty(rejChan)
    for r = 1:length(rejChan)
        if regexp(rejChan{r}, '3')
            replChan{r} = dataCL.label(find(strcmp(rejChan{r}, dataCL.label))-1);
        elseif regexp(rejChan{r}, '2')
            replChan{r} = dataCL.label(find(strcmp(rejChan{r}, dataCL.label))+1);
        elseif regexp(rejChan{r}, '1')
            replChan{r} = dataCL.label(find(strcmp(rejChan{r}, dataCL.label))+10);
        end
    end
    
    for t = 1:length(dataCL.trial)
        for r = 1:length(replChan)
            dataCL.trial{t}(find(strcmp(rejChan{r}, dataCL.label)),:) = dataCL.trial{t}(find(strcmp(replChan{r}, dataCL.label)),:);
        end
    end
    
end
%% ICA
cfg =[];
cfg.method = 'runica';
cfg.numcomponent = 100;
cfg.channel = 'MEG';
cfg.topolabel = dataCL.label(strncmp(dataCL.label, 'MEG',3));
cfg.runica.maxsteps = 100;
dataICA = ft_componentanalysis(cfg, dataCL);

save(fullfile(PATHOUT, [SUBJ{s},'_ICAcomp.mat']), 'dataICA','-v7.3')
end

%


