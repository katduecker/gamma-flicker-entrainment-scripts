%% Entrainment gamma rFT
% Phd project 1
% 
%
% Reject trials identified in a_preprocessing
% reject noisy components
% replace broken sensor with neighbouring sensors
% separate into two conditions:
% entrainment: 7 sec, resonance: 5 sec

% note: resonance: flicker
% entrainment: flicker&gratings

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function c_separate_cond(s, BATCH)

% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
ft_defaults;
PATHRAW = fullfile(MAINPATH,'subjects',['Batch_',num2str(BATCH)]);                        % raw data
PATHFILT = fullfile(MAINPATH,'results','preprocessing','filtered',['Batch_',num2str(BATCH)]); 
PATHICA = fullfile(MAINPATH,'results','preprocessing','ICA',['Batch_',num2str(BATCH)]);     % path ICA components and 'bad components' structure
PATHOUT = fullfile(MAINPATH,'results','preprocessing','conditions',['Batch_',num2str(BATCH)]);

mkdir(PATHOUT)
folds = dir(PATHRAW);                                       % read in subjects

% list subjects
SUBJ = {};
for f = 1:length(folds)
    SUBJ = [SUBJ; folds(f).name];
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];


SUBJPATH = fullfile(PATHRAW, SUBJ{s});                  % subject folder

% all MEG files per subject
files = dir(fullfile(SUBJPATH, 'rawdata'));
for f = 1:length(files)
    testfiles{f} =files(f).name;
end
testfiles(1:2) = [];

% load grad structure
grad = {};
for t = 1:length(testfiles)
    grad = [grad; ft_read_sens(fullfile(SUBJPATH, 'rawdata', testfiles{t}))];
end
% load rejected trials
load(fullfile(PATHFILT, [SUBJ{s},'_trl_chan.mat']))

%% read in data
% read in resonance/flicker condition
for t = 1:length(testfiles)
    cfg = [];
    cfg.dataset = fullfile(SUBJPATH, 'rawdata', testfiles{t});
    cfg.trialdef.eventtype = 'Trigger';
    cfg.trialdef.eventvalue = 8;
    cfg.trialdef.prestim = 0;
    cfg.trialdef.poststim = 5;
    cfg = ft_definetrial(cfg);
    % remove trials with wrong trigger values = button presses
    if rejTRLRes{1}
        cfg.trl(rejTRLRes{t},:) = [];
    end
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

    cfg.detrend = 'yes';
    entrainment{t} = ft_preprocessing(cfg);
end

% average grad positions
mGrad = grad{1};
for g = 2:length(grad)
    mGrad.chanpos = mGrad.chanpos + grad{g}.chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);

ft_plot_sens(grad{1})
% append to one data set
data = ft_appenddata([],resonance{1:end}, entrainment{1:end});

data.grad = mGrad;

clear entrainment resonance
%% Reject trials & replace broken sensor

% reject trials but keep channels
keepTrl = true([1, length(data.time)]);
keepTrl(rejTrl) = 0;

cfg = [];
cfg.trials = keepTrl;
dataCL = ft_selectdata(cfg, data);

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

clear data
%% Reject ICA components

% load ICA components & bad components
load(fullfile(PATHICA,[SUBJ{s},'_ICAcomp.mat']))
load(fullfile(PATHICA,[SUBJ{s},'_badcomp.mat']))

cfg = [];
cfg.component = badComp;
cfg.demean = 'no';                          % don't demean
dataCLEAN = ft_rejectcomponent(cfg, dataICA, dataCL);
clear dataICA dataCL

%% Separate into conditions
% change timevector: stimulus onset = 0
timevec = cellfun(@(x) (x-1),dataCLEAN.time, 'UniformOutput', 0);
dataCLEAN.time = timevec;


trllgth = cell2mat(cellfun(@length, dataCLEAN.trial, 'UniformOutput', 0));

cfg = [];
cfg.trials = trllgth == 5000;
resonance = ft_selectdata(cfg, dataCLEAN);



mkdir(fullfile(PATHOUT, SUBJ{s}))
save(fullfile(PATHOUT, SUBJ{s}, 'resonance.mat'), 'resonance', '-v7.3');
clear resonance
cfg.trials = trllgth == 7000;
entrainment = ft_selectdata(cfg, dataCLEAN);

save(fullfile(PATHOUT, SUBJ{s}, 'entrainment.mat'), 'entrainment', '-v7.3');

end

