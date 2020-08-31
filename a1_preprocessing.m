%% Entrainment gamma rFT
% PhD project 1
%
% preprocessing I
% - check trigger
% - store artefactual trials
% - check artefacts by hand and inspect data

% note: resonance: flicker
% entrainment: flicker&gratings

% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

clear all; close all; clc;

MAINPATH = 'X:\';
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;
addpath(fullfile(MAINPATH,'matlab'))
PATHIN = fullfile(MAINPATH,'subjects');
BATCH = 3;   
PATHIN = fullfile(PATHIN,['Batch_',num2str(BATCH)]);
PATHOUT = fullfile(MAINPATH,'results', 'preprocessing', 'filtered', ['Batch_',num2str(BATCH)]);
mkdir(PATHOUT)

% list subjects (Linux friendly)
folds = dir(PATHIN);
SUBJ = {};
for f = 1:length(folds)
    SUBJ = [SUBJ; folds(f).name];
end

% delete files that are not SUBJ data
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];
                                                % which subject batch?
s = 28;                                                      % start subject
while s <= length(SUBJ)
    SUBJPATH = fullfile(PATHIN, SUBJ{s});
    % all MEG files per subject
    files = dir(fullfile(SUBJPATH, 'rawdata'));
    for f = 1:length(files)
        testfiles{f} =files(f).name;
    end
    testfiles(1:2) = [];
    clear files
    
    % store all headers and events in a cell
    event = {};
    grad = {};
    % load events and grad structure
    for t = 1:length(testfiles)
        event = [event; ft_read_event(fullfile(SUBJPATH, 'rawdata', testfiles{t}))]; % read events
        grad = [grad; ft_read_sens(fullfile(SUBJPATH, 'rawdata', testfiles{t}))];
    end
    
    
    %% read in data
    
    %sanity check: how many trials per condition, have all triggers been sent?
    r = 0;
    TRIGR = cell(1,length(event));
    e = 0;
    TRIGE = cell(1,length(event));
    
    for n = 1:length(event)
        for b = 1:length(event{n})-6
            if strcmp(event{n}(b).type, 'Trigger') && event{n}(b).value == 8
                r = r +1;
                % resonance trigger
                TRIGR{n} = [TRIGR{n}; event{n}(b).value, event{n}(b+3).value, event{n}(b+6).value];
            elseif strcmp(event{n}(b).type, 'Trigger') && event{n}(b).value == 16
                e = n + 1;
                % entrainment trigger
                TRIGE{n} = [TRIGE{n}; event{n}(b).value, event{n}(b+3).value, event{n}(b+6).value, event{n}(b+9).value];
            end
        end
    end
    
    % identify trials with wrong triggers
    rejTRLEnt = cell(1,length(event));
    rejTRLRes = cell(1,length(event));
    % check!
    for t = 1:length(TRIGE)
        rejTRLEnt{t} = unique([find(TRIGE{t}(:,1)~= 16)', find(TRIGE{t}(:,2)~= 32)', find(TRIGE{t}(:,3)~= 128)', find(TRIGE{t}(:,4)~= 64)']);
        rejTRLRes{t} = unique([find(TRIGR{t}(:,1)~= 8)', find(TRIGR{t}(:,2)~= 128)', find(TRIGR{t}(:,3)~= 64)']);
    end
    clear event
    
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
        entrainment{t} = ft_preprocessing(cfg);
    end
    
    %% check diode
    
    t = resonance{1}.time{1}(find(resonance{1}.time{1} == 1) :find(resonance{1}.time{1} == 3));
    diode_res = resonance{1}.trial{1}(find(strcmp(resonance{1}.label,'MISC004')),find(resonance{1}.time{1} == 1) :find(resonance{1}.time{1} == 3));
    
    diode_ent = entrainment{1}.trial{1}(find(strcmp(entrainment{1}.label,'MISC004')),find(entrainment{1}.time{1} == 3) :find(entrainment{1}.time{1} == 5));
    
    plot(t,diode_res)
    hold on
    plot(t,diode_ent)
    clear diode_res diode_ent
    
    %% average grad positions
    mGrad = grad{1};
    for g = 2:length(grad)
        mGrad.chanpos = mGrad.chanpos + grad{g}.chanpos;
    end
    mGrad.chanpos = mGrad.chanpos./length(grad);
    
    % append to one data set
    data = ft_appenddata([],resonance{1:end}, entrainment{1:end});
    
    data.grad = mGrad;
    clear event hdr hdshp grad entrainment resonance 
    %% detrend
    
    cfg = [];
    cfg.detrend = 'yes';                                                   % remove linear trend
    dataF = ft_preprocessing(cfg, data);
    
    clear entrainment resonance TRIGE TRIGR data
    
    %% Artifact rejection
    disp(SUBJ{s})
    % Gradiometers
    cfg=[];
    cfg.channel  ='MEGGRAD';
    MEGgrad      = ft_selectdata(cfg,dataF);
    
    % visually reject using summary
    cfg=[];
    cfg.method  = 'summary';
    cfg.layout  = 'neuromag306planar.lay';
    MEGclgr  = ft_rejectvisual(cfg, MEGgrad);
    
    % broken sensors 
   
    rejTrl = input('noisy trials:');
    rejChan = input('Noisy channels');
    clear MEGgrad MEGclgr
    % magnetometers
    cfg=[];
    cfg.channel  ='MEGMAG';
    MEGmag      = ft_selectdata(cfg,dataF);
    pause
   % visually reject using summary
    cfg=[];
    cfg.method  = 'summary';
    cfg.layout  = 'neuromag306mag.lay';
    MEGmag  = ft_rejectvisual(cfg, MEGmag);
    clear MEGmag
    rejTrlmag = input('noisy trials:');
    rejTrl = unique([rejTrl, rejTrlmag]); 
    
    %% check noisy channels and trials by hand
    
    % align sample info
    sampleinfo = [1 size(dataF.trial{1},2)];
    for c = 2:size(dataF.trial,2)
        sampleinfo = [sampleinfo; [sampleinfo(c-1,1), sampleinfo(c-1,2)]+[size(dataF.trial{c-1},2),size(dataF.trial{c},2)]];
    end
    
    dataF.sampleinfo = sampleinfo;
    
    for f = 1:length(dataF.trial) dataF.trial{f}(find(strncmp(dataF.label,'ECG',3)),:) = dataF.trial{f}(find(strncmp(dataF.label,'ECG',3)),:) ./ 10^11; end
    % look at faulty trials first
       
    cfg = [];
    cfg.channel = 'MEGMAG';
   % cfg.channel = rejChan{1};
    %cfg.trl = dataF.sampleinfo(rejTrl,:);
    
    ft_databrowser(cfg,dataF)

    visrejTrl = input('Reject trials:[]');
    goodnorejTrl = input('Good trials, take off list:[]');
    rejChan = input('Noisy channels');
    rejTrl =  visrejTrl;%unique([rejTrl, visrejTrl]); 
    rejTrl(find(ismember(rejTrl,goodnorejTrl))) = []; 
  
    clear dataF MEGmag mGrad  rejTrlmag replChan keepTrl b e f g n r t visrejTrl goodnorejTrl
    save(fullfile(PATHOUT, [SUBJ{s},'_trl_chan.mat']), 'rejChan', 'rejTrl', 'rejTRLEnt', 'rejTRLRes')
    clear rejChan rejTrl  rejTRLEnt rejTRLRes rej testfiles

     prompt = 'Next subject? y/n';
     nextsubj = input(prompt,'s');
     if strcmp(nextsubj, 'y')
         s = s + 1;
     else
         s = length(SUBJ)+1;
     end
    
end