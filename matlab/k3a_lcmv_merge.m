%% Entrainment gamma rFT
% PhD project 1: entrainment
% 
%
% Concatenate virtual channels to one data set

% Inputs
% s: subject id
% condit: condition (gam = gamma oscillation; fli: flicker response)

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function k4_lcmv_merge(s,condition)

%% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile(MAINPATH, 'matlab','narrowtag'))
addpath(fullfile(MAINPATH, 'matlab','kd fun'))
addpath(fullfile(MAINPATH,'fieldtrip'));

PATHBEAM = fullfile(MAINPATH, 'results', 'beamformer MNI','LCMV');
ft_defaults;

RAWPATH = fullfile(MAINPATH, 'subjects','Batch_3');
PATHVOL = fullfile(MAINPATH, 'results', 'anatomy - vol, leadf, sourcemodel','Batch_3');
PATHCOND = fullfile(MAINPATH, 'results', 'preprocessing', 'conditions','Batch_3');
PATHGAM = fullfile(MAINPATH, 'results','power','gammatron','Batch_3');
PATHENT = fullfile(MAINPATH, 'results','power','entrainment','Batch_3');
PATHRES = fullfile(MAINPATH, 'results','power','resonance','Batch_3');
PATHPLOT = fullfile(MAINPATH, 'results','plots','beamformer');

% read in subj
folds = dir(PATHBEAM);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

% for s = 1:length(SUBJ)
SUBJPATH = fullfile(PATHBEAM, SUBJ{s});

%% Merge power spectra

% list files
folds = dir(SUBJPATH);
for f = 1:length(folds)
    files{f} = folds(f).name;
end
bslfiles = files(find(strncmp(files,['bsl',condition],6)));
stimfiles = files(find(strncmp(files,condition,3)));
% for b = 1:length(bslfiles)
% delete(fullfile(SUBJPATH,bslfiles{b}))
% 
% delete(fullfile(SUBJPATH,stimfiles{b}))
% end
% end
% 
% bslint = [];
% for f = 1:length(bslfiles)
%     beg = strfind(bslfiles{f},'[')+1;
%     mid = strfind(bslfiles{f},' ');
%     strend = strfind(bslfiles{f},']')-1;
%     bslint = [bslint;str2num(bslfiles{f}(beg:mid(2)-1)),str2num(bslfiles{f}(mid(2)+1:strend))];
% end

% load in power values and channel indices
bsl = [];
stim = [];

wb = [];
wg = [];
for f = 1:length(bslfiles)
    load(fullfile(SUBJPATH,bslfiles{f}))
    % if only one channel in data set, something went wrong (only for
    % debugging)
    if length(bslpow) == 1
        wb = [wb;f];
    else
        bsl = [bsl;vc',bslpow];
    end
    clear bslpow vc
    load(fullfile(SUBJPATH,stimfiles{f}))
    % weird stuff?
    if length(stimpow) == 1
        wg = [wg;f];
    else
        stim = [stim;vc',stimpow];
    end

    clear stimpow vc
end

% sort based on channel idx 
[B,I] = sort(bsl);
bsl= bsl(I(:,1),:);
[B,I] = sort(stim);
stim = stim(I(:,1),:);

if strcmp(condition, 'gam')
    condition = 'igf';
end
% compute relative power change
% for IGF
if strcmp(condition,'igf')
    relch = stim;
    relch(:,2) = stim(:,2)./bsl(:,2)-1;
    % for flicker
elseif strcmp(condition,'fli')
    rcprep = stim;
    rcprep = stim(:,2:end)./bsl(:,2:end)-1;
    relch(:,1) = stim(:,1);
    relch(:,2) = mean(rcprep,2);
end

save(fullfile(SUBJPATH,['lcmv_',condition,'.mat']),'bsl','stim','relch')

for f = 1:length(bslfiles)
    delete(fullfile(SUBJPATH,bslfiles{f}))
    delete(fullfile(SUBJPATH,stimfiles{f}))
end


