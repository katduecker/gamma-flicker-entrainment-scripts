%% Entrainment gamma rFT
% PhD project 1: entrainment
% 
%
% Find peak source of oscillatory activity

% Inputs
% s: subject id
% condit: condition (igf = gamma oscillation; fli: flicker response)

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

function k4_lcmv_findpeak(s,condition)

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

SUBJPATH = fullfile(PATHBEAM, SUBJ{s});

%% Interpolate

% load spatial filter & filtered data
load(fullfile(SUBJPATH,'spatfilt.mat'))
load(fullfile(SUBJPATH,['lcmv_',condition,'.mat']))

% load templates
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
% load 8 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d4mm'));
[ftver, ftdir] = ft_version;
% standard MRI
templatedir = fullfile(ftdir, 'external', 'spm8', 'templates');
anatodir = fullfile(ftdir, 'template', 'anatomy');
templmri = ft_read_mri(fullfile(templatedir,'T1.nii'));

% fill in power change values at right virtual channels
srcRC = sourcetime_data;
srcRC.avg.pow(relch(:,1)) = relch(:,2);

% interpolate
srcRC.dim = template.sourcemodel.dim;
srcRC.pos = template.sourcemodel.pos;
cfg              = [];
cfg.interpmethod = 'linear';
cfg.parameter    = 'pow';
srcint = ft_sourceinterpolate(cfg, srcRC, templmri);
clear srcRC

% % plot to double check
% cfg = [];
% cfg.method = 'ortho';
% cfg.nslices = 1;
% 
% cfg.crosshair = 'no';
% cfg.slicedim = 3;
% cfg.funparameter = 'pow';
% cfg.funcolorlim = 'maxabs';
% cfg.maskparameter = 'pow';
% %cfg.atlas = aalAtlas;
% cfg.funcolormap = 'jet';
% ft_sourceplot(cfg,srcint)

% save power (for grandaverage later)
srcrc = srcint.pow;
% find peak
peak_mni = srcint.pos(find(srcint.pow == max(srcint.pow)),:);

save(fullfile(SUBJPATH,['src_peak_',condition,'.mat']),'srcrc','peak_mni')