%% Entrainment gamma rFT
% PhD project 1
% 
%
% Headmodel, sourcemodel, leadfield 

% subject without t1 scan: use template scan

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

% clear all; close all; clc;

MAINPATH = '/rds/projects/2018/jenseno-entrainment';
BATCH = 3;
addpath(fullfile(MAINPATH, 'matlab'))
addpath(fullfile('C:\Users\katha\Documents\MATLAB', 'fieldtrip'))
ft_defaults;
PATHVOL = fullfile(MAINPATH, 'results', 'anatomy - vol, leadf, sourcemodel',['Batch_',num2str(BATCH)]);
RAWPATH = fullfile(MAINPATH, 'subjects',['Batch_',num2str(BATCH)]);

s = 18;
% load templates
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
% load 8 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d4mm'));
[ftver, ftdir] = ft_version;
% standard MRI
templatedir = fullfile(ftdir, 'external', 'spm8', 'templates');
anatodir = fullfile(ftdir, 'template', 'anatomy');
templmri = ft_read_mri(fullfile(anatodir,'single_subj_T1.nii'));


% read in subjects - linux friendly
folds = dir(RAWPATH);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];


% get MEG files
megFiles = dir(fullfile(RAWPATH, SUBJ{s}, 'rawdata'));
megFiles(1:2) = [];

%% Realign headshape and MRI

% read headshape
digshape = ft_read_headshape(fullfile(RAWPATH, SUBJ{s}, 'rawdata', megFiles(1).name));
ft_determine_coordsys(digshape, 'interactive', 'no');
close all


% realign headshape and MRI
cfg = [];
cfg.method = 'headshape';
cfg.headshape = digshape;
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'yes';
cfg.coordsys = 'neuromag';
cfg.paramter = 'anatomy';
cfg.viewesult = 'yes';
mri_realign = ft_volumerealign(cfg, templmri);

cfg.headshape.interactive = 'yes';
mri_realign = ft_volumerealign(cfg, mri_realign);
% resulting realigned mri expressed in neuromag head coordinates
mri_realign.coordsys = 'neuromag';
    
cfg = [];
    cfg.dim = [256 256 256];
    mri_reslice = ft_volumereslice(cfg,mri_realign);


%     clear sens
    %% Headmodel
    
  % segment brain
    cfg = [];
    cfg.output = 'brain';
    segmmri = ft_volumesegment(cfg, mri_reslice);

    % prepare headmodel 
    cfg = [];
    cfg.method = 'singleshell';
    cfg.tissue = 'brain';
    vol = ft_prepare_headmodel(cfg,segmmri);
    vol = ft_convert_units(vol, 'cm');
%     
    % load in sensors
    for m = 1:length(megFiles)
        sens{m} = ft_read_sens(fullfile(RAWPATH, SUBJ{s}, 'rawdata', megFiles(m).name),'senstype','meg');
        sens{m} = ft_convert_units(sens{m}, 'cm');
    end
    
    % average sensor position over all data sets
    mGrad = sens{1};
    for m = 2:length(sens)
        mGrad.chanpos = mGrad.chanpos + sens{m}.chanpos;
    end
    mGrad.chanpos = mGrad.chanpos./length(sens);
    clear sens
    

    
    h =  figure;
    ft_plot_sens(mGrad, 'coilshape', 'square');
    hold on
    ft_plot_vol(vol,'edgecolor', 'brain','vertexcolor', 'none','facecolor','white');
    print(h,'ver_brain','-dpng','-r1000')
    hold on
    ft_plot_mesh(digshape.pos, 'vertexcolor', 'b', 'vertexsize', 8)
    pause
    close all
    
    cfg = [];
    cfg.method = 'singleshell';
    cfg.tissue = 'brain';
    vol = ft_prepare_headmodel(cfg,segmmri);
    vol = ft_convert_units(vol, 'cm');
    
    % location where fieldtrip is installed
    [ftver, ftdir] = ft_version;
    templatedir = fullfile(ftdir, 'template', 'sourcemodel');
    % load 4 mmm sourcemodel
    template = load(fullfile(templatedir, 'standard_sourcemodel3d4mm'));
    % warp template grid tp subject specific coordinates
    cfg = [];
  
    cfg.grid.warpmni    = 'yes';
    cfg.grid.template   = template.sourcemodel;
    cfg.grid.nonlinear  = 'yes';
    cfg.mri             = mri_reslice;   % aligned, resliced mri
    sourcemdl           = ft_prepare_sourcemodel(cfg);
    

    %% leadfield: use averaged sensor positions
    cfg = [];
    cfg.channel = ft_channelselection('meggrad',mGrad);
    cfg.grad = mGrad;
    cfg.grid = sourcemdl;
    cfg.headmodel = vol;
    leadf = ft_prepare_leadfield(cfg);

    save(fullfile(PATHVOL, SUBJ{s}, ['vol_lf_source_',sensortype,'4mm.mat']), 'vol', 'sourcemdl', 'mri_realign','mri_reslice','mGrad','leadf')
