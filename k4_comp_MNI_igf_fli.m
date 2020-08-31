%% Entrainment gamma rFT
% PhD project 1: entrainment
% 
%
% Load in peak coordinates per subject and condition
% plot and prepare for statistical comparison in R

% [c] PGR: K. Duecker
%              k.duecker@bham.ac.uk
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Centre for Human Brain Health

%% settings
clear all; close all; clc
ticksize = 12;
labelsize = 16;
MAINPATH = 'Z:';
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

% read in coordinates
cs= [];
for s = 1:length(SUBJ)
    SUBJPATH = fullfile(PATHBEAM, SUBJ{s});
    % igf
    load(fullfile(SUBJPATH,['src_peak_igf.mat']),'peak_mni')
    ROI(s,1) = {peak_mni};
    clear peak_mni
    % flicker
    load(fullfile(SUBJPATH,['src_peak_fli.mat']),'peak_mni')
    ROI(s,2) = {peak_mni};
    clear peak_mni
    
    % store in cell
    if size(ROI{s,1},1)~=size(ROI{s,2},1)
        cs = [cs;s]
        if intersect(ROI{s,1},ROI{s,2},'row')
            ROI{s,1} = intersect(ROI{s,1},ROI{s,2});
            ROI{s,2} = ROI{s,1};
        elseif size(ROI{s,1},1)<size(ROI{s,2},1)
            ROI{s,2} = ROI{s,2}(1:size(ROI{s,1},1),:);
        elseif size(ROI{s,1},1)>size(ROI{s,2},1)
            ROI{s,1} = ROI{s,1}(1:size(ROI{s,2},1),:);
        end
    end
    
end

ROI_mat = cell2mat(ROI);
subjcol = [1:length(SUBJ)]';

% compute distance between coordinates per participant
DIST = ROI_mat(:,1:3) - ROI_mat(:,4:6);
% mean distance at x,y,z
DIR(1) = mean(DIST(:,1));
DIR(2) = mean(DIST(:,2));
DIR(3) = mean(DIST(:,3));

% standard error
STEIR = std(DIST(:,3))/sqrt(size(ROI_mat,1));

% store in matrix -> excel
ROI_mat = [subjcol,ROI_mat];
ROI = array2table(ROI_mat,'VariableNames',{'subj','x_IGF','y_IGF','z_IGF','x_FL','y_FL','z_FL'});
writetable(ROI,fullfile(PATHBEAM,'ROI_comp4mm.xlsx'))
save(fullfile(PATHBEAM,'ROI_comp4mm.mat'),'ROI_mat','DIR','DIST')

ROI_mat(:,1) = [];
ROI_igf = mean(ROI_mat(:,1:3),1);
ROI_rft = mean(ROI_mat(:,4:6),1);

%% plot coordinates as scatters
fig = figure;
set(gcf, 'Position', [0, 0, 1920/4.5, 1080/2.5],'renderer','painters')
scatter3(ROI_mat(:,1),ROI_mat(:,2),ROI_mat(:,3),50,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.25)
hold on
% mean scatter: larger
scatter3(mean(ROI_mat(:,1)),mean(ROI_mat(:,2)),mean(ROI_mat(:,3)),150,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
hold on
scatter3(ROI_mat(:,4),ROI_mat(:,5),ROI_mat(:,6),50,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.25)
hold on
% mean scatter: larger
scatter3(mean(ROI_mat(:,4)),mean(ROI_mat(:,5)),mean(ROI_mat(:,6)),150,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980])

xlabel('x')
ylabel('y')
zlabel('z')
xticks(-20:20:40)
yticks(-100:40:0)
zticks(-40:20:20)
pbaspect([2 3.33 2])
a = gca;
a.FontName = 'Arial';
a.FontSize = ticksize;
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
print(fig,fullfile(PATHPLOT,'roi_dist4mm_scatter'),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,'roi_dist4mm_scatter'),'-dpng','-r600')

% plot distance (only for z -> significant)
fig = figure;
set(gcf, 'Position', [0, 0, 1920/7, 1080/4],'renderer','painters')
% scatter(repmat(1,size(DIST(:,1),1),1),DIST(:,1),'filled','MarkerFaceColor',[0 00 0])
% hold on
% line([.75:.25:1.25],repmat(DIR(1),1,3), 'Color', 'black','LineWidth',2)
% hold on
% scatter(repmat(2,size(DIST(:,2),1),1),DIST(:,2),'filled','MarkerFaceColor',[0 00 0])
% hold on
% line([0.75:.25:1.25],repmat(DIR(2),1,3), 'Color', 'black','LineWidth',2)
% hold on
scatter(repmat(1,size(DIST(:,3),1),1),DIST(:,3),'filled','MarkerFaceColor',[0 00 0],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.25)
hold on
line([0.75:.25:1.25],repmat(DIR(3),1,3), 'Color', 'black','LineWidth',2)
hold on
line([0.75:.25:1.25],repmat(DIR(3)+STEIR,1,3), 'Color', 'black','LineWidth',2,'LineStyle','-.')
hold on
line([0.75:.25:1.25],repmat(DIR(3)-STEIR,1,3), 'Color', 'black','LineWidth',2,'LineStyle','-.')

xticks(1)
xticklabels({'z'})
ylabel('distance (mm)')
xlabel('coordinates')
 xlim([.5 1.5])
a = gca;
a.FontName = 'Arial';
a.FontSize = ticksize;
a.XLabel.FontSize = labelsize;
a.YLabel.FontSize = labelsize;
print(fig,fullfile(PATHPLOT,'roi_dist4mm_dist'),'-dsvg','-r600')
print(fig,fullfile(PATHPLOT,'roi_dist4mm_dist'),'-dpng','-r600')