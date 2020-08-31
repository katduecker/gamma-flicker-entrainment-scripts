%% Entrainment gamma RFT
% phase slips

% script loads in preprocessed data, separated into conditions
% computes phase similarity between photodiode and meg signal using
% function in oj_phasediff.m (with hanning taper, method used in
% manuscript)
% (also computes shannon entropy, not shown in manus due to redundancy with
% plv)

%
% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health
function n_phase_slips(BATCH,condit,Nbins,method,filterord)

% INPUT: BATCH (sample number)
% condit: 1: flicker&gratings or 2: flicker
% Nbins: number of bins to estimate probability density function using histogram
% method: hann = hanning taper of width 3 (cycles) ...
%         condit 1: 3 cyc of IGF; condit 2: 3 cyc of RFT
% or -   hilb = hilbert transform
% filterord: only define for hilb


%% settings
MAINPATH = '/rds/projects/2018/jenseno-entrainment';
addstr = {'entrainment', 'resonance'};
MATPATH = fullfile(MAINPATH,'matlab');
addpath(fullfile(MATPATH,'narrowtag'));
addpath(fullfile(MATPATH,'kd fun'));
addpath(fullfile(MATPATH,'template function'));
addpath(fullfile(MAINPATH,'fieldtrip'));
ft_defaults;

PATHCOND = fullfile(MAINPATH, 'results','preprocessing','conditions',['Batch_',num2str(BATCH)]);
PATHIN = fullfile(MAINPATH,'results','power',addstr{condit},['Batch_',num2str(BATCH)]);
PATHGAM = fullfile(MAINPATH,'results','power','gammatron',['Batch_',num2str(BATCH)]);
PATHPLOTS = fullfile(MAINPATH,'results','plots','arnold tongue',['Batch_',num2str(BATCH)],'phase slips',addstr{condit});
mkdir(PATHPLOTS)
PATHOUT = fullfile(MAINPATH,'results','phase slips'); 
PATHOUT = fullfile(PATHOUT,method,addstr{condit});

mkdir(PATHOUT)
folds = dir(PATHGAM);
for f = 1:length(folds)
    SUBJ{f} = folds(f).name;
end
SUBJ(find(~strncmp(SUBJ,'201',3))) = [];

for s = 1:length(SUBJ)
    load(fullfile(PATHGAM, SUBJ{s}, 'SOI_freq.mat'))                        % load individual Gamfreq
    IGF(s) = gamFreq;
    SOI_all{s} = SOI;
    numSens(s) = numel(SOI);
    clear SOI gamFreq
end


%test: keep subjects whose IGF is > 56
SUBJ = SUBJ(IGF > 56);
SOI_all = SOI_all(IGF > 56);
IGF = IGF(IGF > 56);
%find min and max IGF
minIGF = min(IGF);
maxIGF = max(IGF);

% upper and lower IGF-+ limit (only for hilbert)
leftLim = 6;
rightLim = 6;

% separate combined SOI into uncombined planar
for h = 1:length(SOI_all)
    w = 1;              % help index
    for l = 1:length(SOI_all{h})
        SOInc{h}{w} = SOI_all{h}{l}(1:strfind(SOI_all{h}{l},'+')-1);         % first sensor
        w = w + 1;
        SOInc{h}{w} = ['MEG',SOI_all{h}{l}(strfind(SOI_all{h}{l},'+')+1:end)];
        w = w + 1;
    end
end


flickfreq = 52:2:90;

n = 10;
m = 0;
ShEnt_IGF = [];                     % dummy shannon entropy (not shown in manus)

for s = 1:length(SUBJ)

disp(['Loading in subject ',num2str(s),' - ',SUBJ{s}])
load(fullfile(PATHCOND, SUBJ{s},[addstr{condit},'.mat']))

%% Frequency conditions

% separate into freq conditons & adjust sample info

if condit == 1
    cfg = [];
    cfg.latency = [2 4];
    entrainment = ft_selectdata(cfg,entrainment);
    [COND, ~] = kd_freqconditions(entrainment, 'MISC004', flickfreq);
elseif condit==2
    cfg = [];
    cfg.latency = [0 2];
    resonance = ft_selectdata(cfg,resonance);
    [COND, ~] = kd_freqconditions(resonance, 'MISC004', flickfreq);
end
if s ==1
    grad = COND{1}.grad;
end

if strcmp(method,'hilb')
    freqvec = IGF(s)-leftLim:2:IGF(s)+leftLim;

    COND = COND(ismember(flickfreq,freqvec));
end
dioind = find(strcmp(COND{1}.label,'MISC004'));

for soi = 1:length(SOInc{s})
    soiind = find(strcmp(COND{1}.label,SOInc{s}{soi}));

    for c = 1:size(COND,2)
        for tr = 1:length(COND{c}.trial)
            
            megsig = COND{c}.trial{tr}(soiind,:);
            diodesig = COND{c}.trial{tr}(dioind,:);
            
            
            % compute phase similarity
            if strcmp(method,'hann')                % hanning taper
                s_diff{soi}{c}{tr} = oj_phasediff(megsig,diodesig,flickfreq(c),flickfreq(c),COND{c}.fsample,0);
                
            elseif strcmp(method,'hilb')        % hilbert transform
                s_diff{soi}{c}{tr} = kd_phasediff(megsig,diodesig,COND{c}.fsample,filterord,[IGF(s)-8 IGF(s)+8]);
            end
            
            mod_s_diff{soi}{c}{tr} = mod(s_diff{soi}{c}{tr},(2*pi));                  % cyclic relative phase
            
            %  Shannon entropy
            % get probability density function
            h = histcounts(mod_s_diff{soi}{c}{tr},Nbins, 'Normalization','probability');         % number of bins length signal/length of cycle
            %
            h(h==0) = eps;
            H = -dot(h,log(h));                     % entropy of probability distribution of mod_s_diff
            ShEnt(c,s,soi,tr) = (log(Nbins)-H)/log(Nbins);
            
        end
    end
    
    
    
end

delete(fullfile(PATHOUT,[SUBJ{s},'mod_s_diff.mat']))
delete(fullfile(PATHOUT,[SUBJ{s},'s_diff.mat']))
delete(fullfile(PATHOUT,[SUBJ{s},'mod_s_diff_new.mat']))

save(fullfile(PATHOUT,[SUBJ{s},'mod_s_diff.mat']),'mod_s_diff','s_diff')

clear COND s_diff mod_s_diff freqvec



end

if strcmp(method,'hann')
    % transform matrix to be fed into ANOVA
    ShEnt_matr = [];
    for f = 1:length(flickfreq)
        for s = 1:length(SUBJ)
            for soi = length(SOInc{s})
                ShEnt_matr = [ShEnt_matr;repmat(flickfreq(f),1,size(ShEnt(f,s,soi,:),4))',...
                    repmat(s,1,size(ShEnt(f,s,soi,:),4))',...
                    repmat(find(strcmp(grad.label,SOInc{s}{soi})),1,size(ShEnt(f,s,soi,:),4))',...
                    squeeze(ShEnt(f,s,soi,:))];
            end
        end
    end
    
    [r,c] = find(isnan(ShEnt_matr));
    
    ShEnt_matr(r,:) = [];
    
    for s = 1:length(SUBJ)
        ShEnt_IGF(:,s,:,:) = ShEnt(flickfreq >= IGF(s)-leftLim&flickfreq <= IGF(s)+rightLim,s,:,:);
    end
elseif strcmp(method,'hilb')
    ShEnt_IGF = ShEnt;
    clear ShEnt
end

freqvec = -leftLim:2:leftLim;
ShEnt_IGF_matr = [];
for f = 1:length(freqvec)
    for s = 1:length(SUBJ)
        for soi = length(SOInc{s})
            ShEnt_IGF_matr = [ShEnt_IGF_matr;repmat(freqvec(f),1,size(ShEnt_IGF(f,s,soi,:),4))',...
                repmat(s,1,size(ShEnt_IGF(f,s,soi,:),4))',...
                repmat(find(strcmp(grad.label,SOInc{s}{soi})),1,size(ShEnt_IGF(f,s,soi,:),4))',...
                squeeze(ShEnt_IGF(f,s,soi,:))];
        end
    end
end

[r,c] = find(isnan(ShEnt_IGF_matr));
ShEnt_IGF_matr(r,:) = [];
if exist(fullfile(PATHOUT,['shannon_ent_',addstr{condit},'_n=',num2str(n),'_m=',num2str(m),'Nbins_',num2str(Nbins),'.mat']))
    delete(fullfile(PATHOUT,['shannon_ent_',addstr{condit},'_n=',num2str(n),'_m=',num2str(m),'Nbins_',num2str(Nbins),'.mat']))
end
if strcmp(method,'hann')
    save(fullfile(PATHOUT,['shannon_ent.mat']),'ShEnt','ShEnt_IGF','ShEnt_matr','ShEnt_IGF_matr')
elseif strcmp(method,'hilb')
    save(fullfile(PATHOUT,['shannon_ent.mat']),'ShEnt_IGF','ShEnt_IGF_matr')
    
end

