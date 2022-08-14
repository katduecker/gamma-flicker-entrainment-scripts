%% Entrainment gamma rFT
% PhD project 1
%
% function separates data into flicker frequency conditions (based on
% photodiode signal) 
% & fixes sample info

% [c] Katharina Duecker
%     k.duecker@bham.ac.uk
%     University of Birmingham, UK
%
% supervisor: O. Jensen
%             University of Birmingham, UK
%             Center for Human Brain Health

function [ COND, freqvec ] = kd_freqconditions( data, misc_chan, flickfreq )

% INPUT: data: preprocessed, epoched
% misc_chan: label of miscellaneous channel (=photodiode)
% flickfreq: flicker frequencies 

%% Frequency vector
% ffT over MISC

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.channel = misc_chan;
cfg.taper = 'hanning';
cfg.foi = 30:2:100;
cfg.tapsmofrq = 2;

% photodiode
MISC = ft_freqanalysis(cfg, data);

% find power peak and save frequencies in frequency vector
freqvec = [];
for p = 1:size(MISC.powspctrm,1)
    [pow pos] = max(squeeze(MISC.powspctrm(p,:,:)));
    freq = round(MISC.freq(pos));
    freqvec = [freqvec freq];
    close all
end


%% Frequency conditions

% categorize into frequency conditions
numfreq = [];                                     % check how many trials per frequency
COND = cell(1, length(flickfreq));             % average resonance to flicker
cfg = [];
cfg.avgoverrpt = 'no';                             % average over trial = avg over freq
for f = 1:length(flickfreq)
    findFRQ = logical(freqvec == flickfreq(f));        % find trials in which this frequency was used
    numfreq = [numfreq; numel(findFRQ)];
    cfg.trials = findFRQ;
    COND{f} = ft_selectdata(cfg, data);
end

% align sample info
for t = 1:length(COND)
    COND{t}.sampleinfo = repmat([1, size(COND{t}.trial{1},2)], size(COND{t}.trial,2),1);
    h = 0;
    for samp = 2:size(COND{t}.sampleinfo)
        h = h+1;
        COND{t}.sampleinfo(samp,:) = COND{t}.sampleinfo(samp,:) + size(COND{t}.trial{1},2)*h;
    end
end
end

