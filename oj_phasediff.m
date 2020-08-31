%% Entrainment gamma RFT
% Phd Project 1

% Function computes phase similarity between photodiode and meg 
% signal using complex hanning taper

% INPUTS:
% s1, s2: meg and photodiode signal
% rftfoi, foi: frequencies of interest (here both are the same)
% Fs: sampling rate (here: 1000)
% make_plot: visualisation to double-check? (0,1)


% [c] O. Jensen & K. Duecker
% University of Birmingham, UK
%             Center for Human Brain Health

function s_diff = oj_phasediff(s1,s2,rftfoi,foi,Fs,make_plot)


width = 3;  % number of cycles for estimating phase
Nfoi = floor(width*Fs/foi);     % length of interval in samples

% complex hanning tapers of length 3 cycles
tap1 = hanning(Nfoi)'.*exp(i*2*pi*foi.*(1:Nfoi)/1000);
tap2 = hanning(Nfoi)'.*exp(i*2*pi*rftfoi.*(1:Nfoi)/1000);

% convolve & take positive frequencies
s1c = conv(s1,tap1);
s1c = s1c(ceil(Nfoi/2):length(s1c)-floor(Nfoi/2));
s2c = conv(s2,tap2);
s2c = s2c(ceil(Nfoi/2):length(s2c)-floor(Nfoi/2));

% phase difference: difference between unwrapped phase angles
s_diff = unwrap(angle(s1c))-unwrap(angle(s2c));

% visualisation
if make_plot
figure; 
subplot(211)
freqlum2 = fft(tap1,Nfoi)/(Nfoi/2);    % two sided spectrum
freqlum = freqlum2(1:Nfoi/2+1);            % single sided spectrum (positive frequencies)
freqvec = (1000/Nfoi)*(1:Nfoi/2+1);
%freqvec = linspace(1000/(N/2),1000/2,N/2+1);
plot(freqvec,real(freqlum).^2)

title(['spectrum wavelet IGF ',num2str(foi)])
subplot(212)
freqlum2 = fft(tap2,Nfoi)/(Nfoi);    % two sided spectrum
freqlum = freqlum2(1:Nfoi/2+1);            % single sided spectrum (positive frequencies)
freqvec = linspace(1000/(Nfoi/2),1000/2,Nfoi/2+1);
plot(freqvec,real(freqlum).^2)
title(['spectrum wavelet RFT ',num2str(rftfoi)])

end
end