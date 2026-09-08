function [ output_args ] = Spectral_Entropy( DSP,Freq_ax );
%Summary of this function goes here
% For Spectral Entropy the steps are:
% 
% 1 Calculate the power spectrum of the signal using FFT command in MATLAB.
% 2 Calculate the Power Spectral Density using the power spectrum or using any other technique.
% 3 Normalize the Power Spectral Density between [0-1], so that it can be treated as a probability density function  pi.
% 3 Calculate the Entropy H(s) = ??(pi)log2(pi)
% Calculate the Entropy sum(P)log2(P), where P = PSD
% 
%   Detailed explanation goes here

fin=max( find(Freq_ax<=30));
inti=1;

DSP = DSP(inti:fin); Freq_ax=Freq_ax(inti:fin);

%    P=sum(abs(fft(data-window)).^2)
   %Normalization
   d=DSP(:);
   d=d/sum(d+ 1e-12);

   %Entropy Calculation
   logd = log2(d + 1e-12);
   
%    Thereafter, the entropy value is normalized 
%    to range between 1 (maximum irregularity) and 0 (complete regularity). 
%    The value is divided by the factor log (N[f1,f2]) 
%    where N[f1,f2] is equal to the total number of frequency components in the range [f1,f2]:

   Entropy = -sum(d.*logd)/log2(length(d));
%    Entropy = -sum(d.*logd);

output_args = Entropy;

end

