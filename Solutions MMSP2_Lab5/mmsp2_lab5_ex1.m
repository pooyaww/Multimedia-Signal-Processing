%% MMSP2 - Lab 5
%  Exercise 1 - Audio coding with DPCM and ADPCM
%  Lucio Bianchi - 23/01/2014

clear
close all
clc

%% 1) Load the file '60.flac', convert to mono and take the segment from 1s to 3s

%% 2) PCM coding. Use an uniform quantizer with R = 2:16 bit. Estimate the SNR for each R.
%%    Plot the PSD of the input signal and of the reconstruction error


%% 3) PCM coding and A-law companding. Repeat point a) and compare SNR. A = 87.7


%% 4) DPCM. Repeat point a) and compare SNR. Predictor is x_hat(n) = alpha*x_tilde(n-1), alpha=0.9.


%% 5) APCM. Use Jayant multipliers for R = 2,3,4 and midrise quantizer.


%% 6) ADPCM. Repeat point a) and compare SNR. Predictor is x_hat(n) = alpha*x_tilde(n-1), alpha=0.9.
%%    Use Jayant multipliers for R = 2,3,4 and midrise quantizer.
