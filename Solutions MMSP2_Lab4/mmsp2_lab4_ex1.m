%% MMSP2 - Lab 4
%  Exercise 4 - Pitch detection and Linear Predictive Coding
%  Lucio Bianchi - 14/01/2014

clear
close all
clc

%% 1) Load the file 'voiced_a.wav' and consider frames of duration 25 ms


% Define frequency bounds for pitch detection

for nn = 1:frame_number
    
    %% 1a) Detect the pitch using zero-crossing-rate on original signal
    
    %% 1aa) Improve pitch detection with zcr by bandpass filtering the input signal
    
    %% 1b) Detect the pitch using autocorrelation on the original signal
    
    
    %% 1c) Detect the pitch using Cepstrum on the orginal signal
    
end

%% 1d) Plot the estimated pitch
% figure(1)
% plot(1:frame_number, pitch_zcr, 'b'), hold on, grid on
% plot(1:frame_number, pitch_zcr_f, 'b--')
% plot(1:frame_number, pitch_acf, 'r')
% plot(1:frame_number, pitch_cepstrum, 'g'), hold off, grid off
% xlabel('Frames')
% ylabel('Pitch [Hz]')
% legend('ZCR', 'ZCR (filtered)', 'ACF', 'Cepstrum');


%% 2) Linear prediction for each frame
for n = 1:frame_number
    
    
    %% 2a) Compute LP coefficients of order 12
    
    
    %% 2b) Plot the prediction error and its magnitude spectrum
    
    
    %% 2c) Build an impulse train with the estimated pitch period
    
    
    
    %% 2d) Consider the impulse train as excitation and build synthetic speech
    
end

%% 2e) Listen to the original and the synthetic speech
