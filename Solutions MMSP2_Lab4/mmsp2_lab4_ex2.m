%% MMSP2 - Lab 4
%  Exercise 4 - Vocoder with voiced/unvoiced classification
%  Lucio Bianchi - 14/01/2014

clear
close all
clc

%% 1) Load the files 'a.wav' and 'shh.wav' and build a single signal


%% 2) Filter the signal with a Butterworth passband filter
%%    pass band: 200 - 3400 Hz
%%    ripple in pass band: 4.5 dB
%%    attenuation in stop band: 10 dB



%% 3) Frame selection and windowing
%%    length = 40 ms
%%    spacing = 10 ms


for n=1:N
    
    
    %% 4) Parameter evaluation
    % Cepstrum intensity below 600 Hz
    
    
    % Zero-crossing rate
    
    
    % Short-time energy
    
end

%% 5) Voiced / Unvoiced classification
% Compute thresholds


% Decision


% plot the parameters versus the thresholds


%% 6) LPC analysis and synthesis


% container for the synthesized speech

for n = 1:N
    
    
    %% 6a) Compute LP coefficients and prediction error
    
    
    %% 6b-1) Voiced segment:
    if voiced(n) == 1
        % Pitch detection
        
        
        % Generate impulse train
        

    %% 6b-2) Unvoiced segment:
    else
        % Generate random noise
        
    end
    
    %% 6c) Normalize the energy of the excitation signal

    
    %% 7) Shaping filter

end

