%% MMSP2 - Lab 4
%  Exercise 4 - Pitch detection and Linear Predictive Coding
%  Lucio Bianchi - 14/01/2014

clear
close all
clc

%% 1) Load the file 'voiced_a.wav' and consider frames of duration 25 ms
[y, Fs] = wavread('voiced_a.wav');

% Fs = 44100;
% freq = 174;
% t = 0:1/Fs:0.5;
% y = sin(2*pi*freq*t);

frame_length = floor(25*1e-3 * Fs);
frame_number = floor(length(y)/frame_length);

% Define frequency bounds for pitch detection
f_low = 60;
f_high = 500;

[ord, Wn] = buttord([f_low f_high]./(Fs/2), [f_low-25 f_high+25]./(Fs/2), 3, 40);
[num, den] = butter(ord, Wn);

pitch_zcr = zeros(frame_number,1);
pitch_zcr_f = zeros(frame_number,1);
pitch_acf = zeros(frame_number,1);
pitch_cepstrum = zeros(frame_number,1);
for nn = 1:frame_number
    frame = y((nn-1)*frame_length+1 : nn*frame_length);

    %% 1a) Detect the pitch using zero-crossing-rate on original signal
    pitch_zcr(nn) = sum(abs(diff(frame>0)))/frame_length * (Fs/2);
    
    %% 1aa) Improve pitch detection with zcr by bandpass filtering the input signal
    frame_f = filter(num,den,frame);
    pitch_zcr_f(nn) = sum(abs(diff(frame_f>0)))/frame_length * (Fs/2);
    
    %% 1b) Detect the pitch using autocorrelation on the original signal
    [acf,lag] = xcorr(frame);
    acf = acf(lag>=0);
    idx_f_low = floor(1/f_low*Fs);
    idx_f_high = floor(1/f_high*Fs);
    [~, idx_max_acf] = max(acf(idx_f_high:idx_f_low));
    idx_max_acf = idx_max_acf + idx_f_high;
    pitch_acf(nn) = 1/(idx_max_acf/Fs);
    
    
    %% 1c) Detect the pitch using Cepstrum on the orginal signal
    C = real(ifft(log(abs(fft(frame)))));
    [~, ind] = max(C(idx_f_high:idx_f_low));
    ind = ind + idx_f_high;
    pitch_cepstrum(nn) = 1/ind * Fs;
end

%% 1d) Plot the estimated pitch
figure(1)
plot(1:frame_number, pitch_zcr, 'b'), hold on, grid on
plot(1:frame_number, pitch_zcr_f, 'b--')
plot(1:frame_number, pitch_acf, 'r')
plot(1:frame_number, pitch_cepstrum, 'g'), hold off, grid off
xlabel('Frames')
ylabel('Pitch [Hz]')
legend('ZCR', 'ZCR (filtered)', 'ACF', 'Cepstrum');


%% 2) Linear prediction for each frame
p = 12; % prediction order
synth = []; % container for the synthetic speech signal
for n = 1:frame_number
    frame = y((nn-1)*frame_length+1 : nn*frame_length);
    
    %% 2a) Compute LP coefficients of order 12
    [r,lag] = xcorr(frame);
    r = r(lag >= 0);
    R = toeplitz(r(1:p));
    phi = r(2:p+1);
    a = R \ phi;
    
    %% 2b) Plot the prediction error and its magnitude spectrum
    if n == 1
        [e, Sf_e] = filter([1; -a], 1, frame);
    else
        [e, Sf_e] = filter([1; -a], 1, frame, Sf_e);
    end
    figure(2)
    subplot(211), plot(e, 'r'), title('Prediction Error'), xlabel('n');
    subplot(212), plot(abs(fft(e)), 'g'), title('Magnitude Spectrum of Prediction Error'), xlabel('f'), ylabel('abs(E(f))');
    
    %% 2c) Build an impulse train with the estimated pitch period
    x = zeros(frame_length, 1);
    ind = round(Fs / pitch_acf(n));
    x(1:ind:end) = 1;
    x = x / sqrt(sum(x.^2)); % energy normalization
    x = x * sqrt(sum(e.^2));
    
    
    %% 2d) Consider the impulse train as excitation and build synthetic speech
    if n == 1
        [frame_synth, Sf_x] = filter(1, [1; -a], x);
    else
        [frame_synth, Sf_x] = filter(1, [1; -a], x, Sf_x);
    end
    close all
    plot(frame,'r--'), hold on, plot(frame_synth,'b'), hold off
    synth = [synth; frame_synth];
end

%% 2e) Listen to the original and the synthetic speech
synth = 0.99 * synth/max(abs(synth));
sound(synth,Fs);