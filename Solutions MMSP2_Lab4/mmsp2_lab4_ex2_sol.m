%% MMSP2 - Lab 4
%  Exercise 4 - Vocoder with voiced/unvoiced classification
%  Lucio Bianchi - 14/01/2014

clear
close all
clc

%% 1) Load the files 'a.wav' and 'shh.wav' and build a single signal
[a, Fs] = wavread('a.wav');
[shh] = wavread('shh.wav');

x = [a; shh];
clear a shh


%% 2) Filter the signal with a Butterworth passband filter
%%    pass band: 200 - 3400 Hz
%%    ripple in pass band: 4.5 dB
%%    attenuation in stop band: 10 dB
f_low = 200;
f_high = 3400;
[ord, Wn] = buttord([f_low f_high]./(Fs/2), [f_low-50 f_high+150]./(Fs/2), 4.5, 10);
[num, den] = butter(ord, Wn);
x_f = filter(num,den,x);


%% 3) Frame selection and windowing
%%    length = 40 ms
%%    spacing = 10 ms
frame_length = floor(40e-3 * Fs);
frame_spacing = floor(10e-3 * Fs);
N = floor((length(x_f) - frame_length)/frame_spacing) + 1; % number of frames

% generate Hamming window
win = hamming(frame_length);

% container for classification parameters
parameter = zeros(N, 3);

for n=1:N
    frame = x_f((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length); % select the frame
    frame = frame .* win; % apply Hamming window
    
    %% 4) Parameter evaluation
    % Cepstrum intensity between 100 Hz and 600 Hz
    f_min = 100;
    f_max = 600;
    ind_min_freq = floor(1 / f_min * Fs);
    ind_max_freq = floor(1 / f_max * Fs);
    C = real(ifft(log(abs(fft(frame)))));
    parameter(n,1) = max(C(ind_max_freq:ind_min_freq));
    
    % Zero-crossing rate
    parameter(n,2) = sum(abs(diff(frame>0)))/frame_length;
    
    % Short-time energy
    parameter(n,3) = sum(frame.^2);
end

%% 5) Voiced / Unvoiced classification
% Compute thresholds
tau_cep = median(parameter(:,1));
tau_zcr = median(parameter(:,2));
tau_ste = median(parameter(:,3));

% Decision
voiced = zeros(N, 1); % voiced(i) = 1 if the i-th frame is voiced, 0 if unvoiced
for n = 1:N
    if parameter(n,1) > tau_cep && parameter(n,2) < tau_zcr && parameter(n,3) > tau_ste
        voiced(n) = 1;
    end
end

% plot the parameters versus the thresholds
threshold = zeros(N, 3);
threshold(:, 1) = tau_cep;
threshold(:, 2) = tau_zcr;
threshold(:, 3) = tau_ste;

subplot(411), plot(parameter(:,1)), title('Cepstrum Energy'), xlabel('Frame number'), ylabel('Peak'), hold on, plot(threshold(:, 1), 'r');
subplot(412), plot(parameter(:,2)), title('Zero Crossing Rate'), xlabel('Frame number'), ylabel('ZCR'), hold on, plot(threshold(:, 2), 'r');
subplot(413), plot(parameter(:,3)), title('Short Time energy'), xlabel('Frame number'), ylabel('STE'), hold on, plot(threshold(:, 3), 'r');
subplot(414), bar(voiced), title('Voiced/Unvoiced Classification'), xlabel('Frame number'), ylabel('V/UV');


%% 6) LPC analysis and synthesis
p = 12; % order of the predictor

% container for the synthesized speech
synth = zeros(length(x_f),1);
for n = 1:N
    frame = x_f((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length); % select the frame
    frame = frame .* win; % apply Hamming window
    
    %% 6a) Compute LP coefficients and prediction error
    [r,lag] = xcorr(frame);
    r = r(lag >= 0);
    R = toeplitz(r(1:p));
    phi = r(2:p+1);
    a = R \ phi;
    if n == 1
        [e, Sf_e] = filter([1; -a], 1, frame);
    else
        [e, Sf_e] = filter([1; -a], 1, frame, Sf_e);
    end
    
    %% 6b-1) Voiced segment:
    if voiced(n) == 1
        % Pitch detection
        f_min = 100;
        f_max = 600;
        ind_min_freq = floor(1 / f_min * Fs);
        ind_max_freq = floor(1 / f_max * Fs);
        [~, ind] = max(r(ind_max_freq:ind_min_freq));
        ind = ind + 1;
        
        % Generate impulse train
        x = zeros(frame_length, 1);
        x(1:ind:end) = 1;

    %% 6b-2) Unvoiced segment:
    else
        % Generate random noise
        x = randn(frame_length, 1);
    end
    
    %% 6c) Normalize the energy of the excitation signal
    x = x / sqrt(sum(x.^2));
    x = x * sqrt(sum(e.^2));
    
    %% 7) Shaping filter
    frame_synth = filter(1, [1; -a], x);
    aux = zeros(length(x_f),1);
    aux((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length) = frame_synth;    
    synth = synth + aux;
end

synth = 0.99 * synth / max(abs(synth));