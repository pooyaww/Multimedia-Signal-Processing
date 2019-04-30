%% MMSP2 - Lab 5
%  Exercise 2 - MDCT
%  Lucio Bianchi - 23/01/2014

clear
close all
clc

%% 1) Load the file '60.flac', convert to mono and take the segment from 1s to 3s
[x_stereo, Fs] = audioread('60.flac');
dur = 3 ;
x_left = x_stereo(:,1);
x_right = x_stereo(:,2);
x = (x_left+x_right)/2;
x = x(Fs:dur*Fs);

%% 2) Generate a sin window of length N = 1024 samples
N = 1024;
n = (1:N)';
s = sin((n-1/2)*pi/N);

%% 3) Zero-pad the signal
N_left_pad = N/2;
N_right_pad = N/2 * ceil(length(x)/(N/2)) - length(x) + N;
x_pad = zeros(N_left_pad+length(x)+N_right_pad, 1);
x_pad(N_left_pad+1:N_left_pad+length(x)) = x;

%% 4) Apply MDCT to the zero-padded signal
n0 = (N/2+1)/2;

y_pad = zeros(length(x_pad),1);
for ii = 1:N/2:(length(x_pad)-N+2)
    % Consider the current frame and window it
    frame = x_pad(ii:ii+N-1);
    frame = frame .* s;
    
    %% 4a) Compute MDCT using FFT
    % Pre-twiddle
    n = (0:N-1)';
    frame_tilde = frame .* exp(-1i*pi*n/N);
    
    % Compute the FFT
    FRAME_tilde = fft(frame_tilde,N);
    
    % Post-twiddle (take only the first half of the frequency domain signal)
    k = (0:N/2-1)';
    FRAME = real(FRAME_tilde(1:N/2) .* exp(-1i*2*pi/N*n0*(k+1/2)));
    
    %% 4b) Compute inverse MDCT using FFT
    % Pre-twiddle
    k = (0:N-1)';
    FRAME_hat = [FRAME; -FRAME(end:-1:1)] .* exp(1i*2*pi/N*k*n0);
    
    % Compute inverse FFT (hint: expolit odd symmetry)
    frame_hat = ifft(FRAME_hat,N);
    
    % Post-twiddle
    frame_bar = real(frame_hat .* exp(1i*pi/N*(n+n0)));    
    
    %% 4c) Multiply by two times the synthesis window and perform overlap-and-add
    frame = 2*s .* frame_bar;
    y_pad(ii:ii+N-1) = y_pad(ii:ii+N-1) + frame;
    
end

%% 5) Compare the reconstructed signal with the original one
e = y_pad-x_pad;
figure
plot(e), title('MDCT (fixed window) - Reconstruction error')
xlabel('Time'), ylabel('Error')
SNR_mdct = db(var(x_pad)/var(x_pad-y_pad));
disp(['MDCT with fixed window - Signal-to-Noise Ratio = ' num2str(SNR_mdct) ' dB'])

%% 6) MDCT with window switching

%% 6a) Define 4 types of windows
% Long window, length = 1024
L = 2048;
n = (1:L)';
long_win = sin((n-1/2)*pi/L);


% Short window, length = 512
S = 1024;
n = (1:S)';
short_win = sin((n-1/2)*pi/S);

% Long-to-short and short-to-long windows
% hint: Edler method
long2short_win = zeros(L,1);
long2short_win(1:L/2) = long_win(1:L/2);
long2short_win(L/2+1 : 3*L/4-S/4) = ones(L/4-S/4, 1);
long2short_win(3*L/4-S/4+1 : 3*L/4+S/4)=short_win(S/2+1:end);

short2long_win = flipud(long2short_win);

%% 6b) Frame switching algorithm
% Compute ZCR for the first frame (HP: long frame)
zcr = sum(abs(diff(x(1:L)>0)))/L;

% First window type selection
if (zcr < 0.3)
    win = long_win;
    win_type = 'long';
    s = long_win;
    N = L;
else
    win = short_win;
    win_type = 'short';
    s = short_win;
    N = S;
end

% Zero-pad the signal
N_left_pad = N/2;
N_right_pad = L/2 * ceil(length(x)/(L/2)) - length(x) + L;
x_pad = zeros(N_left_pad+length(x)+N_right_pad, 1);
x_pad(N_left_pad+1:N_left_pad+length(x)) = x;


nn = 1;
y_pad = zeros(length(x_pad),1);
debug_win = zeros(length(x_pad),1);
while nn < length(x_pad)-L+2
    % Consider the current frame and window it
    frame = x_pad(nn:nn+N-1);
    frame = frame .* s;
    
    %% 4a) Compute MDCT using FFT
    % Pre-twiddle
    n = (0:N-1)';
    frame_tilde = frame .* exp(-1i*pi*n/N);
    
    % Compute the FFT
    FRAME_tilde = fft(frame_tilde,N);
    
    % Post-twiddle (take only the first half of the frequency domain signal)
    k = (0:N/2-1)';
    FRAME = real(FRAME_tilde(1:N/2) .* exp(-1i*2*pi/N*n0*(k+1/2)));
    
    %% 4b) Compute inverse MDCT using FFT
    % Pre-twiddle
    k = (0:N-1)';
    FRAME_hat = [FRAME; -FRAME(end:-1:1)] .* exp(1i*2*pi/N*k*n0);
    
    % Compute inverse FFT (hint: expolit odd symmetry)
    frame_hat = ifft(FRAME_hat,N);
    
    % Post-twiddle
    frame_bar = real(frame_hat .* exp(1i*pi/N*(n+n0)));    
    
    %% 4c) Multiply by two times the synthesis window and perform overlap-and-add
    frame = 2*s .* frame_bar;
    y_pad(nn:nn+N-1) = y_pad(nn:nn+N-1) + frame;
    
    
    %% 7) Next window type selection
    if (nn == 1) % first frame
        Ecurr = sum(x_pad(nn:nn+N).^2);
        nn = nn + N/2;
    elseif strcmp(win_type,'long2short')
        Ecurr = sum(x_pad(nn:nn+N).^2);
        s = short_win;
        N = S;
        win_type = 'short';
        nn = nn + 3/4*L - 1/4*S;
    elseif strcmp(win_type,'short2long')
        Ecurr = sum(x_pad(nn:nn+N).^2);
        s = long_win;
        N = L;
        win_type = 'long';
        nn = nn + 1/2*L;
    else
        Enext = sum(x_pad(nn+N/2:nn+3/2*N).^2);
        if (Enext < 0.7*Ecurr && strcmp(win_type,'long'))
            N = L;
            s = long2short_win;
            win_type = 'long2short';
            nn = nn + N/2;
        elseif (Enext > 1.3*Ecurr && strcmp(win_type,'short'))
            N = L;
            s = short2long_win;
            win_type = 'short2long';
            nn = nn - (1/4*L - 1/4*S) + 1/2*S;
        else
            nn = nn + N/2;
            Ecurr = Enext;
        end
    end
    switch win_type
        case 'long'
            debug_win(nn) = 1;
        case 'short'
            debug_win(nn) = 3;
        case 'long2short'
            debug_win(nn) = 2;
        case 'short2long'
            debug_win(nn) = 4;
    end
end

%% 8) Compare the reconstructed signal with the original one
e = y_pad-x_pad;
figure
plot(e), title('MDCT (switching window) - Reconstruction error')
xlabel('Time'), ylabel('Error')
SNR_mdct = db(var(x_pad)/var(x_pad-y_pad));
disp(['MDCT with switching window - Signal-to-Noise Ratio = ' num2str(SNR_mdct) ' dB'])