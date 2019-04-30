%% MMSP2 - Lab 3
%  Exercise 1 - Predictive coding
%  Lucio Bianchi - 07/01/2014

clear
close all
clc

%% 1) Generate 10000 samples of the random process
%%    x(n) = rho*x(n-1) + z(n), where rho=0.95 and z(n)~N(0,0.1)
N = 10000;
rho = 0.95;
z = sqrt(0.1) * randn(N,1);
x = filter(1,[1 -rho], z);

%% 2) Build a PCM codec. Quantize the signal with an uniform quantizer and 
%%    R bits. Compute the R-D curve for R=1,2,...,8 bits
R = 1:8;

MSE_pcm = zeros(length(R),1);
for ii = 1:length(R)
    Q = (max(x)-min(x)) / (2^(R(ii)));
    x_tilde = Q * floor(x/Q) + Q/2;
    MSE_pcm(ii) = var(x-x_tilde);
end

SNR_pcm = db(var(x)./MSE_pcm);

%% 3) Build a predictive codec in open loop. Use the optimal MMSE predictor.
%%    Use PCM to initialize the codec.

% Remark:
% Given the notation presented in class, the optimal MMSE predictor is
% x_hat(n) = rho*x(n-1), hence the prediction residual is
% d(n) = x(n) - x_hat(n) = x(n) - rho*x(n-1) = z(n).
%
% For the first sample of the process (i.e. n = 1) we would have 
% d(1) = x(1) - rho*x(0), which is impossible to compute. So, in order to
% highlight this fact, we set d(1) = NaN; please notice that the value of d(1)
% will remain unused.
%
% For the second sample of the process (i.e. n = 2) we have 
% d(2) = x(2) - rho*x(1) = z(2).
%
% From that remark follows the definition of the vector d as
d = [NaN; z(2:N)];

MSE_olpc = zeros(length(R),1);
for ii = 1:length(R)
    x_tilde = zeros(N,1);
    
    % first sample: PCM
    Qx = (max(x)-min(x)) / (2^(R(ii)));
    x_tilde(1) = Qx * floor(x(1)/Qx) + Qx/2;
    
    % next samples: OLPC
    Qd = (max(z)-min(z)) / (2^(R(ii)));
    d_tilde = Qd * floor(d(2:N)/Qd) + Qd/2;
    x_tilde(2:N) = filter(1,[1 -rho], d_tilde);
    
    % MSE
    MSE_olpc(ii) = var(x-x_tilde);
end

SNR_olpc = db(var(x)./MSE_olpc);


%% 4) Build a DPCM codec. Use the optimal MMSE predictor.
%%    Use PCM to initialize the codec. 
%%    For each rate R compute G_clp*G_Q.
MSE_dpcm = zeros(length(R),1);
for ii = 1:length(R)
    x_tilde = zeros(N,1);
    
    % first sample: PCM
    Qx = (max(x)-min(x)) / (2^(R(ii)));
    x_tilde(1) = Qx * floor(x(1)/Qx) + Qx/2;
    
    % next samples: DPCM
    Qd = (max(z)-min(z)) / (2^(R(ii)));
    for nn = 2:N
        d = x(nn) - rho*x_tilde(nn-1);
        d_tilde = Qd * floor(d/Qd) + Qd/2;
        x_tilde(nn) = d_tilde + rho*x_tilde(nn-1);
    end
    
    % MSE
    MSE_dpcm(ii) = var(x-x_tilde);
end

SNR_dpcm = db(var(x)./MSE_dpcm);

%% 5) Compare R-D curves for PCM, open-loop DPCM and closed-loop DPCM
figure
plot(R,SNR_pcm,'r'),hold on
plot(R, SNR_olpc,'g')
plot(R, SNR_dpcm, 'b'), hold off
xlabel('Rate [bit]'), ylabel('SNR [dB]')
legend('PCM', 'OLPC', 'DPCM')

%% 6) Build a Delta Modulation codec. Find the
%%    variance of the quantization error and plot both the original signal
%%    and the delta-modulated one. Set Delta = 1 and alpha=0.9
Delta = 1;
alpha = 0.9;

x_tilde = zeros(N,1);

% first sample
x_tilde(1) = Delta * sign(x(1));

% next samples
for nn = 2:N
    d = x(nn) - alpha*x_tilde(nn-1);
    d_tilde = Delta * sign(d);
    x_tilde(nn) = d_tilde + alpha*x_tilde(nn-1);
end

% MSE
MSE_dm = var(x-x_tilde);
SNR_dm = db(var(x)./MSE_dm);

disp(['Delta-modulation: SNR = ' num2str(SNR_dm) ' dB'])

DM_fig = figure;
subplot(2,1,1)
plot(x,'b'), hold on
plot(x_tilde,'r'), hold off
xlabel('Samples'), ylabel('Amplitude')
axis([2500 2800 -3 3])
legend('Original','DM')

%% 7) Build an adaptive Delta Modulation codec. Set epsilon = 1e-2
epsilon = 1e-2;

x_tilde = zeros(N,1);

% first 2 samples: initialization
x_tilde(1) = 0;
Delta = 1;
x_tilde(2) = Delta*sign(x(2)-alpha*x_tilde(1)) + alpha*x_tilde(1);
Delta = sqrt(var(x(1:2)-x_tilde(1:2)));

% next samples
for nn = 3:N
    x_tilde(nn) = Delta*sign(x(nn)-alpha*x_tilde(nn-1)) + alpha*x_tilde(nn-1);
    
    gamma = Delta / sqrt(var(x(1:nn)-x_tilde(1:nn)));
    if gamma < 1-epsilon
        Delta = Delta + 0.05;
        if Delta > 1
            Delta = 1;
        end
    elseif gamma > 1+epsilon
        Delta = Delta - 0.05;
        if Delta < 0.01
            Delta = 0.01;
        end
    end
end

% MSE
MSE_adm = var(x-x_tilde);
SNR_adm = db(var(x)./MSE_adm);

disp(['Adaptive Delta-modulation: SNR = ' num2str(SNR_adm) ' dB'])

figure(DM_fig);
subplot(2,1,2)
plot(x,'b'), hold on
plot(x_tilde,'g'), hold off
xlabel('Samples'), ylabel('Amplitude')
axis([2500 2800 -3 3])
legend('Original','ADM')
