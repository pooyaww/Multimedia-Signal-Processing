%% MMSP2 - Lab 3
%  Exercise 2 - Transform coding
%  Lucio Bianchi - 07/01/2014

clear
close all
clc

%% 1) Load the first 4s of the file 'gb.wav' and quantize it with PCM and R=8 bit.
%%    Compute the MSE and perceptually evaluate the result.
dur = 4;
[x, Fs] = wavread('gb.wav');
x = x(1:dur*Fs);

R = 8;
Q = (max(x)-min(x)) / (2^R);
x_tilde_pcm = Q * floor(x/Q) + Q/2;
MSE_pcm = var(x - x_tilde_pcm);
SNR_pcm = db(var(x)/MSE_pcm);
disp(['PCM coding: SNR = ', num2str(SNR_pcm), ' dB']);

%% 2) Consider groups of 8 symbols and quantize them using an optimal allocation of the 8 bits
X = reshape(x, R, length(x)/R)';

variances = var(X,0,1);
X_tilde_8 = zeros(length(x)/R,R);
r = zeros(R,1);
for ii = 1:R
    r(ii) = R + round(0.5*log2(variances(ii)/geomean(variances)));
    Q = (max(X(:,ii))-min(X(:,ii))) / (2^r(ii));
    X_tilde_8(:,ii) = Q * floor(X(:,ii)/Q) + Q/2;
end
X_tilde_8 = X_tilde_8';
x_tilde_8 = X_tilde_8(:);
MSE_8 = var(x-x_tilde_8);
SNR_8 = db(var(x)/MSE_8);
disp(['8 symbols coding: SNR = ', num2str(SNR_8), ' dB']);

bit_all_Fig = figure;
subplot(3,1,1)
bar(r)
xlabel('Symbol'), ylabel('Bit')
title('Bit allocation - No transform')

%% 3) Consider DCT transformation and repeat step 2 over transformed 
%%    coefficients. Find the distortion and evaluate the perceived quality.
[I, J] = size(X);
N = J;

T = zeros(N,N);

T(1,:) = sqrt(1/N);
for ii = 2:N
    for jj = 1:N
        T(ii,jj) = sqrt(2/N)*cos(pi/(2*N)*(ii-1)*(2*jj-1));
    end
end

% DCT transform
X_dct = zeros(I, J);
for ii = 1:I
    X_dct(ii,:) = T*X(ii,:)';
end

% quantization
variances = var(X_dct,0,1);
X_tilde_dct = zeros(length(x)/R,R);
r = zeros(R,1);
for ii = 1:R
    r(ii) = R + round(0.5*log2(variances(ii)/geomean(variances)));
    Q = (max(X_dct(:,ii))-min(X_dct(:,ii))) / (2^r(ii));
    X_tilde_dct(:,ii) = Q * floor(X_dct(:,ii)/Q) + Q/2;
end

figure(bit_all_Fig);
subplot(3,1,2)
bar(r)
xlabel('Symbol'), ylabel('Bit')
title('Bit allocation - DCT')

% inverse DCT transform
X_tilde_idct = zeros(length(x)/R,R);
for ii = 1:I
    X_tilde_idct(ii,:) = T'*X_tilde_dct(ii,:)';
end

X_tilde_idct = X_tilde_idct';
x_tilde_idct = X_tilde_idct(:);
MSE_dct = var(x-x_tilde_idct);
SNR_dct = db(var(x)/MSE_dct);
disp(['DCT coding: SNR = ', num2str(SNR_dct), ' dB']);


%% 3) Consider a Karhunen-Loève transformation and repeat step 2 over transformed 
%%    coefficients. Find the distortion and evaluate the perceived quality.

RR = zeros(N,N,I);
mean_X = mean(X,1);
for ii = 1:I
    RR(:,:,ii) = (X(ii,:)'-mean_X') * (X(ii,:)-mean_X);
end
RR = mean(RR,3);

[V,D] = eig(RR);
T = fliplr(V)';

% KLT transform
X_klt = zeros(length(x)/R,R);
for ii = 1:I
    X_klt(ii,:) = T*X(ii,:)';
end

% quantization
variances = var(X_klt,0,1);
X_tilde_klt = zeros(length(x)/R,R);
r = zeros(R,1);
for ii = 1:R
    r(ii) = R + round(0.5*log2(variances(ii)/geomean(variances)));
    Q = (max(X_klt(:,ii))-min(X_klt(:,ii))) / (2^r(ii));
    X_tilde_klt(:,ii) = Q * floor(X_klt(:,ii)/Q) + Q/2;
end

figure(bit_all_Fig);
subplot(3,1,3)
bar(r)
xlabel('Symbol'), ylabel('Bit')
title('Bit allocation - KLT')

% inverse KLT transform
X_tilde_iklt = zeros(length(x)/R,R);
for ii = 1:I
    X_tilde_iklt(ii,:) = T'*X_tilde_klt(ii,:)';
end

X_tilde_iklt = X_tilde_iklt';
x_tilde_iklt = X_tilde_iklt(:);
MSE_klt = var(x-x_tilde_iklt);
SNR_klt = db(var(x)/MSE_klt);
disp(['KLT coding: SNR = ', num2str(SNR_klt), ' dB']);