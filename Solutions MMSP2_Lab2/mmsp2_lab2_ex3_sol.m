%% MMSP2 - Lab 2
%  Exercise 3 - Scalar quantization
%  Lucio Bianchi - 17/12/2013

clear
close all
clc

%% 1) Suppose that you are sampling the output of a sensor at 10 kHz. 
%%    Quantize the output with a uniform quantizer at 10 bit per sample.
%%    Assume that the pdf of the signal is gaussian with mean 0 V and variance 4 V^2.
%%    What is the bit rate of the quantized signal?
Fs = 10000;
R = 10;
T = 10;
L = T*Fs;

bitRate = Fs*R;

disp(['Bit rate ' num2str(bitRate) ' bit/sec']);

%% 2) What would be a reasonable choice for the quantization step?
x = sqrt(4) .* randn(L,1);

Q = (max(x)-min(x))/(2^R);

disp(['Quantization step ' num2str(Q) ' V']);

%% 3) What is the MSE?
x_q = Q * floor(x/Q) + Q/2;

% quantization error
e = x_q - x;

% MSE
MSE = mean(e.^2);
MSE = var(e);

disp(['MSE ' num2str(MSE)]);

%% What is the resulting SNR?
SNR = pow2db(var(x)/MSE);

disp(['SNR ' num2str(SNR) ' dB'])