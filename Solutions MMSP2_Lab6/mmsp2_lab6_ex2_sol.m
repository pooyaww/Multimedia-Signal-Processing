%% MMSP2 - Lab 6
%  Exercise 2 - Inter-frame (temporal) predictive encoder
%  Lucio Bianchi - 28/01/2014

clear
close all
clc

%% 1) Load the sequence 'table_tennis.mat' and compute the entropy
%%    hint: assume a DMS
load('table_tennis.mat');
seq = table_tennis;
clear table_tennis
[Nr, Nc, Nf] = size(seq);

h = hist(seq(:),0:255);
p = h / sum(h);
p = p(p~=0);
H = -sum(p.*log2(p));
disp(['Entropy of the DMS: ', num2str(H)]);

%% 2) Lossless open-loop prediction
seq_rec = zeros(Nr, Nc, Nf);
% first frame: I frame
seq_rec(:,:,1) = seq(:,:,1);

% successive frames: P frames
for i = 2:Nf
    seq_rec(:,:,i) = seq(:,:,i) - seq(:,:,i-1);
end

%% 3) Compute the entropy of the prediction error for the P frames
seq_rec = seq_rec(:,:,2:Nf);
h = hist(seq_rec(:),-255:255);
p = h / sum(h);
p = p(p~=0);
H = -sum(p.*log2(p));
disp(['Entropy of the open-loop prediction error: ', num2str(H)])

%% 4) Lossy open-loop prediction 
%% 4a) Quantize the prediction error using 16 levels
M = 16;
Q = (255)/M;
d_tilde = zeros(Nr, Nc, Nf);
seq_tilde = zeros(Nr, Nc, Nf);

% first frame: I frame
seq_tilde(:,:,1) = Q*round(seq(:,:,1)/Q);

% successive frames: P frames
for i = 2:Nf
    d = seq(:,:,i) - seq(:,:,i-1);
    d_tilde(:,:,i) = Q*round(d/Q);
end

%% 4b) Reconstruct the sequence
seq_rec = zeros(Nr, Nc, Nf);
% first frame: I frame
seq_rec(:,:,1) = seq_tilde(:,:,1);

% successive frames: P frames
for i = 2:Nf
    seq_rec(:,:,i) = d_tilde(:,:,i) + seq_rec(:,:,i-1);
end

%% 4c) Verify experimentally the problems caused by drift
%%     Compute the PSNR and plot the original and reconstructed frames
PSNR_O = zeros(Nf,1);
for i = 1:Nf
    MSE = sum(sum((seq_rec(:,:,i) - seq(:,:,i)).^2))/(Nr*Nc);
    PSNR_O(i) = 10*log10(255^2/MSE);
    figure(1)
    subplot(121), imshow(seq(:,:,i),[]), title(['Frame n.', num2str(i), ' of the original sequence'])
    subplot(122), imshow(seq_rec(:,:,i),[]), title(['Frame n.', num2str(i), ' of the open-loop reconstructed sequence'])
    pause
end
figure(2)
plot(PSNR_O,'b'),hold on
title('PSNR of open/closed-loop coding')
xlabel('Frames'), ylabel('PSNR [dB]')
legend('PSNR open-loop prediction')

%% 5) Lossy closed-loop prediction
%% 5a) Quantize the prediction error using 16 levels
M = 16;
Q = (255)/M;
d_tilde = zeros(Nr,Nc,Nf);
seq_tilde = zeros(Nr,Nc,Nf);

% first frame: I frame
seq_tilde(:,:,1) = Q*round(seq(:,:,1)/Q);

% successive frames: P frames
for i = 2:Nf
    d = seq(:,:,i) - seq_tilde(:,:,i-1);
    d_tilde(:,:,i) = Q*round(d/Q);
    seq_tilde(:,:,i) = d_tilde(:,:,i)+seq_tilde(:,:,i-1);
end

%% 5b) Reconstruct the sequence
seq_rec = zeros(Nr, Nc, Nf);

% first frame: I frame
seq_rec(:,:,1) = seq_tilde(:,:,1);

% successive frames: P frames
for i = 2:Nf
    seq_rec(:,:,i) = d_tilde(:,:,i) + seq_rec(:,:,i-1);
end

%% 5c) Verify experimentally the problems caused by drift
%%     Compute the PSNR and plot the original and reconstructed frames
PSNR_C = zeros(Nf,1);
for i = 1:Nf
    MSE = sum(sum((seq_rec(:,:,i) - seq(:,:,i)).^2))/(Nr*Nc);
    PSNR_C(i) = 10*log10(255^2/MSE);
    figure(1)
    subplot(121), imshow(seq(:,:,i),[]), title(['Frame n.', num2str(i), ' of the original sequence'])
    subplot(122), imshow(seq_rec(:,:,i),[]), title(['Frame n.', num2str(i), ' of the DPCM reconstructed sequence'])
    pause
end

figure(2)
plot(PSNR_C, 'r'), hold off
legend('PSNR open-loop prediction', 'PSNR closed-loop DPCM','Location','BestOutside')