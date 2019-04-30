%% MMSP2 - Lab 6
%  Exercise 3 - DCT/KLT comparison
%  Lucio Bianchi - 22/01/2013

clear
close all
clc

%% 1) Load the image 'lena512_color.tiff' and extract the luminance component
Image = imread('lena512color.tiff');
R = double(Image(:,:,1));
G = double(Image(:,:,2));
B = double(Image(:,:,3));

Y = 0.299*R + 0.587*G + 0.114*B;

[Nr,Nc] = size(Y);

figure(1)
subplot(131),imshow(Y,[]), colormap('gray'), title('Original');

%% 2) Consider blocks of dimension 8x8 and estimate the correlation matrix
N = 8;
M = Nr/N * Nc/N;
B = zeros(M,N^2);
ii = 1;
for rr = 1:Nr/N
    for cc = 1:Nc/N
        aux = Y((rr-1)*N+1:rr*N,(cc-1)*N+1:cc*N);
        B(ii,:) = aux(:) - mean(aux(:));
        ii = ii+1;
    end
end

meanY = mean(Y(:));
R = zeros(N^2,N^2,M);
for ii = 1:M
    R(:,:,ii) = (B(ii,:)'-meanY') * (B(ii,:)-meanY);
end
R = mean(R,3);

%% 3) Perform KLT coding
[V, D] = eig(R);
T = fliplr(V)';
Qmatrix = ones(size(N,N));

Y_klt_rec = zeros(size(Y));
for r = 1:Nr/N
    for c = 1:Nc/N
        %% 3a) consider block of size 8x8
        block = Y((r-1)*N+1 : r*N, (c-1)*N+1 : c*N);
        %% 3b) compute the KLT
        block_klt = T*block(:);
        %% 3c) threshold quantization
        block_klt_q = Qmatrix .* round(reshape(block_klt,N,N)./Qmatrix);
        %% 3d) reconstruct the image from quantized coefficients
        Y_klt_rec((r-1)*N+1 : r*N, (c-1)*N+1 : c*N) = reshape(T'*block_klt_q(:),N,N);
    end
end


figure(1)
subplot(132),imshow(Y_klt_rec,[]), colormap('gray'), title('KLT');

%% 4) Compute PSNR for the KLT transform
E = Y_klt_rec - Y;
MSE = sum(sum(E.^2))/(Nr*Nc);
PSNR_klt = 10*log10(255^2/MSE);
disp(['PSNR for the KLT transform: ', num2str(PSNR_klt)]);

%% 5) JPEG baseline coding - estimate the PSNR
Y_dct_rec = zeros(size(Y));
load('Qjpeg.mat'); % quantization matrix
Q = 1; % scaling factor (fine --> coarse)
Qmatrix = QJPEG * Q;

for r = 1:size(Y,1)/N
    for c = 1:size(Y,2)/N
        %% 3a) consider block of size 8x8
        block = Y((r-1)*N+1 : r*N, (c-1)*N+1 : c*N);
        %% 3b) compute the DCT
        block_dct = dct2(block);
        %% 3c) threshold quantization
        block_dct_q = Qmatrix .* round(block_dct./Qmatrix);
        %% 3d) reconstruct the image from quantized coefficients
        Y_dct_rec((r-1)*N+1 : r*N, (c-1)*N+1 : c*N) = idct2(block_dct_q);
    end
end

figure(1)
subplot(133),imshow(Y_dct_rec,[]), colormap('gray'), title('DCT');

E = Y_dct_rec - Y;
MSE = sum(sum(E.^2))/(Nr*Nc);
PSNR_dct = 10*log10(255^2/MSE);
disp(['PSNR for the DCT transform: ', num2str(PSNR_dct)]);