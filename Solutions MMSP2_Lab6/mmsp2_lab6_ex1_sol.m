%% MMSP2 - Lab 6
%  Exercise 1 - JPEG (baseline profile) encoder
%  Lucio Bianchi - 28/01/2014

clear
close all
clc

%% 1) Load the image 'lena512color.tiff' and transform it from RGB color space
%%    to YCbCr color space
Image = double(imread('lena512color.tiff'));
R = Image(:,:,1);
G = Image(:,:,2);
B = Image(:,:,3);

Y = 0.299*R + 0.587*G + 0.114*B;

[Nr, Nc] = size(Y);

figure(1)
subplot(121), imshow(Y,[]),colormap('gray'),title('Original')

%% 2) Consider only the luminance component Y, compute the DCT of the whole image
%%    and show the result
clear R G B Cb Cr
Y_dct = dct2(Y);
% show the DCT of the image
% hint: use log(abs(...))
figure(2)
imshow(log(abs(Y_dct)),[]); colormap('gray')

%% 3) Simulate JPEG baseline profile
N = 8; % NxN is the size of a block
load('Qjpeg.mat'); % quantization matrix
Q = 1; % scaling factor (fine --> coarse)
Qmatrix = QJPEG * Q;
Y_rec = zeros(size(Y));

for r = 1:Nr/N
    for c = 1:Nc/N
        %% 3a) consider block of size 8x8
        block = Y((r-1)*N+1 : r*N, (c-1)*N+1 : c*N);
        %% 3b) compute the DCT
        block_dct = dct2(block);
        %% 3c) threshold quantization
        block_dct_q = Qmatrix .* round(block_dct./Qmatrix);
        %% 3d) reconstruct the image from quantized coefficients
        Y_rec((r-1)*N+1 : r*N, (c-1)*N+1 : c*N) = idct2(block_dct_q);
    end
end

%% 3e) Show the reconstructed image
figure(1)
subplot(122),imshow(Y_rec,[]), colormap('gray'), title('JPEG baseline')

%% 3f) Compute the PSNR
E = Y_rec - Y;
MSE = sum(sum(E.^2))/(Nr*Nc);
PSNR = 10*log10(255^2/MSE);
disp(['PSNR = ', num2str(PSNR), 'dB'])