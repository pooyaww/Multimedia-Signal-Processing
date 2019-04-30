%% MMSP2 - Lab 1
%  Exercise 2 - Image signal encoding
%  Lucio Bianchi - 10/12/2013

clear
close all
clc

%% 1) Load the image 'lena512color.tiff' and display the normalized histogram 
%%    for all (R,G,B) components
%%    hint: 'lena512color.tiff' is an 8-bit RGB image, better to convert into double
%%    hint: express the RGB components as vectors
Image = imread('lena512color.tiff');
R = double(Image(:,:,1));
G = double(Image(:,:,2));
B = double(Image(:,:,3));

R = R(:);
G = G(:);
B = B(:);

alphabet = 0:255;
d_R = hist(R, alphabet);
d_G = hist(G, alphabet);
d_B = hist(B, alphabet);
p_R = d_R/length(R);
p_G = d_G/length(G);
p_B = d_B/length(B);

figure
title('PDF')
subplot(3,1,1),bar(alphabet,p_R),title('Red channel')
subplot(3,1,2),bar(alphabet,p_G),title('Green channel')
subplot(3,1,3),bar(alphabet,p_B),title('Blu channel')

%% 2) Compute the entropy of each channel
p_R = p_R(p_R ~= 0);
p_G = p_G(p_G ~= 0);
p_B = p_B(p_B ~= 0);

H_R = -sum(p_R .* log2(p_R));
H_G = -sum(p_G .* log2(p_G));
H_B = -sum(p_B .* log2(p_B));

disp(['Red channel = ', num2str(H_R)]);
disp(['Green channel = ', num2str(H_G)]);
disp(['Blu channel = ', num2str(H_B)]);

%% 3) Let X represent the source of the red channel and Y the source of the
%%    green channel. Compute and show p(X,Y).
%%    hint: use the function imagesc() to show p(X,Y)
d_joint = hist3([R G], {alphabet, alphabet});
p_joint = d_joint/length(R);

figure
imagesc(p_joint), colorbar, axis xy
xlabel('Y symbols')
ylabel('X symbols')
title('Joint pdf')


%% 4) Compute and display the joint entropy H(X,Y)
p_joint = p_joint(p_joint ~=0);
H_joint = -sum(p_joint .* log2(p_joint));

%% 5) Suppose to encode Y with H(Y) bits and transmit source N = aX+b-Y instead of X
%%    Compute the entropy of N and the conditional entropy H(X|Y)
% Compute coefficients a and b with LS
X_hat = [R ones(length(R),1)];
Y = G;
C = inv(X_hat'*X_hat)*X_hat'*Y;
a = C(1);
b = C(2);

% Compute H(N)
N = round(a*R+b-G);
d_N = hist(N,alphabet);
p_N = d_N/length(N);
p_N = p_N(p_N ~= 0);
H_N = -sum(p_N .* log2(p_N));

% Conditional entropy H(X|Y) = H(X,Y) - H(Y)
H_cond = H_joint - H_G;
