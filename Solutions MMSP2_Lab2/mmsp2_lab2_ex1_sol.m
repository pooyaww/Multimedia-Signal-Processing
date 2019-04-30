%% MMSP2 - Lab 2
%  Exercise 1 - Uniform and optimal quantizer
%  Lucio Bianchi - 17/12/2013

clear
close all
clc

%% 1) Generate a 1000-sample realization of s_g(n)~N(0,2) and s_u(n)~U with variance 2 and mean 0
%%    hint: use the functions randn() and rand()
mean_s = 0;
var_s = 2;
% Normal distribution
s_g = mean_s + sqrt(var_s) .* randn(100000,1);
% Uniform distribution
b = sqrt(3*var_s) + mean_s;
a = 2*mean_s - b;
s_u = a + (b-a) .* rand(100000,1);

% disp(['Mean of gaussian ' num2str(mean(s_g))])
% disp(['Variance of gaussian ' num2str(var(s_g))])
% disp(['Mean of uniform ' num2str(mean(s_u))])
% disp(['Variance of uniform ' num2str(var(s_u))])

%% 2) Quantize s_g(n) and s_u(n) with M=[4,8,16,32,64,128] levels and uniform
%%    quantizer. Plot R-D curve for each number of levels. Compare with Shannon's
%%    lower bound.
%%    hint: use MSE as distortion metric and plot the SNR
% Levels = [4, 8, 16, 32, 64, 128]
R = 2:13;

MSE_g = zeros(length(R),1);
MSE_u = zeros(length(R),1);

% for each level l in levels do
for ii = 1:length(R)
    % compute the quantization step
    Q_g = (max(s_g)-min(s_g)) / (2^R(ii));
    Q_u = (max(s_u)-min(s_u)) / (2^R(ii));
  
    % perform quantization
    y_g = Q_g * floor(s_g/Q_g) + Q_g/2;
    y_u = Q_u * floor(s_u/Q_u) + Q_u/2;
    
    % quantization error
    e_g = y_g - s_g;
    e_u = y_u - s_u;
    
    % MSE
    MSE_g(ii) = mean(e_g.^2);
    MSE_u(ii) = mean(e_u.^2);
end

% compute SNR
SNR_g = pow2db(var(s_g)./MSE_g);
SNR_u = pow2db(var(s_u)./MSE_u);
SNR_s = pow2db(2.^(2*R));

% plot R-D curves
figure
plot(R,SNR_g,'r-*'), hold on
plot(R,SNR_u,'g-*')
plot(R,SNR_s,'b-*'), hold off
xlabel('Rate [bps]'), ylabel('SNR [dB]')
title('Uniform quantizer: R-D curve')
legend('Normal distribution','Uniform distribution','Shannon bound')

%% 3) Design an optimal quantizer for s_g(n) with the same levels defined in the
%%    previous step. Compare the R-D curves obtained with optimal quantizer
%%    with those obtained with uniform quantizer
% Levels = [4, 8, 16, 32, 64, 128]
MSE_g_opt = zeros(length(R),1);
MSE_g_lm = zeros(length(R),1);

% for each level l in levels do
for ii = 1:length(R)
    % compute the quantization step Q
    M = 2^R(ii);
    Q = (max(s_g)-min(s_g))/M;
    
    % build the Q_matrix
    Q_matrix = zeros(M,3);
    
    % first row: overloading zone
    Q_matrix(1,1) = -inf;
    Q_matrix(1,2) = -Q*(M/2-1);
    Q_matrix(1,3) = mean(s_g(s_g <= Q_matrix(1,2)));
    
    % granular zone
    for jj = 2:M-1
        Q_matrix(jj,1) = Q_matrix(jj-1,2);
        Q_matrix(jj,2) = Q_matrix(jj,1) + Q;
        Q_matrix(jj,3) = mean(s_g(s_g > Q_matrix(jj,1) & s_g <= Q_matrix(jj,2)));
    end
    
    % last row: overloading zone
    Q_matrix(M,1) = Q*(M/2-1);
    Q_matrix(M,2) = inf;
    Q_matrix(M,3) = mean(s_g(s_g >= Q_matrix(M,1)));
    
    y_g_opt = zeros(length(s_g),1);
    % for each sample n in s_g(n) do
    for n = 1:length(s_g)
        % for each quantization level
        for mm = 1:M
            % if sample is inside the quantization interval
            if s_g(n) < Q_matrix(mm,2)
                % quantize with the corresponding value
                y_g_opt(n) = Q_matrix(mm,3);
                break
            end
        end
    end
    
    % quantization error
    e_g_opt = y_g_opt - s_g;
    
    % MSE
    MSE_g_opt(ii) = mean(e_g_opt.^2);
    
    % Lloyd-Max quantizer with Matlab implementation
    [partition,codebook] = lloyds(s_g,M);
    [index,quants,distor] = quantiz(s_g,partition,codebook);
    e_g_lm = quants' - s_g;
    MSE_g_lm(ii) = mean(e_g_lm.^2);
end

% compute SNR
SNR_g_opt = pow2db(var(s_g)./MSE_g_opt);
SNR_g_lm = pow2db(var(s_g)./MSE_g_lm);
SNR_s = pow2db(2.^(2*R));

% plot R-D curves
figure
plot(R,SNR_g_opt,'r-*'), hold on
plot(R,SNR_g_lm,'k-*')
plot(R,SNR_g,'g-*')
plot(R,SNR_s,'b-*'), hold off
xlabel('Rate [bps]'), ylabel('SNR [dB]')
title('Quantization of gaussian process: R-D curve')
legend('Optimal quantizer','Lloyd-Max quantizer','Uniform quantizer','Shannon bound')