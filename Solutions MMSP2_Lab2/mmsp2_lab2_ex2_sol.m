%% MMSP2 - Lab 2
%  Exercise 2 - Vector quantization
%  Lucio Bianchi - 17/12/2013

clear
close all
clc

%% 1) Consider a 2D vector quantizer with codebook y1=(1,2), y2=(1,4), y3=(-1,2), y4=(0,-2)
%%    Show optimal assignement regions
y = zeros(2,4);
y(:,1) = [1,2];
y(:,2) = [1,4];
y(:,3) = [-1,2];
y(:,4) = [0,-2];

voronoi(y(1,:),y(2,:));

%% 2) Quantize the sequence x=(-4:5) and compute the MSE.
%%    Assume to quantize groups of 2 samples at time
x = -4:5;

hvq = dsp.VectorQuantizerEncoder(...
    'Codebook', y,...
    'CodewordOutputPort', true,...
    'QuantizationErrorOutputPort', true);

x = reshape(x,2,5);

[idx, x_q, err] = step(hvq, x);

plot(x_q(1,:), x_q(2,:),'ro'), hold on
plot(x(1,:), x(2,:),'bo'), hold off
