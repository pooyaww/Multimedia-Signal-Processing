%% MMSP2 - Lab 1
%  Exercise 4 - Discrete memoryless source coding
%  From 19 June 2006 exam
%  Lucio Bianchi - 04/12/2012

clear
close all
clc

%% Generate N=10000 samples of an AR(1) random process
%% x(n) = rho*x(n-1)+z(n), rho=0.99, z(n) is gaussian with variance=1
%% E[x(n)*z(n)] = 0
N = 10000;
rho = 0.99;
x = zeros(N+1,1);
z = randn(N+1,1);
lambda_2 = 1+rho^2;

for ii = 2:N+1
    x(ii) = rho*x(ii-1)+sqrt(lambda_2)*z(ii);
end
x = x(2:end);
z = z(2:end);


%% Clip sample values in the range [-20,20] and round to nearest integer
x(x<-20) = -20;
x(x>20) = 20;
x = round(x);

%% Compute H(X) assuming that X is a discrete memoryless source
d = hist(x, -20:20);
p = d/N;
p = p(p ~= 0);
H = -sum(p .* log2(p));

M = log2(41);
