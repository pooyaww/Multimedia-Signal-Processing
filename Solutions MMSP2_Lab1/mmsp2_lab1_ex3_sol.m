%% MMSP2 - Lab 1
%  Exercise 3 - Discrete memoryless source coding
%  Lucio Bianchi - 10/12/2013

clear
close all
clc

%% 1) Generate one realization of length 1000000 of the following process:
%%    x(n)=min(max(0,round(rho*x(n-1)+z(n))),15)
%%    where rho=0.95 and z(n) is Gaussian with variance=1
N = 1000000;
rho = 0.95;
z = randn(N+1,1);
x = zeros(N+1,1);

for ii = 2:N+1
    x(ii) = min(max(0,round(rho*x(ii-1)+z(ii))),15);
end
x = x(2:end);
z = z(2:end);


%% 2) Determine the size of the alphabet of the source
%%    hint: use the function unique()
L = length(unique(x));

%% 3) Find H(X) assuming that x is a discrete memoryless source
d = hist(x,L);
p = d/N;
p = p(p ~= 0);
H = -sum(p .* log2(p));

%% 4) Let X=x(n) and Y=rho*x(n-1). Compute p(X,Y) and H(X,Y)
y = round(rho*x(1:end-1));
x = x(2:end);

d_joint = hist3([x y], {1:L, 1:L});
p_joint = d_joint / length(y);
p_joint = p_joint(p_joint ~= 0);
H_joint = -sum(p_joint .* log2(p_joint));

%% 5) Compute the conditional entropy H(X|Y)
H_cond = H_joint - H;