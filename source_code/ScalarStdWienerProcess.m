function [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
% ScalarStdWienerProcess Ns realizations of a scalar std Wiener process
%
% Syntax: [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
% W : Standard Wiener process in [0,T]
% Tw : Time points
% dW : White noise used to generate the Wiener process
%
% T : Final time
% N : Number of intervals
% Ns : Number of realizations
% seed : To set the random number generator (optional)

if nargin == 4
rng(seed);
end
dt = T/N;
dW = sqrt(dt)*randn(Ns,N);
W = [zeros(Ns,1) cumsum(dW,2)];
Tw = 0:dt:T;