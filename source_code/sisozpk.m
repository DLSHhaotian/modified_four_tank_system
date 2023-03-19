function [z,p,k] = sisozpk(A,B,C,D)
% poles
p = eig(A);
% zeros
M = [A B; C D];
N = [eye(size(A)) zeros(size(B));zeros(size(C)) zeros(size(D))];
zz = eig(M,N);
idz = (zz < inf);
z = zz(idz);
% Kzp
s = randn(1);
den = prod(s-p);
num = prod(s-z);
G = C * ( (s*eye(size(A))-A) \ B ) + D;
k = G*den/num;
end