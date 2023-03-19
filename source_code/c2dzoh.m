function [Abar,Bbar]=c2dzoh(A,B,Ts)
[nx,nu]=size(B);
M = [A B; zeros(nu,nx) zeros(nu,nu)];
Phi = expm(M*Ts);
Abar = Phi(1:nx,1:nx);
Bbar = Phi(1:nx,nx+1:nx+nu);