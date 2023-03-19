function [u,I] = MIMOPID(ubar,ybar,y,yold,I,KP,KI,KD,dt,umin,umax,decoup)
% Syntax: [u,I] = MIMOPID(ubar,ybar,y,yold,I,KP,KI,KD,dt,umin,umax,decop)
%         ubar: Reference input
%         ybar: Reference output
%         y: measured output
%         yold: Last measured output
%         I: Integral term
%         KP: Proportional coefficient
%         KI: Integral coefficient
%         KD: Derivative coefficient
%         dt: Unit time
%         u_min: Minimum of input
%         u_max: Maximum of input
%         decop: Decopling matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         u: The input to be executed
%         I: Updated Integral term
e = ybar-y;
e_decop=decoup*e;
P = KP*e_decop;
D = -KD*(y-yold)/dt;
u = ubar + P + I + D;
flag=1;
for n=1:2
if (u(n) >= umax)
u(n) = umax;
flag=0;
elseif (u(n) <= umin)
u(n) = umin;
flag=0;
end
end
if(flag==1)
I = I + KI*e_decop*dt;
end