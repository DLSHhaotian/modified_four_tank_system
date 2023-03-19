function [u,I] = SISOPID(ubar,ybar,y,yold,I,KP,KI,KD,dt,umin,umax)

e = ybar-y;
P = KP*e;
D = -KD*(y-yold)/dt;
u = ubar + P + I + D;

if (u >= umax)
u = umax;
elseif (u <= umin)
u = umin;
else
I = I + KI*e*dt;
end