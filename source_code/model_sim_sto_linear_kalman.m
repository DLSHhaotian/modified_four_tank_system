function [x,x_bar,d_bar,y,y_bar]=model_sim_sto_linear_kalman(u,d,w,v,t,N,x0,xbar0,dbar0,nxbar,A,B,C,D,L_kf,AL,BL,BdL,G_w,CL,DL,p,ys)
% --------------------------------------------------------------
nx = nxbar; nu = 2; ny = 4; nz = 2; nd = 2;
x = zeros(nx,N);
x_bar = zeros(nxbar,N);
d_bar = zeros(nd,N);
y_bar = zeros(ny,N);
y = zeros(ny,N);
y_no_noise = zeros(ny,N);
z = zeros(nz,N);
X = zeros(0,nx);
T = zeros(0,1);
x(:,1) = x0;
x_bar(:,1) =xbar0;
d_bar(:,1) = dbar0*ones(nd,1);
y_bar(:,1)=ys;
for k = 1:N-1

%[Tk,Xk] = ode15s(@ModifiedFourTankSystem,[t(k) t(k+1)],x(:,k),[],...
%u(:,k),d(:,k),p);
x(:,k+1)=AL*x(:,k)+BL*u(:,k)+BdL*d(:,k)+G_w*w(:,k);
y(:,k)=CL*x(:,k)+DL*u(:,k)+v(:,k);
y_no_noise(:,k)=CL*x(:,k)+DL*u(:,k);
z(:,k)=y_no_noise(1:2,k);

x_kf=[x_bar(:,k);d_bar(:,k)];
u_kf=u(:,k);
X_bar=A*x_kf+B*u_kf+A*L_kf*(y(:,k)-(C*x_kf+D*u_kf));
Y_bar=C*x_kf;
x_bar(:,k+1) = X_bar(1:nxbar);
d_bar(:,k+1) = X_bar(nxbar+1:end);
y_bar(:,k+1) = Y_bar;
end
k = N;
y(:,k)=CL*x(:,k)+DL*u(:,k)+v(:,k);
y_no_noise(:,k)=CL*x(:,k)+DL*u(:,k);
z(:,k)=y_no_noise(1:2,k);
%x_kf=[x(:,k);d_bar(:,k)];
%u_kf=u(:,k);
%X_bar=A*x_kf+B*u_kf+L_kf*(z(:,k)-(C*x_kf-D*u_kf));
%x_bar(:,k) = X_bar(1:nx);
%d_bar(:,k) = X_bar(nx+1:end);
end