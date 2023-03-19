function [x,x_bar,d_bar,y,y_bar]=model_sim_sto_kalman(u,d,v,t,N,x0,xbar0,dbar0,nxbar,A,B,C,D,L_kf,p,ys)
% --------------------------------------------------------------
nx = 4; nu = 2; ny = 4; nz = 2; nd = 2;
x = zeros(nx,N);
x_bar = zeros(nxbar,N);
d_bar = zeros(nd,N);
y_bar = zeros(ny,N);
y = zeros(ny,N);
z = zeros(nz,N);
X = zeros(0,nx);
T = zeros(0,1);
x(:,1) = x0;
y_bar(:,1)=ys;
x_bar(:,1) =xbar0;
d_bar(:,1) = dbar0*ones(nd,1);
for k = 1:N-1
%y(:,k) = FourTankSystemSensor(x(:,k),p); % Sensor function
y(:,k) = FourTankSystemSensorNoise(x(:,k),p,v(:,k));
z(:,k) = FourTankSystemOutput(x(:,k),p); % Output function
y_s=y(:,k);
x_kf=[x_bar(:,k);d_bar(:,k)];
u_kf=u(:,k);
X_bar=A*x_kf+B*u_kf+A*L_kf*(y_s-(C*x_kf+D*u_kf));
Y_bar=C*x_kf;
x_bar(:,k+1) = X_bar(1:nxbar);
d_bar(:,k+1) = X_bar(nxbar+1:end);
y_bar(:,k+1) = Y_bar;
[Tk,Xk] = ode15s(@ModifiedFourTankSystem,[t(k) t(k+1)],x(:,k),[],...
u(:,k),d(:,k),p);
x(:,k+1) = Xk(end,:)';
T = [T; Tk];
X = [X; Xk];
end
k = N;
y(:,k) = FourTankSystemSensorNoise(x(:,k),p,v(:,k)); % Sensor function
z(:,k) = FourTankSystemOutput(x(:,k),p); % Output function
%x_kf=[x(:,k);d_bar(:,k)];
%u_kf=u(:,k);
%X_bar=A*x_kf+B*u_kf+L_kf*(z(:,k)-(C*x_kf-D*u_kf));
%x_bar(:,k) = X_bar(1:nx);
%d_bar(:,k) = X_bar(nx+1:end);
end