function u_new=LMPCcompute_inputcons_outputSoftCons(R,X0,D,U,Ad,Bd,Bd_d,Css,u_min,u_max,u_delta_min,u_delta_max,z_min,z_max,N,Q_cof,u_delta_cof,eta_cof1,eta_cof2)
% Syntax: u_new=LMPCcompute_inputcons_outputSoftCons(R,X0,D,U,Ad,Bd,Bd_d,Css,u_min,u_max,u_delta_min,u_delta_max,z_min,z_max,N,Q_cof,u_delta_cof,eta_cof1,eta_cof2)
%         R: Reference throughout the prediction period
%         X0: Observed state x0
%         D: Observed disturbance d0
%         U: Last input U-1
%         Ad,Bd,Bd_d,Czss: Linear model state space matrix
%         u_min: Minimum of input
%         u_max: Maximum of input
%         u_delta_min: Minimum of input rate
%         u_delta_max: Maximum of input rate
%         z_min: Minimum of output
%         z_max: Maximum of input
%         N: prediction length
%         Q_cof: Weight matrix for tracking the reference 
%         u_delta_cof: Weight matrix for input change rate 
%         eta_cof1: Weight matrix for Minimum of output
%         eta_cof2: Weight matrix for Maximum of output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         u_new: The first input to be executed
num_x=size(Ad,1);
num_u=size(Bd,2);
num_d=size(Bd_d,2);
num_y=4;
num_z=2;
num_eta=num_y;
num_u_delta_cons=N*num_u;
num_z_cons=N*num_y;
Phi=zeros(N*num_y,num_x);
Phi_d=zeros(N*num_y,num_d);
Gam=zeros(N*num_y,N*num_u+2*N*num_eta);
%Gam_d=zeros(N*num_y,N*num_d);
%Weight matrix for tracking the reference Qz
Qz=diag(zeros(N*num_y,1));
%Weight matrix for input rate S
Sz=diag(u_delta_cof*ones(num_u,1));
%Quadratic objective term Hs
Hs=zeros(N*num_u+2*N*num_eta,N*num_u+2*N*num_eta);
%Linear objective term Mu_delta
Mu_delta=zeros(N*num_u+2*N*num_eta,num_u);
Mu_delta(1:num_u,:)=-Sz;
%Quadratic objective term H_eta
H_eta1=diag([zeros(N*num_u,1);eta_cof1*ones(N*num_eta,1);zeros(N*num_eta,1)]);
H_eta2=diag([zeros(N*num_u,1);zeros(N*num_eta,1);eta_cof2*ones(N*num_eta,1)]);
U0=zeros(N*num_u+2*N*num_eta,1);
Iu=diag(ones(num_u,1));
Au_delta_cons=zeros(num_u_delta_cons,N*num_u+2*N*num_eta);
for i=1:N
    %fill the Phi
    Phi((i-1)*num_y+1:i*num_y,:)=Css*Ad^i;
    Phi_d((i-1)*num_y+1:i*num_y,:)=Css*Ad^i*Bd_d;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam
    Gam_i=zeros(num_y,N*num_u+2*N*num_eta);%[num_x * num_u*N]
    for j=1:i
        Gam_i(:,(j-1)*num_u+1:j*num_u)=Css*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Gam((i-1)*num_y+1:i*num_y,:)=Gam_i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam_d
%     Gam_di=zeros(num_y,N*num_d);%[num_x * num_u*N]
%     for j=1:i
%         Gam_di(:,(j-1)*num_d+1:j*num_d)=Css*(Ad^(i-j))*Bd_d;%[Hn,Hn-1,Hn-2]
%     end
%     Gam_d((i-1)*num_y+1:i*num_y,:)=Gam_di;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Qz
    qz=diag([Q_cof*ones(num_z,1);zeros(num_y-num_z,1)]);
    Qz((i-1)*num_y+1:i*num_y,(i-1)*num_y+1:i*num_y)=qz;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:i*num_u)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the U0
    U0((i-1)*num_u+1:i*num_u,:)=U;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Au_delta_cons
    if i==1
        Au_delta_cons(1:num_u,1:num_u)=Iu;
    else
        Au_delta_cons((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:i*num_u)=[-Iu,Iu];
    end
end
Mx0=Gam'*Qz*Phi;
Md=Gam'*Qz*Phi_d;
%Linear objective term Mr
Mr=-Gam'*Qz;
%Quadratic objective term Hr
Hr=Gam'*Qz*Gam;
Hu=Hr+Hs+H_eta1+H_eta2;
g=Mx0*X0+Mr*R+Mu_delta*U+Md*D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Az_cons_min=Gam;
Az_cons_min(:,N*num_u+1:N*num_u+N*num_eta)=diag(1*ones(N*num_eta,1));
Az_cons_min(:,N*num_u+N*num_eta+1:end)=diag(zeros(N*num_eta,1));
Az_cons_max=Gam;
Az_cons_max(:,N*num_u+1:N*num_u+N*num_eta)=diag(zeros(N*num_eta,1));
Az_cons_max(:,N*num_u+N*num_eta+1:end)=diag(-1*ones(N*num_eta,1));
Z_min=repmat(z_min,N,1)-Phi*X0-Phi_d*D;
Z_max_fake=5000*ones(num_z_cons,1);
Z_max=repmat(z_max,N,1)-Phi*X0-Phi_d*D;
Z_min_fake=zeros(num_z_cons,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_min=[u_min*ones(N*num_u,1);zeros(2*N*num_eta,1)];
U_max=[u_max*ones(N*num_u,1);5000*ones(2*N*num_eta,1)];
U_delta_min=u_delta_min*ones(num_u_delta_cons,1);
U_delta_max=u_delta_max*ones(num_u_delta_cons,1);
U_delta_min(1)=u_delta_min+U(1);
U_delta_min(2)=u_delta_min+U(2);
U_delta_max(1)=u_delta_max+U(1);
U_delta_max(2)=u_delta_max+U(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_ieq=[Au_delta_cons;Az_cons_min;Az_cons_max];
b_ieq_min=[U_delta_min;Z_min;Z_min_fake];
b_ieq_max=[U_delta_max;Z_max_fake;Z_max];
[u_opt,~] = qpsolver(0.5*(Hu+Hu'),g,U_min,U_max,A_ieq,b_ieq_min,b_ieq_max,U0);
u_new=u_opt(1:num_u);
end