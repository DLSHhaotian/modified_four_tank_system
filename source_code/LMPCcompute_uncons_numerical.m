function u_new=LMPCcompute_uncons_numerical(R,X0,D,U,Ad,Bd,Bd_d,Czss,N,Q_cof,u_delta_cof)
% Syntax: u_new=LMPCcompute_uncons_numerical(R,X0,D,U,Ad,Bd,Bd_d,Czss,N,Q_cof,u_delta_cof)
%         R: Reference throughout the prediction period
%         X0: Observed state x0
%         D: Observed disturbance d0
%         U: Last input U-1
%         Ad,Bd,Bd_d,Czss: Linear model state space matrix
%         N: prediction length
%         Q_cof: Weight matrix for tracking the reference 
%         u_delta_cof: Weight matrix for input change rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         u_new: The first input to be executed
num_x=size(Ad,1);
num_u=size(Bd,2);
num_d=size(Bd_d,2);
num_y=4;
num_z=2;
Phi=zeros(N*num_z,num_x);
Phi_d=zeros(N*num_z,num_d);
Gam=zeros(N*num_z,N*num_u);
%Gam_d=zeros(N*num_z,N*num_d);
%Weight matrix for tracking the reference Qz
Qz=diag(Q_cof*ones(N*num_z,1)) ;
%Weight matrix for input rate S
Sz=diag(u_delta_cof*ones(num_u,1));
%Quadratic objective term Hs
Hs=zeros(N*num_u,N*num_u);
%Linear objective term Mu_delta
Mu_delta=zeros(N*num_u,num_u);
Mu_delta(1:num_u,:)=-Sz;
U0=zeros(N*num_u,1);
for i=1:N
    %fill the Phi
    Phi((i-1)*num_z+1:i*num_z,:)=Czss*Ad^i;
    Phi_d((i-1)*num_z+1:i*num_z,:)=Czss*Ad^i*Bd_d;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam
    Gam_i=zeros(num_z,N*num_u);%[num_x * num_u*N]
    for j=1:i
        Gam_i(:,(j-1)*num_u+1:j*num_u)=Czss*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Gam((i-1)*num_z+1:i*num_z,:)=Gam_i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %fill the Gam_d
%     Gam_di=zeros(num_z,N*num_d);%[num_x * num_u*N]
%     for j=1:i
%         Gam_di(:,(j-1)*num_d+1:j*num_d)=Czss*(Ad^(i-j))*Bd_d;%[Hn,Hn-1,Hn-2]
%     end
%     Gam_d((i-1)*num_z+1:i*num_z,:)=Gam_di;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs(end-num_u+1:end,end-2*num_u+1:end)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the U0
    U0((i-1)*num_u+1:i*num_u,:)=U;
end
Mx0=Gam'*Qz*Phi;
Md=Gam'*Qz*Phi_d;
%Linear objective term Mr
Mr=-Gam'*Qz;
%Quadratic objective term Hr
Hr=Gam'*Qz*Gam;
Hu=Hr+Hs;
g=Mx0*X0+Mr*R+Mu_delta*U+Md*D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u_opt,~] = qpsolver(0.5*(Hu+Hu'),g,[],[],[],[],[],U0);
u_new=u_opt(1:num_u);
end