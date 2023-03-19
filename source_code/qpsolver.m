function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
%A QP solver interface
% min   J=0.5*x'*H*x+g'*x
% s.t.  l<=x<=u
% s.t.  bl<=A*x<=bu
% Syntax: [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
%         info.fval: minimum value
%         info.exitflag: Reason quadprog stopped
%         info.output: optimization process
%         info.lambda: Lagrange multipliers at the solution
Ai=[A;-A];
bi=[bu;-bl];
%options = optimoptions('quadprog','Display','iter','Algorithm',"active-set");
options = optimoptions('quadprog','Algorithm',"active-set");
[x,fval,exitflag,output,lambda] = quadprog(H,g,Ai,bi,[],[],l,u,xinit,options);
info.fval=fval;
info.exitflag=exitflag;
info.output=output;
info.lambda=lambda;
end
