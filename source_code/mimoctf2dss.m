function [Ad,Bd,Cd,Dd,sH]=mimoctf2dss(num,den,tau,Ts,Nmax,tol)
%
% mimoctf2dss   Computes a discrete-time state space model from a
%               continuous-time transfer function model (MIMO).
%
%               Let a MIMO system (p outputs, m inputs) be specified by a continuous-time
%               transfer function model
%
%                       Y(s) = G(s)*U(s)
%
%                       Y(s) = [y1(s) y2(s) ... yp(s)]'
%                       U(s) = [u1(s) u2(s) ... um(s)]'
%
%                              | g11(s) ... g1m(s) |
%                       G(s) = |   .         .     |
%                              | gp1(s) ... gpm(s) |
%
%                The individual transfer functions are of the form
%
%                       yi(s) = gij(s)*uj(s) 
%                             = hij(s)*exp(-tau_ij*s)*uj(s)
%                             = hij(s)*uj(s-tau_ij)
%
%                       hij(s) = num_ij(s)/den_ij(s)
%
%                       den_ij(s) = a_ij(0)*s^na_ij + ... + a_ij(na_ij-1)*s + a_ij(na_ij)
%                       num_ij(s) = b_ij(0)*s^nb_ij + ... + b_ij(nb_ij-1)*s + b_ij(nb_ij)
%
%                 The time delay, tau_ij, must be non-negative and the transfer
%                 function can only be realized if the denominator degree is
%                 at least the same size as the numerator degree.
%
%                 Assuming zero-order-hold with sampling time Ts
%
%                      u(t) = u(k*Ts)       k*Ts <= t < (k+1)*Ts
%               
%                 a discrete-time state space model
%
%                       x[(k+1)*Ts] = Ad*x[k*Ts] + Bd*u[k*Ts]
%                           y[k*Ts] = Cd*x[k*Ts] + Dd*u[k*Ts]
%
%                is computed by mimoctf2dss. 
%
%                The realization is conducted by computing the impulse
%                repsonse of the transfer function and doing a balanced
%                realization from the Hankel matrix of the impulse response
%                matrices. Nmax is the final impulse response matrix
%                computed and used in the construction of the Hankel
%                matrix. If Nmax >= 2*n, in which n is the state order of the
%                true system, then the realization will be exact within the tolerance, tol.
%                tol specifies the cut-off Hankel singular value
%                used to determine the order of the model.
%           
%
% Syntax: [Ad,Bd,Cd,Dd,sH]=mimoctf2dss(num,den,tau,Ts,Nmax,tol)
%
%           num             :   p-times-m cell array with the numerator polynomials
%           den             :   p-times-m cell array with the denominator polynomials
%           tau             :   p-times-m matrix with the time delays
%           Ts              :   Sampling time period
%           Nmax            :   Final sampling time (tfinal = Nmax*Ts)
%           tol             :   Hankel singular value cut off tolerance
%                           (gives the precision of the realization)
%           
%           (Ad,Bd,Cd,Dd)   :   Discrete-time state space matrices
%           sH              :   Hankel singular values
%

H = mimoctf2dimpulse(num,den,tau,Ts,Nmax);
[Ad,Bd,Cd,Dd,sH]=mimodimpulse2dss(H,tol);

