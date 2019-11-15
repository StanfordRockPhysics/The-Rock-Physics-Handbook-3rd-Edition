function [vp,vs]=lm2v(lambda,mu,rho)
% [VP,VS]=LM2V(LAMBDA,MU,RHO)
% convert (LAMBDA,MU,RHO) to VP and VS
%
% input:   LAMBDA, MU, RHO
% output:  VP,VS

% Written by Frank Liu

vp=sqrt((lambda+2.*mu)./rho); vs=sqrt(mu./rho);
