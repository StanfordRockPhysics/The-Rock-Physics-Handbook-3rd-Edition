function [c,rho]=bkusc(f,vp,vs,den)
%        [c,rho]=bkusc(f,vp,vs,den)
%
% Backus average for the effective anisotropic (TI) stiffness, C
% of layered medium made of thin isotropic layers.
%
% input:   f       - volume fraction (layer thickness/total thickness)
%          vp      - p velocity
%          vs      - s velocity
%          den     - density
% output:  c       - 6x6 stiffness matrix as defined in The Rock Physics Handbook
%          rho     - volume average density
%
% See also C2VTI 

% Written by Gary Mavko 9/9/99

% force f to sum to unity, by normalizing
%
f = f/sum(f);
%
% compute isotropic lambda and mu for each thin layer
%
mu  = den.*vs.*vs;
lam = den.*vp.*vp - 2*mu;
%
% perform Backus averages
%
x = sum(f.*mu.*(lam+mu)./(lam+2*mu));
y = sum(f.*mu.*lam./(lam+2*mu));
z = sum(f.*lam./(lam+2*mu));
u = sum(f./(lam+2*mu));
v = sum(f./mu);
w = sum(f.*mu);
%
% assign to 6x6 format
%
A = 4*x + z*z/u;
B = 2*y + z*z/u;
E = 1/u;
F = z/u;
D = 1/v;
M = w;
%
c = zeros(6,6);
c(1,1) = A;
c(1,2) = B;
c(1,3) = F;
c(2,1) = B;
c(2,2) = A;
c(2,3) = F;
c(3,1) = F;
c(3,2) = F;
c(3,3) = E;
c(4,4) = D;
c(5,5) = D;
c(6,6) = M;
%
rho=sum(f.*den);
