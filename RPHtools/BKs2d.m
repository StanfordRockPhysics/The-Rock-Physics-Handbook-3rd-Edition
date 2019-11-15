function Sdry=BKs2d(Ssat,K0,mu0,Kfl,phi)
% function  Sdry=BKs2d(Ssat,K0,mu0,Kfl,phi)
% Brown and Korringa saturated to dry fluid substitution calculation 
% for arbitrary anisotropic medium.  
%
% Inputs:
%    Ssat - 6x6 comliance matrix of the initially saturated rock, as defined 
%           in Mavko, Mukerji, and Dvorkin - Rock Physics Handbook
%    K0  -  Isotropic mineral bulk modulus
%    mu0 -  Isotropic mineral shear modulus
%    Kfl -  Initial fluid bulk modulus
%    phi -  porosity
% Output:
%    Sdry - 6x6 compliance matrix of the dry rock
%	
% see also BKs2s, BKc2c, BKd2s

% written by Gary Mavko 7/9/99
% based on derivation by Tapan Mukerji and Gary Mavko 
%

%
% compute compressibilities of the fluid, mineral, and saturated rock
beta0 = 1/K0;
betasat = sum(sum(Ssat(1:3,1:3)));
%
% make isotropic mineral compliance tensor
%
[S,C] = CSiso(K0,mu0);
%
factor = Kfl/(-Kfl*(betasat - beta0) + phi*(1 - Kfl*beta0));
%
Svect = sum(Ssat(1:3,:)) - sum(S(1:3,:));
Sdry = Ssat + factor*Svect'*Svect;
