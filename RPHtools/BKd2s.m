function Ssat=BKd2s(Sdry,K0,mu0,Kfl,phi)
% function Ssat=BKd2s(Sdry,K0,mu0,Kfl,phi)
% Brown and Korringa dry to saturated fluid substitution calculation 
% for arbitrary anisotropic medium.  
%
% Inputs:
%    Sdry - 6x6 comliance matrix of the initially dry rock, as defined in 
%           Mavko, Mukerji, and Dvorkin - Rock Physics Handbook
%    K0 -   Isotropic mineral bulk modulus
%    mu0 -  Isotropic mineral shear modulus
%    Kfl -  New fluid bulk modulus
%    phi -  porosity
% Output:
%    Ssat - 6x6 compliance matrix of the saturated rock
%	
% see also BKs2s, BKc2c, BKs2d

% written by Gary Mavko 7/9/99
%

%
% compute compressibilities of the fluid, mineral, and dry rock
beta0 = 1/K0;
betadry = sum(sum(Sdry(1:3,1:3)));
%
% make isotropic mineral compliance tensor
%
[S,C] = CSiso(K0,mu0);
%
factor = Kfl/(Kfl*(betadry - beta0) + phi*(1 - Kfl*beta0));
%
Svect = sum(Sdry(1:3,:)) - sum(S(1:3,:));
Ssat = Sdry - factor*Svect'*Svect;
