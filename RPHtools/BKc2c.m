function Csat2=BKc2c(Csat1,K0,mu0,Kfl1,Kfl2,phi)
% function  Csat2=BKc2c(Csat1,K0,mu0,Kfl1,Kfl2,phi)
% Brown and Korringa  fluid substitution calculation for arbitrary 
% anisotroic medium, from arbitrary initial fluid to arbitrary final fluid.  
% Stiffness tensor inputs and outputs.
%
% Inputs:
%    Csat1 - 6x6 stiffness matrix of the initially saturated rock, as defined 
%            in Mavko, Mukerji, and Dvorkin - Rock Physics Handbook
%    K0 -    Isotropic mineral bulk modulus
%    mu0 -  Isotropic mineral shear modulus
%    Kfl1 -  Initial fluid bulk modulus
%    Kfl2 -  New fluid bulk modulus
%    phi  -  porosity
% Output:
%    Csat2 - 6x6 stiffness matrix of the rock saturated with the new fluid
%	
% see also BKs2s, BKs2d, BKd2s

% written by Gary Mavko 7/9/99
% based on derivation by Tapan Mukerji and Gary Mavko
%

Ssat1 = inv(Csat1);
Sdry = BKs2d(Ssat1,K0,mu0,Kfl1,phi);
Ssat2 = BKd2s(Sdry,K0,mu0,Kfl2,phi);
Csat2 = inv(Ssat2);
