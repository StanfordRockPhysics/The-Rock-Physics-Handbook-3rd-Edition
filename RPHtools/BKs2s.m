function Ssat2=BKs2s(Ssat1,K0,mu0,Kfl1,Kfl2,phi)
% function  Ssat2=BKs2s(Ssat1,K0,mu0,Kfl1,Kfl2,phi)
% Brown and Korringa  fluid substitution calculation for arbitrary 
% anisotroic medium, from arbitrary initial fluid to arbitrary final fluid.  
% Compliance tensor inputs and outputs.
%
% Inputs:
%    Ssat1 - 6x6 comliance matrix of the initially saturated rock, as defined 
%            in Mavko, Mukerji, and Dvorkin - Rock Physics Handbook
%    K0 -    Isotropic mineral bulk modulus
%    mu0 -  Isotropic mineral shear modulus
%    Kfl1 -  Initial fluid bulk modulus
%    Kfl2 -  New fluid bulk modulus
%    phi  -  porosity
% Output:
%    Ssat2 - 6x6 compliance matrix of the rock saturated with the new fluid
%	
% see also BKc2c, BKs2d, BKd2s

% written by Gary Mavko 7/9/99
% based on derivation by Tapan Mukerji and Gary Mavko
%
%
Sdry = BKs2d(Ssat1,K0,mu0,Kfl1,phi);
Ssat2 = BKd2s(Sdry,K0,mu0,Kfl2,phi);
