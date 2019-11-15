function [S,C] = CSiso(K,mu)
% function  [S,C] = CSiso(K,mu)
% Creates 6x6 stiffness and compliance tensors for an isotropic rock  
%
% Inputs:
%    K  -   Isotropic rock bulk modulus
%    mu -   Isotropic rock shear modulus
% Output:
%    S  - 6x6 compliance matrix of the isotropic rock
%    C  - 6x6 stiffness  matrix of the isotropic rock
%   both as defined in Mavko, Mukerji, Dvorkin - Rock Physics Handbook
%	

% written by Gary Mavko 7/9/99
%

lambda = K - 2*mu/3;
c11 = lambda + 2*mu;
c12 = lambda;
c44 = mu;
%
C = zeros(6,6);
C(1,1) = c11;
C(1,2) = c12;
C(1,3) = c12;
C(2,1) = c12;
C(2,2) = c11;
C(2,3) = c12;
C(3,1) = c12;
C(3,2) = c12;
C(3,3) = c11;
C(4,4) = c44;
C(5,5) = c44;
C(6,6) = c44;
%
S = inv(C);
