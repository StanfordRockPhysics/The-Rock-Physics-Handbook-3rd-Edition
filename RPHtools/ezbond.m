function [cnew, M] = ezbond(c, theta)
% function [cnew, M] = ezbond(c, theta)
% calculates the bond transformation of the elastic stiffness
% matrix cij, for the simplest possible scenario where only one rotation is
% performed. 
% The old axes are x1 (horizontal), x2(horizontal), x3(vertical)
% the new axes are x1'(horizontal), x2'(horizontal), x3'=x3 (vertical)
% That is, the vertical axis stays the same; while the x1 axis rotates
% counter-clockwise for a certain degree.
%
% Inputs: theta --- the angle required to rotate x1 to x1'. 
%                                the unit is in degree.
% c: the stiffness tensor (cij, i,j=1 to 6) under the old coordinate system
% Output: cnew: the stiffness tensor under the new coordinate system.
%         M: (optional output) the Bond transformation matrix
% This transformation may be applied to TIH media, e.g., that caused by
% vertically aligned fractures. The angle theta, then, may be assigned to 
% be the angle between the fracture normal and a seismic line.
% To test the validity, assign to c an arbitrary 6X6 matrix, apply ezbond twice,
% first by degree phi, then by degree -phi, to see whether the original
% matrix is recovered.

% written by Haibin Xu, 2002
% modified by G. Mavko, 2003, reversing sign on theta, so sense of rotation
% is positive when theta > 0
% Reference: Rock physics handbook, chapter 1.4, "coordinate transformatins"

% calculate the cosines
theta = theta*pi/180;
b11 = cos(theta); b22 = cos(theta); b33 = 1;
b21 = cos(theta + pi/2); b12 = cos(pi/2 - theta);
b13 = 0; b31 = 0;
b23 = 0; b32 = 0;

% calculate the bond matrix
M1 = [b11^2, b12^2, b13^2; b21^2, b22^2, b23^2; b31^2, b32^2, b33^2];
M2 = [b12*b13, b13*b11, b11*b12; ...
      b22*b23, b23*b21, b21*b22;  ...
      b32*b33, b33*b31, b31*b32];
M3 = [b21*b31, b22*b32, b23*b33; ...
      b31*b11, b32*b12, b33*b13; ...
      b11*b21, b12*b22, b13*b23];
M4 = [b22*b33+b23*b32, b21*b33+b23*b31, b22*b31+b21*b32; ...
      b12*b33+b13*b32, b11*b33+b13*b31, b11*b32+b12*b31; ...
      b22*b13+b12*b23, b11*b23+b13*b21, b22*b11+b12*b21];
M = [M1, 2*M2; M3, M4];

% calculate the new stiffness tensor
cnew = M*c*M';

