function [c, rho, por] = bkuslog(depth, vp, vs, den, phi)
% function  [C, denav, phiav] = bkuslog(depth, vp, vs, den, phi)
% Computes Backus Average from log data, assuming that each depth point
% corresponds to a thin isotropic layer.  Average is over the entire log
% given as input.
% Simple volume averages of density and porosity are also returned.
%
% Inputs:
%    depth  -   depth, in any units
%    vp, vs -   P and S wave velocities
%	 den    -	densities
%	 phi    -   porosities
% depth, vp, vs, den, phi are vectors of same size
% Output:
%    C      - 6x6 stiffness  matrix of the Backus Average
%          matrix as defined in Mavko, Mukerji, Dvorkin - Rock Physics Handbook
%    rho    - volume average density
%    por    - volume average porosity
%	
% Assumes that each depth point is the top of the thin layer.  
% Extrapolates another thin layer at the bottom with thickness of last layer.

% written by Gary Mavko 9/9/99
%

% compute vector of layer thicknesses.  Assume that each depth point is the
% top of the thin layer.  We extrapolate another thin layer at the bottom.
%

thick = diff(depth);
thick(end+1)=thick(end);
%
% divide by total thickness to get volume fractions.  This also
% returns the values back to positive numbers, in case depths were input
% in reverse order
%
tot = sum(thick);
thick = abs(thick/tot);
%
% call Backus Average routine
%
[c,rho]=bkusc(thick,vp,vs,den);
%
por = sum(thick.*phi);
