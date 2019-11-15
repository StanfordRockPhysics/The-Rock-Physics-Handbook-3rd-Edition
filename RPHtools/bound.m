function [k_u,k_l,u_u,u_l,ka,ua]=bound(ib,f,k,u)
% [K_U,K_L,U_U,U_L,KA,UA]=BOUND(IB,F,K,U)
% calculate the elastic bounds (upper & lower) of an aggregate.
%
% input:    IB - type of bound (ib=0: Voigt-Reuss; ib=1: Hashin-Shtrikman);
%           F - volume fractions (<=1); K - bulk moduli; U - shear moduli.
% output:   K_U, K_L, U_U, U_L - elastic bounds of the aggregate.
%           KA, UA - arithmetic average of the upper and lower bounds
%                    (equals the Hill average for Hashin-Shtrikman bounds)
% 
% note:     1. VR bounds are the simplest;
%           2. HS bounds are the narrowest possible;
%           3. assumption: rock is isotropic.

% source:   Berryman, J.G., 1993, Mixture theories for rock properties
%           Mavko, G., 1993, Rock physics formulas
%  

% Written by Xingzhou 'Frank' Liu
% Modified by Isao Takahashi 4/27/99

lf=length(f);lk=length(k);lu=length(u);
if lf~=lk|lf~=lu,
	error('Input f, k, and u must have the same length')
end
if sum(f)~=1,
	error('F must sum up to 1')
end

if ib==0		% use Voigt-Reuss bounds

k_u=sum(f.*k);			% Voigt bound
k_l=1/sum(f./k);		% Reuss bound

u_u=sum(f.*u);			% Voigt bound
u_l=1/sum(f./u);		% Reuss bound

ka=(k_u+k_l)/2;			% Hill average
ua=(u_u+u_l)/2;

elseif ib==1		% use Hashin-Shtrikman bounds

c=4/3;

kmx=max(k);
kmn=min(k);
umx=max(u);
umn=min(u);

k_u=1/sum(f./(k+c*umx))-c*umx;	% HS upper bound
k_l=1/sum(f./(k+c*umn))-c*umn;	% HS lower bound

etamx=umx*(9*kmx+8*umx)/(kmx+2*umx)/6;
etamn=umn*(9*kmn+8*umn)/(kmn+2*umn)/6;

u_u=1/sum(f./(u+etamx))-etamx;	% HS upper bound
u_l=1/sum(f./(u+etamn))-etamn;	% HS lower bound

ka=(k_u+k_l)/2;			% simple arithmetic average
ua=(u_u+u_l)/2;

end

if nargout==0, 
disp('k: '),disp([k_u,k_l]),disp('u: '),disp([u_u,u_l]),
disp('ave: '),disp([ka,ua]), end;
