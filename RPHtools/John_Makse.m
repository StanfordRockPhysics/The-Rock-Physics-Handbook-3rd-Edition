function [Vp1c,Vp3c,sigma1c,sigma3c,N,C]=John_Makse(mu,poisson,grainsize,phi,epsilon,e3,rho,Cn)
%function [Vp1c,Vp3c,sigma1c,sigma3c,N,C]=John_Makse(mu,poisson,grainsize,phi,epsilon,e3,rho,Cn)
%
% Norris and Johnson's model (1997) for uniaxial strain in a random sphere pack
% using eqs. of Johnson et al. 1998 and also using the Makse et al.(1999) 
%correction  (i.e. coordination dependent with pressure)
%
% Inputs: 
%			mu=shear modulus of the grains
%			poisson=poisson of the grains
%			grainsize=diameter of the grains
%			n = number of contacts (number of coordination)
%			phi = initial porosity
%			epsilon = hydrostatic strain (assuming negative compression)
%			e3 = axial strain
%			rho = density
%           Cn=normal stiffness in the spheres as an adjustment parameter
%           originally must be Cn= (4* mu)/(1-poisson)
% Outputs:
%           Vp1c=Vp perpendicular to the applied stress direction
%           Vp3c=Vp in the direction of the applied stress
%           sigma1c=induced stress perpendicular to the applied stress
%           sigma3c=applied stress
%           N=coordination number as function of stress
%           C=Cijkl anisotropic stiffness tensor
%
% See also Johnson
%
%references: 
%Norris, A. N., and Johnson, D. L., 1997, Nonlinear elasticity of granular
%media: ASME Journal of Applied Mechanics, 64, 39-49.
%
%Johnson, D. L., Schwartz, L. M., Elata, D., Berryman, J.G., Hornby, B., and
%Norris, A. N., 1998, Linear and nonlinear elasticity of granular media:
%stress-induced anisotropy of a random sphere pack: Transactions of the ASME,
%65, 380-388.
%
%Makse, H. A., Gland, N., Johnson, D. L., Schwartz, L. M., 1999, Why
%effective medium theory fails in granular materials: Physical Review Letter,
%83, 5070-5073.


% Written by S. Vega, 2001

%
% stress determination 
%
% Lamme constant
lamda = mu * ( (2*poisson)/(1-(2*poisson)));
B = (1/(4*pi))*( (1/mu) + (1/(lamda+mu)));
C = (1/(4*pi))* ( (1/mu) - (1/(lamda+mu)));
sigma3c = - (((-e3).^(3/2))*(1-phi)*Z*(3*B+C))/((6*pi^2)*B*(2*B+C)) ; 
sigma1c = - (((-e3).^(3/2))*(1-phi)*Z*C)/((24*pi^2)*B*(2*B+C)) ; 
%
precision=1e50;
tolerance=1e-10;
Z=6;
% Coordination number estimation
while precision > tolerance
n = Z + (((-sigma3c-2*sigma1c)/3)/6e4).^(1/3);
nn=n;
sigma3c = - (((-e3).^(3/2)).*n*(1-phi)*(3*B+C))/((6*pi^2)*B*(2*B+C)) ;
sigma1c = - (((-e3).^(3/2)).*n*(1-phi)*C)/((24*pi^2)*B*(2*B+C)) ; 
before3=sigma3c; before1=sigma1c;
n = Z + (((-sigma3c-2*sigma1c)/3)/6e4).^(1/3);
sigma3c = - (((-e3).^(3/2)).*n*(1-phi)*(3*B+C))/((6*pi^2)*B*(2*B+C)) ;
sigma1c = - (((-e3).^(3/2)).*n*(1-phi)*C)/((24*pi^2)*B*(2*B+C)) ; 
precision = abs(before3 - sigma3c);
end
%
N=n;
% grain radius estimation
radius=grainsize/2;
%
% normal and tangential stiffness in the spheres (grains)
%Cn= (4* mu)/(1-poisson);
%Cn=9.5e10;
Ct= (8* mu)/(2-poisson);
%
% transversely isotropic case
%
% stiffness tensor Cijkl
%
% some arregment terms
%gamma = ( 3*n*(1-phi)*((-epsilon)^(1/2)) )/ (4*(pi^2)*Bw*(2*Bw+Cw));
gamma = (3/32)*n.*Cn*Ct*(1-phi)*(-epsilon)^(1/2);
alfa = sqrt(epsilon./e3);
Bw = 2/(pi*Cn);
Cw = (4/pi)*( (1/Ct) - (1/Cn) );
Io = (1/2)*( sqrt(1+alfa.^2) + (alfa.^2).*log( (1+sqrt(1+alfa.^2))./alfa) );
I2 = (1/4)*( ((1+(alfa.^2)).^(3/2)) - (alfa.^2).*Io);
I4 = (1/6)*( ((1+(alfa.^2)).^(3/2)) - 3*(alfa.^2).*Io);
%
% C11 = C1111
C11 = (gamma./alfa).*( (2*Bw*(Io-I2)) + ((3*Cw/4).*(Io-2*I2+I4)) );
%
% C13 = C1133
C13 = (gamma./alfa).* (Cw*(I2-I4));
%
% C33 = C3333
C33 = (gamma./alfa).* (4*Bw*I2 + 2*Cw*I4);
%
% C44 = C2323
C44 = (gamma./alfa).* ( ((Bw/2)*(Io+I2)) + Cw*(I2-I4));
%
% C66 = C1212
C66 = (gamma./alfa).*(Bw*(Io-I2) + ((Cw/4)*(Io-2*I2+I4))) ;
%
C=zeros(6,6,length(C11));
C(1,1,:)=C11; C(2,2,:)=C11; C(3,3,:)=C33;
C(1,2,:)=C12; C(1,3,:)=C13; C(2,3,:)=C13;
C(4,4,:)=C44; C(5,5,:)=C44; C(6,6,:)=C66;

% calculattion of the corresponding velocities
%
Vp3c = (sqrt(2*C33))/((2*rho)^(1/2))
% Vp1 = Vp2
Vp1c = (sqrt(2*C11))/((2*rho)^(1/2))
%

