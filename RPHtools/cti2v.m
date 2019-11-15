function [vps,vss,vpf,vsf,e,g,d]=cti2v(C,ro)
%function [vps,vss,vpf,vsf,e,g,d]=cti2v(C,ro)
% vps, vss   - slow P and S wave velocities (along the crack normals)
% vpf, vsf   - fast P and S wave velocities (orthogonal to the cracks normals)
% e          - Thomsen's parameter epsilon
% g          - Thomsen's parameter gamma
% d          - Thomsen's parameter delta
% input arguments:
% ro         - density of the rock
% C          - Stiffness matrix

% written by Tapan Mukerji & Diana Sava; Jul 2000;


C33=squeeze(C(3,3,:));
C11=squeeze(C(1,1,:));
C44=squeeze(C(4,4,:));
C66=squeeze(C(6,6,:));
C13=squeeze(C(1,3,:));

C33=C33(:);
C11=C11(:);
C44=C44(:);
C66=C66(:);
C13=C13(:);
ro  =ro(:);

c11=max(C11,C33);
c33=min(C11,C33);
c66=max(C44,C66);
c44=min(C44,C66);

vpf=sqrt(c11./ro);
vps=sqrt(c33./ro);
vsf=sqrt(c66./ro);
vss=sqrt(c44./ro);

e=(c11-c33)./(2*c33);
g=(c66-c44)./(2*c44);
d=((C13+c44).^2-(c33-c44).^2)./(2*c33.*(c33-c44));


