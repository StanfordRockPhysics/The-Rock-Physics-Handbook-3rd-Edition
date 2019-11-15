function [C,den]=hudsonF(cd,ar,Kfl,rofl,K,G,ro,s)
%function [C,den]=hudsonF(cd,ar,Kfl,rofl,K,G,ro,s)
%hudsonF - calculates the anisotropic elastic parameters for cracked rocks
%with a Fisher distribution of crack orientations around 3-axis.
%K,G bulk and Shear moduli
%ro density of the rock
%Kfl,rofl bulk modulus and density of the fluid inside the cracks
%cd crack density 
%ar aspect ratio
%s  standard deviation of Fisher distribution 
%C stiffness matrix
%den density of the cracked rock

%References: Hudson,1990, Geophys. J. Int. v102, 465-469; Rock Physics Handbook 

%written by Diana Sava; Jul 2000

lam=K-2/3.*G;
mu=G;
kapa=Kfl.*(lam+2.*mu)./(pi.*ar.*mu.*(lam+mu));
u3=4/3*(lam+2*mu)./((lam+mu).*(1+kapa));
u1=16/3.*(lam+2*mu)./(3.*lam+4.*mu);

ex=exp(1/(s.^2));
e11  =(-1+2.*s^2.*ex-2*s.^4.*(ex-1))./(2*(ex-1));
e1111=3/8*(-1+4*s.^4.*(2*ex+1)-24*s.^6.*ex+24*s.^8.*(ex-1))./(ex-1);

e12=0;e23=0;e13=0;
e22=e11;
e33=1-2*e11;
e2222=e1111;
e1122=e1111/3;
e1212=e1111/3;
e3333=8/3*e1111-4*e11+1;
e1133=e11-4/3*e1111;
e2233=e1133;
e1313=e1133;
e2323=e1133;

c1111=4*cd.*mu.*u1.*(e1111-e11)-cd./mu.*u3.*(lam.^2+4*lam.*mu.*e11 +4*mu.^2.*e1111);
c1122=4*cd.*mu.*u1.*e1122-cd./mu.*u3.*(lam.^2+2.*lam.*mu.*(e22+e11)+4*mu.^2.*e1212);
c1133=4*cd.*mu.*u1.*e1133-cd./mu.*u3.*(lam.^2+2.*lam.*mu.*(e33+e11)+4*mu.^2.*e1133);
c3333=4*cd.*mu.*u1.*(e3333-e33)-cd./mu.*u3.*(lam.^2+4*lam.*mu.*e33 +4*mu.^2.*e3333);
c2323=  cd.*mu.*u1.*(4*e2323-e22-e33)-4*cd./mu.*u3.*e2323;
c1313=  cd.*mu.*u1.*(4*e1313-e11-e33)-4*cd./mu.*u3.*e1313;
c1212=  cd.*mu.*u1.*(4*e1212-e11-e22)-4*cd./mu.*u3.*e1212;
c2222=c1111;
c2233=c1133;

%C=zeros(6);
C(1,1,:)=lam+2*mu+c1111;
C(1,2,:)=lam+c1122;
C(1,3,:)=lam+c1133;
C(2,2,:)=lam+2*mu+c2222;
C(2,3,:)=lam+c2233;
C(3,3,:)=lam+2*mu+c3333;
C(4,4,:)=mu+c2323;
C(5,5,:)=mu+c1313;
C(6,6,:)=mu+c1212;
C(2,1,:)=C(1,2);
C(3,1,:)=C(1,3);
C(3,2,:)=C(2,3);
phicrack=4*pi*ar./(3*cd);
den=(1-phicrack).*ro+phicrack.*rofl;
