function [vpcr,vscr,rocr,mcr,kcr,mucr]=critpor(vp1,vs1,ro1,vp2,vs2,ro2,phicr)
%[vpcr,vscr,rocr,mcr,kcr,mucr]=critpor(vp1,vs1,ro1,vp2,vs2,ro2,phicr)
%Velocities, density and moduli at critical porosity PHICR. 
%Porosity is volume fraction of material 2.

%Written by T. Mukerji

m1=ro1.*vp1.^2; m2=ro2.*vp2.^2; mu1=ro1.*vs1.^2; mu2=ro2.*vs2.^2;
k1=m1-(4/3)*mu1; k2=m2-(4/3)*mu2;

mcr=(m1.*m2)./((1-phicr).*m2+phicr.*m1);
mucr=(mu1.*mu2)./((1-phicr).*mu2+phicr.*mu1);
kcr=(k1.*k2)./((1-phicr).*k2+phicr.*k1);
rocr=(1-phicr).*ro1+phicr.*ro2; vscr=sqrt(mucr./rocr);
vpcr=sqrt((kcr+(4/3)*mucr)./rocr);

