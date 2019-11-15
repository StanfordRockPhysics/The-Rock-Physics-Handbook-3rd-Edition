function [vpu,vpl,vsu,vsl,por]=hashv(vp1,vs1,ro1,vp2,vs2,ro2)
%function [vpu,vpl,vsu,vsl,por]=hashv(vp1,vs1,ro1,vp2,vs2,ro2)
%HASHV - Hashin-Shtrikman upper and lower bound curves for velocities
%	 as a function of fraction of material 2.
%VP1, VS1, RO1: P and S velocity and density of material 1.
%VP2, VS2, RO2: P and S velocity and density of material 2.
%VPU, VPL, VSU, VSL: Upper and lower bounds on P and S velocities.
%POR: Porosity vector (volume fraction of material 2).
%Assumes material 1 has higher velocity than material 2. If not, then
%upper and lower bounds should be interchanged in the output.
%With no output arguments HASH plots the bounds as a function of porosity
%or fraction of phase 2 material. 
%
%See also HASH, BOUND

%Written by T. Mukerji

%********* HS upper and lower bounds **************
%
mu1=ro1.*vs1.^2; mu2=ro2.*vs2.^2;
k1=ro1.*vp1.^2-(4/3).*mu1; k2=ro2.*vp2.^2-(4/3).*mu2;
     	por=[0.0:0.01:1];por(1)=1e-7;
	ku=k2+(1.-por)*(k1-k2)*(k2+4.*mu1/3.)./(k2+4.*(mu1/3.)+por*(k1-k2));
	kl=k2+(1.-por)*(k1-k2)*(k2+4.*mu2/3.)./(k2+4.*(mu2/3.)+por*(k1-k2));
	fgu=mu1*(9.*k1+8.*mu1)/(6.*(k1+2.*mu1));
	fgl=mu2*(9.*k2+8.*mu2)/(6.*(k2+2.*mu2));
	gu=mu2+(mu1-mu2)*(1.-por)*(mu2+fgu)./(mu2+fgu+por*(mu1-mu2));
	gl=mu2+(mu1-mu2)*(1.-por)*(mu2+fgl)./(mu2+fgl+por*(mu1-mu2));

%	gl=mu1*(1.-por)/(1.+por*mu1);
	ro=(1.-por)*ro1+por*ro2;
	vpu=sqrt((ku+(4.*gu/3.))./ro);
	vpl=sqrt((kl+(4.*gl/3.))./ro);
	vsu=sqrt(gu./ro);
	vsl=sqrt(gl./ro);

if nargout==0
plot(por,vpu,'-g',por,vpl,'-g',por,vsu,'--c',por,vsl,'--c','linewidth',1)
end;
