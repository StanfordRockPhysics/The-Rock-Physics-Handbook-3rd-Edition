function [ku,kl,gu,gl,por]=hash(k1,mu1,k2,mu2)
%HASH - Hashin-Shtrikman upper and lower bound curves
%
%[KU,KL,GU,GL,POR]=HASH(K1,MU1,K2,MU2)
%	K1,MU1, K2,MU2:	Bulk and shear moduli of the two constituents
%	KU,KL, GU,GL: Upper and lower bounds on bulk and shear moduli
%       POR: volume fraction of material 2
%
%With no output arguments HASH plots the bounds as a function of porosity
%or fraction of phase 2 material. 
%
%See also BOUND

%Written by T. Mukerji

%********* HS upper and lower bounds **************
%
     	por=[0:0.01:1]; por(1)=1e-7;
	ku=k2+(1.-por)*(k1-k2)*(k2+4.*mu1/3.)./(k2+4.*(mu1/3.)+por*(k1-k2));
	kl=k2+(1.-por)*(k1-k2)*(k2+4.*mu2/3.)./(k2+4.*(mu2/3.)+por*(k1-k2));
	fgu=mu1*(9.*k1+8.*mu1)/(6.*(k1+2.*mu1));
	fgl=mu2*(9.*k2+8.*mu2)/(6.*(k2+2.*mu2));
	gu=mu2+(mu1-mu2)*(1.-por)*(mu2+fgu)./(mu2+fgu+por*(mu1-mu2));
	gl=mu2+(mu1-mu2)*(1.-por)*(mu2+fgl)./(mu2+fgl+por*(mu1-mu2));
%por=[0,por]; ku=[k1,ku]; kl=[k1,kl]; gu=[mu1,gu]; gl=[mu1,gu];
if nargout==0
plot(por,ku,'-g',por,kl,'-g',por,gu,'--c',por,gl,'--c','linewidth',1)
end;
%	gl=mu1*(1.-por)/(1.+por*mu1);
%	ro=(1.-por)*ro1+por*ro2;
%	vpu=sqrt((ku+(4.*gu/3.))./ro);
%	vpl=sqrt((kl+(4.*gl/3.))./ro);
%	vsu=sqrt(gu./ro);
%	vsl=sqrt(gl./ro);

