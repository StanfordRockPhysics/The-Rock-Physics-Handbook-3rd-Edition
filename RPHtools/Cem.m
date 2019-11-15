function result = Cem(varargin)
%function result = Cem(Phic,C,Gs,nus,Gc,nuC,Kf,scheme)
%Calculates the water-saturated elastic moduli versus porosity lines
%using Dvorkin's contact cement model
%Inputs are by input dialog box, if called without input arguments. 
%     PhiC: critical porosity
%     C:    coordination number
%     Gs:   grain shear modulus
%     nus:  grain Poisson's ratio
%     Gc:   cement shear modulus
%     nuC:  cement Poisson's ratio
%     Kf:   fluid bulk modulus
%     scheme: (1 or 2) Cementation scheme 1: cement at contact, 2: on surface
%Outputs returned in result matrix.
%result=[porosity, M-modulus, shear-modulus];
%Plots effective moduli vs. porosity when called with no output arguments.


%Written by Jack Dvorkin
%I/O modifications, Tapan Mukerji, 6/1999.
%ref. Dvorkin & Nur, 1996, Geophysics, 61, 1363-1370
%     Rock Physics Handbook, section 5.2

prompt={'PhiC','Coord.#','Gs (GPa)','Nus','Gc (GPa)','Nuc','Kf (GPa)','Cem. Scheme (1 or 2)' };
defans={'.38','8.5','45','.064','45','0.064','0.','2'};

if nargin==0
getpar=inputdlg(prompt,'Contact Cement Model',1,defans);
for k=1:length(getpar), param(k)=str2num(getpar{k}); end;
PhiC=param(1); C=param(2); G=param(3); nu=param(4); Gc=param(5); nuC=param(6); 
Kf=param(7); schopt=param(8);
else
PhiC=varargin{1}; C=varargin{2}; G=varargin{3}; nu=varargin{4};
Gc=varargin{5}; nuC=varargin{6}; Kf=varargin{7}; schopt=varargin{8};
end;
format short

K = G.*2.*(1.+nu)./(3.*(1.-2.*nu));
Kc = Gc.*2.*(1.+nuC)./(3.*(1.-2.*nuC));
%Porosity loop Cement
i = (1:100)';
      Phi0 = PhiC-(i-1).*(PhiC-0.15)./100;
%Fraction of cement in the rock 
      fc = PhiC-Phi0;
%Fraction of grain in the solid 
      fgs = (1-PhiC)./(1-Phi0);
%Fraction of cement in the solid
      fcs = (PhiC-Phi0)./(1-Phi0); 
%Bulk modulus of the solid
      Ks = (fgs.*K+fcs.*Kc+1./(fgs./K+fcs./Kc))./2;  
%Shear modulus of the solid 
      Gs = (fgs.*G+fcs.*Gc+1./(fgs./G+fcs./Gc))./2;  
%M-modulus of the solid 
      Ms = Ks+4.*Gs./3.; 
%Kframe and Gframe at Phi=Phi0, Cementation Scheme 1
      if schopt==1
      a = 2.*(((PhiC-Phi0)./(3.*C.*(1.-PhiC))).^0.25);
      end;
%Cementation Scheme 2
      if schopt==2
      a = sqrt((2.*(PhiC-Phi0))./(3.*(1.-PhiC)));
      end;
%Capital Lambdas
      alam = (2./3.14).*(Gc./G).*(1.-nu).*(1.-nuC)./(1.-2.*nuC);
      alamtau = (1./3.14).*(Gc./G);
%Effective bulk modulus
      r1 = Kc+4.*Gc./3;
      r2 = C.*(1.-PhiC)./6;
      r3 = -0.024153.*(alam.^(-1.3646)).*(a.^2)+0.20405.*(alam.^(-0.89008)).*a+0.00024649.*(alam.^(-1.9864));
      Kframe = r1.*r2.*r3;
%Effective shear modulus
      r1 = Gc;
      r2 = 3.*C.*(1.-PhiC)./20;
      a1t = -0.01.*(2.2606.*nu.*nu+2.0696.*nu+2.2952);
      a2t = 0.079011.*nu.*nu+0.17539.*nu-1.3418;
      b1t = 0.05728.*nu.*nu+0.09367.*nu+0.20162;
      b2t = 0.027425.*nu.*nu+0.052859.*nu-0.87653;
      c1t = 0.0001.*(9.6544.*nu.*nu+4.9445.*nu+3.1008);
      c2t = 0.018667.*nu.*nu+0.4011.*nu-1.8186;
      r3 = (a1t.*(alamtau.^a2t)).*a.*a+(b1t.*(alamtau.^b2t)).*a+c1t.*(alamtau.^c2t);
      Gframe = 0.6.*Kframe + r1.*r2.*r3;
      Mframe = Kframe+4.*Gframe./3;
%Gassmann
      Ksat = Ks.*(Phi0.*Kframe-(1+Phi0).*Kf.*Kframe./Ks+Kf)./((1-Phi0).*Kf+Phi0.*Ks-Kf.*Kframe./Ks);   
      Msat = Ksat+4.*Gframe./3;
	  nuSat = 0.5.*(Msat./Gframe-2)./(Msat./Gframe-1);
	  VpVs = (Msat./Gframe).^0.5;
      Msat2 = Ms.*(Phi0.*Mframe-(1+Phi0).*Kf.*Mframe./Ms+Kf)./((1-Phi0).*Kf+Phi0.*Ms-Kf.*Mframe./Ms);

if nargout==0
subplot(1,2,1)
plot(Phi0,Msat,'r-')
axis([0.1 0.4 0 45])
set(gca,'fontname','bookman','fontsize',9)
xlabel('Porosity','fontname','bookman','fontsize',11)
ylabel('M-Modulus (GPa)','fontname','bookman','fontsize',11)
hold on
subplot(1,2,2)
plot(Phi0,Gframe,'r-')
axis([0.1 0.4 0 20])
set(gca,'fontname','bookman','fontsize',9)
xlabel('Porosity','fontname','bookman','fontsize',11)
ylabel('G-Modulus (GPa)','fontname','bookman','fontsize',11)
hold on
end;
result=[Phi0 Msat Gframe];

