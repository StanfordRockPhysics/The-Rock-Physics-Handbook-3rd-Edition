function e = PandaLake(Phi,Tau,S,Dpm,Cdp) 
%function e = PandaLake(Phi,Tau,S,Dpm,Cdp) 
%
%The Manmath N. Panda and Larry W. Lake relation
%provides a quantitative, theoritical relationship 
%between single phase permeability and porosity(Phi),
%tortuosity (Tau),and statistics of a particle size
%distribution using a Modified Carman-Kozeny equation.
%
%K=((((s*Cdp^3)+(3*Cdp^2)+1)^2)/(1+Cdp^2)^2).*((
%   (Dpm^2).*Phi.^3)/(72*tau.*(1-Phi).^2))
%Where: s  = skewness of particles size distribution (psd),
%		 Cdp= coefficient of variation of psd,
%       Dpm= mean particles size (micron-m),
%       tau= tortuosity, scalar or vestor, and
%       Phi= porosity, scalar or vector.
%
%The output is two column matrix, porosity in fraction
%in first column, and permeability (in md) in second
%column.
%
%ASSUMPTIONS:Permeability is most sensitive to porosity 
%            followed by mean grain size, and sorting; Variation 
%            in mean grain size contributes more to permeability
%            variance (heterogeneity); F constant, lack of 
%            sensitivity of pore throat, pore body dimension.
%            It would be possible to combine particle size 
%            distribution information from drill cutting, cores,
%            or thin section in conjuction with porosity 
%            obtained from wireline logs to predict permeability.
%LIMITATIONS:Homogeneous, isotropic in a sens of representative 
%            elementary. Kozeny-Carmen corrected for sorting and 
%            skewness. 
%VALIDITY: The model is not valid in tigh sand.
%REFERENCE:"Estimation of Single-Phase permeability from 
%           Parameters of Particle-Size Distribution," 
%           Manmath N. Panda and Larry W. Lake, AAPG 1994.
%
% See also BERNABE, BLOCH, COATDUM, COATES, FREDRICH,
%          KOZCARM, MODKOZCARM, PANDALAKEKC, TIMUR, TIXIER,
%          OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Tortuosity (Tau)','Skewness of PSD (S)','Mean Grain Size (Dpm(, micron-m','Coefficient of Variation of PSD (Cdp)'};
   defans={'[0.05;0.1;0.2;0.25]','2','0.25','650','0.4'};
   getpar=inputdlg(prompt,'Panda & Lake Model Default Input Parameters',1,defans);
   if isempty(getpar)
	return;
   end
   for k=1:length(getpar),
       if isempty(str2num(getpar{k}))==1
          param{k} = evalin('base',getpar{k});
       else
          param{k} = str2num(getpar{k});
       end
   end
   e = PL(param{1},param{2},param{3},param{4},param{5});  
else
e = PL(Phi,Tau,S,Dpm,Cdp);
end
   
   function e = PL(Phi,Tau,S,Dpm,Cdp)
   
% Permeability calculation

K=((((S.*Cdp.^3)+(3.*Cdp.^2)+1).^2)./(1+Cdp.^2).^2).*(((Dpm.^2).*Phi.^3)./(72.*Tau.*(1-Phi).^2));

semilogy(Phi,K,'m-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Panda & Lake Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

