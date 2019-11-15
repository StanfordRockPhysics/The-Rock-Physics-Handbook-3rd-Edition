function e = FredrichE(Phi,d)
% function e = FredrichE(Phi,d)
%
% Quantitative prediction of permeability.
% Formation factor, F
%		F = Tau./Phi = (2.5)./Phi;
% or   F = satuarted rock resistivity/fluid resistivity.
% Permeability, K
%		K = (1./(b.*F)).*(Phi./Sv).^2;
%substitute Sv with 6.*(1-Phi)./d
%      K = (1000/450).*(d.^2).*(Phi.^3)./(1-Phi).^2;
% Where: Phi = porosity, scalar or vector, 
%        Tau = tortuosity, Tau^2 = 2.5 , 
%		  Sv  = pore surface area per unit sample volume,and
%        b   = shape factor, scalar or vector. 
%
% The output is two column matrix, porosity in fraction
% in first column, and permeability (in md) in second
% column.
%
% ASSUMPTIONS:A key assumption of the model is that
%             the flow paths for electrical current
%             and fluid are identical.The shape factor 
%             b is equal to 2 for circular tubes
%             and equal to 3 for cracks. 
% LIMITATIONS: Porosity higher than 10%.
% REFERENCE:"Pore Geometry and Transport Properties of 
%            Fontainebleau Sandstone," J.T. Fredrich,K.H.
%            Greaves, J.W. Martin,Int. J. Rock Mech. Min. 
%            Sci & Geomech. Abstr. Vol.30, No.7.  
%
% See also BERNABE, BLOCH,COATDUM, COATES, KOZCARM, 
%          MODKOZCARM, PANDALAKE,PANDALAKEKC, TIMUR,
%          TIXIER, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Pore Diameter (d), micronmeter'};
   defans={'[0.01:0.01:0.35]','100'};
   getpar=inputdlg(prompt,'Fredrich Model Default Input Parameters',1,defans);
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
   e = Fred(param{1},param{2});  
else
e = Fred(Phi,d);
end
   
function e = Fred(Phi,d)    
% Formation factor, F
F = (2.5)./Phi;
% Permeability, K
K = (1000/450).*(d.^2).*(Phi.^3)./(1-Phi).^2;

NK=(10^6).*K./(d.^2);
semilogy(Phi,NK,'m-');
xlabel('Porosity','fontsize',9);
ylabel('K(md)/D (mm)','fontsize',9);
title('Fredrich K(md)/D(mm) vs. Phi','fontsize',10);


%semilogy(Phi,,'m-')
%semilogy(Phi,K,'m-');
%set(gca,'fontsize',9);
%xlabel('Porosity','fontsize',9);
%ylabel('Permeability (md)','fontsize',9);
%title('Fredrich permeability versus Porosity','fontsize',10);

grid on;
hold on;
e = [Phi K];
