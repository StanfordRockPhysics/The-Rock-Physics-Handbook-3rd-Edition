function e = KozCarmE(Phi,d) 
% function e = KozCarmE(Phi,d) 
%
% The original Kozeny-Carman relation provides a way
% to estimate the permeability from porosity (Phi), 
% pore diameter (d),geometric factor (b).
%
% K=(1/72*Tau^2).*((d^2).*Phi.^3)/(1-Phi).^2)
% or
% K= (1000/450).*(d.^2).*((Phi.^3)./(1-Phi).^2));
%
% Where: Phi = porosity, fraction, 
%        Tau = tortuosity, Tau = 2.5, and 
%		  d   = Pore Diameter, micronmeter.
%
% The output is two column matrix, porosity in fraction
% in first column, and permeability (in md) in second
% column.
%
% ASSUMPTIONS: Based on flow through circular cross-section pipe,
%              and spherical grains.
%              Tau=2.5 
% LIMITATION: Homogenous and uniform system, unconsolidated.
% REFERENCE: "L'Ecoulement des Gaz a Travers les Milieux Poreux,
%             Bibliotheque des Science et Techniques Nucleaires,"
%             Carman, P. C., 1996, Presses Universitaires de
%             France, Paris, 198pp.
% See also BERNABE, BLOCH,COATDUM, COATES, FREDRICH,
%          MODKOZCARM, PANDALAKE,PANDALAKEKC,TIMUR,
%          TIXIER, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Pore Diameter (d), micronmeter'};
   defans={'[0.01:0.01:0.35]','250'};
   getpar=inputdlg(prompt,'Kozeny-Carman Model Default Input Parameters',1,defans);
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
   e = KC(param{1},param{2});  
else
e = KC(Phi,d);
end

function e = KC(Phi,d)
% permeability


K= (1000/450).*(d.^2).*((Phi.^3)./(1-Phi).^2);

%NK=(10^6).*K./(d.^2);
%semilogy(Phi,K,'r-');
%xlabel('Porosity','fontsize',9);
%ylabel('K(md)/D(mm)','fontsize',9);
%title('Kozeny-Carman K(md)/D(mm) vs. Phi','fontsize',10);


semilogy(Phi,K,'r-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Original Kozeny-Carman Permeability versus Porosity','fontsize',10);

grid on;
hold on;
e = [Phi K];
