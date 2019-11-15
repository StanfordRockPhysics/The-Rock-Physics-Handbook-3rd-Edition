function e = Coates(Phi,Swr) 
% function e = Coates(Phi,Swr) 
% 
% Coates's relation provides a way to estimate the
% permeability of sandstone from in situ measurements 
% of porosity(Phi)and residual fluid saturation (Swr). 
%
%  K=(A*Phi^B)/Swr^C
% 
% where: K is the Permeability,
%        Phi is the porosity, scalar or vector, and 
%        Swr is the irreducible water saturation, scalar or vector.
%
% Coates, Tixier, Coates-Dumanoir, and Timur have the same form
% as above, but the coefficient A, B, C are different in each case.
%
% The output is two column matrix, porosity in fraction in first 
% column, and permeability (in md) in second column.
%
% ASSUMPTIONS: A,B,C are determined empirically 
% VALIDITY: Unconsolidated Sands
% REFERENCE:"Permeability Estimation: The Various Sources and Their
%            Interrelationships" Coates et al.,1991
%            JPT, vol.43, n. 5, p. 578-587.
%           "An Empirical Expression for Permeability
%            in Unconsolidated Sands of Eastern Niger
%            Delta,"Journal Of petroleum Geology, 
%            Vol. 17 (1),January 1994.
%
% See also BERNABE, BLOCH, COATDUM, FREDRICH, KOZCARM,
%          MODKOZCARM,PANDALAKE,PANDALAKEKC, TIMUR,
%          TIXIER, OWOLABI, REVIL,WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Irreducible Water Saturation (Swr)'};
   defans={'[0.01:0.01:0.35]','0.15'};
   getpar=inputdlg(prompt,'Coates Model Default Input Parameters',1,defans);
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
   e = Coat(param{1},param{2});  
else
e = Coat(Phi,Swr);
end
   
function e = Coat(Phi,Swr)    

K = (10000.*(Phi.^4).*(1-Swr).^2)./(Swr.^2);

semilogy(Phi,K,'r-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Coates Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

