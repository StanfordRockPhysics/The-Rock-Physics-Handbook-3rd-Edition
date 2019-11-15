function e=CoatDum(Phi,Swr) 
% function e=CoatDum(Phi,Swr)
%
% Coates-Dumanoir's relation provides a way to estimate 
% permeability of sandstone from in situ measurements 
% of porosity(Phi)and residual fluid saturation (Swr). 
%
%  K=(A*Phi^B)/Swr^C
%where: K is the Permeability,
%       Phi is the porosity, scalar or vector, and 
%       Swr is the residual water saturation, scalar or vector.
%
% Coates-Dumanoir, Tixier, Coates, and Timur have the same form
% as above, but the coefficient A, B, C are different in each case.
%
% The output is two column matrix, porosity in fraction in first 
% column, and permeability (in md) in second column.
%
% ASSUMPTIONS: A,B,C are determined empirically 
% VALIDITY: Unconsolidated Sands
% REFERENCE:"A New Approach to Improved Log-Derived Permeability,"
%            by George R. Coates and J. L. Dumanoir, 1974
%            Log Analyst, Vol. 15, n. 1, p. 17-31, 
%            
%
% See also BERNABE, BLOCH, COATES, FREDRICH, KOZCARM, 
%          MODKOZCARM,PANDALAKE,PANDALAKEKC, TIMUR, TIXIER,
%          OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Irreducible Water Saturation (Swr)'};
   defans={'[0.01:0.01:0.35]','0.15'};
   getpar=inputdlg(prompt,'Coates-Dumanoir Model Default Input Parameters',1,defans);
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
   e = CoDu(param{1},param{2});  
else
e = CoDu(Phi,Swr);
end
   
function e = CoDu(Phi,Swr)    

K = 352.*(Phi.^4)./(Swr.^4);

semilogy(Phi,K,'b-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Coates-Dumanoir Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

