function e = Owolabi(Phi,Swi) 
%function e = Owolabi(Phi,Swi) 
%
%Owolabi et al model provides a way to estimate 
%the permeability in uncosonlidated sands
%(of Pleistocene to Oligocene age)in Eastern 
%Niger Delta from log derived porosity (Phi)
%and irreducible water saturation (Swi). 
%
%For Oil Sand
%Koil = 307 + (26552.*(Phi.^2))-(34540.*(Phi.*Swi).^2)
%For Gas Sand
%Kgas = 30.7 + (2655.*(Phi.^2))-(3454.*(Phi.*Swi).^2);
%
%The output is three column matrix, porosity (fraction)
%in first column, permeability (md) for oil reservoir
%in second column, and permeability (md)for gas reservoir .
%in third column
%
%ASSUMPTIONS:Strongly dependent on the local lithology,
%            and the properties and distribution of 
%            reservoir fluids in the wells.
%VALIDITY:Permeability ranges from 203 to 4530md, 
%         porosity from 0.039 to 0.323, and water 
%         saturation from 0.113 to 0.824.
%REFERENCES:"An Empirical Expression for Permeability 
%           in Unconsolidated Sands of the Eastern Niger
%           Delta," O. O. Owolabi, T.F. LongJohn, and J.A.
%           Ajienka; Journal of Petroleum Geology, 
%           Vol.17(1), January 1994.
%       
%See also BERNABE, BLOCH, COATES, COATDUM, FREDRICH,
%         KOZCARM, MODKOZCARM, PANDALAKE,PANDALAKEKC,
%         TIMUR, TIXIER, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Irreducible Water Saturation (Swi)'};
   defans={'[0.1;0.15;0.2;0.25;0.3]','[1;0.9;0.85;0.8;0.75]'};
   getpar=inputdlg(prompt,'Owolabi Model Default Input Parameters',1,defans);
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
   e = Owo(param{1},param{2});  
else
e = Owo(Phi,Swr);
end
   
function e = Owo(Phi,Swi)    
% For Oil Sand
Koil = 307 + (26552.*(Phi.^2))-(34540.*(Phi.*Swi).^2);
% For Gas Sand
Kgas = 30.7 + (2655.*(Phi.^2))-(3454.*(Phi.*Swi).^2);

semilogy(Phi,Koil,'r-');
hold on
semilogy(Phi,Kgas,'b-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Owolabi Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi Koil Kgas];

