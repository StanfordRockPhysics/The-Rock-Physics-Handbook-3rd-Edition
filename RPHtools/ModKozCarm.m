function e = ModKozCarm(Phi,d,B,Phic) 
% function e = ModKozCarm(Phi,d,B,Phic) 
%
% The Modified-Kozeny-Carman relation provides 
% a way to estimate the permeability of a porous 
% medium in terms of generalized parameters such
% as porosity (Phi), pore diameter (d),geometric 
% factor (B),and percolation porosity (Phic).
%
% K = (B*d.^2).*((Phi-Phic).^3)./(1+Phic-Phi).^2]
% where: Phi = porosity, scalar or vector,
%        Phic= percolation porosity, and
%        B   = geometric factor,scalar or vector.
%
% The output is two column matrix, porosity in fraction
% in first column, and permeability (in md) in second
% column.
%
% ASSUMPTIONS:Based on flow through circular cross-section
%             pipe, 2<Tau<3, Kozeny-Carman constant=5.
%             Phi=pi*(R^2)/A, S=2*pi*R/A
% LIMITATION: Homogenous and uniform medium; distribution 
%             of pore size is small; saturated rock;small 
%             length scale compared to the average pore size;
%             unconsolidated. 
%REFERENCE:  "Kozeny-Carman Relation for Flow", Mavko,Mukerji,
%             Dvorkin,The Rock Physics Handbook,1998
%
% See also BERNABE, BLOCH,COATDUM, COATES, FREDRICH,
%          KOZCARM, PANDALAKE,PANDALAKEKC,TIMUR, TIXIER, 
%          OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Pore Diameter (d), micronmeter','Geometric Factor (B)','Percolation Porosity (Phic)'};
   defans={'[0.05;0.1;0.2;0.25]','60','2','0.02'};
   getpar=inputdlg(prompt,'Modified Kozeny-Carman Default Input Parameters',1,defans);
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
   e = MKC(param{1},param{2},param{3},param{4});  
else
e = MKC(Phi,d,B,Phic);
end
   
function e = MKC(Phi,d,B,Phic) 
% connected porosity
Phix = Phi - Phic; 
% permeability

K = (B*d.^2).*((Phix.^3)./(1-Phix).^2);

semilogy(Phi,K,'b-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Modified-Kozeny-Carman Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

