function e = PandaLakeKCE(Phi,Dp) 
%function e = PandaLakeKCE(Phi,Dp) 
%
%Panda & Lake provide a mathematical model to quantify
%the effects of cements on the single phase 
%permeability estimate of clastic rocks. 
%
% K = (Phi.^3)/(2*Tau*(av^2).*(1-Phi).^2)
%
% Introducing the statistical parameters of the 
% Particles size distribution (psd):
%  
% K = (B./114).*(Dp.^2).*(Phi.^3)/(1-Phi).^2;
%
%where: Phi = porosity, scalar or vector,
%       Tau = tortuosity Tau^2=2.5, and
%		   av  = specific surface area,
%       B   = particle size distribution parameter, and 
%       Dp  = mean particle size.
%
%The output is two column matrix, porosity in fraction
%in first column, and permeability (in md) in second
%column.
%
%VALIDITY: Homogeneous, Isotropic, Unconsolidated medium. 
%REFERENCE:"Basic Kozeny-Carman by Manmath N. Panda 
%          and Larry W. Lake," AAPG 1994.
%
%See also BERNABE, BLOCH,COATDUM, COATES, FREDRICH,
%         KOZCARM, MODKOZCARM,PANDALAKE, TIMUR,
%         TIXIER, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Mean Particle Size (d), micronmeter'}; 
   defans={'[0.01:0.01:0.35]','250'};
   getpar=inputdlg(prompt,'PandaLakeKC Model Default Input Parameters',1,defans);
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
   e = PLKC(param{1},param{2});  
else
e = PLKC(Phi,Dp);
end
   
function e = PLKC(Phi,Dp)

%Permeability of an Unconsolidated Medium 

K =3.34.*(Dp.^2).*((Phi.^3)./(1-Phi).^2);

semilogy(Phi,K,'b-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Basic Kozeny Carman by Panda & Lake Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

