function e = BernabeE(Phi,CrackFraction,w,r)
%function BernabeE(Phi,CrackFraction,w,r)
%	
%Bernabe models permit to compute the permeability and porosity of strongly 
%pressure dependent pores such as cracks and approximately constant pores  
%associated with tubes and nodal pores.
%
% Crack Permeability
% Kcrack = (w.^2).*Phicrack./(12.*Taucrack.^2);
%        = (w^2).*(CrackFraction.*phi)./30       
% tube Permeability 
% Ktube = (r.^2).*Phitube./(8.*Tautube.^2);
%       = (r.^2).*((1-CrackFraction).*Phi)./20;
% K = Ktubes + Kcracks;
%
% where: w=width or aperture of the equivalent crack (micrometer),
%        Taucrack = crack tortuosity, Taucrack.^2 = 2.5, 
%        Phicrack = porosity of wetted cracks, a scalar or a vector,
%        r = radius of the tube (micrometer),
%        Tautube = tube tortuosity, Tautube.^2 = 2.5, and
%        Phitube = porosity of wetted tube, a scalar or a vector.
% 
%The output is four column matrix [Phicrack,Kcrack,Phitube,Ktube]
%
%LIMITATION: Clean sandstones, inappropriate in rocks with complex 
%            microstructure. 
%REFERENCE:"Pore geometry and pressure dependence of the transport properties 
%           in sandstones," Y Bernabe, Geophysics, Vol.56,No.4(April 1991).
%
% See also BLOCH,COATDUM, COATES, FREDRICH, KOZCARM, MODKOZCARM, PANDALAKE
%          PANDALAKEKC,TIMUR,TIXIER, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==5
   Ber1(CrackFraction,w,r); 
else    
   prompt={'Porosity (Phi)','Crack Volume Fraction in Pore Space','Crack Width (w), micrometer','Tube Radius (r), micrometer'};
   defans={'[0.01:0.01:0.35]','0.8','200','150'};
   getpar=inputdlg(prompt,'Bernabe1 Model Default Input Parameters',1,defans);
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
   Ber1(param{1},param{2},param{3},param{4});  
 end

function e = Ber1(Phi,CrackFraction,w,r)
%Porosity Vector
%Phi = [0:0.01:0.4];
%Crack Porosity
Phicrack = Phi.*CrackFraction;
% Crack Permeability
Kcrack = (w.^2).*Phicrack./30;
%Tube Porosity
Phitube = Phi - Phicrack;
% tube Permeability 
Ktube = (r.^2).*Phitube./20;
%Total perm is sum of kcrack and ktube
K = Kcrack+Ktube;

if nargout==0
   semilogy(Phi,K,'b-');
   grid on;
   %semilogy(Phi,Ktube,'g-', Phi, Kcrack,'b-');
   set(gca,'fontsize',9);
   xlabel('Porosity','fontsize',9);
   ylabel('Permeability (md)','fontsize',9);
   title('Bernabe Model:  Permeability vs Total Porosity', 'fontsize',10);
   hold on;
end   
e = [Phicrack,Kcrack,Phitube,Ktube];

