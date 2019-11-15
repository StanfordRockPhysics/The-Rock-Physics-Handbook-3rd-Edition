function e = Bloch(Size,Sort,Content) 
% function Bloch(Size,Sort,Content) 
%
% Bloch models permit to predict porosity and permeability in 
% sandstones prior to drilling
%
% Phi(%) = -6.1 + (9.8./Sort) + (0.17.*Content);
% K(md)=10.^(-4.67 + (1.34.*Size)+(4.08./Sort)+3.42.*(Content./100));
% where: Size = Grain Size (mm),
%        Sort = Trask Sorting Coefficient(computed using sorting 
%               coefficient equation developed by Trask (1932)), and 
%        Content = Rigid Grain Content(%).
%
% The output is two column matrix, porosity (in percent) in first 
% column, and permeability (in md) in second column. 
%
% ASSUMPTIONS:Cement dissolution was not the dominant control of 
%             the reservoir quality. Rather, reservoir quality 
%             was controlled by decrease of porosity and permeability
%             with increasing time and temperature as indentified 
%             from the calibration data set.
% LIMITATION:Successful predictions clearly require understanding of
%            the key factors controlling reservoir quality.
%            Sensitive to the limits imposed by the calibration dataset.
%            Phi and K predictions by statistical model are limited to 
%            samples containing less than 5-10% pore-filling cements.
% REFERENCE:"Empirical Prediction of Porosity and Permeability in
%            sandstones, S Bloch," AAPG Bulletin, Vol 75, No 7 (July 1991).
%
% See also BERNABE, COATDUM, COATES, FREDRICH, KOZCARM, MODKOZCARM, 
%          PANDALAKE, PANDALAKEKC,TIMUR, TIXIER, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Grain Size (mm)','Trask Sorting Coefficient','Rigid Grain Content(%)'};
   defans={'[1.2;0.80;0.20;0.1]','[2;1.7;1.6;1.9]','[10;50;90;130]'};
   getpar=inputdlg(prompt,'Bloch Model Default Input Parameters',1,defans);
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
   e = Blo(param{1},param{2},param{3});  
else
e = Blo(Size,Sort,Content);
end

function e = Blo(Size,Sort,Content) 
% Porosity Calculation
Phi = -6.1 + (9.8./Sort) + (0.17.*Content);
% Permeability Calculation
K = 10.^(-4.67 + (1.34.*Size) + (4.08./Sort)+ 3.42.*(Content./100));

semilogy(Phi./100,K,'c-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Bloch Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];
