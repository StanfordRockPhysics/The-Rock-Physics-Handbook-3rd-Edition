function flpropui
%FLPROPUI Reservoir fluid properties from Batzle-Wang relations. 
%Inputs and outputs are from interactive dialog boxes.
%
%See also: FLPROP

%Written by L. Teng, T. Mukerji, 1997

%[Kreuss,rhoeff,Kvoigt,vpb,rhob,Kb,vpo,rhoo,Ko,vpg,rhog,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,P,T,So,Sg)
% batzle - Batzle-Wang relations
% vpb,vpo,vpg: Vp of brine, oil, and gas in km/s
% rhob, rhoo, rhog: Density of brine, oil, and gas in g/cm3
% Kb, Ko, Kg: Bulk Moduli of brine, oil, and gas in GPa
% rhoeff: effective density of mixed fluid
% Kreuss: Reuss bound of mixed fluid's bulk modulus, homogeneous saturation
% Kvoigt: Voigt bount of mixed fluid's bulk modulus, patchy saturation
% method: 1, use Gas Index in Oil
%         other numbers, use GOR (L/L)
% sal: NaCl Sal (ppm)
% gg: Gas Gravity (Specific Gravity)
% og: Oil Gravity (API number)
% gor: Gas Oil Ratio (L/L)
% giib: Gas Index in Brine
% giio: Gas Index in Oil
% P: Pore Pressure (MPa)
% T: Rock Temperature (oC)
% So: Saturation of oil
% Sg: Saturation of gas
% Sb = 1-So-Sg: Saturation of water

% The results have been compared with Petrotools results.
% The errors are within 0.2%, and maybe caused by the precision in Petrotools.

prmt={'Method: 0 = use GOR (L/L), 1 = use Gas Index in Oil', ...
      'Brine Salinity (Nacl ppm)', ...
      'Oil Gravity (API)', ...
      'Gas Gravity  (Specific Gravity)',...
      'GOR  Gas Oil Ratio (L/L)', ...
      'Gas Index in Brine',...
      'Gas Index in Oil',...
      'Pore Pressure (MPa)',...
      'Temperature (Celsius)',...
      'S_o Oil Saturation',...
      'S_g Gas Saturation'};
defans={0, 30000, 25, 0.8, 300, 0, 0.1, 20, 80, 1, 0};
for k=1:length(defans),defans{k}=num2str(defans{k});end;

button=1;
outh=[];
while 1
 switch button
 case 0
%if ishandle(outh), delete(outh); end;
break;

otherwise

bwinpt=inputdlg(prmt,'Batzle-Wang inputs',1,defans);
defans=bwinpt;
button=prod(size(bwinpt));

if button ~= 0,
for k=1:length(bwinpt), bwinptnum(k)=str2num(bwinpt{k}); end;
method=bwinptnum(1); sal=bwinptnum(2); og=bwinptnum(3); gg=bwinptnum(4);
gor=bwinptnum(5); giib=bwinptnum(6); giio=bwinptnum(7); 
P=bwinptnum(8); T=bwinptnum(9); So=bwinptnum(10); Sg=bwinptnum(11);

% salinity is in ppm divided by 1e6
sal=sal./1e6;

% ideal gas constant
R=8.31441;
%R=8.31;

% gas density=============================================================
Pr=P./(4.892-0.4048.*gg);
Tr=(T+273.15)./(94.72+170.75.*gg);
E=0.109.*(3.85-Tr).^2.*exp(-(0.45+(8.*(0.56-1./Tr).^2)).*(Pr.^1.2./Tr));
Z=(0.03+0.00527.*(3.5-Tr).^3).*Pr+(0.642.*Tr-0.007.*Tr.^4-0.52)+E;
rhog=28.8.*gg.*P./(Z.*R.*(T+273.15)); 

% gas adiabatic bulk modulus==============================================
gamma=0.85+5.6./(Pr+2)+27.1./(Pr+3.5).^2-8.7.*exp(-0.65.*(Pr+1));
f=E.*1.2.*(-(0.45+8.*(0.56-1./Tr).^2).*Pr.^0.2./Tr)+(0.03+0.00527.*(3.5-Tr).^3);
Kg=P.*gamma./(1-Pr./Z.*f)./1000;
vpg=sqrt(Kg./rhog);


% oil density=============================================================
rho0=141.5./(og+131.5);

% calculate GOR for method 1
if (method ==1)
  gormax=2.03.*gg.*(P.*exp(0.02878.*og-0.00377.*T)).^1.205;
  gor=gormax.*giio;
end

% dead oil vs. live oil
if (gor==0)
  rhoog=rho0;
  rhop=rhoog+(0.00277.*P-1.71e-7.*P.^3).*(rhoog-1.15).^2+3.49e-4.*P;
  rhoo=rhop./(0.972+3.81e-4.*(T+17.78).^1.175);
else
  B0=0.972+0.00038.*(2.4.*gor.*sqrt(gg./rho0)+T+17.8).^1.175;
  rhoog=(rho0+0.0012.*gg.*gor)./B0;
  rhoo=rhoog+(0.00277.*P-1.71e-7.*P.^3).*(rhoog-1.15).^2+3.49e-4.*P;
end

% oil velocity============================================================
% live oil use pseudo density
if (gor~=0)
  rho0=rho0./B0./(1+0.001.*gor);
end
% the following formula is for dead oil only
%vpo=15450./sqrt(77.1+og)-3.7.*T+4.64.*P+0.0115.*(0.36.*sqrt(og)-1).*T.*P;
% the following formula is for dead oil and live oil
vpo=2096.*sqrt(rho0./(2.6-rho0))-3.7.*T+4.64.*P+0.0115.* ...
(4.12.*sqrt(1.08./rho0-1)-1).*T.*P;
vpo=vpo./1000;
Ko=vpo.*vpo.*rhoo;


% brine density===========================================================
rhow=1+1e-6.*(-80.*T-3.3.*T.^2+0.00175.*T.^3+489.*P-2.*T.*P+ ...
0.016.*T.^2.*P-1.3e-5.*T.^3.*P-0.333.*P.^2-0.002.*T.*P.^2);
rhob=rhow+sal.*(0.668+0.44.*sal+1e-6.*(300.*P-2400.*P.*sal+ ...
T.*(80+3.*T-3300.*sal-13.*P+47.*P.*sal)));

% brine velocity==========================================================
matrixw=[1402.85 4.871 -0.04783 1.487e-4 -2.197e-7
1.524 -0.0111 2.747e-4 -6.503e-7 7.987e-10
3.437e-3 1.739e-4 -2.135e-6 -1.455e-8 5.230e-11
-1.197e-5 -1.628e-6 1.237e-8 1.327e-10 -4.614e-13]';

% water velocity
velw=0;
for i=1:5
for j=1:4
  velw=velw+matrixw(i,j).*T.^(i-1).*P.^(j-1);
end
end

% gas water ratio
gwrmax=10.^(log10(0.712.*P.*(abs(T-76.71)).^1.5+3676.*P.^0.64)-4- ...
7.786.*sal.*(T+17.78).^(-0.306));
gwr=gwrmax.*giib;

% gas-free brine and brine with gas
vpb0=velw+sal.*(1170-9.6.*T+0.055.*T.^2-8.5e-5.*T.^3+2.6.*P- ...
0.0029.*T.*P-0.0476.*P.^2)+sal.^1.5.*(780-10.*P+0.16.*P.^2)-1820.*sal.^2;
vpb=vpb0./(sqrt(1+0.0494.*gwr));
vpb=vpb./1000;
Kb=vpb.*vpb.*rhob;

% fluid mixer=============================================================

% Effective density
Sb=1-So-Sg;
rhoeff = Sb.*rhob+So.*rhoo+Sg.*rhog;

% Reuss Average
if ( Kb.*Ko.*Kg~=0 )
  Kreuss=1./(Sb./Kb+So./Ko+Sg./Kg);
elseif (Sb==0)
  Kreuss=1./(So./Ko+Sg./Kg);
elseif (So==0)
  Kreuss=1./(Sb./Kb+Sg./Kg);
elseif (Sg==0)
  Kreuss=1./(Sb./Kb+So./Ko);
else
  Kreuss=0;
end

% Voigt Average
Kvoigt=Sb.*Kb+So.*Ko+Sg.*Kg;

outprmt={'Effective Fluid Bulk Modulus (Kreuss, GPa)',...
         'Effective Fluid density (g/cc)', ...
         'Patchy Fluid Bulk Modulus (Kvoigt, GPa)',...
         'Brine Vp (km/s)',...
         'Brine Density (g/cc)',...
         'Brine Bulk Modulus (GPa)',...
         'Oil Vp (km/s)',...
         'Oil Density (g/cc)',...
         'Oil Bulk Modulus (GPa)',...
         'Gas Vp (km/s)',...
         'Gas Density (g/cc)',...
         'Gas Bulk Modulus (GPa)',...
         'GOR (L/L)'};
defout={Kreuss,rhoeff,Kvoigt,vpb,rhob,Kb,vpo,rhoo,Ko,vpg,rhog,Kg,gor};
if isempty(outh),
outh=bwoutdlg(outprmt,'Batzle-Wang outputs',1,defout);
else
 TempHide=get(0,'ShowHiddenHandles');
  set(0,'ShowHiddenHandles','on');
txth=findobj(outh,'style','edit');
for j=1:length(txth), set(txth(j),'string',defout{end-j+1}); end;
set(0,'ShowHiddenHandles',TempHide);
end;
end; %end if button~=0
          
end; % end switch
end; % end while 1

