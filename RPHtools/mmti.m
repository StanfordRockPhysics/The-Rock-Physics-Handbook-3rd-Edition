function [ssuf]=mmti(ssd,ssdp)
% [SSUF]=MMTI(SSD,SSDP)
% calculate wet unrelaxed frame compliances from dry frame compliances for
% anisotropic (TI) rocks using the squirt model (Mukerji and Mavko, 1994)
%
% input:    SSD  - dry rock compliances at any pressures
%           SSDP - dry rock compliances at high pressure
%           (SSD=[s11 s12 s13 s33 s44]), (SSDP=[s11 s12 s13 s33 s44])
% output:   SSUF - unrelaxed wet frame compliances at any pressure
% note:     check applicability before use (remove linear trend
%           if necessary)

% Written by Frank Liu

for i=1:5,
dss(:,i)=ssd(:,i)-ssdp(i);				% DSijkl
end

dsaabb=2*(dss(:,1)+dss(:,2)+2*dss(:,3))+dss(:,4);	% DSaabb (1 column)
dsabab=2*dss(:,1)+dss(:,4)+4*dss(:,5)+4*(dss(:,1)-dss(:,2));
a=(dsabab./dsaabb-1)/4;					% alpha=(DSaabb-1)/2

for i=1:5,
tdss(:,i)=dss(:,i)./dsaabb;				% D~Sijkl
end

b=1-4*a;
g1=tdss(:,1)-(4*a./b).*(tdss(:,2)+tdss(:,3));
g2=tdss(:,2)./b;
g3=tdss(:,3)./b;
g4=tdss(:,4)-(8*a./b).*tdss(:,3);
g5=tdss(:,5)./b-(tdss(:,1)+tdss(:,4))./b/4+(g1+g4)/4;

gg=[g1 g2 g3 g4 g5];

for i=1:5,
ssuf(:,i)=ssd(:,i)-dsaabb.*gg(:,i);
end

