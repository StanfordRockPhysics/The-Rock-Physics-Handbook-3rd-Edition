function [anis]=c2anis(cc)
% [ANIS]=C2ANIS(CC)
% get Thomsen's anisotropy parameters (e,r,d,dsv) from elastic constants
% for a TI medium 
%
% input:     CC=[c11 c33 c44 c66 c13]; 
% output:    ANIS=[epsilon gamma delta deltasv]; 

% Written by Frank Liu

c11=cc(:,1);
c33=cc(:,2);
c44=cc(:,3);
c66=cc(:,4);
c13=cc(:,5);

e=(c11-c33)./c33/2;
g=(c66-c44)./c44/2;
c3=c13+c44;
c4=c33-c44;
c5=c11-c44;
d=(c3.*c3-c4.*c4)./c33./c4/2;
dsv=(c5.*c4-c3.*c3)./c44./c4/2;

anis=[e g d dsv];

