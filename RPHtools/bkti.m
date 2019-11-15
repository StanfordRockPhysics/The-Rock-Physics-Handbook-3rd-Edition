function [ssw]=bkti(phi,bfl,sso,ssd)
% [SSW]=BKTI(PHI,BFL,SSO,SSD)
% Calculate saturated anisotropic rock compliances from dry compliances using 
% Brown and Korringa relations
% note:   symmetry of mineral and dry rock TI or isotropic
%
% input:  fluid     - BFL (compressibility,  = 1/K = 1/(rVp^2) )
%         mineral   - SSO (compliance, =[s11o s12o s13o s33o s44o] )
%         dry rock  - SSD (compliance, =[s11d s12d s13d s33d s44d] )
%         porosity  - PHI
% output: sat. rock - SSW (compliance, =[s11w s12w s13w s33w s44w] )
%

% Written by Xingzhou 'Frank' Liu 1994


% mineral

s11o=sso(:,1);
s12o=sso(:,2);
s13o=sso(:,3);
s33o=sso(:,4);
s44o=sso(:,5);

s1ao=s11o+s12o+s13o;
s2ao=s1ao;
s3ao=2*s13o+s33o;
bo=2*(s11o+s12o+2*s13o)+s33o;
sso=bo;

% dry rock

s11d=ssd(:,1);
s12d=ssd(:,2);
s13d=ssd(:,3);
s33d=ssd(:,4);
s44d=ssd(:,5);

s1ad=s11d+s12d+s13d;
s2ad=s1ad;
s3ad=2*s13d+s33d;
ssd=2*(s11d+s12d+2*s13d)+s33d;

ssb=ssd-sso+(bfl-bo)*phi;

s11s=s11d-(s1ad-s1ao).*(s1ad-s1ao) ./ssb;
s12s=s12d-(s1ad-s1ao).*(s2ad-s2ao) ./ssb;
s13s=s13d-(s1ad-s1ao).*(s3ad-s3ao) ./ssb;
s33s=s33d-(s3ad-s3ao).*(s3ad-s3ao) ./ssb;
s44s=s44d;

ssw=[s11s s12s s13s s33s s44s];

