function [ippn, ipsn, ispn, ipp, ips, isp]=eimp2(vp,vs,ro,theta_deg,K)
%[ippn, ipsn, ispn, ipp, ips, isp]=eimp2(vp,vs,ro,theta_deg,K)
%EIMP2 computes elastic (far offset) impedances for PP and PS waves.
%VP, VS, RO: inputs, scalar or vectors; P-velocity, S-velocity, and bulk density.
%THETA_DEG:   scalar angle(degrees) for EI. can be a column vector of length N
%          for varying angle with depth. This is the incident P-wave angle
%          for P-to-S wave, and is not the same as the reflected S-wave
%          angle. For P-to-P, incidence and reflection angles are the same.  
%K:  (optional) constant K (vs/vp)
%IPPN, IPSN, ISPN: output normalized and scaled elastic impedances for P-to-P, P-to-S and S-to-P waves.
%              Calculation uses mean VS/VP ratio. Normalization is done by dividing VP, VS, RO by
%             their means, and scaling the result by the zero offset impedance (Whitcombe, 2002).
%IPP, IPS, ISP:  optional outputs are the raw, non-normalized far-offset impedances.
%see also EIMP

%T. Mukerji, 1999; modified to include normalization 2002
%modified by Ezequiel Gonzalez, Nov. 2002

if nargin<5
    vsvp=mean(vs./vp);  %K=vs/vp ->constant value
    vpvs=mean(vp./vs);  %K=vs/vp ->constant value
else
    vsvp=K;
    vpvs=1/K;
end

theta=(pi/180).*theta_deg;
ip=ro.*vp; 
is=ro.*vs;

vpn=vp./mean(vp); vsn=vs./mean(vs); ron=ro./mean(ro);
ipn=mean(vp)*mean(ro);          %>>>included
isn=mean(vs)*mean(ro);          %>>>included

%vsvpsin2=(vs./vp).^2.*sin(theta).^2; 
vsvpsin2=vsvp.^2.*sin(theta).^2;        %>>>modified: vsvp

%IPP calculation<<<<<
x1=1+tan(theta).^2; x2=1-4*vsvpsin2; x3=-8*vsvpsin2;
ipp=vp.^x1 .* ro.^x2 .* vs.^x3;
%ippn= (ip).*(vpn.^x1 .* ron.^x2 .* vsn.^x3);
ippn= (ipn).*(vpn.^x1 .* ron.^x2 .* vsn.^x3);        %>>>modified: ipn

%IPS calculation<<<<<  theta is the INCIDENT ave (P-wave)
facroo=sqrt(vpvs.^2-sin(theta).^2);
facall=sin(theta)./(vpvs.*facroo);
a=facall.*(2.*sin(theta).^2 - vpvs.^2 - 2.*cos(theta).*facroo);
b=4.*facall.*(sin(theta).^2 - cos(theta).*facroo);
ips=ro.^a.*vs.^b;
%ipsn = (is).*(ron.^a.*vsn.^b);
ipsn = (ipn).*(ron.^a.*vsn.^b);        %>>>modified: is -> ipn

%ISP calculation<<<<<
x1 = 2*vsvpsin2 -1 -2*vsvp.*cos(theta).*sqrt(1 -vsvpsin2);
a=vsvp.*tan(theta).*x1;
x2 = vsvpsin2 - vsvp.*cos(theta).*sqrt(1 -vsvpsin2);
b=4*vsvp.*tan(theta).*x2;
isp=ro.^a.*vs.^b;
%ispn=(is).*(ron.^a.*vsn.^b);
ispn=(isn).*(ron.^a.*vsn.^b);        %>>>modified: isn

%PLOTS<<<<<
if nargout==0
lbl=['(' num2str(theta_deg) '^o)'];
    subplot(2,2,1), plot(is,ip,'o'); 
xlabel('I_S'); ylabel('I_P'); axis tight;
    subplot(2,2,2), plot(ippn,ip,'o'); 
xlabel(['I_{PP}' lbl ]); ylabel(['I_P']); axis tight;
    subplot(2,2,3), plot(ipsn,ip,'o'); 
xlabel(['I_{PS}' lbl ]); ylabel(['I_P']); axis tight;
    subplot(2,2,4), scatter(ipsn,ippn,'o'); 
xlabel(['I_{PS}' lbl]); ylabel(['I_{PP}' lbl]); axis tight;
end;
