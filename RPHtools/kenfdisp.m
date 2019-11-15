function [f,vel]=kenfdisp(lyr,f)
%function [f,vel]=kenfdisp(lyr,f)
%[F,VEL] = KENFDISP(LYR, F)
%
%Velocity dispersion for 1-D layered media, normal incidence plane wave
%using the Kennett-Frazer algorithm. This is the scattering dispersion,
%not intrinsic, viscoelastic dispersion.
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	F	frequency range (vector) e.g. obtained from LOGSPACE
%	VEL	velocity (vector of length = length of F).
%
%Without any output arguments, plots velocity versus frequency. 

%Written by T. Mukerji

v=lyr(:,1); ro=lyr(:,2); d=lyr(:,3); m=ro.*v.^2; dfr=d./sum(d);
n=length(v);
v=[v;v(n)]; ro=[ro;ro(n)]; d=[d;d(n)];
s=1./v;
rt=cumsum(d./v);
i=sqrt(-1); w=2*pi*f;
rru=0; l=0; tt=[]; rttt=0; xsum=0; emtt=[];
for k=1:n
  den=ro(k)*v(k) + ro(k+1)*v(k+1);
  rd=(ro(k+1)*v(k+1)-ro(k)*v(k))/den; ru=-rd;
  td=2*sqrt(ro(k)*v(k)*ro(k+1)*v(k+1))/den; tu=td;
  theta=exp((i*d(k)/v(k)).*w);
  t1=td./(1-rd.*rru.*theta.^2);
  rru= (ru + rru.*theta.^2)./(1-rd*rru.*theta.^2);
  t2=td./(1-rd.*rru.*theta.^2);
  l=l+d(k);
  xsum=xsum+  (1./(i.*w)).*(log(t1)) ;
  xtt=real(xsum);
  memt=sum(dfr(1:k).*(1./m(1:k))); memt=1./memt;
  vemt=sqrt(memt./sum(dfr(1:k).*ro(1:k)));
  emtt=[emtt;l/vemt];
  tt=[tt;rt(k)+xtt];
end
vel=l./tt(n,:); 
%vrt=l./rt(n), vemt=l./emtt(n),
if nargout==0
vrt=l./rt(n), vemt=l./emtt(n),
semilogx(f,vel); 
end
