function [tt,rt,emtt]=kenfrtt(lyr,f)
%[TT,RT,EMTT]=KENFRTT(LYR,F)
%
%Traveltimes for normally incident wave propagation in layered media, taking
%into account all multiples and thin layer effects. 
%	LYR: [velocity, density, thickness], Nx3 matrix for N layers
%	F:   frequency (Hz)
%
%	TT: The exact traveltimes
%	RT: The ray theory traveltimes 
%	EMTT: The effective medium traveltimes
%
%Based on the Kennet-Frazer algorithm.
%
%See also KENNET

%Written by T. Mukerji

v=lyr(:,1); ro=lyr(:,2); d=lyr(:,3); m=ro.*v.^2;
n=length(v);
rt=cumsum(d./v);
v=[v;v(n)]; ro=[ro;ro(n)]; d=[d;d(n)];
s=1./v;
i=sqrt(-1); w=2*pi*f;
rru=0; l=0; tt=[]; rttt=0; sum=0; emtt=[];
for k=1:n
  den=ro(k)*v(k) + ro(k+1)*v(k+1);
  rd=(ro(k+1)*v(k+1)-ro(k)*v(k))/den; ru=-rd;
  td=2*sqrt(ro(k)*v(k)*ro(k+1)*v(k+1))/den; tu=td;
  theta=exp(i*w*d(k)/v(k));
  t1=td/(1-rru*theta^2*rd);
  rru=ru + t1*rru*theta^2*tu;
  sum=sum+ log(t1);
  l=l+d(k);
  sst=(1/(i*w*l))*sum;
  xtt=l*real(sst);
  tt=[tt;rt(k)+xtt];
  vemt=sqrt(harmmean(m(1:k))/mean(ro(1:k)));
  emtt=[emtt;l/vemt];
end

if nargout==0
 plot(1:n,rt,1:n,emtt,1:n,tt,'--','linewidth',1);
end
