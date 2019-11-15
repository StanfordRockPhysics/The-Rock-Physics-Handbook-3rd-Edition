function [pz,wz]=pgator(lyr,wvlt,dt,alpha)
% PGATOR - Synthetic seismograms for plane wave, normal incidence propagation 
% in 1-D layered media using propagator matrix method.
%
%[PZ,WZ]=PGATOR(LYR,WVLT,DT,ALPHA)
%
%	PZ:   seismogram at bottom of layers (pz=pressure)
%	WZ:   seismogram at the top (wz=particle velocity)
%	LYR:  [velocity, density, thickness] of layers
%	      LYR is a matrix of 3 columns with 
%             number of rows = number of layers
%	WVLT: input source wavelet vector (pressure)
%	      use [] to specify default wavelet (see also SOURCEWVLT)
%	DT:   time sampling interval of wavelet
%	ALPHA:attenuation coefficient, optional parameter (default alpha=0)
%
%Without any output arguments, PGATOR plots the transmitted and
%reflected seismograms

%Written by T. Mukerji

i=sqrt(-1);
if nargin<4, alpha=0; end;
if (length(wvlt)==0) wvlt=sourcewvlt; end; wvlt=[wvlt(:)]';
[nlr,nc]=size(lyr); ro=lyr(:,2); v=lyr(:,1); d=lyr(:,3);
n=length(wvlt); om=(2*pi/(n*dt))*[0:1:(n/2 -1)];
p0=ifft(wvlt); w0=(1/(ro(1)*v(1)))*p0;
a11=ones(size(om)); a21=zeros(size(om)); a12=a21; a22=a11;
     for j=1:nlr
        k=om./v(j); ck=k+i*alpha*k; wdv=d(j)*ck;
        c11=cos(wdv); c12=i*ro(j)*v(j)*sin(wdv);
        c21=(i/(ro(j)*v(j)))*sin(wdv); c22=c11;
        
        b11=c11.*a11 + c12.*a21; b12=c11.*a12 + c12.*a22;
        b21=c21.*a11 + c22.*a21; b22=c21.*a12 + c22.*a22;
        
        a11=b11; a12=b12; a21=b21; a22=b22;
    end
   rzvz=ro(nlr)*v(nlr); 
  pz=(rzvz*(a12.*a21-a11.*a22)./(a12-rzvz*a22)).*p0(1:length(om));
  wz=((rzvz*a21-a11)./(a12-rzvz*a22)).*p0(1:length(om));
     
  fltr=hanning(n/2); fltr=[ones(n/4,1);fltr(n/4+1:n/2)]'; 
  pz=pz.*fltr; wz=wz.*fltr;
  pz=[pz(1),pz(2:length(pz)),0,conj(fliplr(pz(2:length(pz))))];
  wz=[0,wz(2:length(wz)),0,conj(fliplr(wz(2:length(wz))))];
  
pz=real(fft(pz)); wz=real(fft(wz));
figure(1),plot([0:dt:dt*(length(pz)-1)],(pz)),title('transmitted seismogram');
figure(2),plot([0:dt:dt*(length(wz)-1)],wz),title('reflected seismogram');
