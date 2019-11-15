function [Iamp,Iphi,Ifreq]=iatrib(XX);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Iamp,Iphi,Ifreq]=iatrib(XX);
%
% This function calculates the instantaneous attributes of the image XX 
%Outputs: Iamp- Instantaneuos amplitude
% Iphi- Instantaneuos phase
% Ifreq- Instantaneuos frequency
%Input:  XX-  Seismic section [nt,nx]=size(xx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Written by R. Bachrach

ut=hilbert(XX);  %Matlab 'hilbert(x)' returns x+i*hilbert_trans_of_x
%ut=XX+i*hilbert(XX);
Iamp=abs(ut);
Iphi=imag(log(ut));
Ifreq=diff(Iphi);

if nargout==0
subplot(4,1,1), imagesc(XX), title('Input Section');
subplot(4,1,2), imagesc(Iamp),title('Instantaneous Amplitude');
subplot(4,1,3), imagesc(Iphi), title('Instantaneous Phase');
subplot(4,1,4), imagesc(Ifreq), title('Instantaneous Frequency');
colormap pink;
end;

