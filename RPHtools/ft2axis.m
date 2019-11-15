function cx=ft2axis(adj,sig,cx);
% Fourier Transform in space axis(row-wise).
% fft shift corrected.
% Reference : J. Claerbout, Basic Earth Imaging

% 1999.12. Youngseuk Keehm

if(adj==0)
   cx(:,2:2:end) = -cx(:,2:2:end);
   if(sig==1) 	cx=ifft(cx,[],2);
   else 			cx=fft(cx,[],2);
   end
else
   if(sig==1) 	cx=fft(cx,[],2);
   else			cx=ifft(cx,[],2);
   end
   cx(:,2:2:end) = -cx(:,2:2:end);
end
