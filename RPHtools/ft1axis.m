function cx=ft1axis(adj,sig,cx);
% Fourier Transform in time axis(column-wise).
% fft shift corrected.
% Reference : J. Claerbout, Basic Earth Imaging

% 1999.12. Youngseuk Keehm

if(adj==0)
   cx(2:2:end,:) = -cx(2:2:end,:);
   if(sig==1) 	cx=ifft(cx);
   else 			cx=fft(cx);
   end
else
   if(sig==1) 	cx=fft(cx);
   else			cx=ifft(cx);
   end
   cx(2:2:end,:) = -cx(2:2:end,:);
end

      
