function [amp,phase,ds]=fftplot(data,dt)
%FFTPLOT plot amplitude and phase spectrum of time series.
%  Given real valued time series, data, with sampling interval,
%  dt, FFTPLOT plots amplitude and phase spectrum. FFTPLOT draws
%  only positve side since the spectra are even functions.
%
%  data:	Input real valued time series vector.
%  dt:		Sampling interval in seconds.
%  amp:		Amplitude spectrum (Magnitude of complex fft).
%  phase:	Phase spectrum in radian (Phase angle of complex fft).
%  ds:		Frequency step vector corresponding to AMP and PHASE.
%
%  See also FFT
%

%  written by Isao Takahashi, 1/20/1999.
%
nsample=length(data);
spectrum=fft(data);
sc=.5/dt;
ds=2*sc/nsample;
sindex=[0:ds:sc];
amp=abs(spectrum(1:length(sindex)));
phase=angle(spectrum(1:length(sindex)));
subplot(2,1,1)
plot(sindex,amp)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
subplot(2,1,2)
plot(sindex,phase)
xlabel('Frequency (Hz)')
ylabel('Phase (radian)')
