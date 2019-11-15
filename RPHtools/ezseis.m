function [ss,ttax,b]=ezseis(vel,dens,varargin);
%[ss,ttax]=ezseis(vel,dens,dx,dz,top,freq,snr);
%
% Quick and easy normal incidence synthetic seismic section.
% Calculates and displays wiggle trace and color plot of seismograms. 
% Inputs (reguired):  VEL, DENS, velocity and density.
% These can be vectors (e.g. well logs) or 2-D matrices of the same size.  
% Optional input parameters: DX, DZ, horizontal and vertical grid spacing,
% TOP, depth to top of velocity image (or log), FREQ, seismic frequency,
% SNR, Noise to Signal energy ratio (usually < 1) [Note: not S/N].
% When not specified on input, the default values are: DX, DZ =1,
% TOP = 0, FREQ = 25 (Hz), SNR = 0 (no random noise added).
% Outputs (optional): SS synthetic seismogram or seismic section, TTAX time axis.
% If VEL, DENS are vectors, the single seismogram is repeated 25 times 
% before plotting the section.
%                 
% The algorithm consists of calculation of reflectivity from impedance,
% depth to time conversion, low pass filtering the reflectivity sequence
% based on the value of FREQ, and horizontal averaging over a Fresnel zone.
%
% See also: KENNET, KENFRTT, PSPEC2DSH, BORN2D, BORNFILT, EZSEIS2,
%           QUICKSEIS
% For plotting seismic sections see SEISPLOT, SEISRWB

%Written by T. Mukerji, 1998
%Modified by G. Mavko, 2000

prompt={'dx', 'dz', 'depth to top', 'frequency', 'noise to signal ratio'};
def={'1', '1', '0','25','0'};

switch length(varargin)
case 0, answer=inputdlg(prompt, 'SRB EZSeis', 1, def);
	    dx=str2num(answer{1});dz=str2num(answer{2});top=str2num(answer{3});
		freq=str2num(answer{4});snr=str2num(answer{5});
case 1, dx=varargin{1};
		answer=inputdlg(prompt(2:end), 'SRB EZSeis', 1, def(2:end));
	    dz=str2num(answer{1});top=str2num(answer{2});
		freq=str2num(answer{3});snr=str2num(answer{4});
case 2, dx=varargin{1};dz=varargin{2};
		answer=inputdlg(prompt(3:end), 'SRB EZSeis', 1, def(3:end));
	    top=str2num(answer{1});
		freq=str2num(answer{2});snr=str2num(answer{3});
case 3, dx=varargin{1};dz=varargin{2};top=varargin{3};
		answer=inputdlg(prompt(4:end), 'SRB EZSeis', 1, def(4:end));
		freq=str2num(answer{1});snr=str2num(answer{2});
case 4, dx=varargin{1};dz=varargin{2};top=varargin{3};freq=varargin{4};
		answer=inputdlg(prompt(end), 'SRB EZSeis', 1, def(end));
		snr=str2num(answer{1});
case 5, dx=varargin{1};dz=varargin{2};top=varargin{3};freq=varargin{4};
        snr=varargin{5};
otherwise, error('Incorrect number of input arguments'); end;

if size(vel,1)==1, vel=vel(:); end;
if size(dens,1)==1, dens=dens(:); end;

imped=vel.*dens; 

% set the time to the top, using the top depth and mean of upper layer velocities
% time to bottom is then computed from the depth to time conversion

ttop=2*top./(mean(vel(1,:))); 
tt=cumsum(2*dz./vel)+ttop*ones(size(vel));

% time step is half the smallest raw time interval
dttu=0.5*min(min(diff(tt)));

% create equally spaced time steps
ttu=[min(tt(1,:)):dttu:max(tt(end,:))].';
% add extra time steps before and after to aid interpolation
tt=[(min(ttu)-dttu)*ones(1,size(tt,2));tt;(max(ttu)+dttu)*ones(1,size(tt,2))];

imped=[imped(1,:);imped;imped(end,:)];
impedtt=zeros(length(ttu),size(imped,2));

% sample the impedance at equally spaced times, and take difference
for k=1:size(impedtt,2), impedtt(:,k)=interp1q(tt(:,k),imped(:,k),ttu); end;
rz2tt=0.5*diff(log(impedtt));

fn=0.5/(ttu(2)-ttu(1));  %actual Nyquist from input data
fnn = 5.*freq;   % To minimize dead band, choose smaller Nyquist 
R = floor(fn/fnn);
if R>1, 
	for k = 1:size(rz2tt,2),rz2tto(:,k) = decimate(rz2tt(:,k), R, 'fir');  end;
fc=min(R*freq/fn,0.99); 
%%%%%%%%%ttu=[min(ttu):R*(ttu(2)-ttu(1)):max(ttu)]'; to be consistent with
%%%%%%%%%quickseism.m
else, rz2tto=rz2tt;
end;

b=fir1(9,fc);
ss=filtfilt(b,1,rz2tto);

dst=top+dz*0.5*size(vel,1); lmda=mean2(vel)/freq; frnlz=sqrt(dst*lmda);
boxn=max(2,floor(frnlz/dx)); bb=boxcar(boxn); bb=bb./sum(bb);

if size(ss,2)>1, ss=conv2(1,bb,ss,'same'); end;
%ss=filtfilt(bb,1,ss.').'; 
if size(ss,2)==1, ss=ss(:,ones(25,1)); end;

if snr~=0, ss=ss+(sqrt(snr))*std(ss(:))*randn(size(ss)); end;

ttax=[min(ttu):R*(ttu(2)-ttu(1)):min(ttu)+(length(ss(:,1))-1)*R*(ttu(2)-ttu(1))]';
xax = [0:dx:dx*(size(ss,2)-1)];
imagesc(xax, ttax, ss);
%set(gca,'xticklabel',num2str(dx*(str2num(get(gca,'xticklabel')))) );

