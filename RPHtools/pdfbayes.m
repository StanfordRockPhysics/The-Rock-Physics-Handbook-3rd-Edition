function [bayesresult,smcpdf,vec,error,handle]=pdfbayes(indata,code,nbin,nfilt,obs,plotop,varargin);
%[bayesresult,smcpdf,vec,error,handle]=pdfbayes(indata,code,nbin,nfilt,obs,plotop,varargin)
%
%PDFBAYES: creates pdf plot and computes bayesian 
%	classification statistics from discrete data.
%  INPUTS:
%    INDATA: input data in a matrix with either 1,2, or 3 columns.
%        Columns for different types of data (attributes).
%    CODE: 1st column:facies (groups or categories) codes for INDATA. 
%        2nd column(optional):weight of each data point (usually 1).	
%    NBIN(optional): number of bins (scalar) for histogram, 
%        default=[10]; 
%        can be cell array with vectors defining the desired binning.
%        Example: {[X0:dX:Xn],[Y0:dY:Yn],[Z0:dZ:Zn]}, 
%        vectors for different axis can be of different lengths.
%    NFILT(optional): size of smoothing filter, default=[3]; 
%        usually nfilt << nbin.
%    OBS (optional):	scale factor for computing std. dev. of smoothing kernel
%        for each axis, default=[ 0.1]; 
%        std. dev. of smoothing kernel=(obs*std(indata)).
%    PLOTOP (optional): 0 or 1, choice for plotting outputs, 
%        plotop=1(default) draws plots and panels.		
%    VARARGIN (optional): arguments to define plot properties.
%
%  OUTPUTS (optional):
%    BAYESRESULT: structure with the following fields:
%	  total: total Bayes classification error.
%	  successrate: classification success rate.
%	  prob: joint probability of predicted and true facies, 
%	        p{predicted,true}.  
%	  condprob: conditional probability of true facies given 
%	        predicted facies, p{true|predicted}.
%	  entro: unconditional Shannon's entropy of facies, 
%	        H(facies).
%	  condentro: conditional Shannon's entropy of facies, 
%	        H(facies|attribute).
%	  condinfo: conditional Shannon's information of facies, 
%	        I(facies|attribute).
%	  condentronorm: differential entropy, H(facies|attribute)/H(facies).
%	  Each field contains cell arrays with dimensions {i}{j}, 
%	  depending on the number of attributes and their combinations.
%	  For 3 attributes (3 column indata):
%	  i=1, j=1,2,3; results for single attributes taken one at a time.
%	  i=2, j=1,2,3; results for attribute pairs: 1&2, 2&3, 1&3 respectively.
%	  i=3, j=1; results for all three attributes taken together.
%	  For 2 attributes (2 column indata):
%	  i=1, j=1,2; results for single attributes taken one at a time.
%         i=2, j=1;   results for the attribute pair taken together.
%	  For 1 attribute (1 column indata):
%	  i=1, j=1; results for the single attribute.
%    SMCPDF: estimated non-parametric multivariate probability density  
%        function; 4-D numerical array, [nXbin, nYbin, nZbin, nfacies+1]
%        when data attributes<3 the corresponding dimension is 
%        a singleton dimension.
%    VEC: cell array of the vectors defining the axes of SMCPDF.
%    ERROR: std. dev. of the smoothing kernel which can be 
%        also interpreted as the measurement error of the data.
%    HANDLE: structure with handles to plot objects.
%    tempstat.txt: output ascii file containig the computation results
%        written to disk when plotop=1(default).
%
%  EXAMPLES:
%     pdfbayes(indata,code,15,5); computes non-parametric pdf with 15 bins along 
%     each axis and smoothing filter size of 5 along each axis. Outputs are
%     displayed in plots and graphical interfaces.
%
%     [bayesresult,smcpdf]=pdfbayes(indata,code,[],[],[],0); computes 
%     non-parametric pdf and bayesian classification statistics using default
%     values of NBIN, NFILT, OBS, but does not make any output plot.
%     The results are returned in bayesresult and smcpdf as defined above.
%
% See also ENTROPDF, HIST2D. 

%
% Written by Isao Takahashi 1999
% Modifications and graphical interfaces 
% by C. Arroyo-Garcia, T. Mukerji, 7/2000
% Plotting (marker, linecolor) options in control_group changed to take
% care of conflicts with countourgroup handle in matlab7, Mukerji, 2005



if nargin==2
   obs=.1; nbin=10; nfilt=3; plotop=1;
elseif nargin==3
   nfilt=3; obs=.1; plotop=1;
elseif nargin==4
   obs=.1; plotop=1;
elseif nargin==5
   plotop=1;
end;
if isempty(obs), obs=.1; end;
if isempty(nbin), nbin=10; end;
if isempty(nfilt), nfilt=3; end;
if isempty(plotop), plotop=1; end;

%       Default parameter definition for scatter plot.
%%%%%%%%%%%%%%SCAT_PROPERTY=struct('linestyle','-','marker','none');
%%%%%%%%%%%%%%%CONTROL_PROPERTY=struct('marker','auto','color','auto','type','pdf');
%%%%%%%2 lines above changed to take care of conflicts with new contourgroup proprties in matlab7
SCAT_PROPERTY=struct('linestyle','-');
CONTROL_PROPERTY=struct('color','auto','type','pdf');

%       Save input plot properties in cell format in a structure array.
%       Overwrite default properties.
if ~isempty(varargin)
	for i=2:2:length(varargin)
        SCAT_PROPERTY=setfield(SCAT_PROPERTY,varargin{i-1:i});
	end
end;

%       Define control properties, which show if color and marker are 
%       given in input arguments or not.
if (isfield(SCAT_PROPERTY,'color')|isfield(SCAT_PROPERTY,'edgecolor')), CONTROL_PROPERTY=setfield(CONTROL_PROPERTY,'color','user'); end;
if isfield(SCAT_PROPERTY,'marker'), CONTROL_PROPERTY=setfield(CONTROL_PROPERTY,'marker','user'); end;

if isempty(code), code=indata(:,1)*0+1; end;
if size(code,2)==1, code(:,2)=1; end;

[ucode,i,j]=unique(code(:,1));
nlist=length(ucode);
for k=1:length(ucode); pfac(k)=sum(code(code(:,1)==ucode(k),2))/sum(code(:,2));
end

[smcpdf,vec,handle,error]=pdfgendraw(indata,code,obs,nbin,nfilt,plotop,CONTROL_PROPERTY,SCAT_PROPERTY);
[bayesresult]=pdfstat(smcpdf,ucode,pfac,[1:size(indata,2)],error);

ndim=size(indata,2);

%		  If plot option is activated the program calls here the
%		  different programs that create the results windows.
if plotop==1
	type tempstat.txt
	Std= std(indata); Mean=mean(indata);  
	ErrorPercen=error(2,:)*100;
	if ndim==1 
	  figure_1d(Std,Mean,ErrorPercen,bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
	elseif ndim==2
	  figure2d(Std,Mean,ErrorPercen,bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
    figure3d_1(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
    figure3d_2(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
    figure_12d(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
  elseif ndim==3
     figure3d_1(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
     figure3d_2(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
     figure3d_3(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
     figure_12d(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
	   figure_13d(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
	   figure_23d(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
	   figure_123(bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);
	   figure3d(Std,Mean,ErrorPercen,bayesresult.total,bayesresult.successrate,bayesresult.condprob,bayesresult.entro,bayesresult.condentro,bayesresult.condinfo,bayesresult.condentronorm);   
	end;   
end;

