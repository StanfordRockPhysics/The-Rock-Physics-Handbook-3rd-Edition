function simparam= monteccdf(Param,n)
%function simparam= monteccdf(Param,N)
%Monte-Carlo simulation of a single or multiple linearly correlated
%variables based on non-parametric cdf of primary variable.
%Input: Param: matrix of size [number_of_data x number_of_variables]
%              The primary variable is in the first column and all others
%              are simulated by conditioning to the draws of the first variable.
%              The secondary variables are drawn from the conditional 
%              non-parametric pdf. 
%              Param can be a single column vector.
%       N: number of monte-carlo draws to be simulated
%Output: simparam: [N x number_of_variables] matrix of simulated values.
%With no output argument, plots histogram and scatterplot of inputs
%and simulated data, using PLOTMATRIX.

%written by T. Mukerji, March 2003
%Based on earlier versions by P. Avseth and E. Gonzalez

if min(size(Param))==1, Param=Param(:); end;
[ndat npar]=size(Param);

simparam=zeros(n,npar);
r=rand(n,1);
[CDFpar, XX1]=makecdf(Param(:,1));
simparam(:,1)=interp1(CDFpar,XX1,r);

for k=2:npar
[ccdf, binpar1]=makeccdf(Param(:,1), Param(:,k));
simparam(:,k)=drawccdf(ccdf,binpar1, simparam(:,1));
end; %% end for k=2:npar

if nargout==0
figure; plotmatrix(Param); title('Input data');
figure; plotmatrix(simparam); title('Simulated data');
%figure; plot(XX1,CDFpar,'linewidth',3);
%figure;
%for k=1:length(binpar1), plot(ccdf(k).bin, ccdf(k).cdf), hold on; end;
%title('CCDFs')
end;

function [CDF, bins]=makecdf(x,bins)
%Makes cdf of input x vector

if nargin<2, bins=sort(unique(x)); end;
[N1,bins]=hist(x,bins);
CDF=cumsum(N1)/sum(N1);

indi=find(diff(CDF)==0);
CDF(indi)=[]; bins(indi)=[];

if CDF(1)~=0,
CDF=[0;CDF(:)];
bins=[bins(1)-0.01*abs(bins(1)); bins(:)];
end

function [ccdf, bin1]=makeccdf(x,y);
%make structure of conditional cdf for two random variables in vectors x and y

bin1=prctile(x,[5:5:100]);
[ccdf(1).cdf, ccdf(1).bin]=makecdf(y(x<bin1(1)));
for k=2:length(bin1)
[ccdf(k).cdf,ccdf(k).bin]=makecdf(y(x<bin1(k) & x>=bin1(k-1)));
end;

function simdat=drawccdf(ccdf,bin1,conddat)
% draw from ccdf based on conditioning data conddat

conddat=conddat(:);
n=length(conddat);
simdat=zeros(n,1);
r=rand(n,1);

ndx=conddat<bin1(1);
simdat(ndx)=interp1(ccdf(1).cdf,ccdf(1).bin,r(ndx));

for k=2:length(bin1), 
ndx=conddat<bin1(k) & conddat>=bin1(k-1); 
simdat(ndx)=interp1(ccdf(k).cdf,ccdf(k).bin,r(ndx));
end;
ndx=conddat>=bin1(end);
if sum(ndx)~=0 
simdat(ndx)=interp1(ccdf(end).cdf,ccdf(end).bin,r(ndx));
end

