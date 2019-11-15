function simparam= monte(Param,n)
%function simparam= monte(Param,N)
%Monte-Carlo simulation of a single or multiple linearly correlated
%variables based on non-parametric cdf of primary variable.
%Input: Param: matrix of size [number_of_data x number_of_variables]
%              The primary variable is in the first column and all others
%              are simulated by linear correlation to the first variable
%              with the appropriate Gaussian error around the correlation.
%              Param can be a single column vector.
%       N: number of monte-carlo draws to be simulated
%Output: simparam: [N x number_of_variables] matrix of simulated values.
%With no output argument, plots histogram and scatterplot of inputs
%and simulated data, using PLOTMATRIX.

%written by T. Mukerji, July 2000
%Based on earlier versions by P. Avseth and E. Gonzalez

if nargin<=2, play=0; end;
if min(size(Param))==1, Param=Param(:); end;
[ndat npar]=size(Param);

%[N1,XX1]=hist(Param(:,1),100); %%slightly slower
XX1=sort(unique(Param(:,1)));
[N1,XX1]=hist(Param(:,1),XX1);
CDFpar=cumsum(N1)/sum(N1);

indi=find(diff(CDFpar)==0);
CDFpar(indi)=[]; XX1(indi)=[];

if CDFpar(1)~=0,
CDFpar=[0;CDFpar(:)];
XX1=[.99*XX1(1); XX1(:)];
end

%plot(XX1,CDFpar,'linewidth',3), pause
simparam=zeros(n,npar);
r=rand(n,1);
simparam(:,1)=interp1(CDFpar,XX1,r);

for k=2:npar
[Pa2,Sa2]=polyfit(Param(:,1),Param(:,k),1);	
par2pr2=Pa2(1,1)*Param(:,1)+Pa2(1,2);	
diffPar2=(Param(:,k)-par2pr2);	
errorPar2=std(diffPar2);

%r=rand(n,1);
simparam(:,k)=Pa2(1,1)*simparam(:,1)+Pa2(1,2)+errorPar2*randn(n,1);
end; %% end for k=2:npar

if nargout==0
figure; plotmatrix(Param); title('Input data');
figure; plotmatrix(simparam); title('Simulated data');
figure; plot(XX1,CDFpar,'linewidth',3);
end;

