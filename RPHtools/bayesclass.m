function [code,cprobs,maxprobs]=bayesclass(indata,smcpdf,vec);
%BAYESCLASS bayes classification
% INDATA: input data in a matrix (ndata x ncol) with one, two or three 
%        columns corresponding to the number of attributes.
% SMCPDF:estimated non-parametric multivariate probability density  
%        function; 4-D numerical array, [nXbin, nYbin, nZbin, nfacies+1].
%        When data attributes<3 the corresponding dimension is 
%        a singleton dimension. SMCPDF can be obtained from pdfbayes.m.
% VEC: cell array of the vectors defining the axes of SMCPDF.
% CODE: vector of facies code (ndata x 1).
% CPROBS: conditional probability for each facies (ndata x nfacies).
% MAXPROBS: maximum conditional probability (ndata x 1).
%
%see also PDFBAYES

%Written by C. Arroyo-Garcia, 8/2000

[nrow,ncol]=size(indata);
n=ones(nrow,3);
for i=1:nrow
   for j=1:ncol
      xx2=vec{j};
		xx2=xx2(:)'; 
      binwidth = [diff(xx2) 0];
      xx2 = [xx2(1)-binwidth(1)/2 xx2+binwidth/2];
      nbin2=length(xx2)-1;
   n(i,j)=sum(indata(i,j)>=xx2)+(sum(indata(i,j)>=xx2)==0) - (sum(indata(i,j)>=xx2)>nbin2);  
   cprobs(i,:)=squeeze(smcpdf(n(i,1),n(i,2),n(i,3),1:end-1))';
	end
end
[maxprobs,code]=max(cprobs,[],2);
