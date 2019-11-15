function [code,cprobs,maxprobs]=bayesclass(indata,smcpdf,vec);
%
%
[nrow,ncol]=size(indata);
n=ones(nrow,3);
for i=1:nrow
	for j=1:ncol
   n(i,j)=sum(indata(i,j)>vec{j})+(sum(indata(i,j)>vec{j})==0);  
   cprobs(i,:)=squeeze(smcpdf(n(i,1),n(i,2),n(i,3),1:end-1))';
	end
end
[maxprobs,code]=max(cprobs,[],2);
