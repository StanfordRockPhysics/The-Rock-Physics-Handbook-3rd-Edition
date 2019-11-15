function berror=bayes(smcpdf,l,pfac);
%berror=bayes(smcpdf,pfac);
% Bayes' prediction and error estimation.
% smcpdf:	

nlist=length(pfac);
epsilon=1e-10;
sizesm=size(smcpdf);
ndim=ndims(smcpdf);

juncpfac=repmat(reshape(pfac,[1 1 1 nlist]),sizesm(1:ndim-1));
% pdf with prior included.
wpdf=smcpdf(:,:,:,1:nlist).*juncpfac;
% probability of predicting the mosth likely facies at each attribute.
prediction=max(wpdf,[],ndim);
junkmaxpdf=repmat(prediction,[1 1 1 nlist]);
flagpdf=(junkmaxpdf==wpdf).*(junkmaxpdf~=0);
nummaxprob=sum(flagpdf,ndim);

% For the case when two or more facies give the same maximum values.
while any(nummaxprob(:)>1)
[a1,a2]=find(nummaxprob>1);
for i=1:length(a1);
  for j=1:nlist
  tempwpdf=wpdf(:,:,:,j);
  tempwpdf(a1(i),a2(i))=tempwpdf(a1(i),a2(i))+epsilon*rand(1,1);
  wpdf(:,:,:,j)=tempwpdf;
  end;
end;
prediction=max(wpdf,[],ndim);
junkmaxpdf=repmat(prediction,[1 1 1 nlist]);
flagpdf=(junkmaxpdf==wpdf).*(junkmaxpdf~=0);
nummaxprob=sum(flagpdf,ndim);
end
junkl=repmat(reshape(l,1,1,1,nlist),sizesm(1:ndim-1));
% Predicted facies numnber
predictedfacies=sum(junkl.*flagpdf,ndim);
junkpredicted=repmat(predictedfacies,[1 1 1 nlist]);

%berror=1-sum(prediction(:));

for i=1:nlist
	for j=1:nlist
		tempwpdf=wpdf(:,:,:,j);
		berror(i,j)=sum(tempwpdf(predictedfacies==l(i)));
	end
end
