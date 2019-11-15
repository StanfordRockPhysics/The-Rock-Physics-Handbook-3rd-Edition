function [outcpdf] = cpdf(data,vec,nfac,prob,pfac)

ndim=length(vec);
for i=1:nfac;
[outhist,junkvec]=histnd(data(prob(:,1)==i,1:ndim),vec,prob(prob(:,1)==i,2));
outcpdf(:,:,:,i)=outhist/sum(outhist(:));
end;
rep=ones(1,3);
for i=1:ndim;
rep(i)=length(vec{i});
end
tempmat=repmat(reshape(pfac,1,1,1,nfac),rep);
outcpdf(:,:,:,nfac+1)=sum(outcpdf(:,:,:,1:nfac).*tempmat,4);
