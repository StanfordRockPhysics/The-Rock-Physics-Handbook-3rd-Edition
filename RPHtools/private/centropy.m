function entro=centropy(smcpdf,pfac);

siz=size(smcpdf);
nfac=siz(end)-1;
if nargin==1, pfac=ones(nfac,1)./nfac; end
tempstring=['smcpdf(',repmat(':,',[1,ndims(smcpdf)-1]),'i)'];
for i=1:length(pfac)+1;
	ent(i)=entropdf(eval(tempstring));
end
ent=ent(:);

entro=sum(ent(1:end-1)'.*pfac);		%       H(X1,X2|Y) or H(X1|Y)
entro(2)=ent(end);			%       H(X1,X2)   or H(X1)
entro(3)=entro(2)-entro(1);	  	%       I(X1,X2|Y) or I(X1|Y)
entro(4)=entropdf(pfac);			%       H(Y)
entro(5)=entro(4)-entro(3);       	%       H(Y|X1,X2) or H(Y|X1)
if entro(4)==0, entro(6)=entro(5);
else, entro(6)=entro(5)/entro(4);       %       H(Y|X1,X2)/H(Y)
end;
entro(7)=1-entro(6);
entro=entro(:);

