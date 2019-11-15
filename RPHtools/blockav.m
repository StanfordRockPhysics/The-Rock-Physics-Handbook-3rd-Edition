function [logb]=blockav(log, nb)
%logb = blockav(log,nb)
%
%Block average of input LOG over NB points. LOG can be a single column
%vector or multiple column matrix with Nan as missing values. 
%LOG is padded by the last row if the total number of rows in LOG is 
%not an integer multiple of NB. Output is the blocked log LOGB with
%the same number of rows as LOG.

% Written by T. Mukerji, 1998

if size(log,1)==1,log=log(:); end; [nr,nc]=size(log);

npad=nb-rem(nr,nb);
if npad~=0, lrow=log(end,:); logpad=[log;lrow(ones(npad,1),:)]; end;

logpad=reshape(logpad,nb,(nr+npad)/nb,nc); logpad=nanmean(logpad);
logpad=reshape(logpad(ones(nb,1),:,:),nr+npad,nc);
logb=logpad(1:nr,:);
