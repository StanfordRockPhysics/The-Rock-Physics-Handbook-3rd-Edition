function [no,xo1,xo2,xo3] = hist3d(y,x1,x2,x3,weight)
%HIST3D  3 Dimensional Histogram.
%   N = HIST3D(Y) bins the elements of Y into 15 equally spaced bin
%   and returns the number of elements in each bin.  
%   N = HIST3D(Y, X1, X2, X3) uses the vectors X1, X2, X3 to
%   define the bin centers along the three axes.
%   
%   Y must be a three-column matrix to conduct three dimensional binning.
%
%   If y is a vector, [N,X1,X2,X3] = HIST3D(Y,X1,X2,X3) does the same
%   operation as [N,X1,X1,X1] = HIST(Y,X1).
%
%   If y is a two-column matrix, [N,X1,X2,X3] = HIST3D(Y,X1,X2,X3) does
%   the same operation as [N,X1,X2,X1] = HIST2D(Y,X1,X2).
%
%   N = HIST3D(Y,X1,X2,X3), where X1,X2, and X3 are scalars,
%   uses X1*X2*X3 bins.
%
%   N = HIST3D(Y,X1,X2,X3), where X1,X2, and X3 are vectors, returns
%   the distribution of Y among bins with centers specified by X1,X2, and X3.
%
%   N = HIST3D(Y,X1) does the same operation as N = HIST3D(Y,X1,X1,X1)
%
%   N = HIST3D(Y,X1,X2,X3,W), where W is a vector with the same length as Y,
%	returns N as the sum of the weights (W) in each bin.
%
%   [N,X1,X2,X3] = HIST3D(...) also returns the position of the bin centers
%   in X1,X2, and X3.
%
%   HIST3D(...) without output arguments produces a slice plot of
%   the result.
%
%   See also HIST,HIST2D.

%   J.N. Little 2-06-86
%   Revised 10-29-87, 12-29-88 LS
%   Revised 8-13-91 by cmt, 2-3-92 by ls.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.13 $  $Date: 1997/12/02 19:27:08 $

%   Created by Isao.Takahashi 1998/7/30

if nargin == 0
   error('Requires one or two input arguments.')
end


[m,n]=size(y);

if (nargin == 1), x1 = 15; x2 = 15; x3 = 15; end
if (nargin == 2)
	if length(x1)==1
		x2 = x1; x3 = x1; 
	elseif (length(x1)==2) & (n==2)
		temp=x1;
		x1 = temp(1); x2 = temp(2);
	elseif (length(x1)==3) & (n==3)
		temp=x1;
		x1 = temp(1); x2 = temp(2); x3 = temp(3);
	end
end
if (nargin <=4), weight = 1; end;
if n==1
	if nargout == 0
   		hist1d(y,x1,weight);
	else
   		[no, xo1] = hist1d(y,x1,weight);
		xo1 = xo1'; xo2 = xo1; xo3 = xo1;
	end

elseif n==2
	if nargout == 0
		hist2d(y,x1,x2,weight);
	else
		[no, xo1, xo2] = hist2d(y,x1,x2,weight);
		xo3 = xo2;
	end

elseif n==3
if length(x1) == 1
	nbin1 = x1;
   	miny1 = min(y(:,1)); maxy1 = max(y(:,1));
   	if miny1 == maxy1,
        	miny1 = miny1 - floor(x1/2) - 0.5; 
        	maxy1 = maxy1 + ceil(x1/2) - 0.5;
   	end
   	binwidth1 = (maxy1 - miny1) ./ x1;
   	xx1 = miny1 + binwidth1*(0:x1);
   	xx1(length(xx1)) = maxy1;
   	x1 = [xx1(1:length(xx1)-1) + binwidth1/2]';

	y1 = y(:,1);
	y1 = ceil((y1 - miny1)/binwidth1);
	y1 = y1 + (y1==0) - (y1>nbin1);
else
   	xx1 = x1(:)';
   	miny1 = min(y(:,1)); maxy1 = max(y(:,1));
	binwidth1 = [diff(xx1) 0];
   	xx1 = [xx1(1)-binwidth1(1)/2 xx1+binwidth1/2];
	xx1(1) = miny1; xx1(length(xx1)) = maxy1;

	nbin1 = length(xx1)-1;
	y1 = y(:,1);
	for k=1:length(y1);
		y1(k) = sum(y1(k)>=xx1);
	end
	y1 = y1 + (y1==0) - (y1>nbin1);
end

if length(x2) == 1
	nbin2 = x2;
	miny2 = min(y(:,2)); maxy2 = max(y(:,2));
   	if miny2 == maxy2,
        	miny2 = miny2 - floor(x2/2) - 0.5; 
        	maxy2 = maxy2 + ceil(x2/2) - 0.5;
   	end
	binwidth2 = (maxy2 - miny2) ./ x2;
   	xx2 = miny2 + binwidth2*(0:x2);
   	xx2(length(xx2)) = maxy2;
   	x2 = [xx2(1:length(xx2)-1) + binwidth2/2]';

	y2 = y(:,2);
	y2 = ceil((y2 - miny2)/binwidth2);
	y2 = y2 + (y2==0) - (y2>nbin2);
else
   	xx2 = x2(:)';
   	miny2 = min(y(:,2)); maxy2 = max(y(:,2));
   	binwidth2 = [diff(xx2) 0];
   	xx2 = [xx2(1)-binwidth2(1)/2 xx2+binwidth2/2];
   	xx2(1) = miny2; xx2(length(xx2)) = maxy2;

	nbin2 = length(xx2)-1;
	y2 = y(:,2);
	for k=1:length(y2);
		y2(k) = sum(y2(k)>=xx2);
	end
	y2 = y2 + (y2==0) - (y2>nbin2);
end

if length(x3) == 1
	nbin3 = x3;
	miny3 = min(y(:,3)); maxy3 = max(y(:,3));
   	if miny3 == maxy3,
        	miny3 = miny3 - floor(x3/2) - 0.5; 
        	maxy3 = maxy3 + ceil(x3/2) - 0.5;
   	end
	binwidth3 = (maxy3 - miny3) ./ x3;
   	xx3 = miny3 + binwidth3*(0:x3);
   	xx3(length(xx3)) = maxy3;
   	x3 = [xx3(1:length(xx3)-1) + binwidth3/2]';

	y3 = y(:,3);
	y3 = ceil((y3 - miny3)/binwidth3);
	y3 = y3 + (y3==0) - (y3>nbin3);
else
   	xx3 = x3(:)';
   	miny3 = min(y(:,3)); maxy3 = max(y(:,3));
   	binwidth3 = [diff(xx3) 0];
   	xx3 = [xx3(1)-binwidth3(1)/2 xx3+binwidth3/2];
   	xx3(1) = miny3; xx3(length(xx3)) = maxy3;

	nbin3 = length(xx3)-1;
	y3 = y(:,3);
	for k=1:length(y3);
		y3(k) = sum(y3(k)>=xx3);
	end
	y3 = y3 + (y3==0) - (y3>nbin3);
end

ytemp = (y3-1)*nbin1*nbin2 + (y2-1)*nbin1 + y1;
S = sparse(ytemp,1:length(ytemp),weight,nbin1*nbin2*nbin3,length(ytemp));
nn = reshape(full(sum(S')),[nbin1 nbin2 nbin3]);

th=90;
if nargout == 0
	[xx,yy,zz]=ndgrid(x1,x2,x3);
        handle=patch(isosurface(yy,xx,zz,nn,prctile(nn(:),th)));
        handle2=patch(isocaps(yy,xx,zz,nn,prctile(nn(:),th)));
	set(handle,'facecolor','r','edgecolor','none')
	set(handle2,'facecolor','interp','edgecolor','none')
        isonormals(yy,xx,zz,nn,handle);
	%daspect([1 1 .5]);
	view(3)
	camlight(0,0);
	lighting phong
else
   	no = nn; xo1 = x1; xo2 = x2; xo3 = x3;
end

elseif n~=3
	error('Requires input vector or matrix with two or three columns')

elseif isstr(y) | isstr(x1) | isstr(x2) | isstr(x3)
   	error('Input arguments must be numeric.')
end
