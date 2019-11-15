function multiscat_property=colormarkerset(CONTROL_PROPERTY,SCAT_PROPERTY,nlist,skip)
%colormarkerset(CONTROL_PROPERTY,SCAT_PROPERTY,nlist,skip)
%		CONTROL_PROPERTY=
%		SCAT_PROPERTY=
%		nlist=
%		skip=

if nargin==3, skip=1; end;
faciescolor={'b','r','g','k','c','m','g','y'}';
faciesmarker={'h','+','o','*','s','V','d','-'}';

j=0;
for i=1:nlist;
j=j+1;
	if strcmp(CONTROL_PROPERTY.color,'auto');
		if (nlist==1)|((nlist==2)&(strcmp(CONTROL_PROPERTY.type,'line')))
		col='b';
     elseif (nlist>8), cmat=colormap;
		cind=round((i-1)/(nlist-1)*63+1); col=cmat(cind,:); 
     else 
		col=faciescolor{i};
		end;
		switch CONTROL_PROPERTY.type
			case {'line','scat'}
			field='color';
			case 'pdf'
			field='edgecolor';
		end
     SCAT_PROPERTY=setfield(SCAT_PROPERTY,field,col);
	end;
 if isfield(CONTROL_PROPERTY, 'marker')  %%% added 2005 to take care of pdfbayes contour color problems
	if strcmp(CONTROL_PROPERTY.marker,'auto');
 		if (nlist==1)|(nlist>5), marker='x';
		else    marker=faciesmarker{i};
        end
     SCAT_PROPERTY=setfield(SCAT_PROPERTY,'marker',marker);
   end;
 end
multiscat_property{i}=str2cell(SCAT_PROPERTY);
end

