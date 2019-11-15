function [datastr,data,colnames,header]=loadlas(file)
%[DATASTR,DATA,COLNAMES,HEADER]=LOADLAS(FILE)
% Load LAS (Log Ascii Standard) formatted well log file
%
% FILE (optional): Input filename of file in LAS format.
%       FILE should be character array within single quotes 'file'.
% Without input arguments displays dialogbox for file name.
% DATA: the numeric data matrix.
% COLNAMES: names of the columns as specified in file (cell array)
% HEADER: Header of the file as string matrix.
% DATASTR: Structure array containing HEADER and the COLNAMES as its fields.

% Written by Isao Takahshi, Tapan Mukerji 12/4/1999

if nargin==0, [file,pathname]=uigetfile('*.las');  else pathname=''; end;
if ischar(file)
fid=fopen([pathname,file]);
k=0; header=''; colblock=0; datablock=0; 

while ~datablock
jline=fgetl(fid);
if ~ischar(jline), break, end;
    header=strvcat(header,jline); 
if any(strmatch('~A',upper(jline))), colblock=0; datablock=1; end;
if colblock&any(strmatch('~',upper(jline))), colblock=0; end;
      if colblock & ~strcmp(jline(1),'#'),
         k=k+1; colnames{k}=strtok(jline);
      end;
if any(strmatch('~CURVE',upper(jline))), colblock=1; datablock=0; end;
end;
disp(' ');
disp(['Reading data for ' num2str(k) ' log curves from ' file]) 
data=fscanf(fid,'%f');
data=reshape(data,[k length(data)/k])';
fclose(fid);
data(data<=-999)=nan;

colnames=regexprep(colnames,'\W','_');
%colnames=strrep(colnames,'.','_');
%colnames=strrep(colnames,'-','_'); 
colnames=lower(colnames);

datastr.header=header;
for k=1:length(colnames)
datastr = setfield(datastr,colnames{k},data(:,k));
end;
end;
