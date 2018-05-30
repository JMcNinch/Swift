function Iout=filt2d(Iin,type,r,c)
%FILT2D Filter a 2d matrix using min,max,sum,or avg
%   Iout = filt2d(Iin,type,r,c) uses a tile the size of r x c
%
%   Iout = filt2d(Iin,type,tile) uses a tile of ones and zeros input by the
%   user.
%
%   type options are
%   'min'
%   'max'
%   'avg'
%   'sum'

if nargin>3
    h=ones(r,c);%user input the number of rows and columns for the matrix
elseif nargin ==3
    h=r; %user input a different ones matrix
else
    error('NOT ENOUGH INPUT ARGUMENTS');
end

indx=isnan(Iin);
if strcmp(type,'max')
    Iin(indx)=-inf;
    h=strel(h);
    Iout=imdilate(Iin,h);
elseif strcmp(type,'min')
    Iin(indx)=inf;
    h=strel(h);
    Iout=imerode(Iin,h);
elseif strcmp(type,'minnans')
    Iin(indx)=inf;
    h=strel(h);
    Iout=imerode(Iin,h);
    indx=[];
elseif strcmp(type,'sum')
    Iin(indx)=0;
    Iout=imfilter(Iin,h);
elseif strcmp(type,'avg') || strcmp(type,'mean')
    k=imfilter(double(~indx),h);
    Iin(indx)=0;
    ISum=imfilter(Iin,h);
    Iout=ISum./k;
elseif strcmp(type,'avgnans')
    k=imfilter(double(~indx),h);
    Iin(indx)=0;
    ISum=imfilter(Iin,h);
    Iout=ISum./k;
    indx=[];
else
    error('INVALID FILTER TYPE: options are min,max,sum,avg,avgnans');
end
Iout(indx)=nan;
end
