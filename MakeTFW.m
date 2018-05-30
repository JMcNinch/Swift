%%%%%MAKE A TFW file%%%%%%%%%%
%%%%%%%%STEPS%%%
%%%1) save image from GE with know positions of the upper left and lower right
%%%2) import image in irfanview and trim to exact corner edge plus record
%%%the number of pixels in x and y directions
%%%3) save trimmed image as a *.tif
%%%4) convert geographic corners to rectilinear coordinate (e.g. UTM)
%%http://www.rcn.montana.edu/resources/converter.aspx
%%%5) calculate the difference in easting and northing between corners
%%%6) determine pixel size in x and y directions by dividing delta distance
%%%by the number of pixels
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%TEMPLATE%%%%%%%%%%
% % LatDMS=[ULdeg ULmin ULsec; LRdeg LRmin LRsec];
% % LonDMS=[ULdeg ULmin ULsec; LRdeg LRmin LRsec];
% 
%%%%%%%%%%%%%%%%%USER INPUT%%%
%%%upper left corner%%
ULlatdeg=36;
ULlatmin=55;
ULlatsec=33.79;
ULlondeg=75;
ULlonmin=59;
ULlonsec=56.13;
%%%%%%%%%%%Lower Right corner
LRlatdeg=36;
LRlatmin=55;
LRlatsec=0.92;
LRlondeg=75;
LRlonmin=58;
LRlonsec=54.83;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoXpixels=651;
NoYpixels=426;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%END USER INPUT%%%%%%%%%%%%%%%

LatDD=[ULlatdeg+(ULlatmin+ULlatsec/60)/60; LRlatdeg+(LRlatmin+LRlatsec/60)/60];
LonDD=-[ULlondeg+(ULlonmin+ULlonsec/60)/60; LRlondeg+(LRlonmin+LRlonsec/60)/60];

[x,y,~] = deg2utm(LatDD,LonDD);


A=abs(diff(x))./NoXpixels;%%pixel size in x direction (delta easting/#pixels in x direction)
B=0; %rotation about the y axis
C=0; %rotation about the x axis
D=abs(diff(y)./NoYpixels).*-1; %pixel size in y direction typically negative (delta northing/#pixels in y direction)
E=x(1);% x coordinate in rectilinear coordinate system
F=y(1);% y coordinate in rectilinear coordinate system

BaseTFW=[A; B; C; D; E; F];

fileID = fopen('Base.tfw','w');
fprintf(fileID,'%3.8f\r\n %1.1f\r\n %1.1f\r\n %3.8f\r\n %10.2f\r\n %10.2f\r\n', BaseTFW);
fclose(fileID);


function  [x,y,utmzone] = deg2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
% Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.
%
% Inputs:
%    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
%    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
%
% Outputs:
%    x, y , utmzone.   See example
%
% Example 1:
%    Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
%    Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
%    [x,y,utmzone] = deg2utm(Lat,Lon);
%    fprintf('%7.0f ',x)
%       458731  407653  239027  230253  343898  362850
%    fprintf('%7.0f ',y)
%      4462881 5126290 4163083 3171843 4302285 2772478
%    utmzone =
%       30 T
%       32 T
%       11 S
%       28 R
%       15 S
%       51 R
%
% Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
%    LatDMS=[40 18 55.56; 46 17 2.04];
%    LonDMS=[-3 29  8.58;  7 48 4.44];
%    Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
%    Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
%    [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06, Aug/06
% Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern 
%    hemisphere coordinates. 
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------

% Argument checking
%
error(nargchk(2, 2, nargin));  %2 arguments required
n1=length(Lat);
n2=length(Lon);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end


% Memory pre-allocation
%
x=zeros(n1,1);
y=zeros(n1,1);
utmzone(n1,:)='60 X';

% Main Loop
%
for i=1:n1
   la=Lat(i);
   lo=Lon(i);

   sa = 6378137.000000 ; sb = 6356752.314245;
         
   %e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   %alpha = ( sa - sb ) / sa;             %f
   %ablandamiento = 1 / alpha;   % 1/f

   lat = la * ( pi / 180 );
   lon = lo * ( pi / 180 );

   Huso = fix( ( lo / 6 ) + 31);
   S = ( ( Huso * 6 ) - 183 );
   deltaS = lon - ( S * ( pi / 180 ) );

   if (la<-72), Letra='C';
   elseif (la<-64), Letra='D';
   elseif (la<-56), Letra='E';
   elseif (la<-48), Letra='F';
   elseif (la<-40), Letra='G';
   elseif (la<-32), Letra='H';
   elseif (la<-24), Letra='J';
   elseif (la<-16), Letra='K';
   elseif (la<-8), Letra='L';
   elseif (la<0), Letra='M';
   elseif (la<8), Letra='N';
   elseif (la<16), Letra='P';
   elseif (la<24), Letra='Q';
   elseif (la<32), Letra='R';
   elseif (la<40), Letra='S';
   elseif (la<48), Letra='T';
   elseif (la<56), Letra='U';
   elseif (la<64), Letra='V';
   elseif (la<72), Letra='W';
   else Letra='X';
   end

   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log( ( 1 +  a) / ( 1 - a ) );
   nu = atan( tan(lat) / cos(deltaS) ) - lat;
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   ta = ( e2cuadrada / 2 ) * epsilon ^ 2 * ( cos(lat) ) ^ 2;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000;
   yy = nu * v * ( 1 + ta ) + Bm;

   if (yy<0)
       yy=9999999+yy;
   end

   x(i)=xx;
   y(i)=yy;
   utmzone(i,:)=sprintf('%02d %c',Huso,Letra);
end
end

