function [X, Y] = ll2MCBCP(ALat, ALon)
%  function [X, Y] = f2ll(ALat, ALon)
%  X = MCBCP cross-shore (m)
%  Y = MCBCP longshore (m)
%  ALat = latitude (decimal degrees)
%  ALon = longitude (decimal degrees, positive)
%
%   **********   ll2f.pl  **************
%   14 Feb 2011 - does lat/lon to FRF, vice versa
%      Used angle computed from C and D, and
%         origin of FRF coordinates at 36 10' 39.369" N and
%         75 44' 58.882"
%         Origin position and angle were computed from an average
%         determined by PIERORIENT2.FOR.  Values were found between
%         C and "west rail" (WR), C and "east rail" (ER), D and WR,
%         and D and ER.  These values needed to be recomputed since
%         the previous positions of C and D, given to us by the 
%         Wilmington Dist were incorrect.
%
%         This also gave a new position for the origin and angle:
%   Uses southwest corner of local airstrip at MCBCP
%               ANGLE = 215d math coordinates
%               LAT   = 33 17' 09.13"N   
%               LONG  = 117 27  32.0"W   
%
%  
%  Rev 31 Jan 2012 - used Bill Birkemeier's adjustments

%  ************************************************************************

r2d = 180.0 / pi;

% Eom=901951.70333;                % E Origin State Plane
% Nom=274093.13928;                % N Origin State Plane
ALat0=17.1521666;             % Origin Lat minutes
ALon0=27.5333;             % Origin Lon minutes
DegLat = 110909.5166;             % m/deg at 33.2849d NLat
DegLon = 93151.67698;        % m/deg long
GridAngle=215./r2d;

LatDeg = 33;
LonDeg = 117;

LatDeg=floor(ALat);              % degrees
ALat = (ALat-LatDeg)*60;      % minutes
if (ALon < 0); ALon = -ALon; end
ALatLeng = (ALat - ALat0) * DegLat/60.0;

%LonDeg=floor(ALon);          % degrees
if (ALon < 117)
%	LonDeg=LonDeg-1;
	ALon = (ALon-LonDeg)*60;       % minutes
	ALonLeng = -(ALon - ALon0) * DegLon/60.0;
else
	ALon = (ALon-LonDeg)*60;       % minutes
	ALonLeng = -(ALon - ALon0) * DegLon/60.0;
end

R = sqrt(ALatLeng.^2 + ALonLeng.^2 );
Ang1 = atan2(ALonLeng, ALatLeng);
Ang2 = Ang1 + GridAngle;
X = ((R .* sin(Ang2)).*-1)+200;%+300;%%shifts zero to upper beach at Red Beach
Y = (R .* cos(Ang2)).*-1;%%switches convention so alongshore values increase northward on west coast
return;
