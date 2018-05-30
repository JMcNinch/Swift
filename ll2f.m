function [X, Y] = ll2f(ALat, ALon)
%  function [X, Y] = f2ll(ALat, ALon)
%  X = FRF cross-shore (m)
%  Y = FRF longshore (m)
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
%
%               ANGLE = 18.14182902
%               LAT   = 36 10' 39.368875"   (10.6561479238')
%               LONG  = 75 44  58.882497"   (44.9813749534')
%
%  
%  Rev 31 Jan 2012 - used Bill Birkemeier's adjustments

%  ************************************************************************

r2d = 180.0 / pi;

Eom=901951.70333;                % E Origin State Plane
Nom=274093.13928;                % N Origin State Plane
ALat0=10.65583950;             % Origin Lat minutes
ALon0=44.9811435;             % Origin Lon minutes
DegLat = 110962.3573;             % m/deg Lat
DegLon = 89955.86413;            % m/deg long
GridAngle=18.1391./r2d;

LatDeg = 36;
LonDeg = 75;

LatDeg=floor(ALat);              % degrees
ALat = (ALat-LatDeg)*60;      % minutes
if (ALon < 0); ALon = -ALon; end;
ALatLeng = (ALat - ALat0) * DegLat/60.0;

%LonDeg=floor(ALon);          % degrees
if (ALon < 75)
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
X = R .* sin(Ang2);
Y = R .* cos(Ang2);
return;
