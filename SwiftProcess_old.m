function [XGrid, YGrid, IGrid, Time_total, r, c, p]=SwiftProcess;


%%%%%%%%%%%%% USE for multiple collections %%%% INSURE just binary files
%%%%%%%%%%%%% needed to process are in RAW directory
grdimage=imread('SwiftLogo.jpg');
%%%%%%%%%%CREATE INSTRUCTION MESSAGE BOX%%%%%%%%%%%%%%%%
f = msgbox({'Select 1 for boom or 2 for 4G dome Antenna Type';'Suggest 1m grid node for 400m (1/8nm Navico)';'Suggest 2m grid node for 800m range (1/4nm Navico)';'Suggest 3m grid node for 1600m range (1/2nm Navico)';'Suggest 5m grid node for 2500m (3/4nm Navico)';'Compass offset rotate +CW & -CCW to align shoreline y-axis'},'Guidance','custom',grdimage);
set(f, 'position', [50 300 300 100]); %makes box bigger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  USER INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Antenna Type:','Grid node size (m):','Tide elevation:','Wet/dry elevation:', 'CompassOffset:' };
dlg_title = 'Swift Processing Input';
defaultans = {'1','1','0','0.5','12'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
%%%%%%%%%%%%%%%%%%%%%%%%
AntennaType=str2num(answer{1}); %%%1 for JRC or other boom antenna; 2 for 4G dome antenna

GridNode=str2num(answer{2});%% m grid node spacing in cross and longshore (good for 4G dome at quarter nm setting or 0.67m/collect over 688m

TideLevel=str2num(answer{3});%%tidal water level at time of collect (+values above 0; - values below)

WetDryZ=str2num(answer{4});

CompassOffset=str2num(answer{5}); % This will need to be calculated when the compass is mounted + rotates CW
%%%%%%%%%%%% Compass offset must include deviation plus degrees difference
%%%%%%%%%%%% between shoreline orientation and 0deg Magnetic
%%%%% For example, a shoreline orientation of 325dM (35d from 360) plus 10d
%%%%% deviation would result in a total CompassOffset = 45d


%%pre-allocate matrices to speed processing MAY NEED TO ALTER %%%%
s=1;
% IGrid=NaN(800,800,1000);%%large size for mobile collects
% XGrid=NaN(800,800,10);
% YGrid=NaN(800,800,10);
IGrid=NaN(500,500,500);
XGrid=NaN(500,500,10);
YGrid=NaN(500,500,10);

RawFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\';
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';

%% Change the filename below
if nargin==0
%     BRIOSfiles=dir('*.bin');
    BRIOSfiles=dir([RawFold '*.bin']);
    
    for p=1:length(BRIOSfiles)
    filename=BRIOSfiles(p,1).name;
    filedate=BRIOSfiles(p,1).date;
    disp(filedate);

    %% NEW BINARY FORMAT%%% CALLS THIS FUNCTION AT END OF SCRIPT FILE%%%%
    [data, shotTime, ENH]=BinaryFMCW2doubles_WRT(filename,filedate,RawFold);

%% Read in/Convert Data
    Azimuth=data(:,15)/4095*360; %Convert 12bit Azimuth to Geographic Degrees
    Compass=data(:,16)/4095*360; %Convert 12bit Azimuth to Geographic Degrees
    scannum=1+[0; cumsum(diff(Azimuth)<0)]; % calculates the scan number for each shot.

% A new shot starts when the Azimuth goes back to 0
    disp('max scan number =')
    disp(max(scannum))
%%Time keeper
    cntr=1:1:length(Azimuth);
    cntr=cntr';
%%Compass collection using KVH is 1Hz; zeros are filled between
    tempgood=find(Compass>0);
    [tcon,Ccon] = consolidator(cntr(tempgood),Compass(tempgood));
    CompassFill=interp1(tcon,Ccon,cntr);
    CompassFill(1:tempgood(1))=Compass(tempgood(1));%%fill the first part with first good value
    CompassFill(tempgood(end):end)=Compass(tempgood(end));%%fill the last part with last good value
% % % %%%debugging
% % % % CompassFill=(CompassFill.*0)+358;
%%create a smoothing filter for compass
    b=ones(1,8000)./8000; %~5s
   CompassSmooth=filtfilt(b,1,CompassFill);
% CompassSmooth=CompassFill;


    Azimuth=Azimuth + CompassSmooth + CompassOffset; %Geographic with Compass Smoothing
%     Azimuth=Azimuth + CompassFill + CompassOffset; %Geographic with Compass Compensation

    Azimuth=90-Azimuth; %Geographic to Math

    disp('collection duration =')
    disp(datestr(shotTime(end)-shotTime(1),13));

    Rind=1-0.5:data(1,11)-0.5; % (-0.5) map the pixel to the center coordinate

    % R=(rangeCellSize_mm * 2 * rangeCellsDiv2) / nOfSample (FROM DOCUMENTATION)
    Rshot=Rind.*(data(1,13).*2.*data(1,12))./data(1,11) / 1000; %/1000 conver from mm to meters
    % Rshot=Rind.*(data(1,13).*data(1,12))./data(1,11) / 1000; %/1000 conver from mm to meters

    disp('max range and range interval')
    disp(max(Rshot))
    disp(Rshot(10)-Rshot(9));
    
    disp('processing Collection number')
    disp(p)

%     figure
%     plot(shotTime,CompassFill,'b.'); hold on
%     plot(shotTime,CompassSmooth,'r');

    
    intensityData=data(:,19:end); %reshape so it's just intensity data
%     [m,n]=size(intensityData);


    %%extract geographic and transform to rectilinear and remove bogus values
    tempLat=floor(data(:,5));
    tempLon=floor(data(:,6));
    %%DEFINE REGIONAL GEOGRAPHIC COORDINATES
    Latitude=mode(tempLat);
    Longitude=mode(tempLon);

    ind=find(abs(tempLat-Latitude)<2&abs(tempLon-Longitude)<2);%%%finds good values excludes zeros and spikes
%     indBad=find(abs(tempLat+tempLon)-abs(Latitude+Longitude)>2);% %%finds bad values excludes zeros and spikes
    Lat=data(ind,5);
%     Lat(indBad)=mean(data(ind,5));
    Lon=data(ind,6);
%     Lon(indBad)=mean(data(indBad,6));
   
    %%%%%%%%%%%If sitting in static position -- mean the gps
    %%%%%%%%%%%signal%%%%%%%%%%%%%
    L=length(Lat);
    LatMean=repmat(nanmean(Lat),[L 1]);
    LonMean=repmat(nanmean(Lon),[L 1]);
   
%     %%smoothing position data
%     LatSmooth=filtfilt(b,1,Lat);
%     LonSmooth=filtfilt(b,1,Lon);
        %%smoothing position data from meaned data
    LatSmooth=filtfilt(b,1,LatMean);
   LonSmooth=filtfilt(b,1,LonMean);
%     LatSmooth=LatMean;
%     LonSmooth=LonMean;
  
    %%%use data with good positions only
    intensityData=intensityData(ind,:);
    Azimuth=Azimuth(ind);
    [m,n]=size(intensityData);
    shotTime=shotTime(ind);
    scannum=scannum(ind,:);

   
    %%%%%%%%%%%%%%%%%%% insert for static tests only %%%%%%%%%
    % easting=(Lat_dspk*0)+1;
    % northing=(Lon_dspk*0)+1;

    [easting,northing,~]=deg2utm(LatSmooth,LonSmooth); %%transforms to UTM rectilinear
%       [easting,northing]=ll2f(Lat,-Lon); %%transforms to local FRF rectilinear
%      [easting,northing]=ll2f(LatSmooth,-LonSmooth); %%transforms Smooth geo to local FRF rectilinear
%      [easting,northing]=ll2MCBCP(LatSmooth,-LonSmooth); %%transforms to local FRF rectilinear
EastingCenter=nanmean(easting);
NorthingCenter=nanmean(northing);

%       figure
%     plot(easting,northing,'+'); hold on
    

    %% make a Range, Azimuth, and Scan Number for each point in intensityData
    % This speeds up the code by avoiding for loops
     R=repmat(Rshot,[m 1]); 
    Az=repmat(Azimuth,[1 n]);
    X=repmat(easting,[1 n]);%%%just fills matrix with Easting position of antenna at each scan
    Y=repmat(northing,[1 n]);%%just fills matrix with Northing position of antenna at each scan

    scannum=repmat(scannum,[1 n]);
    
    %%%create time string that starts at each new scan rotation
    Time=shotTime(2);
    temp=find(diff(scannum(:,100))>0);
    Time(2:length(temp)+1)=shotTime(temp+1);

    %% Convert from polar coordinates to Cartesian coordinates relative to (0,0)

    [x,y]=pol2cart(Az.*pi/180,R);%%%transforms to xy using rotated angle so waves are on X-axis
%%%%%%%%%%%%%%%%%%%%%leave coordinate system relative to antenna and
%%%%%%%%%%%%%%%%%%%%%transform later after wave analysis
%     x=x+X; 
%     y=y+Y;

      %% Grid the data
    GridCellSize=GridNode;%% m grid node spacing in cross and longshore;
    % ** Current increment is set to 3 meters... this will need to be adapted for
    % longer range or mobile collects... the files will get too big.  Will have to 
    % do some temp small files like in the CLARIS code 
    % [xgrid, ygrid]=meshgrid(min(x(:)/1.2):GridCellSize:max(x(:)/1.2),min(y(:)/1.2):GridCellSize:max(y(:)/1.2));
%     [xgrid, ygrid]=meshgrid(0:GridCellSize:max(x(:)/1.0),min(y(:)/1):GridCellSize:max(y(:)/1)); %grid boundaries for FRF north
%     [xgrid, ygrid]=meshgrid(min(x(:)):GridCellSize:max(x(:)/1.0),min(y(:)/1):GridCellSize:max(y(:)/1)); %grid boundaries for entire coverage

%      %%%gridding approach for only narrow cross-shore section%%%232m range
        %%%%%%%%%%%%%%%%%%%%%%%% FROM BEACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-50),min(northing)-100:GridCellSize:max(northing)+100); %grid boundaries for FRF north
%         %%%%%%%%%%%%%%%%%%%%%%%% FROM BEACH at MCBCP%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [xgrid, ygrid]=meshgrid(min(x(:)+125):GridCellSize:max(easting)+100,min(northing)-150:GridCellSize:max(northing)+150); %grid boundaries for FRF north
        
%%%%%%%%%%%%%%%%%%%%%%%% FROM FRF South BEACH TOWER%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [xgrid, ygrid]=meshgrid(min(easting)-20:GridCellSize:max(x(:)-100),min(northing)-100:GridCellSize:max(northing)+100); %grid boundaries for FRF north

%%%%%%%%%%%%%%%%%%%%%%%% FROM FRF 4G Antenna%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xgrid, ygrid]=meshgrid(min(easting)-25:GridCellSize:max(x(:)-180),min(northing)-250:GridCellSize:max(northing)+250); %higher resolution for drifters

%%%%%%%%%%%%%%%%%%%%%%%% FROM FRF JRC Antenna SET AT 200m Range and 0.2m/sample (navico 1/16mi)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%HIGH RESOLUTION FOR SWASH AND WET/DRY line%%%%%%%%%%%%%
% [xgrid, ygrid]=meshgrid(min(easting)-20:GridCellSize:max(x(:)),min(northing)-150:GridCellSize:max(northing)+150); %higher resolution for drifters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%AUTOMATED GRIDDING OPTIONS%%%%%%%%%
if AntennaType==1 %%%JRC or other boom antenna
    if GridNode< 1.2
    % %%%%%%%%%%%%%Range settings for JRC and 4G at 400m setting (1/8miNavico or <) and 1m grid node
    [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-400:GridCellSize:max(northing)+400); %trim for short JRC antenna
    elseif GridNode > 1.2 | GridNode < 2.2
    %%%%%%%%%%%%%grid settings for JRC at 800 setting (1/4miNavico) and 2m grid node
    [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-600:GridCellSize:max(northing)+600); %trim for short JRC antenna
    elseif GridNode > 2.2 | GridNode < 3.2
    %%%%%%%%%%%%%grid settings for JRC at 1.6km setting (1/2miNavico) and 3m grid node
    [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-1400:GridCellSize:max(northing)+1400); %trim for short JRC antenna
    elseif GridNode > 3.2
    %%%%%%%%%%%%%grid settings for JRC at 2.5km setting (3/4miNavico) and 5m grid node
    [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-2300:GridCellSize:max(northing)+2300); %trim for short JRC antenna
    else
    end
elseif AntennaType==2%%%4G or other dome antenna
    if GridNode< 1.2
    % %%%%%%%%%%%%%Range settings for JRC and 4G at <800m setting (1/8miNavico or <) and 1m grid node
%     [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-150),min(northing)-200:GridCellSize:max(northing)+200); %trim for short JRC antenna
    [xgrid, ygrid]=meshgrid(-50:GridCellSize:max(x(:)-100),-200:GridCellSize:200); %trim for short JRC antenna

    elseif GridNode > 1.2 | GridNode < 2.2
    %%%%%%%%%%%%%grid settings for JRC at 800 setting (1/4miNavico) and 2m grid node
    [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-350:GridCellSize:max(northing)+350); %trim for short JRC antenna
%     elseif GridNode > 2.2 | GridNode < 3.2
%     %%%%%%%%%%%%%grid settings for JRC at 1.6km setting (1/2miNavico) and 3m grid node
%     [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-300:GridCellSize:max(northing)+300); %trim for short JRC antenna
%     elseif GridNode > 3.2
%     %%%%%%%%%%%%%grid settings for JRC at 2.5km setting (3/4miNavico) and 5m grid node
%     [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-300:GridCellSize:max(northing)+300); %trim for short JRC antenna
    else
    end
else
end


%%%%%%%%%%%%%%%%%%%USING HALO6 & JRC antennas in 1/2nm setting (actual 1.66km range)
% [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)),min(y(:)):GridCellSize:max(y(:))); %higher alongshore resolution for Halo6
%%%%%%%%%%%%%grid settings for JRC at 1.66km setting and 3m grid node
% [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-200),min(northing)-1000:GridCellSize:max(northing)+1000); %trim for short JRC antenna
% %%%%%%%%%%%%%grid settings for JRC at 0.8km setting (1/4miNavico) and 2m grid node
% [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-600:GridCellSize:max(northing)+600); %trim for short JRC antenna
%%%%%%%%%%%%%grid settings for JRC at 1.6km setting (1/2miNavico) and 2m grid node
% [xgrid, ygrid]=meshgrid(min(easting)-50:GridCellSize:max(x(:)-100),min(northing)-600:GridCellSize:max(northing)+600); %trim for short JRC antenna


    
    %     %      %%%gridding approach for only narrow cross-shore section%%%232m range
%         %%%%%%%%%%%%%%%%%%%%%%%% FROM OFFSHORE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [xgrid, ygrid]=meshgrid(min(x(:)+50):GridCellSize:max(easting)+20,min(northing)-100:GridCellSize:max(northing)+100); %grid boundaries for FRF north

         %%%gridding approach for only SWASH cross-shore section%%%232m range
%     [xgrid, ygrid]=meshgrid(min(easting)-20:GridCellSize:max(easting)+150,min(northing)-50:GridCellSize:max(northing)+50); %grid boundaries for FRF north

    %preallocate grid
    igrid=nan([size(xgrid) max(scannum(:))]);
    % Grid each scannum into igrid (THIS WILL INTERPOLATE OVER NANS!)
    fprintf('Gridding Data...using GridRoundAvg...so 1st and last scan could look funny\n');

    doRound=1; %experimental round gridding algorithm

        for i=2:1:max(scannum(:)-1);%%skip the first and last rotation
        tic
        fprintf('Processing\t%i/%i\t',i,max(scannum(:)));
        ind=scannum==i;
        if doRound==1
        igrid(:,:,i)=roundgridavg(x(ind),y(ind),intensityData(ind),xgrid,ygrid);
        else
        F=TriScatteredInterp(x(ind),y(ind),intensityData(ind));
        igrid(:,:,i)=F(xgrid,ygrid);
        end
        fprintf('%.2f sec\n',toc);
        end
    imean=nanmean(igrid,3);%calculates the mean intensity

%     close all
    % Code to show the movie
        figure
    [~,~,P]=size(igrid);
        for i=1:P
        pcolor(xgrid,ygrid,igrid(:,:,i));shading flat
        title(num2str(i));
        colorbar
        pause(0.25);
        end
clf
figure
pcolor(xgrid,ygrid,imean); shading flat
    
    [r, c]=size(xgrid);
    e=length(Time)+s;
    IGrid(1:r,1:c,s:s+size(igrid,3)-1)=igrid;
    Time_total(s:e-1)=Time;
    XGrid(1:r,1:c,p)=xgrid;
    YGrid(1:r,1:c,p)=ygrid;
    s=length(Time_total)+1;
    end
    pause(5)
% stop
% close all
end
IGrid=IGrid(1:r,1:c,1:length(Time_total));
XGrid=XGrid(1:r,1:c,1:p);
YGrid=YGrid(1:r,1:c,1:p);
% save RutHut_31Jan2018 IGrid XGrid YGrid Time_total r c p
save([ProcFold 'ProcessedFile'],'IGrid', 'XGrid', 'YGrid', 'Time_total', 'r', 'c', 'p', 'GridNode', 'TideLevel','WetDryZ', 'EastingCenter', 'NorthingCenter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% CHOICE TO EXAMINE SURFACE CURRENTS%%%%%%%%%%
% %%%%%%%%%create gui to continue processing or stop%%%%%
choice=questdlg('Would you like examine Swift data for currents?',...
                'Yes', 'No');
switch choice
    case 'Yes'
        disp([choice '....As you wish'])
        products=1;
        ExamineSurfaceCurrents_Swift;%%%execute product processing script file
    case 'No'
        disp([choice '....Continue with bathy processing'])
        products=2;
        clear
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%create gui to continue processing or stop%%%%%
choice=questdlg('Would you like to create Swift products now?',...
                'Yes', 'No');
switch choice
    case 'Yes'
        disp([choice '....As you wish'])
        products=1;
        XtraSwiftProductProcessing_Mar2018;%%%execute product processing script file
    case 'No'
        disp([choice '....Go have a beer'])
        products=2;
        clear
end

end




function [data, shotTime, ENH]=BinaryFMCW2doubles_WRT(filename,filedate,RawFold)
% %% File Format %%%old format from Brian/GIA
% %4 UTCHH,UTCMM,UTXSS,UTChh
% %4 Lat,Lon,Track,BearingError
% %4 BitsPerSamp,CompassInv,Samples,CellsDiv2
% %4 CellSize,Sequence,Azimuth,Compass
% %3 SpokeLength,TrueNorth,Data

%%%%%%%binary format saved using WRT Swift collection%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Format
%4 UTCHH,UTCMM,UTXSS,UTChh (int8)
%3 Lat,Lon,Track, (float)
%3 BearingError,BitsPerSamp,CompassInv, (int8)
%8 Samples,CellsDiv2,  CellSize,Sequence,Azimuth,Compass, SpokeLength,TrueNorth (int 16)
%1024 Data (int8)
%1 isgood (int8)

%numElements = 1043;

% % %     fwrite(fid,data(n,1:4),'int8');
% % %     fwrite(fid,data(n,5:7),'float');
% % %     fwrite(fid,data(n,8:10),'int8');
% % %     fwrite(fid,data(n,11:18),'int16');
% % %     fwrite(fid,data(n,19:(18+1024)),'int8');
% % %     fwrite(fid,data(n,1043),'int8');
% %      fwrite(fid,[UTCHH UTCMM UTCSS UTChh],'int8');
% %     fwrite(fid,[lat lon track],'float');
% %     fwrite(fid,[bearingZeroError bitPerSize compassInvalid],'int8');
% %     fwrite(fid,[spokeSamples rangeCellsDiv2 rangeCellSize_mm spokeSequenceNumber spokeAzimuth spokeCompass ...
% %                 spokeLengthBytes trueNorth],'int16');
% %     fwrite(fid,intensity,'int8');
% %     fwrite(fid,isgood,'int8');
%%%readmultBin(fid,nread,str,byteperline,startbyte)
% fid=fopen(filename);
fid=fopen(fullfile(RawFold,filename));

numlinebytes=1060;
% numlinebytes=1043;%WRT

%Convert HMS to time... Year Month Day is hardcoded
Year=year(datetime(filedate));
Mo=month(datetime(filedate));
Day=day(datetime(filedate));
HMShh=readmultBin(fid,4,'*int8',numlinebytes,0);
shotTime=datenum(Year,Mo,Day,HMShh(1,:),HMShh(2,:),HMShh(3,:)+HMShh(4,:)./100);

%Read Easting Northing Heading
ENH=readmultBin(fid,3,'*float',numlinebytes,4);

%
dat1=readmultBin(fid,3,'*int8',numlinebytes,16);
% dat1=readmultBin(fid,3,'*int8',numlinebytes,7);%%WRT 
% BitZeroError=dat1(1,:);
% BitsPerSample=dat1(2,:);
% CompassInvalid=dat1(3,:);

% 
dat2=readmultBin(fid,8,'*int16',numlinebytes,19);
% dat2=readmultBin(fid,8,'*int16',numlinebytes,10);%%WRT
% NumSamples=dat2(1,:);
% RangeCellDiv2=dat2(2,:);
% RangeCellSizemm=dat2(3,:);
% SequenceNum=dat2(4,:);
% SpokeAz=dat2(5,:);
% SpokeCompass=dat2(6,:);
% SpokeLengthBytes=dat2(7,:);
% TrueNoth=dat2(8,:);

idata=readmultBin(fid,1024,'*int8',numlinebytes,35);
% idata=readmultBin(fid,1024,'*int8',numlinebytes,18);%%WRT
% idata=readmultBin(fid,1024,'*int8',numlinebytes,27);%%WRT

isgood=readmultBin(fid,1,'*int8',numlinebytes,1059);
% isgood=readmultBin(fid,1,'*int8',numlinebytes,1042);

fclose(fid);

data=[HMShh; ENH; dat1; dat2;idata]';

end

function data=readmultBin(fid,nread,str,byteperline,startbyte)
fseek(fid,startbyte,'bof');
if strcmp('*int8',str)
   nbytes=1;
elseif strcmp('*float',str)
   nbytes=4; 
elseif strcmp('*int16',str)
   nbytes=2;
elseif strcmp('*int64',str)
   nbytes=8;
end
data=fread(fid,inf,[num2str(nread) str],byteperline-nread*nbytes);

data=reshape(data,nread,numel(data)/nread);

end

