
%%Initialize
close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%PATHWAYS%%%%%%%%%%%%%%%
RawFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\';
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';
SupportFold ='C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\SupportParms\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grdimage=imread('SwiftLogo.jpg');
%%%%%%%%%%%%%%%%%%%  USER INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter Comm port # for GGA string at 9600baud' };
dlg_title = 'Swift GPS Input';
defaultans = {'COM1'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
%%%%%%%%%%%%%%%%%%%%%%%%
comPort=string(answer{1});
clear answer
%%%%%%%%%%CREATE INSTRUCTION MESSAGE BOX%%%%%%%%%%%%%%%%
f = msgbox({'Halo6: port 6132 & 236.6.7.100';'4G: port 7507 & 236.6.10.19 ';'JRC&4G SN696: port 6678 & 236.6.7.8 '},'Ethernet guide','custom',grdimage);
set(f, 'position', [50 300 200 100]); %makes box bigger
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT FOR IP RADAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANTENNA%%CONNECTION%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter radar antenna destination port', 'IP address' };
dlg_title = 'Swift Ethernet Input';
defaultans = {'6678','236.6.7.8'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
socketport=str2double(answer{1});%
ip=string(answer{2});
clear answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT RANGE SELECTED IN NAVICO%%%%
%%%%%%%%%%CREATE INSTRUCTION MESSAGE BOX%%%%%%%%%%%%%%%%
f = msgbox({'1/8nm Navico=~400m Range';'1/4nm Navico=~800m Range';'1/2nm Navico=~1600m Range';'3/4nm Navico=~2500m Range'},'Range Settings','custom',grdimage);
set(f, 'position', [50 100 200 100]); %makes box bigger
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT FOR IP RADAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANTENNA%%CONNECTION%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter range of Swift collect'};
dlg_title = 'Swift Range';
defaultans = {'400'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
radarCoverageRange=str2double(answer{1});%
clear answer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%file for GeoTiff overlays%%%%%%%%%%%%%%
%%%%%%%%%%CREATE INSTRUCTION MESSAGE BOX%%%%%%%%%%%%%%%%
f = msgbox({'FortStory_150.tif';'FRF_SouthPier.tif'},'Ethernet guide','custom',grdimage);
set(f, 'position', [50 450 200 100]); %makes box bigger
%%%%%%%%%%%%%%%%5555
prompt = {'Enter name of geotiff file:' };
dlg_title = 'Swift base map';
defaultans = {'Jmap_FStory.tif'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
%%%%%%%%%%%%%%%%%%%%%%%%
parms.tifFile=char(answer{1});
clear answer

%%%%%%%%%%%%% needed to process are in RAW directory
grdimage=imread('SwiftLogo.jpg');
%%%%%%%%%%CREATE INSTRUCTION MESSAGE BOX%%%%%%%%%%%%%%%%
f = msgbox({'Input degrees between alongshore orientation and true north +CCW'},'Guidance','custom',grdimage);
set(f, 'position', [50 200 300 50]); %makes box bigger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  USER INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Rotation angle between shoreline and T North:'};
dlg_title = 'Swift Processing Input';
defaultans = {'40'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
%%%%%%%%%%%%%%%%%%%%%%%%
RotationTN=str2num(answer{1}); %%%1 for JRC or other boom antenna; 2 for 4G dome antenna

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GPS STUFF%%%%%%%%%%%%%%
% comPort = 'COM1';
cnt = 1;
disableCnt = 1;
incMsg_GPS = [];
%% Init Parameters
parms.gpsBuffSz = 30;
parms.gpsYExtent = 500;
parms.fontSize = 12;
parms.textColor = 'k';
parms.cmap = jet(16);
parms.cmap(1,:) = 0;
parms.amap = linspace(0,1,256);
parms.plotUpdateDistance = .25;
% parms.tifFile = 'FRF_SouthPier.tif';

radarCoverageHeading = 0; %degrees
radarCoverageExtent = 180; %degrees

moveDist = 5;  %Draw a new radar patch everytime this distance is moved

[gpsE,gpsN,phi,refX,refY,speed,prevGpsE,prevGpsN,bearingPhi,scanCnt,...
    gpsBuff] = initSwiftNavProc(parms);

initFlag = 1;
[parms,handles] = setDisplayLayout(parms); %Go in here to change location of processing window
handles.hPatch = zeros(100000,1);

patchCnt = 1;
ignoreFirstPatch = 1;
currPos = [0 0];
refreshCnt = 1;
gpsRefreshCnt = 0;

currGPSData.t = ' ';
currGPSData.lat = ' ';
currGPSData.lon = ' ';
currGPSData.Q = ' ';
currGPSData.numSats = ' ';
currGPSData.HDOP = ' ';
currGPSData.utmE = ' ';
currGPSData.utmN = ' ';
currGPSData.utmZone = ' ';
currGPSData.compassHeading = ' ';
currGPSData.track = ' ';
ggaValidReadFlag = 0;
interrupt.QA_QC = 0;
track = 0;

spokeCnt = 1;
updateRadarCompass = tic;
updateRadarImage = tic;
tmpCnt = [];

isgood = 1;
maxInt = 0;

mask1 = uint16(bin2dec('0000111111111111'));
mask2 = uint16(bin2dec('0001111111111111'));
mask3 = uint16(bin2dec('0011111111111111'));

intIdx = 1:1024;
while handles.isFig
    if initFlag
        [gpsE,gpsN,phi,refX,refY,speed,prevGpsE,prevGpsN,bearingPhi,scanCnt,...
            gpsBuff] = initSwiftNavProc(parms);
        
        %% Open Ports
        inst = instrfind;
        delete(inst);
        clear inst;
        clear portObj;
        try
%                         ;portObj_radar = tcpip('169.254.122.56',10001,'timeout',5,'inputBufferSize',50000)
            %%%%%%%%%%%%%%%%%%OPENS DESTINATION PORT%%%%%%%%%%%%%
%                         fopen(portObj_radar);
% %         %    socket = java.net.MulticastSocket(6132);%halo6 antenna
% %        %  socket = java.net.MulticastSocket(7507);%4Gsn139antenna
% %         % socket = java.net.MulticastSocket(6678);%6kw JRC
   socket = java.net.MulticastSocket(socketport);%6kw JRC

% %          %   socket.joinGroup(java.net.InetAddress.getByName('236.6.10.19'));%4Gsn139           
% %        %   socket.joinGroup(java.net.InetAddress.getByName('236.6.7.100'));%%Halo6
   socket.joinGroup(java.net.InetAddress.getByName(ip));%JRC 6kw
            socket.setReuseAddress(1);
            socket.setSoTimeout(1);
            
            packet = java.net.DatagramPacket(zeros(1, 17160, 'int8'), 17160);
            initFlag = 0;
        catch
            warndlg('UDP Port Opening Error','Comms Error')
            break
        end
        gpsConnected = true;
        try
            portObj_GPS = serial(comPort,'BaudRate',9600,'Timeout',5,'inputBufferSize',50000,'byteOrder','bigEndian');
            fopen(portObj_GPS);
        catch
            warndlg('Serial Port Opening Error','Comms Error')
            gpsConnected = false;
        end
    end
    %% GPS Data
    %Always read GPS data
    if gpsConnected
        [ggaData, numScansCurr, incMsg_GPS, ggaValidReadFlag, rawMsg, vtgValidReads,vtgData] = readPosData(portObj_GPS, incMsg_GPS);
        
        
        if length(ggaData) > 0
            [posTime, lat, latDir, lon, lonDir, gpsQuality, numPosReads, numSats, hdop, alt, heightAboveGeoid] = convertGGAData(ggaData);
            [lat,lon] = convertNMEA2DD(lat,lon,latDir,lonDir);
            [gpsE,gpsN,utmzone] = deg2utm(lat,lon);
            currGPSData.t = num2str(posTime(end),'%10.2f');
            currGPSData.lat = num2str(lat(end),'%3.10f');
            currGPSData.lon = num2str(lon(end),'%3.10f');
            currGPSData.Q = num2str(gpsQuality(end));
            currGPSData.numSats = num2str(numSats(end));
            currGPSData.utmE = num2str(gpsE(end),'%10.2f');
            currGPSData.utmN = num2str(gpsN(end),'%10.2f');
            currGPSData.utmZone = utmzone(end,:);
            currGPSData.HDOP = num2str(hdop(end),'%5.2f');
        end
        
        [track, magneticTrack,speed_kts,speed_kph] = convertVTGData(vtgData);
        if length(vtgData) > 0
            currGPSData.track = num2str(track(end),'%5.2f');
        end
        
        
        gpsRefreshCnt = gpsRefreshCnt + length(ggaData);
    end
    
    %try
    %% Radar data
    radarReadCnt = 1;
    hdr = int8([]);
    dataMat = int8([]);
    while 1%radarReadCnt < 20
        try
            socket.receive(packet);
            %len1(cnt) = packet.getLength;
            spokeFrame = packet.getData;
            spokeFrame_noFrameHeader = spokeFrame(9:17160);
            
            spokeMat = reshape(spokeFrame_noFrameHeader,[size(spokeFrame_noFrameHeader,1)/32 32]);
            hdr(radarReadCnt,:,:)  = spokeMat(1:24,:);
            dataMat(radarReadCnt,:,:) = int8(spokeMat(25:end,:));
            disp(radarReadCnt);
            radarReadCnt = radarReadCnt + 1;
            
        catch
            %disp(size(dataMat,1))
            intensity = zeros(size(dataMat,1),1024,32,'uint8');
            spokeLengthBytes = uint16([]);
            spokeSequenceNumber = uint16([]);
            spokeSamples = uint16([]);
            bitsPerSize = uint8([]);
            rangeCellSize_mm = uint16([]);
            spokeAzimuth = uint16([]);
            bearingZeroError = int8([]);
            spokeCompass = uint16([]);
            trueNorth = int8([]);
            compassInvalid = int8([]);
            rangeCellsDiv2 = uint16([]);
            
            for radarCnt = 1:size(dataMat,1)
                for datNum = 1:size(dataMat,3)
                    
                    %intensity(1:2:end,datNum) = bitand(uint8(15),dataMat(:,datNum));
                    intensity(radarCnt,1:2:end,datNum) = typecast(bitand(int8(15),dataMat(radarCnt,:,datNum)),'uint8');
                    intensity(radarCnt,2:2:end,datNum) = bitshift(typecast(dataMat(radarCnt,:,datNum),'uint8'),-4);
                    %         maxInt = max(maxInt,max(max(intensity(2:2:end,datNum))));
                    %         disp(maxInt)
                    spokeLengthBytes(radarCnt,datNum) = typecast(hdr(radarCnt,1:2,datNum),'uint16');
                    spokeSequenceNumber(radarCnt,datNum) = typecast(hdr(radarCnt,3:4,datNum),'uint16');
                    spokeSamples(radarCnt,datNum) = bitand(typecast(hdr(radarCnt,[5 6],datNum),'uint16'),mask1);
                    bitsPerSize(radarCnt,datNum) = bitand(typecast(hdr(radarCnt,6,datNum),'uint8'),uint8(15));
                    rangeCellSize_mm(radarCnt,datNum) = typecast(hdr(radarCnt,7:8,datNum),'uint16');
                    spokeAzimuth(radarCnt,datNum) = bitand(typecast(hdr(radarCnt,[9 10],datNum),'uint16'),mask2);
                    spokeAzimuth(radarCnt,datNum) = min(spokeAzimuth(radarCnt,datNum),4095);
                    bearingZeroError(radarCnt,datNum) = bitget(hdr(radarCnt,10,datNum),7);
                    spokeCompass(radarCnt,datNum) = bitand(typecast(hdr(radarCnt,[11 12],datNum),'uint16'),mask3);
                    trueNorth(radarCnt,datNum) = bitget(hdr(radarCnt,12,datNum),7);
                    compassInvalid(radarCnt,datNum) = bitget(hdr(radarCnt,12,datNum),8);
                    rangeCellsDiv2(radarCnt,datNum) = typecast(hdr(radarCnt,13:14,datNum),'uint16');
                    if handles.logging
                        %Interpolate GPS data to radar data and write
                        %% File Format
                        %4 UTCHH,UTCMM,UTXSS,UTChh (int8)
                        %3 Lat,Lon,Track, (float)
                        %3 BearingError,BitsPerSamp,CompassInv, (int8)
                        %8 Samples,CellsDiv2,  CellSize,Sequence,Azimuth,Compass, SpokeLength,TrueNorth (int 16)
                        %1024 Data (int8)
                        %1 isgood (int8)
                        utc_time = java.lang.System.currentTimeMillis;
                        dv = datevec(datenum(1970,1,1)+(utc_time/1000)/86400);
                        UTCHH = dv(4);UTCMM = dv(5); UTCSS = floor(dv(6));UTChh = rem(dv(6),1);
                        %            writeData = [UTCHH UTCMM UTCSS UTChh lat(end) lon(end) track(end),...
                        %                         bearingZeroError(datNum) bitsPerSize(datNum) compassInvalid(datNum) ,...
                        %                         spokeSamples(datNum) rangeCellsDiv2(datNum) rangeCellSize_mm(datNum) ,...
                        %                         spokeSequenceNumber(datNum) spokeAzimuth(datNum) spokeCompass(datNum) ,...
                        %                         spokeLengthBytes(datNum) trueNorth(datNum) intensity(:,datNum)' isgood];
                        writeSwiftPacket(handles.saveFileId,UTCHH,UTCMM,UTCSS,UTChh,lat(end),lon(end),track(end),...
                            bearingZeroError(radarCnt,datNum), bitsPerSize(radarCnt,datNum), compassInvalid(radarCnt,datNum),...
                            spokeSamples(radarCnt,datNum), rangeCellsDiv2(radarCnt,datNum), rangeCellSize_mm(radarCnt,datNum) ,...
                            spokeSequenceNumber(radarCnt,datNum), spokeAzimuth(radarCnt,datNum), spokeCompass(radarCnt,datNum) ,...
                            spokeLengthBytes(radarCnt,datNum), trueNorth(radarCnt,datNum), intensity(radarCnt,:,datNum)', isgood);
                    end
                end
            end
            break
        end
    end
    
    
    %Counter to start new file
    if handles.logging
        if toc(handles.fileTime) > handles.LoggingTime
            fclose(handles.saveFileId);
            handles.fileNum = handles.fileNum + 1;
            handles.fileName = [get(handles.hPrefix,'string') '_' num2str(handles.fileNum,'%04.0f')];
            set(handles.hFileName,'string',handles.fileName);
            handles.saveFileId = fopen([handles.fileName '.bin'],'w');
            handles.fileTime = tic;
        end
    end
    
    %tmpCnt = [tmpCnt spokeSequenceNumber];
    
    %Set new cdata based on radar data read
    for iii = 1:size(intensity,1)
        for intNum = 1:size(intensity,3)
            %currSpokeNum = mod(spokeCnt,2048) + 1;
            %Not including heading
            %currSpokeNum = round((spokeAzimuth(intNum)+1)/2);
            %Including heading
            currSpokeNum = round((rem(spokeAzimuth(iii,intNum) + spokeCompass(iii,intNum),4095) + 1) / 2);
            %disp(spokeAzimuth(intNum));
            handles.zMat(:,currSpokeNum) = intensity(iii,:,intNum);
            spokeCnt = spokeCnt + 1;
        end
        currGPSData.compassHeading = num2str(double(spokeCompass(iii,end))/4096*360,'%5.2f');
    end
    
    
    %Plot Gps plus radar compass and radar intensity data on screen
    if handles.isFig
      %%%%%%%%%%%comment the if loop lines 249-254 %% to prioritize writing binary file%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %Update Radar Compass data every second
       if toc(updateRadarCompass) > 1
%           set(handles.radarImage,'Cdata',handles.zMat);
            %%%%%%Heading data comes with radar data
           set(handles.hCompass,'string',currGPSData.compassHeading);
           updateRadarCompass = tic;
       end
 %%%%%%%%%%%%%%%%%%%trial method for displaying radar image for short, defined period%%%      
       % %%Update Radar Image for 20s before stopping
       if toc(updateRadarImage) < 20
          set(handles.radarImage,'Cdata',handles.zMat);
            %%%%%%Heading data comes with radar data
%            set(handles.hCompass,'string',currGPSData.compassHeading);
%            updateRadarImage = tic;
       end
       
        %keyboard
        if gpsConnected & gpsRefreshCnt > 0
            
            %% Update GPS Data box
            %Update GPS Text boxes
            set(handles.hGpsTime,'string',currGPSData.t);
            set(handles.hLat,'string',currGPSData.lat);
            set(handles.hLon,'string',currGPSData.lon);
            set(handles.hGpsQuality,'string',currGPSData.Q);
            set(handles.hNumSats,'string',currGPSData.numSats);
            set(handles.hUtmE,'string',currGPSData.utmE);
            set(handles.hUtmN,'string',currGPSData.utmN);
            set(handles.hUtmZone,'string',currGPSData.utmZone);
            set(handles.hHDOP,'string',currGPSData.HDOP);
            set(handles.hTrack,'string',currGPSData.track);
            %Change currPos is gps has moved
            if sqrt( (currPos(1)-gpsE(end)).^2 + (currPos(2)-gpsN(end)).^2) > moveDist;
                currPos = [gpsE(end) gpsN(end)];
                %Plot radar coverage patch if the radar has moved
%                 if handles.logging
                    [handles, patchCnt] = updateCoverageMap(radarCoverageRange,radarCoverageHeading,radarCoverageExtent,currPos(1)-handles.refX,currPos(2)-handles.refY,handles,patchCnt);
%                 end
            end
            
            %Plot radar location
            if min([currPos(1)-handles.refX]) < parms.currXLim(1) || max([currPos(1)-handles.refX]) > parms.currXLim(2) || ...
                    min([currPos(2)-handles.refY]) < parms.currYLim(1) || max([currPos(2)-handles.refY])> parms.currYLim(2)
                halfX = (parms.currXLim(2)-parms.currXLim(1))/2;
                halfY = (parms.currYLim(2)-parms.currYLim(1))/2;
                parms.currXLim = [min(currPos(1))-handles.refX-halfX min(currPos(1))-handles.refX+halfX];
                parms.currYLim = [min(currPos(2))-handles.refY-halfY min(currPos(2))-handles.refY+halfY];
                set(handles.gpsAx,'xlim',parms.currXLim,'ylim',parms.currYLim);
            end
            
            set(handles.hRadarLoc,'xdata',currPos(1)-handles.refX,'ydata',currPos(2)-handles.refY);
            
            
            
            gpsRefreshCnt = 0;
            
            
        end
        
        
    end
    %Increment Scan Count
    scanCnt = scanCnt + 1;
    refreshCnt = refreshCnt + 1;
    if refreshCnt > 9
        drawnow;
        refreshCnt = 0;
    end
    cnt = cnt + 1;
    if interrupt.QA_QC
        [XGrid, YGrid, IGrid, Time_total, r, c, p]=qaqc_script;
        flushinput(portObj_GPS);
        incMsg_GPS = [];
        interrupt.QA_QC = false;
    end
    
    
end
%Cleanup if necessary
socket.close()
if exist('portObj_GPS')
    if strcmp(get(portObj_GPS,'status'),'open')
        fclose(portObj_GPS);
    end
    delete(portObj_GPS);
    clear portObj_GPS;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save variables for%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%%%%%
save([SupportFold 'SupportParms'],'comPort', 'socketport', 'ip', 'radarCoverageRange', 'parms', 'RotationTN');
movefile('*.bin','C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\');


fclose all;
close all;