function [parms,handles] = setDisplayLayout(parms)

%% Radar and Logging Figure
handles.hFig = figure;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(gcf,'javaframe');
jIcon=javax.swing.ImageIcon('WRTfavicon2.jpg');
jframe.setFigureIcon(jIcon);

handles.isFig = 1;

screenSize = get(0,'screensize');
%% Pixel based
set(handles.hFig,'name','SWIFT','units','pixels','position',[round(screenSize(3).*.39) round(screenSize(4).*.05)  round(screenSize(3).*.6)  round(screenSize(4).*.675)],'numbertitle','off',...
    'DeleteFcn','handles.isFig = 0;','render','OpenGL','colormap',parms.cmap,'alphamap',parms.amap,'toolbar','none','menubar','none');

% %% Logging and File Saving
handles.hLogControl = uipanel('position',[0.05    0.85    0.9    0.125],'backgroundcolor',[0.8 0.8 0.8],...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold');
%Logging Button
handles.logging = 0;
handles.hLogging = uicontrol('style','pushbutton','units','normalized','string','Log',...
    'Position', [0.2   0.305   0.1   0.6],'callback','setLogging','backgroundcolor','green',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hLogControl,'value',handles.logging,'enable','on');
handles.hLoggingIndicatorAx = axes('position',[0 0 1 .2],'parent',handles.hLogControl,'visible','off');
handles.hLoggingIndicator = text(0.22,.85,'Logging','parent',handles.hLoggingIndicatorAx,'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','color','r','visible','off');

%File Saving
axes('position',[0    0   1    1],'visible','off','parent',handles.hLogControl)
%Prefix
text(0.48,.75,'Prefix / Time (s)','fontsize',parms.fontSize,'fontname','arial','fontweight','bold');
handles.prefix = 'swift';
handles.hPrefix = uicontrol('style','edit','units','normalized','string',handles.prefix,'callback','setPrefix;',...
    'Position', [0.65    0.55    0.2    0.39],'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hLogControl);

%Logging Time
handles.LoggingTime = 120;
handles.hLoggingTime = uicontrol('style','edit','units','normalized','string',num2str(handles.LoggingTime),...
    'Position', [0.85    0.55    0.1    0.39],'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hLogControl,'callback','setLoggingTime;');

%Filename
handles.fileNum = 1;
text(0.525,0.3,'Filename','fontsize',parms.fontSize,'fontname','arial','fontweight','bold');
handles.fileName = [handles.prefix '_' num2str(handles.fileNum)];
handles.hFileName = uicontrol('style','edit','units','normalized','string',handles.fileName,...
    'Position', [0.65    0.1000    0.300    0.39],'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hLogControl,'enable','off');

% Pull down menu for replay and simulation
menuItems = ...
    {1,'QA/QC','interrupt.QA_QC = 1;'};     

 menuItems = reshape(menuItems, 3, length(menuItems)/3)';

pad = 3.5; % approximated additional padding added around labels for menu items
handles.pullDownhandles = zeros(max([menuItems{:,1}])+1,1);
handles.pullDownhandles(1) = handles.hFig;
menuwidth = pad;

for i=1:size(menuItems,1)
    hm = uimenu(handles.pullDownhandles(menuItems{i,1}), 'Label', menuItems{i,2} );
    handles.pullDownhandles(i+1) = hm;
    
    set(hm,'Callback',menuItems{i,3});
    
    if menuItems{i,1}==1
        menuwidth = menuwidth + length(get(hm,'Label')) + pad;
    end
end

handles.radarAx = axes('parent',handles.hFig,'units','normalized','position',[0.1 0.1 0.8 0.7]);

R = linspace(0,464,1024);  % Range 1xN
theta = linspace(0,360,2048);   % Azimuth 1xM
handles.zMat = zeros(length(R),length(theta),'uint8');
handles.radarImage = polarPcolor(R,theta,handles.zMat,'Ncircles',5,'Nspokes',9);

%% GPS figure
handles.hGpsFig = figure;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(gcf,'javaframe');
jIcon=javax.swing.ImageIcon('WRTfavicon2.jpg');
jframe.setFigureIcon(jIcon);

%Pixel based
set(handles.hGpsFig,'name','GPS Data','units','pixels','position',[round(screenSize(3).*.05) round(screenSize(4).*.8)  round(screenSize(3).*.9)    round(screenSize(4).*.175)],'numbertitle','off',...
    'render','painters','colormap',parms.cmap,'alphamap',parms.amap,'toolbar','none','menubar','none');

handles.gpsDataAx = axes('parent',handles.hGpsFig,'units','normalized','position',[0 0 1 1],'visible','off');

%set(gcf,'toolbar','figure')
%GPS Time
text(0.09,0.15,'Time','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hGpsTime = uicontrol('style','edit','string','','units','normalized','position',[.01 .35 .15 .45],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%Lat
text(0.2273,0.15,'Lat','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hLat = uicontrol('style','edit','string','','units','normalized','position',[0.1829    0.35    0.09    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');
%Lon
text(0.347,0.15,'Lon','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hLon = uicontrol('style','edit','string','','units','normalized','position',[ 0.3    0.35    0.09    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');
%Quality
text(0.4277,0.15,'Q','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hGpsQuality = uicontrol('style','edit','string','','units','normalized','position',[0.4111    0.3500    0.0310    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');
%Num Sats
text(0.4706,0.15,'Sats','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hNumSats = uicontrol('style','edit','string','','units','normalized','position',[0.4567    0.3500    0.0310    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%HDOP
text(0.5329,0.15,'HDOP','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hHDOP = uicontrol('style','edit','string','','units','normalized','position',[0.5022    0.3500    0.06    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%UTME
text(0.625,0.15,'UTM E','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hUtmE = uicontrol('style','edit','string','','units','normalized','position',[0.5753    0.3500    0.1    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');
%UTMN
text(0.74,0.15,'UTM N','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hUtmN = uicontrol('style','edit','string','','units','normalized','position',[0.69    0.3500    0.1    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%UTMZone
text(0.829,0.15,'Zone','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hUtmZone = uicontrol('style','edit','string','','units','normalized','position',[0.805    0.3500    0.0469    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%Compass
text(0.886,0.15,'Compass','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hCompass = uicontrol('style','edit','string','','units','normalized','position',[0.864    0.3500    0.0469    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%Track
text(0.957,0.15,'Track','parent',handles.gpsDataAx,'fontsize',8,'fontname','arial','fontweight','bold','horizontalalignment','center');
handles.hTrack = uicontrol('style','edit','string','','units','normalized','position',[0.932    0.3500    0.0469    0.4500],...
      'fontweight','bold','fontsize',8,'fontname','arial');

%% Coverage Map Figure
handles.hCoverage = figure;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(gcf,'javaframe');
jIcon=javax.swing.ImageIcon('WRTfavicon2.jpg');
jframe.setFigureIcon(jIcon);

%Pixel based
set(handles.hCoverage,'name','Coverage','units','pixels','position',[round(screenSize(3).*.05) round(screenSize(4).*.15)  round(screenSize(3).*.45)    round(screenSize(4).*.45)],'numbertitle','off',...
    'render','painters','colormap',parms.cmap,'alphamap',parms.amap,'toolbar','none','menubar','none');

% place push button in figure to clear patches
clear_button = uicontrol('Style', 'pushbutton', 'String', 'Clear','units','normalized',...
    'Position', [0.45 0.9 0.1 0.05],'FontSize',14,...
    'Callback','clearPatches;'); 
%Coverage
handles.gpsAx = axes;

%Load tif 
[latVecTif,lonVecTif,tifData] = readTifAndTfw(parms.tifFile);
if size(tifData,3) == 4
    tifData = tifData(:,:,1:3);
end

lonStart = lonVecTif(1);lonStop = lonVecTif(end);lonMean = (lonStart+lonStop)/2;
latStart = latVecTif(1);latStop = latVecTif(end);latMean = (latStart+latStop)/2;

lonAvg = lonVecTif - mean(lonVecTif);
latAvg = latVecTif - mean(latVecTif);

handles.hImTif = imagesc(lonAvg,latAvg,tifData);
set(gca,'ydir','normal','fontname','arial','fontweight','bold','fontsize',16);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Position',[0.1010 0.1224 0.7750 0.8150]);

xlabel(['Easting - ' num2str(mean([lonStart lonStop]),'%10.2f') ' (m)']);
ylabel(['Northing - ' num2str(mean([latStart latStop]),'%10.2f') ' (m)']);

handles.refY = mean([latStart latStop]);
handles.refX = mean([lonStart lonStop]);

set(handles.gpsAx,'position',[0.1    0.15    0.8    0.75],'fontsize',parms.fontSize,'fontname','arial','fontweight','bold',...
    'ylim',[-parms.gpsYExtent parms.gpsYExtent],'xcolor',parms.textColor,'ycolor',parms.textColor,'clim',[1 256],'parent',handles.hCoverage);
grid on;box on;
hold on;
axis(handles.gpsAx,'equal');
drawnow;
set(handles.gpsAx,'xlimmode','manual','ylimmode','manual','drawmode','fast');
xTmp = get(handles.gpsAx,'xlim');
xlim([-(xTmp(2)-xTmp(1))/2 (xTmp(2)-xTmp(1))/2])
parms.currYLim = get(handles.gpsAx,'ylim');
parms.currXLim = get(handles.gpsAx,'xlim');

handles.hRadarLoc = line(0,0,'marker','o','color','k','markerfacecolor','k','markersize',10);

%Zoom in and out of gps Axis
handles.hGpsAxZoomIn = uicontrol('style','togglebutton','units','normalized','position',[0.8006    0.1783    0.0400    0.0517],'String','+',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxZoomIn','value',0);
handles.hGpsAxZoomOut = uicontrol('style','togglebutton','units','normalized','position',[0.8006    0.1300    0.0400    0.0517],'String','-',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxZoomOut','value',0);
handles.hGpsAxPanRight = uicontrol('style','togglebutton','units','normalized','position',[0.8398    0.1300    0.0205    0.1000],'String','>',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxPanRight','value',0);
handles.hGpsAxPanLeft = uicontrol('style','togglebutton','units','normalized','position',[0.7813    0.1300    0.0205    0.1000],'String','<',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxPanLeft','value',0);
handles.hGpsAxPanUp = uicontrol('style','togglebutton','units','normalized','position',[0.7813    0.2283    0.0791    0.0400],'String','^',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxPanUp','value',0);
handles.hGpsAxPanDown = uicontrol('style','togglebutton','units','normalized','position',[0.7813    0.0900    0.0781    0.0400],'String','v',...
    'fontsize',parms.fontSize,'fontname','arial','fontweight','bold','parent',handles.hCoverage,'callback','setGpsAxPanDown','value',0);

drawnow;

%Load last file location
figPos = struct;
if exist('figurePosition.mat','file') == 2
   figPos = load('figurePosition');
end
if isfield(figPos,'mainFigPos')
    set(handles.hFig,'position',figPos.mainFigPos)
end  
if isfield(figPos,'gpsFigPos')
    set(handles.hGpsFig,'position',figPos.gpsFigPos)
end
if isfield(figPos,'coverageFigPos')
    set(handles.hCoverage,'position',figPos.coverageFigPos)
end

drawnow;
set(handles.hFig,'SizeChangedFcn','saveFigPos_Main;')
set(handles.hGpsFig,'SizeChangedFcn','saveFigPos_GPS;')
set(handles.hCoverage,'SizeChangedFcn','saveFigPos_Coverage;')

