
%%%%%%%%%%%%%%%%%%%%%%%%%PATHWAYS%%%%%%%%%%%%%%%
RawFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\';
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';
SupportFold ='C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\SupportParms\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([ProcFold 'CurrentProducts']);
        load([SupportFold 'SupportParms']);
        load([ProcFold 'BathyProducts']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screenSize = get(0,'screensize');
% parms.tifFile = 'FRF_SouthPier.tif';%%already loaded%%%
% parms.tifFile = 'FortStory_150.tif';

parms.cmap = jet(16);
parms.cmap(1,:) = 0;
parms.amap = linspace(0,1,256);
parms.fontSize = 12;
parms.gpsYExtent = 500;
parms.textColor = 'k';

handles.hCoverage = figure;

%Pixel based
set(handles.hCoverage,'name','Coverage','units','pixels','position',[round(screenSize(3).*.05) round(screenSize(4).*.15)  round(screenSize(3).*.45)    round(screenSize(4).*.45)],'numbertitle','off',...
    'render','painters','colormap',parms.cmap,'alphamap',parms.amap,'toolbar','none','menubar','none');

%Coverage
handles.gpsAx = axes;
%Load tif 
[latVecTif,lonVecTif,tifData] = readTifAndTfw(parms.tifFile);
if size(tifData,3) == 4
    tifData = tifData(:,:,1:3);
end

lonStart = lonVecTif(1);lonStop = lonVecTif(end);lonMean = (lonStart+lonStop)/2;
latStart = latVecTif(1);latStop = latVecTif(end);latMean = (latStart+latStop)/2;

% lonAvg = lonVecTif - mean(lonVecTif);%%%reduces xy to just distance
                                        % cross-shore and longshore
% latAvg = latVecTif - mean(latVecTif);
lonAvg = lonVecTif;% - mean(lonVecTif);
latAvg = latVecTif;% - mean(latVecTif);

%%%%%%%%%ROTATE GEOTIFF BASE IMAGE BY COMPASS OFFSET
%+is CCW
% RotateBaseTiff=imrotate(tifData,-16);%%Pass Compass Offset in future
% RotateBaseTiff=imrotate(tifData,0);%%Pass Compass Offset in future

handles.hImTif = imagesc(lonAvg,latAvg,tifData);
% handles.hImTif = imagesc(lonAvg,latAvg,RotateBaseTiff);
set(gca,'ydir','normal','fontname','arial','fontweight','bold','fontsize',16);
set(gca,'XMinorTick','on','YMinorTick','on');
% set(gca,'Position',[0.1010 0.1224 0.7750 0.8150]);

xlabel(['Easting - ' num2str(mean([lonStart lonStop]),'%10.2f') ' (m)']);
ylabel(['Northing - ' num2str(mean([latStart latStop]),'%10.2f') ' (m)']);

%%%%%%%%%%%%%%%%%%%%% NOW DRAPE RESULTS OVER%%%%%%%%%%
hold on

% % load([ProcFold 'CurrentProducts']);
pcolor(rectX,rectY,RateMean);shading flat
    title('Surface current rate (m/s)')
%     title('Mean bedform migration rate (m/hr)')
%     caxis([0 50])%m/day
%     caxis([0 2])%m/s
    set(gca, 'Ydir', 'normal')
    axis equal;axis tight;
    colormap(bone)
%     cmocean('thermal')
    colorbar
    hold on
%     quiver(XGrid(1:5:end,1:5:end),YGrid(1:5:end,1:5:end),VxFiltMean(1:5:end,1:5:end),VyFiltMean(1:5:end,1:5:end),10,'g');
    quiver(rectX(1:10:end,1:10:end),rectY(1:10:end,1:10:end),VxFiltMean(1:10:end,1:10:end),VyFiltMean(1:10:end,1:10:end),5,'g');
    set(gca, 'Ydir', 'normal')
    axis equal;axis tight;
%     ylabel(datestr(time(1)))
%     xlabel(datestr(time(end)))
 

plot(rectX(alg_y,Swash_start:SwashStop),rectY(alg_y,Swash_start:SwashStop),'r','LineWidth',3); hold on
plot(rectX(alg_y,Xshore_start:Xshore_stop),rectY(alg_y,Xshore_start:Xshore_stop),'m','LineWidth',3); hold off

%%%%%%%try writing as geotiff

saveas(gcf,[ProcFold 'Currents_Overlay.tif']);




function [latVecTif,lonVecTif,tifData] = readTifAndTfw(tifFileName)
tifData = 0;
lonVecTif = 0;
latVecTif = 0;
try
    tifData = imread(tifFileName);
catch
    msgbox('No TIF File Found','TIF Plot Error','error');
end
tfwFileName = tifFileName;
tfwFileName(end-2:end) = 'tfw';
try
    tfwData = load(tfwFileName);
    lonVecTif = tfwData(5):tfwData(1):tfwData(5)+(size(tifData,2)-1)*tfwData(1).';
    latVecTif = (tfwData(6):tfwData(4):tfwData(6)+(size(tifData,1)-1)*tfwData(4));
catch
    msgbox('No TFW File Found','TIF Plot Error','error');
end
end