%%%%examine surface currents%%%
%%%%%%%%%%%CODE FOR NEXT LEVEL OF PROCESSING SWIFT PRODUCTS%%%%%%%%%%%%
%%%%%%%CALLED IN SWIFTPROCESS CODE%%%%%%%%%%%%%%
%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%PATHWAYS%%%%%%%%%%%%%%%
RawFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\';
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';
SupportFold ='C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\SupportParms\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([ProcFold 'ProcessedFile']);
load([SupportFold 'SupportParms']);

R=r;
 % %%%%%%%%%create gui to continue processing or stop%%%%%
% %%%%%%%%%%%%%%%%%%%  USER INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Length of time to average (2-3X wave period):' };
dlg_title = 'Swift Processing Input';
defaultans = {'25'};
N=50; %this will control the width of the inputdlg
answer = inputdlg(prompt,dlg_title,[1 N],defaultans);
%%%%%%%%%%%%%%%%%%%%%%%%
ts=str2num(answer{1}); %%%1 for JRC or other boom antenna; 2 for 4G dome antenna

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=length(Time_total);
SamplePeriod=mean(diff(Time_total)).*86400;%%mean of time steps
seconds=diff(Time_total).*86400;%%explicit measure of each rotation
M=length(XGrid(:,1,1));
N=length(XGrid(1,:,1));
%%%%%for now while STATIC%%%
XGrid=XGrid(:,:,1);
YGrid=YGrid(:,:,1);
NoCollects=p;%%total number of collects 

%%%%%%%%%%%%pre-scan Swift collections for blanks%%%
for i=1:L;
    tmp(i)=nansum(nansum(IGrid(:,:,i)));
end
r=find(tmp>0);%%finds intensity blanks
%%%recreate variables sans blanks and preserve naming
Time_tmp=Time_total(r);
IGrid_tmp=IGrid(:,:,r);
clear L IGrid Time_total
IGrid=IGrid_tmp; 
Time_total=Time_tmp;
seconds=diff(Time_total).*86400;%%explicit measure of each rotation
L=length(Time_total);
clear IGrid_tmp Time_tmp
 
[rows,cols,layers] = size(IGrid);
perN=1;
IGridNormalized=ones(rows,cols,ceil(layers/perN))*nan;
id2=0;

%%%%%%%%%%%%%%%
%%Fill grid with normalized intensities
for i=1:perN:layers
    id2=id2+1;
    tmp=double(IGrid(:,:,i))./double(nanmax(nanmax(IGrid(:,:,i))));%%normalize by dividing by max
%     tmp(voidpts)=nan;tmp(voidland)=nan;
% %     tmp=allmeans(:,:,i);%removes normalizing switch
    tmp2=tmp-nanmean(tmp(:));%remove mean to normalize
    IGridNormalized(:,:,id2) = tmp2;
    end
%%%create time filtered data to remove waves 
% ts=25;
IGN_ts=tsmooth2(IGridNormalized, ts );%%removing higher period incident band

%%%%%%%%%%%%%%%emamine motions of time series with waves removed
clear tmp2
id2=0;
opticFlow = opticalFlowFarneback;
figure
% for i=1:perN:10
for i=1:perN:layers
    id2=id2+1;
    tmp2=IGN_ts(:,:,i)-nanmean(IGridNormalized,3);%%remove mean to highlight only infragravity range
    tmp2(tmp2<0)=0;
    %%%%%%section for plotting to make sure all looks good%%%
    pcolor(XGrid,YGrid,tmp2);shading flat
    title(i);
    caxis([0 .25])
    colorbar
%     caxis([-2*Mstd 2*Mstd])
   xlabel(i);
%     set(gca, 'Ydir', 'normal')
%     axis equal;axis tight;
%     cmocean('thermal')
%     drawnow
    pause(0.2)
        flow = estimateFlow(opticFlow,tmp2); 
        Vx(:,:,id2)=flow.Vx;
        Vy(:,:,id2)=flow.Vy;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% stop
%

%%%%%%%%%%%%%%%%%%%%%%%%%%quantify current motions%%%%%%%%
Vx_pixels=double(Vx(:,:,2:end));%%%the first result is whacky-big from comparison
Vy_pixels=double(Vy(:,:,2:end));

VxFilt=Vx_pixels;%(:,:,4:10);
% VxFilt(Vx_pixels<-9|Vx_pixels>-0.9)=NaN;%%%just keep the negative shifts and assigns zeros

VyFilt=Vy_pixels;%(:,:,4:10);
% vaFilt(Vy_pixels<-9|Vy_pixels>9)=NaN;%%%just keep the negative shifts

%%%%%%%%%%%determine range magnitude of movement from filtered u v
%%%%%%%%%%%%%%%%%%%%
[angle,range]=cart2pol(VxFilt,VyFilt);
angleD=rad2deg(angle); %%puts in degrees format
range=range.*5;%converts pixel range to meters in RIOS coordinate world
Rate=(range./(perN.*1.25));%*ts));%%converts distance range of surface current movement to rate in m/sec
VxFiltMean=nanmean(VxFilt,3);
VyFiltMean=nanmean(VyFilt,3);
RangeMean=nanmean(range,3);
RateMean=nanmean(Rate,3);

figure
%     pcolor(XGrid,YGrid,RangeMean);shading flat
%     pcolor(XGrid,YGrid,diffall(:,:,5));shading flat
    pcolor(XGrid,YGrid,RateMean);shading flat
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
    quiver(XGrid(1:10:end,1:10:end),YGrid(1:10:end,1:10:end),VxFiltMean(1:10:end,1:10:end),VyFiltMean(1:10:end,1:10:end),5,'g');
    set(gca, 'Ydir', 'normal')
    axis equal;axis tight;
%     ylabel(datestr(time(1)))
%     xlabel(datestr(time(end)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Reproject chart to real rectified coordinate system from
%%%%%%%%%%%%local%%%
tmpx=reshape(XGrid,R.*c,1);
tmpy=reshape(YGrid,R.*c,1);
points=[tmpx, tmpy];
points=points';
coordinate = [0; 0];%Coordinate of rotation
angleInDegrees = RotationTN ;%rotation angle in "degree"
%//////////////////////////////////%
% original  = [[piont1; 1], [piont2; 1], [piont3; 1], [piont4; 1]] ;%corrdinate of rotation
original  = ([points; points(1,:).*0+1]) ;%corrdinate of rotation

angleInRadians = degtorad(angleInDegrees);
m = coordinate(1, 1) ;
n = coordinate(2, 1) ;
first = [1 0 -m;0 1 -n;0 0 1];
third = [1 0 m; 0 1 n; 0 0 1];
second = [cos(angleInRadians) -sin(angleInRadians) 0; sin(angleInRadians) cos(angleInRadians) 0; 0 0 1];
rotated_point = third* second* first*original;
ShiftX = rotated_point (1,:) ;
ShiftY = rotated_point (2,:) ;
figure
% plot(ShiftX,ShiftY,'r+'); hold on
% plot(original(1,:),original(2,:),'k*'); hold on
%%%%%%%%%%%%%%%%Transform rotated local back to UTM coordinate system
ShiftX=reshape(ShiftX,R,c);
ShiftY=reshape(ShiftY,R,c);

rectX=ShiftX+EastingCenter; 
rectY=ShiftY+NorthingCenter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%save key variables for draping over geotiff
save([ProcFold 'CurrentProducts'], 'XGrid', 'YGrid', 'RateMean', 'VxFiltMean', 'VyFiltMean', 'rectX', 'rectY');
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%create gui to continue processing or stop%%%%%
choice=questdlg('Would you like to make Geotiff of currents?',...
                'Yes', 'No');
switch choice
    case 'Yes'
        disp([choice '....As you wish'])
        products=1;
        CreateGeoTiffOverlay_Currents;%%%execute product processing script file
    case 'No'
        disp([choice '....consider next processing options'])
        products=2;
        clear
end

