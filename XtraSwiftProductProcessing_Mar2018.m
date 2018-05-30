%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Collection parameters%%%%%%%%%%%%
% GridNode=5;%%SET IN RAW Processed script and passed to saved processed file grid node spacing in cross and longshore
% TideLevel=0;%%PASSED FROM USER INPUT AT RAW PROCESSING
L=length(Time_total);
SamplePeriod=mean(diff(Time_total)).*86400;%%mean of time steps
seconds=diff(Time_total).*86400;%%explicit measure of each rotation

% %%%%%%%%%%%USER INPUT TIME STEPS %%%%%%%%%%%%%%%%%%%%%%%
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
% % % %%%%%%%%%%%%%%%%%%insure all blanks removed with movie%%%%
% % % figure
% % % for i=1:L;
% % %     pcolor(XGrid,YGrid,IGrid(:,:,i)); shading flat
% % %     title(datestr(Time_total(i)));
% % %     pause(0.05);
% % % end


%%%%%%%%%%% verify data and location of wave parameter sampling
disp('SELECT LOCATION FOR WAVE ANALYSIS');
figure
IgridMean=nanmean(IGrid,3);
pcolor(IgridMean); shading flat
hold on
[xpos, apos]=ginput(1);
xshore=floor(xpos);%%%%just the matrix position of xshore
ashore=floor(apos);
plot(xshore,ashore,'m*'); 
xlabel('matrix xshore indice')
ylabel('matrix ashore indice')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Define alongshore position for xshore profile
disp('INPUT 4 CLICKS (Dune Toe; Wet/Dry Line; Seaward swash; Seaward extent) FOR XShore TRANSECT -- always START AT Dune Toe then seaward');
[trans_x,trans_y]=ginput(4);

alg_y=floor(trans_y(1));
algY_grid=YGrid(alg_y,1);
Beach_start=floor(trans_x(1));%%%dune toe
Swash_start=floor(trans_x(2));%%%upper swash extent
SwashStop=floor(trans_x(3));%%%upper swash extent

Xshore_stop=floor(trans_x(4));%%matrix indice
Xshore_stopGrid=XGrid(1,Xshore_stop);
plot(Beach_start:Xshore_stop,alg_y,Beach_start:Xshore_stop,alg_y,'g.'); hold on

%%%%FIND LOWER SWASH EDGE %%%% Technique IF user can't make the call%%%%
% % % Intensity_tran=squeeze(IGrid(alg_y,Swash_start:Xshore_stop,:));
% % % tmpmean=nanmean(Intensity_tran,2);
% % % tmp=max(tmpmean);
% % % % tmp=nanstd(Intensity_tran,0,2);
% % % % tmp1=nanmean(tmp);%%used std approach
% % % % tmpr=find(tmp>tmp1); %%find where std reaches mean for constant wave action
% % % tmpr=find(tmpmean==tmp);%%uses position of peak swash intensity for mean waterlevel
% % % SwashStop=tmpr(1)+Swash_start;
% % % clear tmp tmp1 tmpr
%%%%%sstart wave processing for depths etc at lower swash%%%%
Xshore_start=floor((SwashStop+Swash_start)./2);%%%average across swash
% Xshore_start=Swash_start;%%%determine from user input
Xshore_startGrid=XGrid(1,Xshore_start);

plot(Swash_start:SwashStop,alg_y,Swash_start:SwashStop,alg_y,'r.');
plot(Xshore_start:Xshore_stop,alg_y,Xshore_start:Xshore_stop,alg_y,'k.');
hold off
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%look at gradient in xshore direction to find
%%%%%%%%%%%% the leading edge of the waves to measure wavelength %%%%%%%%%%%%%%%%%%%
clear tmp
XGrad=zeros(M,N,L);
% figure
for i=1:L
    [tmp, ~]=gradient(IGrid(:,:,i));
    tmp(tmp<3)=0;
    tmp(tmp>3)=10;

        %%%%viewing for QAQC%%%%
%         pcolor(XGrid,YGrid,tmp); shading flat
%         title(i)
%         pause(0.2)
    XGrad(:,:,i)=tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEASURE WAVE SPEEDS WITH MOTION CODE%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%measure wave speeds with optical video approach%%
%%%%%%%%%%%%%%%%%%%%works best with full wave imaged on radar and no
%%%%%%%%%%%%%%%%%%%%normalization of intensity%%%%%%%%%%%%%%%%%
disp('Measuring speed of wave from radar image time series')
opticFlow = opticalFlowFarneback;
tic
 for i=1:L;   
        fprintf(1,'     ----- %1.0d of %1.0d done...\n',i,L)
        TS(i,1)=Time_total(i);%%tracks time of wave spikes
        TS(i,2)=IGrid(ashore,xshore,i);%%tracks wave spikes
%             tmp=edge(IGrid(:,:,i),'canny');
%             temp=flipud(tmp);
%             framegray = uint8(temp).*100;
        framegray = uint8(IGrid(:,:,i));%.*10;
        flow = estimateFlow(opticFlow,framegray); 
        Vx(:,:,i)=flow.Vx;
        Vy(:,:,i)=flow.Vy;
 end
 toc               
%%%%%%%%%%%%EXAMINE area in radar domain that is decent for wave parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmptime=TS(:,1);
tmpspike=TS(:,2);
xi=tmptime(find(~isnan(tmpspike)));yi=tmpspike(find(~isnan(tmpspike)));
TimeSeries=xi(1):1/(86400*2):xi(end);
ts_spline=interp1(xi,yi,TimeSeries,'linear');

% ts_spline=interp1(TS(:,1),TS(:,2),TS(2,1):1/(86400*2):TS(end-1,1),'linear');%%%create wave time series at uniform time steps
DT=0.5;%period between samples in seconds with ts_spline time series

%%%%%%%%%USE WAVELET Approach to peak period through time %%%%%%%%%
[wave,period,~,~]=wavelet(ts_spline,DT,1,0.05);
amp=abs(wave).^2;%%%extract the power spectra
%%%%find periods between 2 and 22s
tmpS=find(period>3);
periodStart=tmpS(1); clear tmpS
tmpS=find(period>20);
periodStop=tmpS(1);
wave_Tp=period(periodStart:periodStop);%% just look at the portion from ~3-22s

tmp_amp=amp(periodStart:periodStop,:);%%just use region between 3-22s
wave_amp=nansum(tmp_amp,2);%%cumulative sum of each power spectra
clear tmp

% figure
for i=1:length(ts_spline)
%     plot(wave_Tp,tmp_amp(:,i)); hold on
    [~,Tp]=max(tmp_amp(:,i));
    PkTp(i)=wave_Tp(Tp);
%     plot(PkTp(i),max(tmp_amp(:,i)),'r*')
%     pause(0.02)
end

%%%determine running smoothed wave periods
b=ones(1,15)/15; %25 hr filter
PkTp_filt=filtfilt(b,1,PkTp);

%%%extract filtered periods at collect times
Tp_collect=interp1(TimeSeries,PkTp_filt,Time_total,'linear');%%%create wave time series at uniform time steps

Tp_mean=nanmean(PkTp);%%%%mean wave period to use with wavelenth measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%exploring explicit wavelength measures and crest
%%%%%%%%%%%%%%%%movement%%%%%
% % % % %%%%%%%%%%%%%%%%along selected transect%%%%%%%%%%%%%%
waveL_x=squeeze(XGrid(alg_y,Xshore_start:Xshore_stop));
waveL_i=squeeze(XGrad(alg_y,Xshore_start:Xshore_stop,:));

clear tmp tmp2
tmpLimUpper=0;%Tp_mean*1.8;
tmpLimLower=500;%Tp_mean*10;

tmp2=1;
for i=1:L
    [~,LOCS] = findpeaks(waveL_i(:,i),waveL_x);
    tmp=length(LOCS);
    Crest_x(tmp2:tmp2+tmp-2)=LOCS(2:tmp);
    WaveL_dist(tmp2:tmp2+tmp-2)=diff(LOCS);
    WaveL=WaveL_dist(WaveL_dist<tmpLimLower&WaveL_dist>tmpLimUpper);%%remove unrealistic data 20s T moving at 10m/s
    CrestPos=Crest_x(WaveL_dist<tmpLimLower&WaveL_dist>tmpLimUpper);
    tmp2=length(Crest_x)+1;

clear LOCS
end
clear tmp tmp2

[CrestX,I]=sort(CrestPos);
WaveLsort=WaveL(I);
[x,indx]=unique(CrestX);
WaveL_interp=interp1(x,WaveLsort(indx),waveL_x,'linear','extrap');

b=ones(1,20)./20; %5*#steps
WaveL_f=filtfilt(b,1,WaveL_interp);
%%%%%%%%%%QA/QC%%%%%%%%%%%%%%%%%5
figure
plot(CrestX,WaveL,'k.'); hold on
plot(waveL_x,WaveL_f,'r.');
plot(waveL_x,WaveL_interp,'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% EXAMINE WAVELENGTHS ACROSS ENTIRE DOMAIN %%%%%%%%%%%%%
disp('Measuring explicit wavelengths throughout domain');
WLength_dist=NaN(M,N,L);
tic
for a=1:M%%each alongshore step
    waveL_X=squeeze(XGrid(a,:));
    waveL_I=squeeze(XGrad(a,:,:));
    fprintf(1,'     ----- %1.0d of %1.0d done...\n',a,M)

    for i=1:L
    [PKS,LOCS] = findpeaks(waveL_I(:,i),waveL_X);
    [~,ind]=findpeaks(waveL_I(:,i));
        if length(LOCS)>1
            tmp=length(LOCS);
            WLength_dist(a,ind(2:tmp),i)=diff(LOCS);
            else
            %%skip to next row in no peaks found
        end
       
    clear PKS LOCS
    end
end
toc
clear tmp
%%%%%%%%%%%%%%%%%%%%%%%DETERMINE SPEED FROM WAVE LENGTH OVER GRID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WLengthMean=nanmean(WLength_dist,3);
WLengthSpeedMean=WLengthMean./Tp_mean;
for id=1:L
    WLengthSpeed(:,:,id)=WLength_dist(:,:,id)./Tp_collect(id);
end
%%%%%%%%%%%%%%%%%%%step%%
for id=1:L
    [h_wl(:,:,id)]=getdepth(WLengthSpeed(:,:,id),Tp_collect(id),'lineartheory');
end

d=-h_wl;%%
d(d<-20)=NaN;
d=d+TideLevel;%%%adjust depths to tidal elevation

tmpdepthmean=nanmean(d,3);

x=reshape(XGrid,M*N,1);
y=reshape(YGrid,M*N,1);
z=reshape(tmpdepthmean,M*N,1);
%%%%%%%%%%%%%gridding wavelength-based depths%%%%%%%%%%%
[wldepthFZgrid,numpts]=roundgridavg(x(~isnan(z)),y(~isnan(z)),z(~isnan(z)),XGrid,YGrid);
% F=TriScatteredInterp(x(~isnan(z)),y(~isnan(z)),z(~isnan(z)));
% wldepthFZgrid(:,:)=F(XGrid,YGrid);
%%%%%%%%%%%%%%%QA/QC Figures%%%%%%%%%%%%%5
figure
subplot(1,2,1)
pcolor(XGrid,YGrid,tmpdepthmean); shading flat
hold on
contour(XGrid,YGrid,tmpdepthmean);
% caxis([-10 0]);
colorbar
title('Depths from wavelengths')
% subplot(1,2,2)
% pcolor(XGrid,YGrid,wldepthFZgrid); shading flat
% caxis([-10 0]);
% colorbar
% pause(5);
%%%%%%%%%%%APPROACHES FOR GRIDDING%%%%%%%%%%%%%%%%
% % % [ScaledZgrid, numpts]=roundgridavg(ScaledX(r1),ScaledY(r1),ScaledZ(r1),XGrid,YGrid);%%grid soundings to local grid coordinates
% % % F=TriScatteredInterp(ScaledX(r1),ScaledY(r1),ScaledZ(r1));
% % % ScaledZgrid(:,:)=F(XGrid,YGrid);
% % % ScaledZ=reshape(ScaledZgrid,M*N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%% examine results of motion code%%%%%%%%%%%%%%%%%%               
Vx_pixels=double(Vx(:,:,2:end));%%%the first result is whacky-big from comparison
Vy_pixels=double(Vy(:,:,2:end));

VxFilt=Vx_pixels;
VxFilt(Vx_pixels<-9|Vx_pixels>-0.9)=NaN;%%%just keep the negative shifts and assigns zeros

VyFilt=Vy_pixels;
VyFilt(Vy_pixels<-9|Vy_pixels>9)=NaN;%%%just keep the negative shifts

%%%%%%%%%%%determine range magnitude of movement from filtered u v
%%%%%%%%%%%%%%%%%%%%
[angle,range]=cart2pol(VxFilt,VyFilt);
Angle=rad2deg(angle); %%puts in degrees format
MeanWaveAngle=nanmean(nanmean(nanmean(Angle,3)));%%check values and orientation convention

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Examine explicit measures of time and pixel movement with each
%%%%%%%%%%% rotation
%%%solve for each time step and preserve
for id=1:L-1
speed_magStep(:,:,id)=range(:,:,id).*GridNode./seconds(id);
end

%%%%%%%%%%%%%%%%%Wave Speed along transect determined from wave period and
%%%%%%%%%%%%%%%%%explicit measures of wavelength 
WaveSpeed_tran=WaveL_f./Tp_mean;
[DepthWaveLength_tran]=getdepth(WaveSpeed_tran,Tp_mean,'lineartheory');

%%%%%%%%%%%%%%%%%%%%%%QUICK PLOT OF WAVE SPEEDS AND WAVE PERIOD%%%%%
figure
subplot(2,1,1)
pcolor(XGrid,YGrid,nanmean(speed_magStep,3)); shading flat
caxis([0 10])
hold on
colorbar
title('mean of each wave speed')
xlabel('X-shore distance (FRF,m)')
ylabel('A-shore distance (FRF,m)')
plot(XGrid(1,xshore),YGrid(ashore,1),'r*'); hold off
subplot(2,1,2)
plot(PkTp); hold on
plot(PkTp_filt,'r'); hold off
title('running wave period')
ylabel('period (s)')
xlabel('time (s)')
pause(5);

%%%%%%%%%%%%%%%%%%%solve for depth with linear dispersion at each time
%%%%%%%%%%%%%%%%%%%step%%
for id=1:L-1
    [h(:,:,id)]=getdepth(speed_magStep(:,:,id),Tp_collect(id),'lineartheory');
end

depth=-h;%%
depth=depth+TideLevel;%%%adjust depths to tidal elevation
depth(depth<-20)=NaN;%%remove completely erroneous solutions
depthStd=nanstd(depth,1,3);%%determine standard deviation of raw depths

%%%remove high std solutions from depth matrix
depthMean=nanmean(depth,3);
depthMean(depthStd>5.0|depthStd<0.5)=NaN;%%arbitrary threshold for depth solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Grid both depth solutions together%%%
Dcomb=(wldepthFZgrid+depthMean)./2;%%mean the two depth approaches

depthSmooth = filt2d(Dcomb,'avgnans',25,5); % smooth that new matrix *** should be user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%5QA/QC plots%%%%%%%%%%%%%%%%%%%%
%%%%
figure
subplot(1,2,1)
pcolor(XGrid,YGrid,depthStd); shading flat
caxis([0 5])
colorbar
set(gca,'FontSize',12)
title('Std of solutions')
subplot(1,2,2)
pcolor(XGrid,YGrid,depthMean); shading flat
% caxis([-8 0])
hold on
contour(XGrid,YGrid,depthMean,'k')
set(gca,'FontSize',12)
colorbar
title('depths from wave speed')
pause(5)

% close all%%close all QA/QC figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%extract depth profile%%%%
figure
pcolor(XGrid,YGrid,depthSmooth); shading flat
caxis([-8 0])
axis equal
hold on
contour(XGrid,YGrid,depthSmooth,'k')
set(gca,'FontSize',11)
colorbar('EastOutside');%,'Direction','reverse')
xlabel('x-shore distance (FRF,m)')
ylabel('alongshore distance (FRF,m)')
title('Depth Chart')
plot(XGrid(alg_y,Swash_start:SwashStop),YGrid(alg_y,Swash_start:SwashStop),'r','LineWidth',3); hold on
plot(XGrid(alg_y,Xshore_start:Xshore_stop),YGrid(alg_y,Xshore_start:Xshore_stop),'k','LineWidth',3); hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LOOK AT PROFILE RESULTS%%%%%%%%

%%%%Extract results from wave speed solutions %%%%%%
% prof_std=depthStd(alg_y,Xshore_start:Xshore_stop);
prof_x=XGrid(alg_y,Xshore_start:Xshore_stop); %FRF 900
tmp_z=depthSmooth(alg_y,Xshore_start:Xshore_stop); %FRF 900
xi=prof_x(find(~isnan(tmp_z)));yi=tmp_z(find(~isnan(tmp_z)));
prof_z=interp1(xi,yi,prof_x,'linear');

%%%%plot depths of transect using wavelength-based speed solutions%%
%%waveL_x=cross shore distances from Xshore_start to XShore_stop
%%%-DepthWaveLength_tran are depth solutions for each cross shore position

%%%%%%%%%%%% Average the two solutions and filter where they greatly
%%%%%%%%%%%% diverge
tmpz=((-DepthWaveLength_tran+prof_z)./2);
clear tmp xi yi
tmp=abs(-DepthWaveLength_tran-prof_z);
xi=prof_x(1)-1; yi(1)=TideLevel;
xi(2:length(prof_x(tmp<2))+1)=prof_x(tmp<2);yi(2:length(prof_x(tmp<2))+1)=tmpz(tmp<2);
prof_Z=interp1(xi,yi,prof_x,'linear','extrap');
%%%%%%%% using zero-phase filter %%%%%%%%%%
b=ones(1,15)./15; %5*#steps
f=filtfilt(b,1,prof_Z);
f=f(tmp<2);
f_x=prof_x(tmp<2);
%%%%%%%%%%%%%%%%%%Integrate SWASH AND BEACH PROFILE%%%%%%%%%%%%
%%%%%%%%%%%%determine a foreshore and beach profile and Shore SLope%%%%
BeachX=XGrid(alg_y,Beach_start:Swash_start);
BeachZ=(BeachX.*0)+TideLevel+0.5;
SwashX=XGrid(alg_y,Swash_start:SwashStop);
rr=find(prof_x>SwashX(end));
clear tmpz 
tmpx=[SwashX(1),mean(SwashX),SwashX(end)];
tmpz=[WetDryZ,TideLevel,prof_Z(rr(1))];
SwashZ=interp1(tmpx,tmpz,SwashX,'linear');
run=abs(SwashX(1)-SwashX(end));
rise=abs(SwashZ(1)-SwashZ(end));
ShoreSlope_degrees=rad2deg(atan(rise./run));
ShoreSlope_degrees=round(ShoreSlope_degrees,1);
prof_x(1),prof_Z(1)
%%%%%%%%%%%%%%%%%%%Plot profiles
figure
plot(f_x,f,'k','LineWidth',3); hold on
plot(SwashX,SwashZ,'r','LineWidth',3);
plot(BeachX,BeachZ,'g','LineWidth',3);
plot(waveL_x(tmp<2),-DepthWaveLength_tran(tmp<2),'c','LineWidth',0.5);
plot(prof_x(tmp<2),prof_z(tmp<2),'--c','LineWidth',0.5); 
% errorbar(prof_x(1:10:end),f(1:10:end),prof_std(1:10:end),'c');
% xlim([min(XGrid(alg_y,:)) max(XGrid(alg_y,:))])
xlim([min(BeachX)-10,max(f_x)+10]);
ylim([min(prof_Z(tmp<2))-2 WetDryZ+1.5]);
set(gca,'FontSize',11)
ylabel('depth (NAVD88,m)')
xlabel('x-shore distance (FRF,m)')
title('Depth Profile')
legend('Avg. depths','swash','beach','wl method','speed method','Location','SouthOutside');
%%%write foreshore slope on profile figure
strSlope = ['Shore slope (deg) = ',num2str(ShoreSlope_degrees)];
text(SwashX(1),WetDryZ+0.5,strSlope,'HorizontalAlignment','left');
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

%%%%%%%%%%%QA/QC
figure
pcolor(rectX,rectY,depthSmooth); shading flat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\05Dec2017\CRAB_SWIFT_20171205\Crab_prof91');
% plot(Crab_x,Crab_z,'.r');
hold off

%%%%%%%%save key variables for draping over geotiff
save([ProcFold 'BathyProducts'], 'alg_y', 'Swash_start', 'SwashStop', 'Xshore_start', 'Xshore_stop', 'f_x', 'f', 'XGrid', 'YGrid', 'depthSmooth', 'rectX', 'rectY', 'EastingCenter', 'NorthingCenter', 'RotationTN', 'PkTp_filt', 'MeanWaveAngle', 'Tp_mean', 'IgridMean');
clear

% %%%%%%%%%create gui to continue processing or stop%%%%%
choice=questdlg('Would you like to make Geotiff of bathy?',...
                'Yes', 'No');
switch choice
    case 'Yes'
        disp([choice '....As you wish'])
        products=1;
        CreateGeoTiffOverlay_Bathy;%%%execute product processing script file
    case 'No'
        disp([choice '....consider next processing options'])
        products=2;
        clear
end

%%%use estimated depth profile to model wave height%%

% %%%%%%%%%%%%%%%%linear dispersion function uses explicit measure of wavelength%%%%%%%%%%%%
function[h]=getdepth_wavelength(wavelength,period,type)

if strcmp(type,'lineartheory')
    omega=2*pi/period;
    k=2*pi./(wavelength);%%uses explicit measure of wavelength
    h=abs(atanh((omega.^2)./(9.8.*k)))./k;
else
    disp(['you wish i solved ' type])
end

end

% %%%%%%%%%%%%%%%%linear dispersion function %%%%%%%%%%%%
function[h]=getdepth(speed,period,type)

if strcmp(type,'lineartheory')
    speed=abs(speed);
    omega=2*pi/period;
    k=2*pi./(speed.*period);%%assumes shallow-water wave relationship
    h=abs(atanh((omega.^2)./(9.8.*k)))./k;
else
    disp(['you wish i solved ' type])
end

end

  