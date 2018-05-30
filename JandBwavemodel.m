%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INGESTS CROSS-SHORE PROFILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETERS AND EXECUTES JandB WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DISSIPATION MODEL TO OUTPUT WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETERS ACROSS SURF ZONE%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PATHWAYS%%%%%%%%%%%%%%%
RawFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Raw\';
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';
SupportFold ='C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\SupportParms\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ProcFold 'ProcessedFile']);
load([SupportFold 'SupportParms']);
load([ProcFold 'BathyProducts']);

%%%use existing variable name conventions
%%insure there are no gaps in cross-shore profile
tmpxint=f_x(1):GridNode:f_x(end);
tmpyint=interp1(f_x,f,tmpxint,'linear');
x = flip(-tmpxint);%%xshore depth positions
zb = flip(tmpyint);%%%smoothed depth profile

% % % provide offshore water elevation at time of survey (tide and surge))%%
% % %%%e.g. swift collection 5 Dec 21:00UTC or 16:00EST
% % % Year Mo Da Time   Gauge  Tide     Pred    Residual
% % %            EST    ID      m         m        m
% % % 2017 12 05 1600    11   -0.558    -0.736     0.178 
% % %wl = -0.56; pushed from user input prior
wl=TideLevel;
period = nanmean(PkTp_filt); %sec deepwater or use Swift mean
angle1 = MeanWaveAngle;%%work to be done but perhaps mean angle at selected wave spot relative to shoreline

h = wl-zb;

% %%%%independently measured wave heights across surf zone
% load('C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\code\revisedCode_Apr2017\PressureResponseFunction\CrabWaveHeights');

% % %offshore wave
% % % Year Mo Da Time   Gauge  dir        HmO      Period
% % %            EST    ID      N          m          s
% % % 2017 12 05 1602   630  120.000     1.353     6.042 
% % % Hrms1 = 1.35; %m
% % Hrms1 = Hmo_offshore; %pushed from prior user input
% % period = 6.04; %sec %pushed from prior user input
setup1 = 0; %0
angle1 = 0;

% gamma = saturated wave height to depth ratio (typically 0.7)
gam = .7;

%evolve energy and mom across surf
[Hrms,setup,Sxx,Sxy,Hm] = b_and_j(x,h,Hrms1,period,setup1,angle1,gam);
[k,n,c] = dispersion (2*pi/period,h);

% %%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%
%%%%%%%Plot for oceanography and Hazards%%%%
Tp_mean=round(Tp_mean,1);
%%%%%%%
figure
[haxes,hline1,hline2]=plotyy(-x,Hrms,-x,zb); hold on
strTp = ['Wave period = ',num2str(Tp_mean)];
text(mean(-x)+100,max(Hrms)+1.5,strTp,'HorizontalAlignment','left');
title('Wave Height across surf zone')
set(hline1,'LineWidth',3);
set(hline1,'Color','b');
set(hline2,'LineWidth',3);
set(hline2,'Color','k');
set(haxes,{'ycolor'},{'b';'k'})  % Left color red, right color blue...
set(haxes(1),'YLim',[0 max(Hrms)+2])
set(haxes(2),'YLim',[-10 0])
ylabel(haxes(1), 'wave height (m)');
ylabel(haxes(2), 'depths (m)');
xlabel('X-shore distance (FRF,m)'); hold off
%%%%%%%%%%%%%MAP VIEW OF HAZARD PLOT%%%%%%%%%%%%%% 
figure
pcolor(XGrid,YGrid,IgridMean); shading flat
xlabel('X-shore distance (FRF,m)')
ylabel('A-shore distance (FRF,m)')
title('Surfzone Hazards')
hold on
plot(XGrid(alg_y,Swash_start:SwashStop),YGrid(alg_y,Swash_start:SwashStop),'r','LineWidth',3); hold on
plot(XGrid(alg_y,Xshore_start:Xshore_stop),YGrid(alg_y,Xshore_start:Xshore_stop),'g','LineWidth',3); 
colormap(bone)
caxis([0 5])

%%%%%%%%save key variables for draping over geotiff
save([ProcFold 'HazardProducts'], 'XGrid', 'YGrid', 'IgridMean', 'rectX', 'rectY', 'alg_y', 'Swash_start', 'SwashStop', 'Xshore_start', 'Xshore_stop');
clear
% %%%%%%%%%create gui to continue processing or stop%%%%%
choice=questdlg('Would you like to make Geotiff of hazards?',...
                'Yes', 'No');
switch choice
    case 'Yes'
        disp([choice '....As you wish'])
        products=1;
        CreateGeoTiffOverlay_Hazards;%%%execute product processing script file
    case 'No'
        disp([choice '....done, go have a beer'])
        products=2;
        clear
end
