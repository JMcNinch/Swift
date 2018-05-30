%%%script for BathyCompariosn
%%%%%%%%%%%
 clear
 %%%%%%%%%%%%pre-processed file containing gridded Swift data
ProcFold = 'C:\Users\Jesse McNinch\Desktop\projects\StandoffRadar\data\Processed\';
load([ProcFold 'BathyProducts']);
load([ProcFold 'TS17FathometerBathy']);

M=size(depthSmooth,1);
N=size(depthSmooth,2);
%%%%%%%%%%%APPROACHES FOR GRIDDING%%%%%%%%%%%%%%%%
% % % [ScaledZgrid, numpts]=roundgridavg(ScaledX(r1),ScaledY(r1),ScaledZ(r1),XGrid,YGrid);%%grid soundings to local grid coordinates
% % % F=TriScatteredInterp(ScaledX(r1),ScaledY(r1),ScaledZ(r1));
% % % ScaledZgrid(:,:)=F(XGrid,YGrid);
% % % ScaledZ=reshape(ScaledZgrid,M*N,1);

% [TS17zGrid, numpts]=roundgridavg(easting_fath,northing_fath,z_fath,rectX,rectY);%%grid soundings to local grid coordinates
F=TriScatteredInterp(easting_fath,northing_fath,z_fath);
TS17zGrid(:,:)=F(rectX,rectY);
TS17zGrid=reshape(TS17zGrid,M,N);

Bathydif=depthSmooth-TS17zGrid;
MeanDif=nanmean(nanmean(Bathydif));
StdDif=nanstd(nanstd(Bathydif));
%%digitize and make polygon of shared area
% [fathXborder, fathYborder]=ginput;
load('TSfathometerBorder');

figure
subplot(1,3,1)
pcolor(rectX,rectY,TS17zGrid); shading flat
caxis([-6 -1])
colorbar
subplot(1,3,2)
pcolor(rectX,rectY,depthSmooth); shading flat
colorbar
caxis([-6 -1])
hold on
plot(fathXborder,fathYborder,'m','LineWidth',3);
subplot(1,3,3)
pcolor(rectX,rectY,Bathydif); shading flat
colorbar
caxis([-1 1])
