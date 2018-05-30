function [zgrid, numpts]=roundgridavg(x,y,z,xgrid,ygrid,makenans)
%% Needs a constant DX and DY
% calculates the mean of the values closest to each grid node
if nargin==5
   makenans=1; 
end
if sum(isnan(z))>0
   error('NO USING NANS IN Z!!'); 
end

ii=ceil(numel(x)/100000);%find out how many chunks to break the data into for speed performance
%preallocate zgrid and numpts
zgrid=zeros(size(xgrid));
numpts=zeros(size(xgrid));
%for each chunk
%fprintf('gridding...\n');
for i=1:ii
    [zgridii, numptsii]=roundgrid(x(i:ii:end),y(i:ii:end),z(i:ii:end),xgrid,ygrid);
    zgrid=(zgrid.*numpts+zgridii.*numptsii)./(numpts+numptsii);%calculate a running zgrid
    numpts=numpts+numptsii;%calculate a running numpts
    zgrid(numpts==0)=0;%when numpts==0, zgrid/0 = inf in running zgrid, so need to set to 0 again
end

if makenans==1 %change 0s in Zgrid to be nans (Default)
    zgrid(numpts==0)=nan;
end
end

function [zgrid, numpts]=roundgrid(x,y,z,xgrid,ygrid)
x=x(:)';y=y(:)';z=z(:)';

xmax=max(xgrid(:));
ymax=max(ygrid(:));
xmin=min(xgrid(:));
ymin=min(ygrid(:));
%calculate limits of the grid
if diff(xgrid(1,1:2))==0
    dx=diff(xgrid(1:2,1));
    dy=diff(ygrid(1,1:2));
else
    dx=diff(xgrid(1,1:2));
    dy=diff(ygrid(1:2,1));
end

%convert to 2D index space
X=round(x/dx-(xmin/dx-1));
Y=round(y/dy-(ymin/dy-1));
%find points out of grid
[m,n]=size(xgrid);
ind=X<1 | X>n | Y<1 | Y>m;
X(ind)=[];
Y(ind)=[];
z(ind)=[]; 
ind=sub2ind(size(xgrid),Y,X);%convert to 1D index
%sort index and zs by index number
sortedIndZ=sortrows([ind;z]',1);
ind=sortedIndZ(:,1);
z=sortedIndZ(:,2);
%calculate the difference in inde number... diff of 0 means its the same
%index
di=[1; diff(ind)];
%preallocate
zgrid=zeros(size(xgrid));
numpts=zeros(size(xgrid));

while(~isempty(z))%while there are still z values to include
    dinot0=di~=0; %if it's not 0, process it
    idi=ind(dinot0);%index number of those values
    numpts(idi)=numpts(idi)+1;%each index is going to have one more point included in the mean
    
    zgrid(idi)=z(dinot0).*(1./numpts(idi))+zgrid(idi).*((numpts(idi)-1)./numpts(idi));%calculate running mean
    
    z=z(~dinot0);%get rid of points that were just processed
    ind=ind(~dinot0);%get rid of their index too
    
    di=[1; diff(ind)];%calculate the difference again
    
end
end