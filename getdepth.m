function [ h ] = getdepth( speed,period,type )
% %%%%%%%%%%%%%%%%linear dispersion function %%%%%%%%%%%%
if strcmp(type,'lineartheory')
    speed=abs(speed);
    omega=2*pi/period;
    k=2*pi./(speed.*period);
    h=abs(atanh((omega.^2)./(9.8.*k)))./k;
else
    disp(['you wish i solved ' type])
end


end


