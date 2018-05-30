% %%%%%%%%%%%%%%%%linear dispersion function %%%%%%%%%%%%
function[h]=getdepth_wavelength(wavelength,period,type)

if strcmp(type,'lineartheory')
    omega=2*pi/period;
    k=2*pi./(wavelength);%%uses explicit measure of wavelength
    h=abs(atanh((omega.^2)./(9.8.*k)))./k;
else
    disp(['you wish i solved ' type])
end

end
