% function [Hrms,setup,Sxx,Sxy,Hm] = bandj(x,d,Hrms1,period,setup1,angle1,gam)
%
function [Hrms,setup,Sxx,Sxy,Hm] = b_and_j(x,d,Hrms1,period,setup1,angle1,gam)
if size(x,1)~=1;x=x';d=d';end  
%if ~exist('angle1')|isempty(angle1); angle1= 0;end
%if angle1==[]|angle1='';angle1= 0;end
if isempty(angle1); angle1= 0;end
omega = 2*pi/period;
a1 = angle1*pi/180;
[k n] = dispersion(2*pi/period,d);c=omega./k;
if ~exist('gam')|isempty(gam)
%if gam==[]|gam=='';
  cdeep = 9.81/omega;ldeep = cdeep*2*pi/omega;ndeep=.5;
  alphadeep = asin(cdeep *sin(a1)/c(1));
  Hrmsdeep = Hrms1*sqrt(cos(a1)/cos(alphadeep)*c(1)*n(1)/(cdeep*ndeep));
  gam = real(0.5+.4*tanh(33*Hrmsdeep/ldeep));
end
dx = x(2)-x(1);
dxcheck = x(2:length(x))-x(1:length(x)-1);
if max(dxcheck) - min(dxcheck)>.001*dx;
  error('Error: x must be uniformly spaced')
end
if min(size(x)==size(d))==0
  error('Error: x and d must have equal size ')
end
if (max(isnan(d))==1|max(isnan(x))==1|max(isnan(Hrms1))==1)
  error('Error: x or d or Hrms contains NaN ')
end
%if min(d)<0;
%  error('Error: d must be positive throughout the domain')
%end

disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp(['Running wave computation across profile, with'])
disp(['Hrms at boundary  = ',num2str(Hrms1),' [m]'])
disp(['Constant period =   = ',num2str(period),' [s]'])
disp(['setup at boundary = ',num2str(setup1),' [m]'])
disp(['angle at boundary = ',num2str(angle1),' [deg]'])
disp(['Gamma = ',num2str(gam)])
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])


Hrms = zeros(size(x));Sxx=Hrms;Sxy=Hrms;h=Hrms;
Hrms(1) = Hrms1;
setup(1) = setup1;
h(1) = d(1)+setup(1);
Hm(1) = 0.88./k(1).*tanh(gam*k(1).*h(1)/.88);
Q(1)=probbreaking(Hrms(1),Hm(1));
Db(1) = .25*9810*omega/(2*pi)*Q(1)*Hm(1)^2;
E(1) = 1/8*9810*Hrms(1)^2;
Ef(1) = E(1)*c(1)*n(1);
a(1) = a1;
Sxx(1) = E(1)*(n(1)*(cos(a(1))^2+1)-.5);
%hhh=waitbar(0,'Please wait...');
for j = 1:length(x)-1
  Db(j+1) = Db(j);
  a(j+1) = a(j);
  if d(j+1)>0;h(j+1) = d(j+1)+setup(j);else h(j+1) = setup(j);end  
  [k(j+1) n(j+1)] = dispersion(omega,h(j+1));
  c(j+1) = omega/k(j+1);
  Sxx(j+1) = Sxx(j);
  
  for i = 1:2
    Ef(j+1) = Ef(j)*cos(a(j))/cos(a(j+1)) - (dx/(2*cos(a(j+1))))*(Db(j)+Db(j+1));
    setup(j+1) = setup(j) - 2/(9810*(h(j+1)+h(j)))*(Sxx(j+1)-Sxx(j));
    if d(j+1)>0;h(j+1) = d(j+1)+setup(j+1);else h(j+1) = setup(j)+d(j+1)-d(j);end
    [k(j+1) n(j+1)] = dispersion(omega,h(j+1));
    c(j+1) = omega/k(j+1);
    E(j+1)= max(Ef(j+1)/(c(j+1)*n(j+1)),0);
    Hrms(j+1) = sqrt(8/9810*E(j+1));
    Hm(j+1) = 0.88./k(j+1).*tanh(gam*k(j+1).*h(j+1)/.88);
    Q(j+1)=probbreaking(Hrms(j+1),Hm(j+1));
    if Hrms(j+1)>=Hm(j+1);
      Hrms(j+1)=Hm(j+1);
      E(j+1) = 1/8*9810*Hrms(j+1)^2;
      Ef(j+1) = E(j+1)*c(j+1)*n(j+1);
    end
    Db(j+1) = .25*9810*omega/(2*pi)*Q(j+1)*Hm(j+1)^2;
    a(j+1) = asin(sin(a(1))*c(j+1)/c(1));
    Sxx(j+1) = E(j+1)*(n(j+1)*(cos(a(j+1))^2+1)-.5);
  end
end

Sxy = 1/8*9810*Hrms.^2.*n.*cos(a).*sin(a);
%close(hhh)
