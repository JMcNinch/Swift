% [Q] = probbreaking (Hrms,Hm)
function [Q] = probbreaking (Hrms,Hm)
Q=.5*ones(size(Hm));
ind_fullbreak = (Hrms>=Hm);
Q(ind_fullbreak)=1;
ind_nobreak = (Hrms./Hm<.2);
Q(ind_nobreak)=0;
ind = ~ind_fullbreak&~ind_nobreak;
err = zeros(size(Hm));err(ind) = 1;
count = 0;
while max(max(err))>10^-10;
  count = count+1;
  if count >10000;disp(Hrms);error('Error: Stuck in the probbreaking call');end
  err(ind)=((Q(ind)-1)./log(Q(ind)) - Hrms(ind).^2./Hm(ind).^2);
  dQ =  ((Q-1)./log(Q) - Hrms.^2./Hm.^2)./(1./log(Q)-(Q-1)./(log(Q).^2.*Q));
  Q(ind) = Q(ind)-.1*dQ(ind);
end

%disp(['local Q =',num2str(Q),' with Hrms/Hm = ',num2str(Hrms/Hm)]) 
