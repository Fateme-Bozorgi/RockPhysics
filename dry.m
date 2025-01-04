
% Lee
function [Kdry,mu_dry]=dry(Km,mu_m,phi,alpha,Gamma,Ek,Pk,Emu,Pmu,effectives)

Kinf=Km.*((1-phi)./(1+(alpha.*phi)));
Kdry=(Kinf)./(1+(Ek.*exp(-effectives./Pk)));

mu_inf=mu_m.*((1-phi)./(1+(alpha*Gamma.*phi)));
mu_dry=(mu_inf)./(1+(Emu.*exp(-effectives./Pmu)));
end