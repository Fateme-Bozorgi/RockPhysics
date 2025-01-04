

function [Kfluid,Kreuss,rhoeff,vpb,rhob,Kb,rhog,Kg]=fluid(sal,gg,P,T,So,Sb)

% salinity is in ppm divided by 1e6
sal=sal./1e6;


% ideal gas constant
R=8.31441;

% gas density=============================================================

Pr=P./(4.892-0.4048.*gg);
Tr=(T+273.15)./(94.72+170.75.*gg);
E=0.109.*(3.85-Tr).^2.*exp(-(0.45+(8.*(0.56-1./Tr).^2)).*(Pr.^1.2./Tr));
Z=(0.03+0.00527.*(3.5-Tr).^3).*Pr+(0.642.*Tr-0.007.*Tr.^4-0.52)+E;
rhog=(28.8.*gg.*P./(Z.*R.*(T+273.15)))*1000; 

% gas bulk modulus========================================================
gamma=0.85+5.6./(Pr+2)+27.1./(Pr+3.5).^2-8.7.*exp(-0.65.*(Pr+1));
f=E.*1.2.*(-(0.45+8.*(0.56-1./Tr).^2).*Pr.^0.2./Tr)+(0.03+0.00527.*(3.5-Tr).^3);
Kg=(P.*gamma./(1-Pr./Z.*f))*1000000;


% % oil density=============================================================
% rho0=141.5./(og+131.5);
% 
% % dead oil vs. live oil
% if (gor==0)
%   rhoog=rho0;
%   rhop=(rhoog+(0.00277.*P-1.71e-7.*P.^3).*(rhoog-1.15).^2+3.49e-4.*P)*1000;
%   rhoo=rhop./(0.972+3.81e-4.*(T+17.78).^1.175);
% else
%    % B0=0.972+0.00038.*(2.4.*gor.*sqrt(gg./rho0)+T+17.8).^1.175;
%    %in this case B0 is calculated and is: 1.07
%   rhoog=(rho0+0.0012.*gg.*gor)./B0;
%   rhop=(rhoog+(0.00277.*P-1.71e-7.*P.^3).*(rhoog-1.15).^2+3.49e-4.*P)*1000;
%   rhoo=rhop./(0.972+3.81e-4.*(T+17.78).^1.175);
% end
% 
% % oil velocity============================================================
% % live oil use pseudo density
% if (gor~=0)
%   rho0=rho0./B0./(1+0.001.*gor);
% end
% 
% % the following formula is for dead oil and live oil
% vpo=2096.*sqrt(rho0./(2.6-rho0))-3.7.*T+4.64.*P+0.0115.* ...
% (4.12.*sqrt(1.08./rho0-1)-1).*T.*P;
% 
% Ko=vpo.*vpo.*rhoo;

Ko=0;  % WE DONT HAVE OIL IN THIS WELL


% brine density===========================================================
rhow=1+1e-6.*(-80.*T-3.3.*T.^2+0.00175.*T.^3+489.*P-2.*T.*P+ ...
0.016.*T.^2.*P-0.000013.*T.^3.*P-0.333.*P.^2-0.002.*T.*P.^2);
rhob=(rhow+sal.*(0.668+0.44.*sal+1e-6.*(300.*P-2400.*P.*sal+ ...
T.*(80+3.*T-3300.*sal-13.*P+47.*P.*sal))))*1000;

% brine velocity==========================================================
matrixw=[1402.85 4.871 -0.04783 1.487e-4 -2.197e-7
1.524 -0.0111 2.747e-4 -6.503e-7 7.987e-10
3.437e-3 1.739e-4 -2.135e-6 -1.455e-8 5.230e-11
-1.197e-5 -1.628e-6 1.237e-8 1.327e-10 -4.614e-13]';

% water velocity
velw=0;
for i=1:5
for j=1:4
  velw=velw+matrixw(i,j).*T.^(i-1).*P.^(j-1);
end
end

vpb=velw+sal.*(1170-9.6.*T+0.055.*T.^2-8.5e-5.*T.^3+2.6.*P- ...
0.0029.*T.*P-0.0476.*P.^2)+sal.^1.5.*(780-10.*P+0.16.*P.^2)-1820.*sal.^2;
Kb=vpb.*vpb.*rhob;


% fluid mix=============================================================
% Effective density
Sg=1-So-Sb;
% rhoeff = Sb.*rhob+So.*rhoo+Sg.*rhog;  %rho fluid
rhoeff = Sb.*rhob+Sg.*rhog;  %rho fluid

% Reuss Average
if ( Kb.*Ko.*Kg~=0 )
  Kreuss=1./(Sb./Kb+So./Ko+Sg./Kg);
elseif (Sb==0)
  Kreuss=1./(So./Ko+Sg./Kg);
elseif (So==0)
  Kreuss=1./(Sb./Kb+Sg./Kg);
elseif (Sg==0)
  Kreuss=1./(Sb./Kb+So./Ko);
else
  Kreuss=0;
end

% % Voigt Average
% Kvoigt=Sb.*Kb+So.*Ko+Sg.*Kg;
% 
% %Hill average
% Kfluid=(Kreuss+Kvoigt)./2;

 Kfluid=Kreuss; %in this case just Reuss Average is used
end
