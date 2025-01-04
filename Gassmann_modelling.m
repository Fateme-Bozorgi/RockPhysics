clc
clear all
close all
%% data prepration
VOL_ANH=xlsread("spd10-08_multimin",'J48:J2904');
VOL_CAL=xlsread("spd10-08_multimin",'K48:K2904');
VOL_DOL=xlsread("spd10-08_multimin",'L48:L2904');
VOL_WCS=xlsread("spd10-08_multimin",'M48:M2904');
Vp_log=xlsread("SPD10-08_logs_Rest",'R47:R4052');
Vs_log=xlsread("SPD10-08_logs_Rest",'S47:S4052');
DEPTH_log_p=xlsread("SPD10-08_logs_Rest",'B47:B4052');
DEPTH_log_S=xlsread("SPD10-08_logs_Rest",'B47:B4052');


K_ANH=45 *1e+9; % Gpa to Pa
K_CAL=76 *1e+9; 
K_DOL=95 *1e+9; 
K_WCS=21.8 *1e+9; 

MU_ANH= 29 *1e+9; % Gpa to Pa
MU_CAL=32 *1e+9;
MU_DOL=45 *1e+9;
MU_WCS=6.66 *1e+9;

RHO_ANH= 2950;
RHO_CAL=2710;
RHO_DOL=2870;
RHO_WCS=2600;


 for i=1:length(VOL_ANH)
     K_VOIT(i)= (VOL_ANH(i)*K_ANH + VOL_CAL(i)*K_CAL + VOL_DOL(i)*K_DOL + VOL_WCS(i)*K_WCS)/(VOL_ANH(i)+VOL_CAL(i)+VOL_DOL(i)+VOL_WCS(i));
     K_REUSS(i)= (1/(VOL_ANH(i)/K_ANH + VOL_CAL(i)/K_CAL + VOL_DOL(i)/K_DOL + VOL_WCS(i)/K_WCS))/(VOL_ANH(i)+VOL_CAL(i)+VOL_DOL(i)+VOL_WCS(i));
     K_VRH(i)=(K_VOIT(i)+K_REUSS(i))/2;
 end
  for i=1:length(VOL_ANH)
     MU_VOIT(i)= (VOL_ANH(i)*MU_ANH + VOL_CAL(i)*MU_CAL + VOL_DOL(i)*MU_DOL + VOL_WCS(i)*MU_WCS)/(VOL_ANH(i)+VOL_CAL(i)+VOL_DOL(i)+VOL_WCS(i));
     MU_REUSS(i)= (1/(VOL_ANH(i)/MU_ANH + VOL_CAL(i)/MU_CAL + VOL_DOL(i)/MU_DOL + VOL_WCS(i)/MU_WCS))/(VOL_ANH(i)+VOL_CAL(i)+VOL_DOL(i)+VOL_WCS(i));
     MU_VRH(i)=(MU_VOIT(i)+MU_REUSS(i))/2;
  end
   for i=1:length(VOL_ANH)
     RHO_VOIT(i)= (VOL_ANH(i)*RHO_ANH + VOL_CAL(i)*RHO_CAL + VOL_DOL(i)*RHO_DOL + VOL_WCS(i)*RHO_WCS)/(VOL_ANH(i)+VOL_CAL(i)+VOL_DOL(i)+VOL_WCS(i));
%      RHO_REUSS(i)= 1/(VOL_ANH(i)/RHO_ANH + VOL_CAL(i)/RHO_CAL + VOL_DOL(i)/RHO_DOL + VOL_WCS(i)/RHO_WCS);
%      RHO_VRH(i)=(RHO_VOIT(i)+RHO_REUSS(i))/2;
     RHO_VRH(i)=RHO_VOIT(i); % just VOIT Average is used
 end
 
Km=K_VRH';       % bulk modulus_matrix
mu_m=MU_VRH';    % shear modulus_matrix
rhom=RHO_VRH';   % density_matrix
 
DEPTH=xlsread("spd10-08_multimin",'B48:B2904'); % depth
P=5290*0.00689;      %pressure (Mpa)
T=103;               %temprature (c)
PGD=1;               %Pressure gradient with Depth
effectives=((14.7+ PGD.*(DEPTH./0.3048))*0.00689)-P; %effective pressure

phi=xlsread("spd10-08_multimin",'E48:E2904'); % porosity (effective)

alpha=7.2;         %Lee 2005
Gamma=2.2;         %...
sal=290000;        %Salinity
          
Pk=6.32;           %MacBeth (2004)
SK=0.6;
Ek=SK/(1-SK);
Pmu=9.5; 
SMU=0.6;
Emu=SMU/(1-SMU);

gg=0.68;            %Gas Gravity (Specific Gravity)

Sb=xlsread("spd10-08_multimin",'H48:H2904');%Saturation of water
So=0;               %Saturation of oil

%% calculation

%Bulk modulus of dry rock
[Kdry,mu_dry]=dry(Km,mu_m,phi,alpha,Gamma,Ek,Pk,Emu,Pmu,effectives);

%Bulk modulus of mix fluid
[Kfluid,Kreuss,rhoeff,vpb,rhob,Kb,rhog,Kg]=fluid(sal,gg,P,T,So,Sb);

%saturated bulk modulus
Ksat=Kdry+(((1-(Kdry./Km)).^2)./((phi./Kfluid)+((1-phi)./Km)-(Kdry./(1-Km).^2)));

%saturated density / or bulk density
rho_sat=(1-phi).*rhom + phi.*rhoeff;

% velocity=============================================================

Vp=sqrt((Ksat+ 4/3.*mu_dry)./(rho_sat)); %compressional velocity (m/s)
Vs=sqrt(mu_dry./(rho_sat)); % shear velocity (m/s)


%%
% Equalization the DEPTHS for Vp ========================================
[row,column]=find(Vp_log == -999.25);
Vp_log(row)=[];
DEPTH_log_p(row)=[];

R=[];l=0;
for ii=1:size(DEPTH,1)
    for j=1:size(DEPTH_log_p,1)
     e(j)=abs(DEPTH(ii,1)- DEPTH_log_p(j,1));
    end
     [m(ii),n(ii)]=min(e);  
      l=l+1;
      R(l)=n(ii); %DEPTH_log_Vp index   
      
end
R2=R';
Vp_measured=Vp_log(R2);
DEPTH_measured_p=DEPTH_log_p(R2);

% DEPTHS Equalization for Vs ========================================
[row2,column2]=find(Vs_log == -999.25);
Vs_log(row2)=[];
DEPTH_log_S(row2)=[];

RR=[];l2=0;
for iii=1:size(DEPTH,1)
    for jj=1:size(DEPTH_log_S,1)
     e2(jj)=abs(DEPTH(iii,1)- DEPTH_log_S(jj,1));
    end
     [m2(iii),n2(iii)]=min(e2);  
      l2=l2+1;
      RR(l2)=n(iii); %DEPTH_log_Vs index   
      
end
RR2=RR';
Vs_measured=Vs_log(RR2);
DEPTH_measured_S=DEPTH_log_S(RR2);

%%
% figures: Vp & Vs  ====================================================
figure;
plot(Vs,DEPTH)
set(gca,'YDir','reverse')
hold on
plot(Vs_measured,DEPTH_measured_S)
set(gca,'YDir','reverse')
xlabel('Vs(m/s)')
ylabel('Depth(m)')
legend('Vs-predicted','Vs-measured')
title('Shear Velocity Log')

figure;
plot(Vp,DEPTH)
set(gca,'YDir','reverse')
hold on
plot(Vp_measured,DEPTH_measured_p)
set(gca,'YDir','reverse')
xlabel('Vp(m/s)')
ylabel('Depth(m)')
legend('Vp-predicted','Vp-measured')
title('Compressional Velocity Log')
xlim([3500 7500])

%% density validation

RHOB_log=xlsread("SPD10-08_logs_Rest",'O47:O4052');  %(g/cm3)--->its supposed to be(kg/m3)
DEPTH_log_rho=xlsread("SPD10-08_logs_Rest",'B47:B4052'); %(m)

[row3,column3]=find(RHOB_log == -999.25);
RHOB_log(row3)=[];
DEPTH_log_rho(row3)=[];

R3=[];l3=0;
for ii3=1:size(DEPTH,1)
    for j3=1:size(DEPTH_log_rho,1)
     e3(j)=abs(DEPTH(ii3,1)- DEPTH_log_rho(j3,1));
    end
     [m3(ii3),n3(ii3)]=min(e3);  
      l3=l3+1;
      R3(l3)=n3(ii3); %DEPTH_log_Vp index   
      
end
R33=R3';
RHOB_measured=RHOB_log(R33)*1000;
DEPTH_measured_rho=DEPTH_log_rho(R33);


figure;
plot(rho_sat,DEPTH)
set(gca,'YDir','reverse')
hold on
plot(RHOB_measured,DEPTH_measured_rho)
set(gca,'YDir','reverse')
xlabel('RHO(kg/m3)')
ylabel('Depth(m)')
legend('RHOB-pred','RHOB-measured')
title('RHOB prediction')



%% error (meter per s)
MSE_p=mean((Vp-Vp_measured).^2);
RMSE_p=sqrt(MSE_p)%/length(Vp)
RRMSE_P=(RMSE_p ./(max(Vp_measured)-min(Vp_measured)))*100

MSE_s=mean((Vs-Vs_measured).^2);
RMSE_s=sqrt(MSE_s)%/length(Vp)
RRMSE_s=(RMSE_s ./(max(Vs_measured)-min(Vs_measured)))*100

MSE_rhob=mean((rho_sat-RHOB_measured).^2);
RMSE_rhob=sqrt(MSE_rhob)%/length(Vp)

%% validation_figures
% figure;
% plot(Kfluid,DEPTH)
% xlabel('K-fluid(pa)')
% ylabel('depth(m)')
% set(gca,'YDir','reverse')
% 
% figure;
% plot(rhoeff,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('density-fluid(kg/m3)')
% ylabel('depth(m)')
% 
% figure;
% plot(rho_sat,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('density-bulk(kg/m3)')
% ylabel('depth(m)')
% 
% figure;
% plot(Km,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('K-matrix(Pa)')
% ylabel('depth(m)')
% 
% figure;
% plot(Kdry,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('K-dry(Pa)')
% ylabel('depth(m)')
% 
% figure;
% plot(Ksat,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('K-sat(Pa)')
% ylabel('depth(m)')
% 
% figure;
% plot(mu_dry,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('mu-dry(Pa)')
% ylabel('depth(m)')
% 
% figure;
% plot(mu_dry,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('mu-dry & Kdry(Pa)')
% ylabel('depth(m)')
% hold on
% plot(Kdry,DEPTH)
% set(gca,'YDir','reverse')
% legend('mu-dry','Kdry')
% 
% figure;
% plot(Ksat,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('Ksat & Kdry(Pa)')
% ylabel('depth(m)')
% hold on
% plot(Kdry,DEPTH)
% set(gca,'YDir','reverse')
% legend('Ksat','Kdry')
% 
% figure;
% plot(mu_dry,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('mu-dry & Vs')
% ylabel('depth(m)')
% hold on
% plot(Vs,DEPTH)
% set(gca,'YDir','reverse')
% legend('mu-dry','Vs')
% 
% figure;
% plot(Vs,DEPTH)
% set(gca,'YDir','reverse')
% xlabel('Vs')
% ylabel('depth(m)')
