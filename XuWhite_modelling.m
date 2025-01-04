clc
clear all
close all

% using AR log (4 aspect ratio)

%% input    
Vp_log=xlsread("SPD10-08_logs_Rest",'R47:R4052');      %(m/s)
Vs_log=xlsread("SPD10-08_logs_Rest",'S47:S4052');      %(m/s)
DEPTH_log_p=xlsread("SPD10-08_logs_Rest",'B47:B4052');  %(m)
DEPTH_log_S=xlsread("SPD10-08_logs_Rest",'B47:B4052');  %(m)
DEPTH_log_rho=xlsread("SPD10-08_logs_Rest",'B47:B4052'); %(m)
RHOB_log=xlsread("SPD10-08_logs_Rest",'O47:O4052');  %(g/cm3)--->it is supposed to be(kg/m3)

gg=0.68;            %Gas Gravity (Specific Gravity)
Sb=xlsread("spd10-08_multimin",'H48:H2904');%Saturation of water
So=0; 
sal=290000;          %Salinity
P=5290*0.00689;      %pressure (Mpa)
T=103;               %temprature (c)

VOL_ANH=xlsread("spd10-08_multimin",'J48:J2904'); %volume fraction
VOL_CAL=xlsread("spd10-08_multimin",'K48:K2904');
VOL_DOL=xlsread("spd10-08_multimin",'L48:L2904');
VOL_WCS=xlsread("spd10-08_multimin",'M48:M2904');

VOL_ANH=VOL_ANH./(VOL_ANH+VOL_CAL+VOL_DOL+VOL_WCS); %normalized volume fraction
VOL_CAL=VOL_CAL./(VOL_ANH+VOL_CAL+VOL_DOL+VOL_WCS);
VOL_DOL=VOL_DOL./(VOL_ANH+VOL_CAL+VOL_DOL+VOL_WCS);
VOL_WCS=VOL_WCS./(VOL_ANH+VOL_CAL+VOL_DOL+VOL_WCS);

RHO_ANH= 2960;  % density
RHO_CAL=2640;
RHO_DOL=2780;
RHO_WCS=2600;

K_ANH=60.52 *1e+9; % Gpa to Pa
K_CAL=62.71 *1e+9; 
K_DOL=67.58 *1e+9; 
K_WCS=21.8 *1e+9; 

MU_ANH= 37.97 *1e+9; % Gpa to Pa
MU_CAL=26.88 *1e+9;
MU_DOL=33.6 *1e+9;
MU_WCS=6.66 *1e+9;

 for i=1:length(VOL_ANH)
     K_VOIT(i)= (VOL_ANH(i)*K_ANH + VOL_CAL(i)*K_CAL + VOL_DOL(i)*K_DOL + VOL_WCS(i)*K_WCS);
     K_REUSS(i)= (1/(VOL_ANH(i)/K_ANH + VOL_CAL(i)/K_CAL + VOL_DOL(i)/K_DOL + VOL_WCS(i)/K_WCS));
     K_VRH(i)=(K_VOIT(i)+K_REUSS(i))/2;
 end
  for i=1:length(VOL_ANH)
     MU_VOIT(i)= (VOL_ANH(i)*MU_ANH + VOL_CAL(i)*MU_CAL + VOL_DOL(i)*MU_DOL + VOL_WCS(i)*MU_WCS);
     MU_REUSS(i)= (1/(VOL_ANH(i)/MU_ANH + VOL_CAL(i)/MU_CAL + VOL_DOL(i)/MU_DOL + VOL_WCS(i)/MU_WCS));
     MU_VRH(i)=(MU_VOIT(i)+MU_REUSS(i))/2;
  end
  

KmHill=K_VRH';
MumHill=MU_VRH';

%Bulk modulus of mix fluid
[Kfluid,Kreuss,rhoeff,vpb,rhob,Kb,rhog,Kg]=fluid(sal,gg,P,T,So,Sb);

Mu_fluid=0;

% pore aspect ratio
alpha_ANH=0.231;    %aspect ratio
alpha_CAL=0.19;
alpha_DOL=0.184;
alpha_WCS=0.05;      

%  inclusion concentration
X_ANH=VOL_ANH;
X_CAL=VOL_CAL;
X_DOL=VOL_DOL;
X_WCS=VOL_WCS;


% % pore aspect ratio
% alpha_IP=0.18;
% alpha_stiff=0.4;
% alpha_IPstiff=0.25;

DEPTH=xlsread("spd10-08_multimin",'B48:B2904'); % depth
PGD=1;              %Pressure gradient with Depth
effectives=((PGD.*(DEPTH./0.3048))*0.00689)-P; %effective pressure

Pk=6.32;           %MacBeth (2004)
SK=0.6;
Ek=SK/(1-SK);
Pmu=9.5; 
SMU=0.6;
Emu=SMU/(1-SMU);

inclusion=[1,2,3,4]; % 1=ANH  2=CAL  3=DOL  4=WCS 
K_inclusion=0;   % K'
mu_inclusion=0;  % mu'

% AR LOGS
load('ARlog_resistivity.mat') % AR calculated by resistivity
load('ARlog_VDL.mat') % AR calculated by VDL
load('alpha_torDEV.mat') % AR calculated by tortuosity deviation
load('ARlog_Azadpour') % AR calculated from TDL
load('ARlog_Azadpour2') % AR calculated from TDL2

% alpha=[alpha_ANH,alpha_CAL,alpha_DOL,alpha_WCS];
% alpha=0.15;
% alpha=AR_pred;  % from resistivity
% alpha=alpha_log;  %from VDL
% alpha=alpha_tor_deviation; % from tortuosity deviation
alpha=alpha_TDL2; % from TDL_constantTOR
% alpha=alpha_TDL; % from TDL_azadpour


porosity=xlsread("spd10-08_multimin",'E48:E2904'); % porosity (effective)

  
phil=[X_ANH,X_CAL,X_DOL,X_WCS];
% phil=porosity;--->wrong!
% phil=1;
% phil=1-porosity;

%% calculation
phi=(alpha./((1-(alpha).^ 2).^ 1.5)).*(acos(alpha)-alpha.*sqrt((1-(alpha).^ 2)));
g=(((alpha.^2)./(1- alpha.^2))).*((3.* phi)-2);
R=(3.*MumHill)./(3.*KmHill + 4.*MumHill);
B=1/3.*((K_inclusion./KmHill)-(mu_inclusion./MumHill));
A=(mu_inclusion./MumHill)-1;
F9=A.*(g.*(R-1) - R.*phi) + B.*phi.*(3-4.*R);
F8=A.*(1-(2.*R)+(g./2).*(R-1)+(phi./2).*(5.*R -3)) + B.*(1-phi).*(3-(4.*R));
F7=2+ (A./4).*(9.*phi + 3.*g - R.*(5.*phi + 3.*g)) + B.*phi.*(3-(4.*R));
F6=1+A.*(1+g-R.*(g+phi)) +B.*(1-phi).*(3-4.*R);
F5=A.*(R.*(g+phi - 4/3)-g) + B.*phi.*(3-(4.*R));
F4=1+ A./4 .*(3.*phi+g - R.*(g-phi));
F3=1+A./2 .*(R.*(2-phi) + ((1+alpha.^2)./alpha.^2).*(g.*(R-1)));
F2=1+ A.*(1 + 3/2.*(g+phi) - R/2.*(3.*g +5.*phi)) + B.*(3- 4.*R) + A./2 .*(A+ 3.*B).*...
    (3- 4.*R).*(g+phi - R.*(g-phi + 2.*(phi.^2)));
F1=1+ A.*(3/2 .*(g+phi) - R.*((3/2 .*g) + (5/2 .*phi) - 4/3));

%sigma for Kdry
Tiijj=(3.*F1)./F2; 
s=(phil.*Tiijj);
sigma=zeros(length(s),1);
for i1=1:length(s)
    sigma(i1)=sum(s(i1,:));
end

%sigma for Mu_dry
FF=(2./F3)+(1./F4)+(((F4.*F5) +(F6.*F7) +(F8.*F9))./(F2.*F4)); 
s2=(phil.*FF);
sigma2=zeros(length(s2),1);
for i2=1:length(phil)
    sigma2(i2)=sum(s2(i2,:));
end




%kdry----------------------------------------------------------
Kinf=(4*MumHill.*porosity.*(sigma/3).*(K_inclusion-KmHill)+KmHill.*(3*KmHill+4*MumHill))./(3*KmHill+4*MumHill-3*porosity.*(sigma/3).*(K_inclusion-KmHill));

%mu_dry---------------------------------------------------------
Mu_inf=MumHill.*(porosity.*sigma2.*(mu_inclusion-MumHill).*(9*KmHill+8*MumHill)+25*MumHill.*(3*KmHill+4*MumHill))./(25*MumHill.*(3*KmHill+4*MumHill)-6*(KmHill+2*MumHill).*porosity.*sigma2.*(mu_inclusion-MumHill));


%macbeth 2004 (pressure effect)------------------------------------------------
%Kdry:
Kdry=(Kinf)./(1+(Ek.*exp(-effectives./Pk)));

mu_dry=(Mu_inf)./(1+(Emu.*exp(-effectives./Pmu)));


%% gassmann

   for i=1:length(VOL_ANH)
     RHO_VOIT(i)= (VOL_ANH(i)*RHO_ANH + VOL_CAL(i)*RHO_CAL + VOL_DOL(i)*RHO_DOL + VOL_WCS(i)*RHO_WCS);
%      RHO_REUSS(i)= 1/(VOL_ANH(i)/RHO_ANH + VOL_CAL(i)/RHO_CAL + VOL_DOL(i)/RHO_DOL + VOL_WCS(i)/RHO_WCS);
%      RHO_VRH(i)=(RHO_VOIT(i)+RHO_REUSS(i))/2;
     RHO_VRH(i)=RHO_VOIT(i); % just VOIT Average is used
 end
 
rhom=RHO_VRH';   % density_matrix

%saturated bulk modulus
Ksat=Kdry+(((1-(Kdry./KmHill)).^2)./((porosity./Kfluid)+((1-porosity)./KmHill)-(Kdry./(1-KmHill).^2)));

%saturated density / or bulk density
rho_sat=(1-porosity).*rhom + porosity.*rhoeff;

% velocity----------------------------------------------------------------

Vp=sqrt((Ksat+ 4/3.*mu_dry)./(rho_sat)); %compressional velocity (m/s)
Vs=sqrt(mu_dry./(rho_sat)); % shear velocity (m/s)

%% figures
figure;
plot(Kdry,DEPTH)
xlabel('Kdry & km (Pa)')
ylabel('Depth(m)')
set(gca,'YDir','reverse')
hold on
plot(KmHill,DEPTH)
set(gca,'YDir','reverse')
% 
% figure;
% plot(Ksat,DEPTH)
% xlabel('Ksat(Pa)')
% ylabel('Depth(m)')
% set(gca,'YDir','reverse')
% 
% figure;
% plot(Vp,DEPTH)
% xlabel('Vp(m/s)')
% ylabel('Depth(m)')
% set(gca,'YDir','reverse')

%% velocities validation
% Equalization the DEPTHS for Vp ----------------------------------------
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

% Equalization the DEPTHS for Vs -----------------------------------------
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
      RR(l2)=n2(iii); %DEPTH_log_Vs index   
      
end
RR2=RR';
Vs_measured=Vs_log(RR2);
DEPTH_measured_S=DEPTH_log_S(RR2);

%===========================================================================
% figures: Vp & Vs  -----------------------------------------------------
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
xlim([2500 7000])


%% density validation

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

%% final figure
%needed data:

measured_AR=xlsread("AR_azad",'D2:D39');    %experimental AR
valid_depth_AR=xlsread("AR_azad",'E2:E39');
measured_alpha=[]; x=0;
for ii2=1:length(measured_AR)
if measured_AR(ii2) < 0.13
   x=x+1;
   measured_alpha(x)=0.15; %AR for negetive deviation
   elseif 0.2 < measured_AR(ii2)
    x=x+1;
   measured_alpha(x)=0.3; %AR for positive deviation
   elseif 0.13 <= measured_AR(ii2) && measured_AR(ii2) <=0.2 
    x=x+1;
   measured_alpha(x)=0.19; %AR for 0 deviation
end
end

tortuosity=xlsread("AR_azad",'G2:G2858');   %predicted tortuosity factor

measured_tortuosity=xlsread("AR_azad",'H2:H70');   %experimental tortuosity factor
valid_depth_tor=xlsread("AR_azad",'I2:I70'); 

LLD=xlsread("AR_azad",'K2:K2858');    %LLD EDIT (OHMM)
%----------------------------------------
figure;
subplot(1,6,1)
plot(alpha,DEPTH)
xlabel('predicted aspect ratio')
ylabel('Depth(m)')
hold on
plot(measured_alpha,valid_depth_AR,'.r')
set(gca,'YDir','reverse')
xlim([0.1 0.4])

subplot(1,6,2)
plot(rho_sat,DEPTH)
xlabel('predicted RHO(kg/m3)')
set(gca,'YDir','reverse')

subplot(1,6,3)
plot(Vp,DEPTH)
xlabel('predicted Vp(m/s)')
set(gca,'YDir','reverse')

subplot(1,6,4)
plot(Vs,DEPTH)
xlabel('predicted Vs(m/s)')
set(gca,'YDir','reverse')

subplot(1,6,5)
plot(tortuosity,DEPTH)
xlabel('predicted tortuosity')
hold on
plot(measured_tortuosity,valid_depth_tor,'.r')
set(gca,'YDir','reverse')

subplot(1,6,6)
plot(LLD,DEPTH)
xlabel('LLD EDIT (OHMM)')
set(gca,'YDir','reverse')

%% crossplot
% Vp(129)=6000;
% Vp(130)=6000;

figure;
plot(Vp,Vp_measured,'.r','MarkerSize',10)
xlabel("Vp est")
ylabel("Vp measured")

%% error (meter per s)
MSE_p=mean((Vp-Vp_measured).^2);
RMSE_p=sqrt(MSE_p)%/length(Vp)
RRMSE_P=(RMSE_p ./(max(Vp_measured)-min(Vp_measured)))*100

MSE_s=mean((Vs-Vs_measured).^2);
RMSE_s=sqrt(MSE_s)%/length(Vp)
RRMSE_s=(RMSE_s ./(max(Vs_measured)-min(Vs_measured)))*100

MSE_rhob=mean((rho_sat-RHOB_measured).^2);
RMSE_rhob=sqrt(MSE_rhob)%/length(Vp)
