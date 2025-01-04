clc
clear all
close all

%% input

% a=1;  % tortuosity factor
m=2; % cementation factor
n=2;  % saturation factor
Rw=0.013; % water resistivity
load('all_logs.mat')
load('all_logs_name.mat')

Rt=C(:,11); % total resistivity (I used LLD)
porosity=C(:,28); % effective porosity
Sw=C(:,31); % effective water saturation
Vsh=C(:,36); % shale volume
Rsh=6;  % shale resistivity
Rw=0.0138;  % water resistivity

% archie
% a=((Sw.^n).*(porosity.^m).*(Rt))./(Rw); % tortuosity
% % deleting noises
% a(1562)=mean(a);
% a(1561)=mean(a);
% a(428)=mean(a);

%indonesia
S= Sw.^(n/2);
B1=sqrt(1./Rt);
F=((Vsh).^(1-(0.5.*Vsh)))./(Rsh);
D=((B1-(S.*F))./(S)).^2;
a=(porosity.^m)./(Rw.*D); % tortuosity
% deleting noises
a(1562)=mean(a);
a(1561)=mean(a);
a(428)=mean(a);

% archie and indonesia have excactly same results

figure;
plot(a,C(:,1))
xlabel("a")
ylabel("depth (m)")
set(gca,'YDir','reverse')
title('formula tortuosity')


%% Tortuosity log
% a derived from **Decision Tree** method
Pred=xlsread('AZADPOUR_a','A1:A2858');
%---------------------

TOR_pred=Pred;
 figure;
plot(TOR_pred,C(:,1)) 
set(gca,'YDir','reverse')
title("tortuosity NN")
xlim([0 4])

%main formula
TDL=TOR_pred-1; 


figure;
% subplot(1,2,1)
% plot(a1,C(:,1))
% subplot(1,2,2)
plot(TDL,C(:,1),'linewidth',1)
xlabel("a",'fontsize',12)
ylabel("Depth (m)",'fontsize',12)
set(gca,'YDir','reverse')
title('Tortuosity Deviation Log')
xlim([-2 4])

%% estimation of AR log

% for 2 categories
alpha=[]; x=0;
for ii2=1:length(TDL)
if TDL(ii2) <= 0
   x=x+1;
   alpha(x)=0.18; %AR for negetive deviation
   elseif 0 < TDL(ii2)
    x=x+1;
   alpha(x)=0.25; %AR for positive deviation
%    elseif -2 <= TDL(ii2) && TDL(ii2) <=2 
%     x=x+1;
%    alpha(x)=0.19; %AR for 0 deviation
end
end
alpha_TDL2=alpha'; %log of AR

%% validation of AR

% for 2 categories
measured_AR=xlsread("ARmatch",'B2:B205');    %experimental AR
valid_depth_AR=xlsread("ARmatch",'A2:A205');
measured_alpha=[]; x2=0;
for ii22=1:length(measured_AR)
if measured_AR(ii22) < 0.25
   x2=x2+1;
   measured_alpha(x2)=0.19; %AR for interparticles
   elseif 0.25 <= measured_AR(ii22)
    x2=x2+1;
   measured_alpha(x2)=0.3; %AR for stiff pores

end
end

%% figures 
figure;
plot(alpha_TDL2,C(:,1),'linewidth',1)
xlabel(" Aspect Ratio",'fontsize',12)
ylabel("Depth (m)",'fontsize',12)
set(gca,'YDir','reverse')
title('Log of aspect ratio')
xlim([0.15 0.3])
%------------------------------
figure;
hold on;
idx_015 = alpha_TDL2 == 0.15;
idx_019 = alpha_TDL2 == 0.19;
idx_03 = alpha_TDL2 == 0.3;

scatter(alpha_TDL2(idx_015), C(idx_015, 1), 5, 'g', 'filled');
scatter(alpha_TDL2(idx_019), C(idx_019, 1), 5, 'r', 'filled');
scatter(alpha_TDL2(idx_03), C(idx_03, 1), 5, 'b', 'filled');
set(gca,'YDir','reverse')

hold off;

hold on
plot(measured_alpha,valid_depth_AR,'ok')
% xlim([0 0.4])
% set(gca,'YDir','reverse')
legend('0.15', '0.19', '0.3','experimental AR');
xlabel("Aspect Ratio")
ylabel("Depth (m)")
xlim([0.1 0.4])
set(gca,'YDir','reverse')
title('Log of Aspect Ratio')

%-----------------------
figure;
plot(alpha_TDL2,C(:,1))
xlabel('predicted aspect ratio')
ylabel('depth(m)')
set(gca,'YDir','reverse')
xlim([0 0.4])


%-----------------------
figure;
scatter(alpha_TDL2,C(:,1),5,alpha_TDL2,'filled')
colormap([0 0.5 0; 1 0 0])
c = colorbar;
c.Ticks = [0.18 0.25];
c.Label.String = 'Aspect Ratio';
xlabel("Aspect Ratio")
ylabel("Depth (m)")
xlim([0.15 0.3])
% ylim([2980 3000])
ylim([2920 3050])
set(gca,'YDir','reverse')
title('Log of Aspect Ratio')

%% just for checking result with thin section
% the final result of this part is in excel file VDL sheet2
% you can comment this part

% load('experimental_AR.mat')
% exp_AR_depth=A(:,1);
% exp_AR=A(:,2);
% 
% R=[];l=0;
% for ii=1:size(exp_AR_depth)
%     for j=1:size(C(:,1))
%      e(j)=abs(exp_AR_depth(ii,1)- C(j,1));
%     end
%      [m(ii),k(ii)]=min(e);  
%       l=l+1;
%       R(l)=k(ii); %tortuosity deviation_depth index   
% end
% RR=R';
% equal_tor_depth=C(RR,1);
% 
% W=TDL(RR);

