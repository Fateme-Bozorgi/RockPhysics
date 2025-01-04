
%% data Prepration
clear;clc
A=xlsread("spd10-08_wire_final_cor.xlsx",'A58:F122829');

B= ["DEPTH","BS","CALI","CGR","DRHO","DT","DTSM","DT EDIT","GR","LLD","LLD EDIT","LLS","LLS EDIT", ...
    "NPHI","PEF","POTA","RHOB","RHOB EDIT","RLA0","RLA1","RLA2","RXO","SGR","THOR","URAN"];

l=size(A,1)/6;

C=zeros(l,length(B));

for i=1:l
C(i,1)=A((i-1)*6+1,1);
C(i,2:6)=A((i-1)*6+2,2:6);
C(i,7:11)=A((i-1)*6+3,2:6);
C(i,12:16)=A((i-1)*6+4,2:6);
C(i,17:21)=A((i-1)*6+5,2:6);
C(i,22:25)=A((i-1)*6+6,2:5);
end

D=xlsread("spd10-08_multimin.xlsx",'B48:M2904');

b=["DEPTH","T PRED","NPHI COR_PRED","PHIE","PHIT","RHOB COR_PRED","SWE","SWT","VOL ANHYDR", ...
    "VOL CALCITE","VOL DOLOM","VOL WCS"];
C(:,26:36)=-999.25*ones(size(C,1),11);
C(17445:20301,26:36)=D(:,2:end);% All well-logs
B=[B,b(2:end)];% well-log names

Markesr_name=["Ilam","Sarvak","Kazhdumi","Dariyan","Lower Dariyan","Upper Gadvan","Khalij Mbr", ...
"Lower Gadvan","Fahliyan","Hith","Upper Surmeh","Mand Member","Lower Limestone", ...
"Lower Surmeh Shale","S8 lower Sudair","K1","K2","K3","K4","Nar"];
Markesr_depth=[1041.31,1136.3,1194.91,1229,1308.1,1335,1352.95,1367.38,1387.52,1561.03,1627.21, ...
1982.47,2078.78,2192.54,2523.94,2756.33,2873,2920.01,3038.36,3204.58];

Ar= xlsread("archie.xlsx",'B19:D87'); 
Ar_names=["depth","a","m"];

inputs=[4 9 10 11 12 13 15 19 20 21 22 23 24 25 28 29 31 32 33 34 35 36];inputs2=[11,28,32,36];% columns of C
output=[2 ];%columns of Ar
Ar(:,4:4+length(inputs)-1)=zeros(size(Ar,1),length(inputs));

for i=1:size(Ar,1)
[~,k]=min(abs(C(:,1)-Ar(i,1)));
Ar(i,4:4+length(inputs)-1)=C(k,inputs);
end


ind=[];c=0;
for i=1:size(C,1)
if sum(double(C(i,:)==-999.25))>10
c=c+1;
ind(c)=i;
end
end
C(ind,:)=[];

    %% Train NN
hiddenLayerSize =[15 10];
net = fitnet(hiddenLayerSize,'trainscg');
% net = feedforwardnet(hiddenLayerSize);

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 90/100;
net.divideParam.valRatio = 5/100;
net.divideParam.testRatio = 5/100;
% Train the Network
IN=Ar(:,4:end);
% Two classes
OUT=Ar(:,output);

net.trainParam.max_fail=150;
[net,tr] = train(net,IN',OUT');

% Test the Network
outputs = net(IN');
errors = gsubtract(OUT',outputs);
performance = perform(net,OUT',outputs);

% View the Network
view(net)

% Plots
figure, plotperform(tr)
figure, plottrainstate(tr)
figure, plotregression(OUT',outputs)
figure, ploterrhist(errors)
%  save("net1","net", "tr")

%% Apply NN

 load net
Pred=net(C(:,inputs)')';


for i=1:length(inputs2)
subplot(1,length(inputs2)+size(Pred,2)+1,i+1)
plot(C(:,inputs2(i)),C(:,1))
set(gca,'YDir','reverse')
title(B(inputs2(i)))
if i==1; ylabel("depth(m)");end
end


for i=length(inputs2)+1:length(inputs2)+size(Pred,2)
subplot(1,length(inputs2)+size(Pred,2)+1,i+1)
plot(Pred(:,i-length(inputs2)),C(:,1))
set(gca,'YDir','reverse')
title(Ar_names(1+i-length(inputs2)))
hold on
plot(Ar(:,i-length(inputs2)+1),Ar(:,1),'.r','MarkerSize',5)
hold off
end

legend("Pridicted","Measured")

subplot(1,length(inputs2)+size(Pred,2)+1,1)
set(gca,'YDir','reverse')
hold on
title("Markers")
for i=1:length(Markesr_depth)
plot(0:1,[Markesr_depth(i) Markesr_depth(i)],'black')
text(-1,Markesr_depth(i),Markesr_name(i))
end

%% figures
% figure;
% subplot(1,2,1)
% markers_depth_new=[2756.33000000000,2873,2920.01000000000,3038.36000000000];;
% markers_name_new=["K1","K2","K3","K4"];
% for ii=1:length(markers_depth_new)
% plot(0:1,[markers_depth_new(ii) markers_depth_new(ii)],'black')
% text(-1,markers_depth_new(ii),markers_name_new(ii))
% end
figure;
plot(Pred,C(:,1))
set(gca,'YDir','reverse')
xlabel('a')
ylabel('Depth(m)')
title('Tortuosity Factor')
hold on
plot(Ar(:,2),Ar(:,1),'.r','MarkerSize',5)

B2=B(inputs);
