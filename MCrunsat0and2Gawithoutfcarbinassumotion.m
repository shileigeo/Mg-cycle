clear all;
load('RMgPh.mat');
L=100000000;
ageinter=0:0.01:2000;
RMginter1=interp1(AgePh,RMgPh,ageinter,'linear');
for i=1:200001
  
    decay=i:1:200001;
    Rdecay=exp(-0.001.*(decay-i)/100);
    swpreinter(i)=sum(RMginter1(i:end).*Rdecay)/sum(Rdecay);
end
RMg2Ga=RMginter1(end);
%% modern data
% fcarbin=[0 1]; fcarbout=[0 1]; carbin=[-2.1 -1.5]; frsil=[0 0.8];
for i=1:20
fcarbin=rand(1,L);
fcarbout=rand(1,length(fcarbin));
carbin=swpreinter(1)+normrnd(-1.5,0.2,[1 length(fcarbin)]);
silin= normrnd(-0.4,0.1,[1 length(fcarbin)]);
frsil=normrnd(0.8,0.4,[1 length(fcarbin)]);
frcarb=normrnd(-1.5,0.2,[1 length(fcarbin)]);
misfit=fcarbin.*(carbin-silin)-fcarbout.*(frcarb-frsil)-frsil+0.83+silin;
Df=abs(misfit);
sindex=find(Df<0.001);
fcarbinsolved=fcarbin(sindex);fcarboutsolved=fcarbout(sindex);

if i==1
fcarbinsolveds(i,:)=fcarbinsolved(1:50000);
fcarboutsolveds(i,:)=fcarboutsolved(1:50000);

else 
fcarbinsolveds=[fcarbinsolveds fcarbinsolved(1:50000)];
fcarboutsolveds=[fcarboutsolveds fcarboutsolved(1:50000)];

end
i
end
mfcarbinsolved=mean(fcarbinsolveds);
mfcarboutsolved=mean(fcarboutsolveds);

%% 2Ga
% fcarbin=[0 1]; fcarbout=[0 1]; carbin=[-2.1 -1.5]; frsil=[0 0.8];
for i=1:20
fcarbin=rand(1,L);
fcarbout=rand(1,length(fcarbin));
carbin=swpreinter(length(swpreinter))+normrnd(-1.5,0.2,[1 length(fcarbin)]);
silin= normrnd(-0.4,0.1,[1 length(fcarbin)]);
frsil=normrnd(0.8,0.4,[1 length(fcarbin)]);
frcarb=normrnd(-1.5,0.2,[1 length(fcarbin)]);
misfit1=fcarbin.*(carbin-silin)-fcarbout.*(frcarb-frsil)-frsil-RMg2Ga+silin;
Df1=abs(misfit1);
sindex1=find(Df1<0.001);
fcarbinsolved1=fcarbin(sindex1);fcarboutsolved1=fcarbout(sindex1);
if i==1
fcarbinsolveds1(i,:)=fcarbinsolved1(1:50000);
fcarboutsolveds1(i,:)=fcarboutsolved1(1:50000);

else 
fcarbinsolveds1=[fcarbinsolveds1 fcarbinsolved1(1:50000)];
fcarboutsolveds1=[fcarboutsolveds1 fcarboutsolved1(1:50000)];


end 
j=i
end

mfcarbinsolved1=mean(fcarbinsolveds1);
mfcarboutsolved1=mean(fcarboutsolveds1);
fcarbinsolvedall=[fcarbinsolveds fcarbinsolveds1];
fcarboutsolvedall=[fcarboutsolveds fcarboutsolveds1];

% plot(fcarbinsolved,fcarboutsolved,'o','Markersize',1);

%%
X=fcarbinsolvedall;
Y=fcarboutsolvedall;
Xmin=min(X);Xmax=max(X);
Ymin=min(Y);Ymax=max(Y);

Nx=200;
Ny=200;


Xedge=linspace(Xmin,Xmax,Nx);
Yedge=linspace(Ymin,Ymax,Ny);


[N,~,~,binX,binY] = histcounts2(X,Y,[-inf,Xedge(2:end-1),inf],[-inf,Yedge(2:end-1),inf]);

XedgeM=movsum(Xedge,2)/2;
YedgeM=movsum(Yedge,2)/2;

[Xedgemesh,Yedgemesh]=meshgrid(XedgeM(2:end),YedgeM(2:end));


%%

fcarbinsolveds2=fcarbinsolveds(1:100000);fcarboutsolveds2=fcarboutsolveds(1:100000);
fcarbinsolveds3=fcarbinsolveds1(1:100000);fcarboutsolveds3=fcarboutsolveds1(1:100000);
scatter(fcarbinsolveds2,fcarboutsolveds2,0.1,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
hold on;
scatter(fcarbinsolveds3,fcarboutsolveds3,0.1,'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
xlabel('fcarbin');ylabel('fcarbout');

