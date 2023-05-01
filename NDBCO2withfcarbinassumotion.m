clear all;
load('RMgPh.mat');
ageinter1=0:0.01:2000;
RMginter1=interp1(AgePh,RMgPh,ageinter1,'linear');
for i=1:200001
   
    decay=i:1:200001;
    Rdecay=exp(-0.001.*(decay-i)/100);
    swpre(i)=sum(RMginter1(i:end).*Rdecay)/sum(Rdecay);
end
ageinter2=0:0.001:2000;
swpreinter=interp1(ageinter1,swpre,ageinter2,'linear');
fcarbininter=0.42.*(3000-ageinter2)./3000;
RMginter=interp1(AgePh,RMgPh,ageinter2,'linear');
%% pick random timespot for calculation
for i=1:100
tt=rand(1,500000000)*2000;
tt=fix(tt*1000);ttindex=tt+1;
a=1:1:length(tt);
% plot(tt,a,'o','Markersize',1);

% definition of paramaters
% fcarbin=[0 1]; fcarbout=[0 1]; carbin=[-2.1 -1.5]; frsil=[0 0.8];
fcarbin=fcarbininter(ttindex)+normrnd(0,0.04,[1 length(tt)]);
indexss=find(fcarbin>0&fcarbin<1);
fcarbin=fcarbin(indexss);
ttindex=ttindex(indexss);
fcarbout=rand(1,length(fcarbin));
 carbin=swpreinter(ttindex)+normrnd(-1.5,0.2,[1 length(fcarbin)]);silin= normrnd(-0.4,0.1,[1 length(fcarbin)]);
frsil=normrnd(0.8,0.4,[1 length(fcarbin)]);
frcarb=normrnd(-1.5,0.2,[1 length(fcarbin)]);
RMgsolv=RMginter(ttindex);
sw=RMgsolv+normrnd(0,0.1,[1 length(fcarbin)]);
misfit=fcarbin.*(carbin-silin)-fcarbout.*(frcarb-frsil)-frsil-sw+silin;
Df=abs(misfit);
sindex=find(Df<0.001);
fcarbinsolved=fcarbin(sindex);fcarboutsolved=fcarbout(sindex);
carbinsolved=carbin(sindex);frsilsolved=frsil(sindex);
ttsolved=ttindex(sindex)-1;swsolved=sw(sindex);
frcarbsolved=frcarb(sindex);silinsolved=silin(sindex);
R=(fcarboutsolved-fcarbinsolved);
R1=(1-fcarbinsolved)./(1-fcarboutsolved);
if i==1
Rsolveds(i,:)=R(1:250000);
ttsolveds(i,:)=ttsolved(1:250000);
R1solveds(i,:)=R1(1:250000);
fcarbinsolveds(i,:)=fcarbinsolved(1:250000);
fcarboutsolveds(i,:)=fcarboutsolved(1:250000);
else 
Rsolveds=[Rsolveds R(1:250000)];
ttsolveds=[ttsolveds ttsolved(1:250000)];
R1solveds=[R1solveds R1(1:250000)];
fcarbinsolveds=[fcarbinsolveds fcarbinsolved(1:250000)];
fcarboutsolveds=[fcarboutsolveds fcarboutsolved(1:250000)];
end 
i
end

%%

X=ttsolveds/1000;
Y=Rsolveds;
Xmin=min(X);Xmax=max(X);
Ymin=min(Y);Ymax=max(Y);

Nx=300;
Ny=300;

Xedge=linspace(Xmin,Xmax,Nx);
Yedge=linspace(Ymin,Ymax,Ny);

[N,~,~,binX,binY] = histcounts2(X,Y,[-inf,Xedge(2:end-1),inf],[-inf,Yedge(2:end-1),inf]);

XedgeM=movsum(Xedge,2)/2;
YedgeM=movsum(Yedge,2)/2;

[Xedgemesh,Yedgemesh]=meshgrid(XedgeM(2:end),YedgeM(2:end));


XX=ttsolveds/1000;YY=Rsolveds;
XX=round(XX);

j=1;
for ii=1:2001
   mRsort(j)=mean(YY(find(XX==(ii-1))));
  eRsort(j)=std(YY(find(XX==(ii-1))));
   Agesort(j)=(ii-1);
    j=j+1;
end
%%
mmRsort(1)= mRsort(1);Agesort2(1)=0;
eeRsort(1)=eRsort(1);
for iii=1:200
    mmRsort(iii+1)=mRsort(iii*10+1);
    eeRsort(iii+1)=eRsort(iii*10+1);
   Agesort2(iii+1)=iii*10;
end
%% 

figure(1)
pcolor(Xedgemesh,Yedgemesh,N');hold on;
plot(Agesort2,mmRsort);hold on;plot(Agesort2,mmRsort-eeRsort);plot(Agesort2,mmRsort+eeRsort);hold off;
ylim([-0.2 0.35]);   xlim([0 500]);
