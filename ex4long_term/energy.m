clear all; 

xright=2.0; yright=1.0; xleft = 0; yleft = 0; nx=128*2; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
time =4;

 ss=sprintf('./data1/Mene.m'); EnEM1 = load(ss);
%  ss=sprintf('./data1/Oene.m'); EnEO1 = load(ss);
ss=sprintf('./data2/Mene.m'); EnEM2 = load(ss);
% ss=sprintf('./data2/Oene.m'); EnEO2 = load(ss);
ss=sprintf('./data3/Mene.m'); EnEM3 = load(ss);
% ss=sprintf('./data3/Oene.m'); EnEO3 = load(ss);
ss=sprintf('./data4/Mene.m'); EnEM4 = load(ss);
%  ss=sprintf('./data4/Oene.m'); EnEO4 = load(ss);




fig=figure(137);
hold on;
t = linspace(0,80,41);
plot(t(1:2:end),EnEM1(1:2:end),'bo--','markersize',10,'linewidth',1);hold on;
% plot(t,EnEO1,'bo-','markersize',7,'linewidth',0.8);hold on;
plot(t(1:2:end),EnEM2(1:2:end),'r*--','markersize',10,'linewidth',1);hold on;
% plot(t,EnEO2,'b*-','markersize',9,'linewidth',0.8);hold on;
plot(t(1:2:end),EnEM3(1:2:end),'mh--','markersize',10,'linewidth',1);hold on;
% plot(t,EnEO3,'bh-','markersize',10,'linewidth',0.8);hold on;
plot(t(1:2:end),EnEM4(1:2:end),'k<--','markersize',10,'linewidth',1);hold on;
% plot(t,EnEO4,'b<-','markersize',8,'linewidth',0.8);hold on;
legend('Modified energy curves \Delta t =0.1','Modified energy curves \Delta t=0.05','Modified energy curves \Delta t=0.01'......
,'Modified energy curves \Delta t=0.005');
xlabel('Time');
ylabel('Energy');
set(gca,'fontsize',24);
box on;
ss = sprintf('longtermenergy.eps');
print(fig,'-depsc',ss);






