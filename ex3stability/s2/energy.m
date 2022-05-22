clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=128*2; ny = 128*2;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
time =4;

 ss=sprintf('./stability001/outdata/Mene.m'); EnEM12 = load(ss);
 ss=sprintf('./stability001/outdata/Oene.m'); EnEO12 = load(ss);
ss=sprintf('./stability0005/outdata/Mene.m'); EnEM22 = load(ss);
ss=sprintf('./stability0005/outdata/Oene.m'); EnEO22 = load(ss);
ss=sprintf('./stability00025/outdata/Mene.m'); EnEM32 = load(ss);
ss=sprintf('./stability00025/outdata/Oene.m'); EnEO32 = load(ss);
ss=sprintf('./stability0001/outdata/Mene.m'); EnEM42 = load(ss);
 ss=sprintf('./stability0001/outdata/Oene.m'); EnEO42 = load(ss);




fig=figure(137);
clf
hold on;
t = linspace(0,4,21);
plot(t,EnEM12,'ro--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO12,'bo-','markersize',10,'linewidth',1);hold on;
plot(t,EnEM22,'r*--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO22,'b*-','markersize',10,'linewidth',1);hold on;
plot(t,EnEM32,'rh--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO32,'bh-','markersize',10,'linewidth',1);hold on;
plot(t,EnEM42,'r<--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO42,'b<-','markersize',10,'linewidth',1);hold on;
% title('Modified and Original energy curves for S =2 and four different dt')
legend('Modified energy curves :\Delta t =0.01','Original energy curves :\Delta t =0.01','Modified energy curves :\Delta t=0.005','Original energy curves :\Delta t=0.005','Modified energy curves :\Delta t=0.0025'......
,'Original energy curves: \Delta t =0.0025','Modified energy curves :\Delta t=0.001','Original energy curves :\Delta t=0.001');
xlabel('Time');
ylabel('Energy');
set(gca,'fontsize',25);
box on;

ss = sprintf('fig2energy.eps');
print(fig,'-depsc',ss);




