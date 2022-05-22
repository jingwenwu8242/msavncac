clear all; 

xright=4.0; yright=4.0; xleft = 0; yleft = 0; nx=128; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt =h ;max_it=300;
time = dt*max_it;

ss=sprintf('./figa/dataout/Mene.m');EnEMa= load(ss);
ss=sprintf('./figb/dataout/Mene.m');EnEMb= load(ss);
ss=sprintf('./figb/dataout/Mene.m');EnEMc= load(ss);




len=length(EnEMc)

fig=figure(137);
hold on;
t = linspace(0,time,len);
plot(t,EnEMa,'ro-','markersize',10,'linewidth',1);hold on;
plot(t,EnEMb,'b*-','markersize',10,'linewidth',1);hold on;
plot(t,EnEMc,'mh-','markersize',10,'linewidth',1);hold on;
xlabel('Time');
ylabel('Energy');
xlim([0,time]);
%title('Modified energy curves');
legend('m=(1/4,1/4,1/4)','m=(1/5,1/5,1/5) ','m=(1/6,1/6,1/6)');
box on;
set(gca,'fontsize',25);
ss = sprintf('fig4energy.eps');
print(fig,'-depsc',ss);
