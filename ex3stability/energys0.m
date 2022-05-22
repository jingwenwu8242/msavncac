clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=128*2; ny = 128*2;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
time =4;



ss=sprintf('./s2/stability0001/outdata/Mene.m'); EnEM42 = load(ss);
 ss=sprintf('./s2/stability0001/outdata/Oene.m'); EnEO42 = load(ss);
ss=sprintf('./s0/stability0001/outdata/Mene.m'); EnEM40 = load(ss);
 ss=sprintf('./s0/stability0001/outdata/Oene.m'); EnEO40 = load(ss);





fig=figure(137);
hold on;
t = linspace(0,4,21);
plot(t,EnEM42,'ro--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO42,'bo-','markersize',10,'linewidth',1);hold on;

plot(t,EnEM40,'r<--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO40,'b<-','markersize',10,'linewidth',1);hold on;


%title('Modified and Original energy curves for S =2 and four different dt')
legend('Modified energy curves :S_{\mu}=2','Original energy curves: S_{\mu}=2','Modified energy curves :S_{\mu}=0','Original energy curves:S_{\mu}=0');
xlabel('Time');
ylabel('Energy');
set(gca,'fontsize',25);

axes('Position',[0.2517 0.1785 0.24 0.24])
plot(t,EnEM42,'ro--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO42,'bo-','markersize',10,'linewidth',1);hold on;

plot(t,EnEM40,'r<--','markersize',10,'linewidth',1);hold on;
plot(t,EnEO40,'b<-','markersize',10,'linewidth',1);hold on;

box on;

ss = sprintf('fig2compare.eps');
print(fig,'-depsc',ss);


