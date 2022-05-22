clear all;
close all;
nx=64; ny=64; nz=64;
xright = 1.0*pi; h = xright/nx;
yright = 1.0*pi; zright = 1.0*pi;


x=linspace(0.5*h,xright-0.5*h,nx);
y=linspace(0.5*h,yright-0.5*h,ny);
z=linspace(0.5*h,zright-0.5*h,nz);

[yy,xx,zz]=meshgrid(y,x,z);

dt = 0.1*h;

max_it =600;
ns =max_it/30;
time= dt*max_it;



ss=sprintf('./a/dataout/Mene.m'); EnEM1= load(ss);
 ss=sprintf('./a/dataout/Oene.m'); EnEO1 = load(ss);
 ss=sprintf('./b/dataout/Mene.m'); EnEM2= load(ss);
 ss=sprintf('./b/dataout/Oene.m'); EnEO2 = load(ss);
 ss=sprintf('./c/dataout/Mene.m'); EnEM3= load(ss);
 ss=sprintf('./c/dataout/Oene.m'); EnEO3 = load(ss);
len =length(EnEM1)
fig=figure(137);
clf;
hold on;
t = linspace(0,time,len);
plot(t,EnEM1,'ro-.','markersize',10,'linewidth',1);hold on;
plot(t,EnEO1,'bo-','markersize',10,'linewidth',1);hold on;
plot(t,EnEM2,'r*-.','markersize',10,'linewidth',1);hold on;
plot(t,EnEO2,'b*-','markersize',10,'linewidth',1);hold on;
plot(t,EnEM3,'rh-.','markersize',10,'linewidth',1);hold on;
plot(t,EnEO3,'bh-','markersize',10,'linewidth',1);hold on;
legend('Modified energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/4,1/4,1/4,1/4)','Original energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/4,1/4,1/4,1/4) ',...
'Modified energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/5,1/5,1/5,2/5)','Original energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/5,1/5,1/5,2/5) ',...
'Modified energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/6,1/6,1/6,3/6)','Original energy curves for (m_{1},m_{2},m_{3},m_{4})=(1/6,1/6,1/6,3/6)')
set(gca,'fontsize',25);
xlabel('Time');
xlim([0,time]);
ylabel('Energy');
box on;
ss = sprintf('fig43Denergy.eps');
print(fig,'-depsc',ss);
