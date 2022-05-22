clear all;
close all;


xright=2.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
max_it =300;
dt =h;
time = dt*max_it;
ns = max_it/30;
count =0;


ss=sprintf('./data1/datac.m'); phia1 = load(ss);
ss=sprintf('./data1/datac2.m'); phia2 = load(ss);




ss=sprintf('./data2/datac.m'); phib1 = load(ss);
ss=sprintf('./data2/datac2.m'); phib2 = load(ss);



ss=sprintf('./data3/datac.m'); phic1 = load(ss);
ss=sprintf('./data3/datac2.m'); phic2 = load(ss);


ss=sprintf('./data4/datac.m'); phid1 = load(ss);
ss=sprintf('./data4/datac2.m'); phid2 = load(ss);


for i=20
%     i = [1,15,20,25,30,40]%i=1:30
    count=count+1;
    
fig= figure(i);
 clf;
 
 t = (i-1)*ns*dt;

Aa = phia1((i-1)*nx+1:i*nx,:);
Ba = phia2((i-1)*nx+1:i*nx,:);
Ca = 1-Aa-Ba;


Ab = phib1((i-1)*nx+1:i*nx,:);
Bb = phib2((i-1)*nx+1:i*nx,:);
Cb = 1-Ab-Bb;



Ac = phic1((i-1)*nx+1:i*nx,:);
Bc= phic2((i-1)*nx+1:i*nx,:);
Cc = 1-Ac-Bc;


Ad = phid1((i-1)*nx+1:i*nx,:);
Bd= phid2((i-1)*nx+1:i*nx,:);
Cd = 1-Ad-Bd;

hh1a=contourf(xx,yy,Aa',[0.5 0.5],'facecolor','None'); 
hh2a=contourf(xx,yy,Ba',[0.5 0.5],'facecolor','None'); 
hh3a=contourf(xx,yy,Ca',[0.5 0.5],'facecolor','None'); 

hh1b=contourf(xx,yy,Ab',[0.5 0.5],'facecolor','None'); 
hh2b=contourf(xx,yy,Bb',[0.5 0.5],'facecolor','None'); 
hh3b=contourf(xx,yy,Cb',[0.5 0.5],'facecolor','None'); 


hh1c=contourf(xx,yy,Ac',[0.5 0.5],'facecolor','None'); 
hh2c=contourf(xx,yy,Bc',[0.5 0.5],'facecolor','None'); 
hh3c=contourf(xx,yy,Cc',[0.5 0.5],'facecolor','None'); 


hh1d=contourf(xx,yy,Ad',[0.5 0.5],'facecolor','None'); 
hh2d=contourf(xx,yy,Bd',[0.5 0.5],'facecolor','None'); 
hh3d=contourf(xx,yy,Cd',[0.5 0.5],'facecolor','None'); 

p1=plot(hh1a(1,2:end),hh1a(2,2:end),'b-','markersize',12,'linewidth',1.2);hold on;
p2=plot(hh1b(1,2:end),hh1b(2,2:end),'r-','markersize',12,'linewidth',1.2);hold on;
p3=plot(hh1c(1,2:end),hh1c(2,2:end),'g-','markersize',12,'linewidth',1.2);hold on;
p4=plot(hh1d(1,2:end),hh1d(2,2:end),'k-','markersize',12,'linewidth',1.2);hold on;

plot(hh2a(1,2:end),hh2a(2,2:end),'b-','markersize',12,'linewidth',1.2);hold on;
plot(hh2b(1,2:end),hh2b(2,2:end),'r-','markersize',12,'linewidth',1.2);hold on;
plot(hh2c(1,2:end),hh2c(2,2:end),'g-','markersize',12,'linewidth',1.2);hold on;
plot(hh2d(1,2:end),hh2d(2,2:end),'k-','markersize',12,'linewidth',1.2);hold on;

leg=legend([p1,p2,p3,p4],'$$\epsilon = h/[2\sqrt{2}\tanh^{-1}(0.9)]$$','$$\epsilon = 2h/[2\sqrt{2}\tanh^{-1}(0.9)]$$','$$\epsilon = 3h/[2\sqrt{2}\tanh^{-1}(0.9)]$$','$$\epsilon = 4h/[2\sqrt{2}\tanh^{-1}(0.9)]$$');
set(leg, 'Interpreter','latex');









xlabel('x');
ylabel('y');

axis image; axis([0 xright 0 yright]); 
set(gca,'fontsize',25)
box on
axes('Position',[0.633928571428572 0.227227227227227 0.19047619047619 0.198198198198198]);
plot(hh1a(1,2:end),hh1a(2,2:end),'b-','markersize',12,'linewidth',1.2);hold on;
plot(hh2a(1,2:end),hh2a(2,2:end),'b-','markersize',12,'linewidth',1.2);hold on;
% plot(hh3a(1,2:end),hh3a(2,2:end),'b-','markersize',12,'linewidth',3);hold on;


plot(hh1b(1,2:end),hh1b(2,2:end),'r-','markersize',12,'linewidth',1.2);hold on;
plot(hh2b(1,2:end),hh2b(2,2:end),'r-','markersize',12,'linewidth',1.2);hold on;
%  plot(hh3b(1,2:end),hh3b(2,2:end),'k-','markersize',12,'linewidth',3);hold on;
% 

plot(hh1c(1,2:end),hh1c(2,2:end),'g-','markersize',12,'linewidth',1.2);hold on;
plot(hh2c(1,2:end),hh2c(2,2:end),'g-','markersize',12,'linewidth',1.2);hold on;
% plot(hh3c(1,2:end),hh3c(2,2:end),'r-','markersize',12,'linewidth',3);hold on;


plot(hh1d(1,2:end),hh1d(2,2:end),'k-','markersize',12,'linewidth',1.2);hold on;
plot(hh2d(1,2:end),hh2d(2,2:end),'k-','markersize',12,'linewidth',1.2);hold on;
% plot(hh3d(1,2:end),hh3d(2,2:end),'g-','markersize',12,'linewidth',3);hold on;
set(gca,'fontsize',25)


pause(0.6);

 ss = sprintf('figepboundary%d.eps',count);
print(fig,'-depsc',ss);


end
% len=length(EnEM);
% 
% fig=figure(137);
% hold on;
% t = linspace(0,time,len);
% plot(t(1:end),EnEM(1:end),'ro-','linewidth',0.8);hold on;
% plot(t(1:end),EnEO(1:end),'bx-','linewidth',0.8);hold on;
% % title('modified and original energy curves')
% legend('modified energy curves','original  energy curves ')
% set(gca,'fontsize',25)
% %% 
% box on
%  ss = sprintf('ep1energy%.eps');
% print(fig,'-depsc',ss);
% 
% 

% figure(138);
% hold on;
% t = linspace(0,time,30);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% %plot(t,MMtotal,'b--','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')