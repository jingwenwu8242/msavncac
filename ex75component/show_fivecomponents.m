clear all;
close all;

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=128; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt =0.2*h;
max_it = 900;
time = dt*max_it;
ns = max_it/30;


ss=sprintf('./dateout/datac.m'); phi = load(ss);
ss=sprintf('./dateout/datac2.m'); phi2 = load(ss);
ss=sprintf('./dateout/datac3.m'); phi3 = load(ss);
ss=sprintf('./dateout/datac4.m'); phi4 = load(ss);
ss=sprintf('./dateout/datac5.m'); phi5 = load(ss);
ss=sprintf('./dateout/Mene.m');EnEM= load(ss);
ss= sprintf('./dateout/Oene.m');EnEO=load(ss);



for i= [10,15,20,25]
    
 figure(i);
t = (i-1)*ns*dt
clf;

A = phi((i-1)*nx+1:i*nx,:);
B = phi2((i-1)*nx+1:i*nx,:);
C= phi3((i-1)*nx+1:i*nx,:);
D= phi4((i-1)*nx+1:i*nx,:);
E = 1-A-B-C-D;


 
% massA = 0; massB = 0;massC =0;
% for ii = 1:nx
%     for jj = 1:ny
%         massA = massA + A(ii,jj);
%         massB = massB + B(ii,jj);
%         massC = massC +C(ii,jj);
%     end
% end
% MMA(i) = massA*h*h;
% MMB(i) = massB*h*h;
% MMC(i) = massC*h*h;

% for i = 1:nx
%     for j = 1:ny
%         
%         if( A(i,j) > 1.0)
%             A(i,j) = 1.0;
%         elseif( A(i,j) < -1.0)
%             A(i,j) = -1.0;
%         else
%             A(i,j) = A(i,j);
%         end
%     end
% end

hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','k','edgecolor','k'); hold on;
hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','b','edgecolor','k'); hold on;
hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','r','edgecolor','k'); hold on;
hh4=contourf(xx,yy,D',[0.5 0.5],'facecolor','g','edgecolor','k'); hold on;
hh5=contourf(xx,yy,E',[0.5 0.5],'facecolor','y','edgecolor','k'); hold on;
xlabel('x');
ylabel('y');



% view(0,90);
% surf(xx,yy,0.5*A'+B');
% shading interp;
% colormap jet;


axis image; axis([0 xright 0 yright]); 
set(gca,'fontsize',25)
box on

pause(0.2);

end

len=length(EnEM);
fig=figure(137);
clf
hold on;
t = linspace(0,time,31);
plot(t,EnEM,'ro-','markersize',10,'linewidth',1);hold on;
plot(t,EnEO,'bx-','markersize',10,'linewidth',1);hold on;
% title('modified  energy curves')
xlabel('Time');
ylabel('Energy');
box on;
xlim([0,time]);
legend('modified energy curves','original energy curves ')
set(gca,'fontsize',25);
ss = sprintf('fig6energy.eps');
print(fig,'-depsc',ss);
% figure(138);
% hold on;
% t = linspace(1,time,30);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')
