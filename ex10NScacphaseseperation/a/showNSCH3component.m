clear all; 
close all;

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=128; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt = 0.001; max_it=1000;
time = dt*max_it;
ns = max_it/50;

ss=sprintf('./dataout/datac1.m'); c1 = load(ss);
ss=sprintf('./dataout/datac2.m'); c2 = load(ss);
ss=sprintf('./dataout/datac3.m'); c3 = load(ss);
ss=sprintf('./dataout/auq.m'); auq = load(ss);
ss=sprintf('./dataout/datau.m'); u = load(ss);
ss=sprintf('./dataout/datav.m'); v = load(ss);
ss=sprintf('./dataout/Mene.m'); EnEM = load(ss);
ss=sprintf('./dataout/Oene.m'); EnEO = load(ss);


for kk=1:50
    %kk = [10,20,30,40,50]
    %

%     
% figure(kk);
%      clf;
% t = (kk-1)*ns*dt;

A = c1((kk-1)*nx+1:kk*nx,:);
A = A';
%mesh(xx,yy,A);hold on
B = c2((kk-1)*nx+1:kk*nx,:);
B =B';
%mesh(xx,yy,B);hold on
C = c3((kk-1)*nx+1:kk*nx,:);
C = C';
% mesh(xx,yy,C);hold on

% 
% U = u((kk-1)*nx+1:kk*nx,:);
% V = v((kk-1)*nx+1:kk*nx,:);
% U = U';
% V = V';
% u_max(kk)=max(max(U));


% P = p((kk-1)*nx+1:kk*nx,:);
% P= P';

% surf(xx,yy,A);
% shading interp;
% colormap gray
% colormap jet
% 
% surf(xx,yy,B);
% shading interp;
% colormap gray
% colormap jet
% 
% surf(xx,yy,C);
% shading interp;
% colormap gray
% colormap jet


% 
% contourf(xx,yy,A,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,B,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,C,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
 
%  contourf(xx,yy,A,[0.5 0.5],'facecolor','r','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,B,[0.5 0.5],'facecolor','g','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,C,[0.5 0.5],'facecolor','y','edgecolor','k','linewidth',1.3);hold on;
%  xlabel('x');
%  ylabel('y')






% vector
%         kk1=4; kk2 =4; h = quiver(xx(2:kk1:end,2:kk2:end),yy(2:kk1:end,2:kk2:end),...
%              0.1*U(2:kk1:end,2:kk2:end),0.1*V(2:kk1:end,2:kk2:end),0,'b');hold on
% 

% 
massA = 0; massB = 0;massC =0;
for ii = 1:ny
    for jj = 1:nx
        massA = massA + A(ii,jj);
        massB = massB + B(ii,jj);
        massC = massC +C(ii,jj);
    end
end
MMA(kk) = massA*h*h;
MMB(kk) = massB*h*h;
MMC(kk) = massC*h*h;


%   axis image; axis([0 xright 0 yright]); 
%  set(gca,'fontsize',25)
% box on
% 
% 
%  view(0,90);


end
fig=figure(137);
clf
hold on;
len=length(EnEM);
t = linspace(0,time,len);
plot(t,EnEM,'ro-','markersize',10,'linewidth',1);hold on;
plot(t,EnEO,'bx-','markersize',10,'linewidth',1);hold on;
% title('modified and original energy curves')
legend('modified energy curves','original energy curves ');
set(gca,'fontsize',20);
xlabel('Time');
ylabel('Energy');
xlim([0,time]);
box on;ss = sprintf('fig8energy.eps');
print(fig,'-depsc',ss);



fig=figure(200);
clf
hold on;
len=length(auq);
t = linspace(0,time,len);
plot(t,auq,'-','markersize',10,'linewidth',1);hold on;
xlabel('Time')  ;
xlim([0,time]);
set(gca,'fontsize',23);
box on;
ss = sprintf('fig8q.eps');
print(fig,'-depsc',ss);

% 
fig=figure(138);
clf
hold on;
len=length(MMA);
t = linspace(0,time,len);
plot(t,MMA,'rx-','markersize',10,'linewidth',1);hold on;
plot(t,MMB,'bo-','markersize',10,'linewidth',1);hold on;
plot(t,MMC,'g*-','markersize',10,'linewidth',1);hold on;
%plot(t,MMtotal,'b--','markersize',10,'linewidth',1);hold on;
% title('mass of A,B,C')
legend('MMA','MMB ','MMC')
xlabel('Time') ;
ylabel('Mass');
set(gca,'fontsize',20);
box on;
xlim([0,time]);
ss = sprintf('fig8mass.eps');
print(fig,'-depsc',ss);
