clear all; 
close all;
xright=1.0; yright=2; xleft = 0; yleft = 0; nx=128; ny = 256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt = 0.001; max_it=4500;
time = dt*max_it;
ns = max_it/30;
count =0;

ss=sprintf('./dataout/datac1.m'); c1 = load(ss);
ss=sprintf('./dataout/datac2.m'); c2 = load(ss);
ss=sprintf('./dataout/datac3.m'); c3 = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example6NScac3componentbubble/drop1/dataout/auq.m'); auq = load(ss);
ss=sprintf('./dataout/datau.m'); u = load(ss);
ss=sprintf('./dataout/datav.m'); v = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example6NScac3componentbubble/drop1/dataout/Mene.m'); EnEM = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example6NScac3componentbubble/drop1/dataout/Oene.m'); EnEO = load(ss);


for kk =[2,4,6,8,12,14,16,18,22,24,26,28]
    %[2,4,6,8,12,14,16,18,22,24,26,28]
    count =count+1;
    
    
 fig=figure(kk);
     clf;
t = (kk-1)*ns*dt

A = c1((kk-1)*nx+1:kk*nx,:);
A = A';
%mesh(xx,yy,A);hold on
B = c2((kk-1)*nx+1:kk*nx,:);
B =B';
%mesh(xx,yy,B);hold on
C = c3((kk-1)*nx+1:kk*nx,:);
C = C';
% mesh(xx,yy,C);hold on


% U = u((kk-1)*nx+1:kk*nx,:);
% V = v((kk-1)*nx+1:kk*nx,:);
% U = U';
% V = V';
% u_max(kk)=max(max(U));
% 

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



 contourf(xx,yy,A,[0.5 0.5],'facecolor','r','edgecolor','k','linewidth',1.3);hold on;
 contourf(xx,yy,B,[0.5 0.5],'facecolor','g','edgecolor','k','linewidth',1.3);hold on;
 contourf(xx,yy,C,[0.5 0.5],'facecolor','y','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,A,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,B,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
%  contourf(xx,yy,C,[0.5 0.5],'facecolor','None','edgecolor','k','linewidth',1.3);hold on;
 xlabel('x');
 ylabel('y')






%vector
%         kk1=6; kk2 =6; s=0.5;h = quiver(xx(2:kk1:end,2:kk2:end),yy(2:kk1:end,2:kk2:end),...
%              s*U(2:kk1:end,2:kk2:end),s*V(2:kk1:end,2:kk2:end),0,'b');hold on
% 


% massA = 0; massB = 0;massC =0;
% for ii = 1:ny
%     for jj = 1:nx
%         massA = massA + A(ii,jj);
%         massB = massB + B(ii,jj);
%         massC = massC +C(ii,jj);
%     end
% end
% MMA(kk) = massA*h*h;
% MMB(kk) = massB*h*h;
% MMC(kk) = massC*h*h;


  axis image; axis([0 xright 0 yright]); 
 

 view(0,90);
  set(gca,'fontsize',25.0);
  pause(0.2); box on;
  ss = sprintf('fig7b%d.png',count);
print(fig,'-dpng',ss);

end
% figure(137);
% hold on;
% t = linspace(0,time,31);
% plot(t,EnEM,'rx-','linewidth',0.8);hold on;
% plot(t,EnEO,'bx-','linewidth',0.8);hold on;
% title('modified and original energy curves')
% legend('modified energy curves','original energy curves ')
% figure(200);
% hold on;
% t = linspace(0,time,31);
% plot(t,auq,'-','linewidth',0.8);hold on;
% xlabel('Time')  ;
% title('the value of auxiliary q');
% box on;


% figure(138);
% hold on;
% t = linspace(0,time,30);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo-','linewidth',0.8);hold on;
% plot(t,MMC,'g*-','linewidth',0.8);hold on;
% %plot(t,MMtotal,'b--','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')
