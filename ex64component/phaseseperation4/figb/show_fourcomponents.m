clear all; 

xright=4.0; yright=4.0; xleft = 0; yleft = 0; nx=128; ny = 128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
max_it =300;
dt =h;
time = dt*max_it;
ns = max_it/30;


ss=sprintf('./dataout/datac.m'); phi = load(ss);
ss=sprintf('./dataout/datac2.m'); phi2 = load(ss);
ss=sprintf('./dataout/datac3.m'); phi3 = load(ss);
ss=sprintf('./dataout/datac4.m'); phi4 = load(ss);
ss=sprintf('./dataout/Mene.m');EnEM= load(ss);
ss= sprintf('./dataout/Oene.m');EnEO=load(ss);
count=0;



for i = [15,20,25,30]%i=1:30
    count=count+1;
fig=figure(i);
 clf;
 
 t = (i-1)*ns*dt;

A = phi((i-1)*nx+1:i*nx,:);
B = phi2((i-1)*nx+1:i*nx,:);
C= phi3((i-1)*nx+1:i*nx,:);
D = 1-A-B-C;

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


% massA = 0; massB = 0;massC =0;
% for ii = 1:ny
%     for jj = 1:nx
%         massA = massA + A(ii,jj);
%         massB = massB + B(ii,jj);
%         massC = massC +C(ii,jj);
%     end
% end
% MMA(i) = massA*h*h;
% MMB(i) = massB*h*h;
% MMC(i) = massC*h*h;


hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','y','edgecolor','k'); hold on;
hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','g','edgecolor','k'); hold on;
hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','b','edgecolor','k'); hold on;
hh4=contourf(xx,yy,D',[0.5 0.5],'facecolor','r','edgecolor','k'); hold on;
xlabel('x');
ylabel('y');
% view(0,90);
% surf(xx,yy,0.5*A'+B');
% shading interp;
% colormap jet;
axis image; axis([0 xright 0 yright]); 
set(gca,'fontsize',25)
box on
ss = sprintf('fig4b%d.eps',count);
print(fig,'-depsc',ss);
pause(0.6);

end

% 
% figure(137);
% hold on;
% t = linspace(0,time,31);
% plot(t(1:end),EnEM(1:end),'ro-','linewidth',0.8);hold on;
% plot(t(1:end),EnEO(1:end),'bx-','linewidth',0.8);hold on;
% title('modified and original energy curves')
% legend('modified energy curves','original  energy curves ')
% set(gca,'fontsize',10)
% box on



% figure(138);
% hold on;
% t = linspace(0,time,30);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% %plot(t,MMtotal,'b--','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')