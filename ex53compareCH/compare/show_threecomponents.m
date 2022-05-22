clear all; 
close all;

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny =256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt = 0.004*h; max_it=1200;
time = dt*max_it;
ns = max_it/40;

ss=sprintf('./dataout/3componentCAC/datac.m'); phi = load(ss);
ss=sprintf('./dataout/3componentCAC/datac2.m'); phi2 = load(ss);
ss=sprintf('./dataout/3componentCAC/datac3.m'); phi3 = load(ss);
% ss=sprintf('./dataout/3componentCAC/Mene.m');EnEM1= load(ss);
% ss= sprintf('./dataout/3componentCAC/Oene.m');EnEO1=load(ss);

ss=sprintf('./dataout/3componentCH/datac.m'); phich = load(ss);
ss=sprintf('./dataout/3componentCH/datac2.m'); phi2ch = load(ss);
ss=sprintf('./dataout/3componentCH/datac3.m'); phi3ch = load(ss);
% ss=sprintf('./dataout/3componentCH/Mene.m');EnEM2= load(ss);
% ss= sprintf('./dataout/3componentCH/Oene.m');EnEO2=load(ss);

for i= [1,5,10,15,20]
    
figure1=figure(i);
 t= (i-1)*ns*dt
 clf;
A = phi((i-1)*nx+1:i*nx,:);
B = phi2((i-1)*nx+1:i*nx,:);
C = 1-A-B;
Ach = phich((i-1)*nx+1:i*nx,:);
Bch = phi2ch((i-1)*nx+1:i*nx,:);
Cch = 1-Ach-Bch;

% % for i = 1:nx
% %     for j = 1:ny
% %         
% %         if( A(i,j) > 1.0)
% %             A(i,j) = 1.0;
% %         elseif( A(i,j) < -1.0)
% %             A(i,j) = -1.0;
% %         else
% %             A(i,j) = A(i,j);
% %         end
% %     end
% % end
% 
% hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','none','edgecolor','r'); hold on;
% hh=contourf(xx,yy,B',[0.5 0.5],'facecolor','none','edgecolor','r'); hold on;
% hh=contourf(xx,yy,C',[0.5 0.5],'facecolor','none','edgecolor','r'); hold on;

%mark contourf
hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','none');
hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','none');
hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','none');

hhs=contourf(xx,yy,Ach',[0.5 0.5],'facecolor','none');
hh2s=contourf(xx,yy,Bch',[0.5 0.5],'facecolor','none');
hh3s=contourf(xx,yy,Cch',[0.5 0.5],'facecolor','none');


plot(hh(1,1:10:end),hh(2,1:10:end),'bo','markersize',8,'linewidth',1);hold on;
plot(hh2(1,1:10:end),hh2(2,1:10:end),'bo','markersize',8,'linewidth',1);hold on;
plot(hh3(1,1:10:end),hh3(2,1:10:end),'bo','markersize',8,'linewidth',1);hold on;

plot(hhs(1,1:10:end),hhs(2,1:10:end),'r+','markersize',8,'linewidth',1);hold on;
plot(hh2s(1,1:10:end),hh2s(2,1:10:end),'r+','markersize',8,'linewidth',1);hold on;
plot(hh3s(1,1:10:end),hh3s(2,1:10:end),'r+','markersize',8,'linewidth',1);hold on;
xlim([0,1]);ylim([0,1]);
xlabel('x');
ylabel('y');

annotation(figure1,'textbox',...
    [0.380357142857143 0.703761904761905 0.0575892857142857 0.0547619047619048],...
    'String',{'$\Omega_{1}$'}, 'FitBoxToText','off','FontSize',25,'LineStyle','none','Interpreter','latex');
annotation(figure1,'textbox',...
    [0.486714285714285 0.526190476190478 0.0575892857142857 0.0547619047619048],...
    'String',{'$\Omega_{2}$'}, 'FitBoxToText','off','FontSize',25,'LineStyle','none','Interpreter','latex');
annotation(figure1,'textbox',...
    [0.532142857142857 0.279952380952381 0.069642857142857 0.0547619047619048],...
    'String',{'$\Omega_{3}$'},...
    'FitBoxToText','off','FontSize',25,'LineStyle','none','Interpreter','latex');


set(gca,'fontsize',25)

% 
% % 
% % % view(0,90);
% % % surf(xx,yy,0.5*A'+B');
% % % shading interp;
% % % colormap jet;
% 
% % massA = 0; massB = 0;massC =0;
% % for ii = 1:nx
% %     for jj = 1:ny
% %         massA = massA + A(ii,jj);
% %         massB = massB + B(ii,jj);
% %         massC = massC +C(ii,jj);
% %     end
% % end
% % MMA(i) = massA*h*h;
% % MMB(i) = massB*h*h;
% % MMC(i) = massC*h*h;
% 
% 
axis image; axis([0 xright 0 yright]); 
box on

pause(0.1);

end
% 
% 
% figure(137);
% hold on;
% t = linspace(0,time,41);
% plot(t,EnEM1,'r*-.','markersize',10,'linewidth',0.8);hold on;
% plot(t,EnEO1,'bo-','markersize',10,'linewidth',0.8);hold on;
% plot(t,EnEM2,'r*-.','markersize',10,'linewidth',0.8);hold on;
% plot(t,EnEO2,'bx-','markersize',10,'linewidth',0.8);hold on;
% %title('modified and original energy curves')
% xlabel('Time');
% ylabel('Energy');
% %  legend('Modified energy curves for CAC','Original energy curves for CAC');
% legend('Modified energy curves for CAC','Original energy curves for CAC','Modified energy curves for CH','Original energy curves for CH ')
% box on;
% set(gca,'fontsize',10);
% 
% figure(138);
% hold on;
% t = linspace(1,16,16);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% plot(t,MMtotal,'b--','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC','MMtotal')





