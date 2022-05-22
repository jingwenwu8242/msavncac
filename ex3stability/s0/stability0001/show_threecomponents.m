clear all; 

xright=2.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny =128;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt = 0.01; max_it=300;
time = dt*max_it;
ns = max_it/30;

ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/cacexample3compareSAVs/compare/dataout/3componentcompare/datac.m'); phi = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/cacexample3compareSAVs/compare/dataout/3componentcompare/datac2.m'); phi2 = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/cacexample3compareSAVs/compare/dataout/3componentcompare/datac3.m'); phi3 = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/cacexample3compareSAVs/compare/dataout/3componentcompare/Mene.m');EnEM1= load(ss);
ss= sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/cacexample3compareSAVs/compare/dataout/3componentcompare/Oene.m');EnEO1=load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example3compareSAVs/compare/datdaout3componentcompare/3componentSAV/datac.m'); phisav = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example3compareSAVs/compare/datdaout3componentcompare/3componentSAV/datac2.m'); phi2sav = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example3compareSAVs/compare/datdaout3componentcompare/3componentSAV/datac3.m'); phi3sav = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example3compareSAVs/compare/datdaout3componentcompare/3componentSAV/Mene.m');EnEM2= load(ss);
% ss= sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example3compareSAVs/compare/datdaout3componentcompare/3componentSAV/Oene.m');EnEO2=load(ss);

for i= 1:30
    
%  figure(i);
%   t= (i-1)*ns*dt
 clf;
A = phi((i-1)*nx+1:i*nx,:);
B = phi2((i-1)*nx+1:i*nx,:);
C = 1-A-B;
% As = phisav((i-1)*nx+1:i*nx,:);
% Bs = phi2sav((i-1)*nx+1:i*nx,:);
% Cs = 1-As-Bs;

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
hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','r','edgecolor','k'); hold on;
hh=contourf(xx,yy,B',[0.5 0.5],'facecolor','y','edgecolor','k'); hold on;
hh=contourf(xx,yy,C',[0.5 0.5],'facecolor','b','edgecolor','k'); hold on;
% hh=mesh(xx,yy,A');hold on;
% hh=mesh(xx,yy,B');hold on;
% hh=mesh(xx,yy,C');hold on;


%mark contourf
% hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','none');
% hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','none');
% hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','none');
% hhs=contourf(xx,yy,As',[0.5 0.5],'facecolor','none');
% hh2s=contourf(xx,yy,Bs',[0.5 0.5],'facecolor','none');
% hh3s=contourf(xx,yy,Cs',[0.5 0.5],'facecolor','none');


% plot(hh(1,1:10:end),hh(2,1:10:end),'k-');hold on;
% plot(hh2(1,1:10:end),hh2(2,1:10:end),'k-');hold on;
% plot(hh3(1,1:10:end),hh3(2,1:10:end),'k-');hold on;
% plot(hhs(1,1:10:end),hhs(2,1:10:end),'r+');hold on;
% plot(hh2s(1,1:10:end),hh2s(2,1:10:end),'r+');hold on;
% plot(hh3s(1,1:10:end),hh3s(2,1:10:end),'r+');hold on;
xlabel('x');
ylabel('y');
% text('Interpreter','latex','String','$c_1$','Position',[0,0],'FontSize',16);
% text('Interpreter','latex','String','$c_2$','Position',[0,0],'FontSize',16);
% text('Interpreter','latex','String','$c_3$','Position',[0,0],'FontSize',16);
set(gca,'fontsize',15)

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

pause(0.5);

end
% 
% 
figure(136);
hold on;
t = linspace(0,time,31);
plot(t,EnEM1,'b*-.','markersize',7,'linewidth',0.8);hold on;
plot(t,EnEO1,'ro-','markersize',7,'linewidth',0.8);hold on;
% plot(t,EnEM2,'rx-.','markersize',10,'linewidth',0.8);hold on;
% plot(t,EnEO2,'bx-','markersize',10,'linewidth',0.8);hold on;
%title('modified and original energy curves')
xlabel('Time');
ylabel('Energy');
% legend('Modified energy curves for mSAV','Original energy curves for mSAV','Modified energy curves for SAV','Original energy curves for SAV ')
legend('Modified energy curves','Original energy curves');
box on;
set(gca,'fontsize',10);
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





