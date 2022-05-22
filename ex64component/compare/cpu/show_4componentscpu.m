clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny =256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
dt =0.04*h;max_it=400;
time = max_it*dt;
[xx,yy]=meshgrid(x,y);


% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentcompare/datac.m'); phi = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentcompare/datac2.m'); phi2 = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentcompare/datac3.m'); phi3 = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentcompare/Mene.m');EnEM= load(ss);
% ss= sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentcompare/Oene.m');EnEO=load(ss);
ss=sprintf('./dataoutput/cpuch.m'); cpu1 = load(ss);

% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentSAV/datac.m'); phisav = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentSAV/datac2.m'); phi2sav = load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/data_output/3componentSAV/datac3.m'); phi3sav = load(ss);
ss=sprintf('./dataoutput/cpucac.m'); cpu2 = load(ss);
len=length(cpu2)
fig=figure(1);
hold on;
t = linspace(0,time,len);
plot(t,cpu1,'ro-','linewidth',0.8);hold on;
plot(t,cpu2,'bx-','linewidth',0.8);hold on;
%title('Cpu time')
h=legend('CH model','CAC model ');
xlabel('Time');
ylabel('CPU time');
box on;
set(gca,'fontsize',20);
ss = sprintf('fig5cpu.eps');
print(fig,'-depsc',ss);


difference =abs((cpu1-cpu2)./cpu2)*100


% for i= [1,8,12,20]
%     
%  figure(i);
% % clf;
% 
% A = phi((i-1)*nx+1:i*nx,:);
% B = phi2((i-1)*nx+1:i*nx,:);
% C = 1-A-B;
% As = phisav((i-1)*nx+1:i*nx,:);
% Bs = phi2sav((i-1)*nx+1:i*nx,:);
% Cs = 1-As-Bs;
% 
% hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','none');
% hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','none');
% hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','none');
% hhs=contourf(xx,yy,As',[0.5 0.5],'facecolor','none');
% hh2s=contourf(xx,yy,Bs',[0.5 0.5],'facecolor','none');
% hh3s=contourf(xx,yy,Cs',[0.5 0.5],'facecolor','none');
% 
% 
% plot(hh(1,1:10:end),hh(2,1:10:end),'yo');hold on;
% plot(hh2(1,1:10:end),hh2(2,1:10:end),'yo');hold on;
% plot(hh3(1,1:10:end),hh3(2,1:10:end),'yo');hold on;
% plot(hhs(1,1:10:end),hhs(2,1:10:end),'r+');hold on;
% plot(hh2s(1,1:10:end),hh2s(2,1:10:end),'r+');hold on;
% plot(hh3s(1,1:10:end),hh3s(2,1:10:end),'r+');hold on;
% 
% % 
% % % 
% % % % view(0,90);
% % % % surf(xx,yy,0.5*A'+B');
% % % % shading interp;
% % % % colormap jet;
% % 
% % % massA = 0; massB = 0;massC =0;
% % % for ii = 1:nx
% % %     for jj = 1:ny
% % %         massA = massA + A(ii,jj);
% % %         massB = massB + B(ii,jj);
% % %         massC = massC +C(ii,jj);
% % %     end
% % % end
% % % MMA(i) = massA*h*h;
% % % MMB(i) = massB*h*h;
% % % MMC(i) = massC*h*h;
% % 
% % 
% axis image; axis([0 xright 0 yright]); 
% box on
% 
% pause(0.1);
% 
% end
% 
% 
% figure(137);
% hold on;
% t = linspace(0,10,61);
% plot(t,EnEM,'ro-','linewidth',0.8);hold on;
% plot(t,EnEO,'bx-','linewidth',0.8);hold on;
% title('modified and original energy curves')
% legend('modified energy curves','original energy curves ')

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





