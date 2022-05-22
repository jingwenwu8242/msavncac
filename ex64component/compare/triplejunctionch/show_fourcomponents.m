clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=256; ny = 256;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
dt =h;max_it=1000;
time = max_it*dt;
ns=max_it/50;

ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/datac.m'); phi = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/datac2.m'); phi2 = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/datac3.m'); phi3 = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/datac4.m'); phi4 = load(ss);
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/Mene.m');EnEM= load(ss);
% ss= sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example44component/triplejunction/dataout/Oene.m');EnEO=load(ss);




for i = [1,3,5,7]
    
 figure(i);
 clf;
t = (i-1)*ns*dt
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

hh=contourf(xx,yy,A',[0.5 0.5],'facecolor','r','edgecolor','k'); hold on;
hh2=contourf(xx,yy,B',[0.5 0.5],'facecolor','g','edgecolor','k'); hold on;
hh3=contourf(xx,yy,C',[0.5 0.5],'facecolor','b','edgecolor','k'); hold on;
hh4=contourf(xx,yy,D',[0.5 0.5],'facecolor','none','edgecolor','k'); hold on;
xlabel('x');
ylabel('y');
text('Interpreter','latex','String','$c_1$','Position',[0,0],'FontSize',16);




% view(0,90);
% surf(xx,yy,0.5*A'+B');
% shading interp;
% colormap jet;
axis image; axis([0 xright 0 yright]); 
set(gca,'fontsize',15);
box on;

pause(0.6);

end


figure(137);
hold on;
t = linspace(0,time,51);
plot(t(1:2:51),EnEM(1:2:51),'ro-','linewidth',0.8);hold on;
% plot(t(1:end),EnEO(1:end),'bx-','linewidth',0.8);hold on;
%title('odified  energy curves');
xlabel('Time');
ylabel('Modified energy');
set(gca,'fontsize',10);
box on;
% legend('modified energy curves','original  energy curves ')