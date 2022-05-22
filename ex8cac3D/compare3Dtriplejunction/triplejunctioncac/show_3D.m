clear all;
nx=128; ny=128; nz=128;
xright = 1.0; h = xright/nx;
yright = 1.0; 
zright = 1.0;


x=linspace(0.5*h,xright-0.5*h,nx);
y=linspace(0.5*h,yright-0.5*h,ny);
z=linspace(0.5*h,zright-0.5*h,nz);

[yy,xx,zz]=meshgrid(y,x,z);

dt = 0.04*h;

max_it =400;
ns =max_it/40;
time= dt*max_it;


for kk = [10,20,30,40]
    
t= (kk-1)*dt*ns;
    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/datac1%d.m',kk); phi1=load(ss);    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/datac2%d.m',kk); phi2=load(ss);    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/datac3%d.m',kk); phi3=load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/datac4%d.m',kk); phi4=load(ss);  
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/Mene.m'); EnEM1= load(ss);
% ss=sprintf('/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/compare3Dtriplejunction/triplejunctioncac/dataout/Oene.m'); EnEO1 = load(ss);
    
figure(kk);
clf;
for k=1:nz
    for j=1:ny
        for i=1:nx
            A(i,j,k)=phi1(nx*ny*(k-1)+nx*(j-1)+i);
            B(i,j,k)=phi2(nx*ny*(k-1)+nx*(j-1)+i);
            C(i,j,k)=phi3(nx*ny*(k-1)+nx*(j-1)+i);
        end
    end
end

% massA = 0; massB = 0;massC=0;
% 
% 
% 
% for ii = 1:nx
%     for jj = 1:ny
%          for ll = 1:nz
%         massA = massA + A(ii,jj,ll);
%         massB = massB + B(ii,jj,ll);
%         massC = massC +C(ii,jj,ll);
%          end
%        
%     end
% end
% MMA(kk) = massA*h*h;
% MMB(kk) = massB*h*h;
% MMC(kk) = massC*h*h;


      
     p3=patch(isosurface(yy,xx,zz,C,0.5)); 
     set(p3,'FaceColor','g','EdgeColor','none'); daspect([1 1 1]);hold on;
      alpha(p3,1);   
      view(52,37); 

      p2=patch(isosurface(yy,xx,zz,B,0.5)); 
      set(p2,'FaceColor','y','EdgeColor','none'); daspect([1 1 1]);hold on;
       alpha(p2,1);   
       view(52,37); 
      
       p=patch(isosurface(yy,xx,zz,A,0.5)); 
      set(p,'FaceColor','b','EdgeColor','none'); daspect([1 1 1]);hold on;    
      alpha(p,1);   
      view(52,37);  
      
      xlabel('x');
      ylabel('y');
      zlabel('z');
      set(gca,'fontsize',15)
      
      
      
      
      
   
      camlight; lighting phong; 







axis on;
line([0 0],[0 0],[0 zright],'color','k');
line([xright xright],[0 0],[0 zright],'color','k');
line([0 0],[xright xright],[0 zright],'color','k');
line([xright xright],[xright xright],[0 zright],'color','k');
line([0 0],[0 xright],[0 0],'color','k');
line([0 xright],[0 0],[0 0],'color','k');
line([xright xright],[0 xright],[0 0],'color','k');
line([0 xright],[xright xright],[0 0],'color','k');
line([0 0],[0 xright],[zright zright],'color','k');
line([0 xright],[0 0],[zright zright],'color','k');
line([xright xright],[0 xright],[zright zright],'color','k');
line([0 xright],[xright xright],[zright zright],'color','k');


axis image;
axis([0 xright 0 yright 0 zright]);box on;

% view(154,27);
% view(-78,24);
% view(-70,15);
pause(0.02);
end 
% figure(137);
% hold on;
% t = linspace(0,time,31);
% plot(t,EnEM1,'ro-.','markersize',7,'linewidth',0.8);hold on;
% plot(t,EnEO1,'bo-','markersize',7,'linewidth',0.8);hold on;
% title('modified and original energy curves')
% legend('modified energy curves','original  energy curves ')
% set(gca,'fontsize',10);
% xlabel('Time');
% ylabel('Energy');
% box on;

% figure(138);
% hold on;
% t = linspace(0,time,40);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')
