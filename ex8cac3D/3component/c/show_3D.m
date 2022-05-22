clear all;
close all;
nx=64; ny=64; nz=64;
xright = 1.0*pi; h = xright/nx;
yright = 1.0*pi; 
zright = 1.0*pi;


x=linspace(0.5*h,xright-0.5*h,nx);
y=linspace(0.5*h,yright-0.5*h,ny);
z=linspace(0.5*h,zright-0.5*h,nz);

[yy,xx,zz]=meshgrid(y,x,z);

dt = 0.1;

max_it =900;
ns =max_it/30;
time= dt*max_it;
count =0;


for kk = [5,10,15,25]
    
t= (kk-1)*dt*ns;
    count=count+1;
ss=sprintf('./dataout/datac1%d.m',kk); phi1=load(ss);    
ss=sprintf('./dataout/datac2%d.m',kk); phi2=load(ss);    
ss=sprintf('./dataout/datac3%d.m',kk); phi3=load(ss);  
ss=sprintf('./dataout/Mene.m'); EnEM1= load(ss);
 ss=sprintf('./dataout/Oene.m'); EnEO1 = load(ss);
    
fig=figure(kk);
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
massA = 0; massB = 0;massC=0;



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
% MMA(kk) = massA*h*h*h;
% MMB(kk) = massB*h*h*h;
% MMC(kk) = massC*h*h*h;



% %because of ny=nz=nx,we use two iteration.in genergal,we use three
% %iteration for i,j,k respectviely.
d=0;

    for j=1:ny
        for i= 1:nx
            if A(i,j,nz) > 0.5 
               A(i,j,nz)=d;
            end
            
            if A(i,j,1) > 0.5
               A(i,j,1)=d;
            end
            
            if  A(nx,i,j) > 0.5 
               A(nx, i,j)=d;
            end
            
            if A(1,i,j) > 0.5 
               A(1,i,j)=d;
            end
            
            if A(i,ny,j) > 0.5 
               A(i,ny,j)=d;
            end
            
            if A(i,1,j) > 0.5 
               A(i,1,j)=d;
            end
            
         end
    end
    
    
   d=0;
    for j=1:ny
        for i= 1:nx
            if B(i,j,nz) > 0.5 
               B(i,j,nz)=d;
            end
            
            if B(i,j,1) > 0.5
               B(i,j,1)=d;
            end
            
            
           if B(nx,i,j) > 0.5 
               B(nx, i,j)=d;
            end
            
            if B(1,i,j) > 0.5 
               B(1,i,j)=d;
            end
            
            if B(i,ny,j) > 0.5 
               B(i,ny,j)=d;
            end
            
            if B(i,1,j) > 0.5 
               B(i,1,j)=d;
            end
            
         end
    end 
d=0;

    for j=1:ny
        for i= 1:nx
            if C(i,j,nz) > 0.5 
               C(i,j,nz)=d;
            end
            
            if C(i,j,1) > 0.5
               C(i,j,1)=d;
            end
            
            if C(nx,i,j) > 0.5              
               C(nx, i,j)=d;
            end
            
            if C(1,i,j) > 0.5 
               C(1,i,j)=d;
            end
            
            if C(i,ny,j) > 0.5 
               C(i,ny,j)=d;
            end
            
            if C(i,1,j) > 0.5 
               C(i,1,j)=d;
            end
            
         end
    end
     p3=patch(isosurface(yy,xx,zz,C,0.5)); 
     set(p3,'FaceColor','r','EdgeColor','none'); daspect([1 1 1]);hold on;
      alpha(p3,1);   
      view(52,37); 

      p2=patch(isosurface(yy,xx,zz,B,0.5)); 
      set(p2,'FaceColor','y','EdgeColor','none'); daspect([1 1 1]);hold on;
       alpha(p2,1);   
       view(52,37); 
      
        p=patch(isosurface(yy,xx,zz,A,0.5)); 
      set(p,'FaceColor','g','EdgeColor','none'); daspect([1 1 1]);hold on;    
      alpha(p,1);    view(52,37);  
      
      xlabel('x');
      ylabel('y');
      zlabel('z');
      set(gca,'fontsize',25)
      
      
      
      
      
   
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

xlim([0,pi]);
ylim([0,pi]);
zlim([0,pi]);

axis image;
axis([0 xright 0 yright 0 zright]);box on;
  ss = sprintf('33Dc%d.eps',count);
print(fig,'-depsc',ss);

% view(154,27);
% view(-78,24);
% view(-70,15);
pause(0.06);
end 
% figure(137);
% hold on;
% t = linspace(0,time,31);
% plot(t,EnEM1,'ro-.','markersize',7,'linewidth',0.8);hold on;
% plot(t,EnEO1,'bo-','markersize',7,'linewidth',0.8);hold on;
% title('modified and original energy curves')
% legend('modified energy curves','original  energy curves ')
% set(gca,'fontsize',25);
% xlabel('Time');
% ylabel('Energy');
% box on;
% figure(138);
% hold on;
% maxM=max(MMA)
% minM=min(MMA)
% t = linspace(0,time,30);
% plot(t,MMA,'rx-','linewidth',0.8);hold on;
% plot(t,MMB,'bo','linewidth',0.8);hold on;
% plot(t,MMC,'g*','linewidth',0.8);hold on;
% title('mass of A,B,C')
% legend('MMA','MMB ','MMC')
% set(gca,'YLim',[8 12]);
