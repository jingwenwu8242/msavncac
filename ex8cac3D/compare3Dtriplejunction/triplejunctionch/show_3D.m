clear all;
nx=64; ny=64; nz=64;
xright = 1.0*pi; h = xright/nx;
yright = 1.0*pi; zright = 1.0*pi;


x=linspace(0.5*h,xright-0.5*h,nx);
y=linspace(0.5*h,yright-0.5*h,ny);
z=linspace(0.5*h,zright-0.5*h,nz);

[yy,xx,zz]=meshgrid(y,x,z);

dt = 0.1;

max_it =500;
ns =max_it/20;
time= dt*max_it;


for kk = [5,10,12,15]
    
t= (kk-1)*dt*ns
    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example93D/a/dataout/datac1%d.m',kk); phi1=load(ss);    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example93D/a/dataout/datac2%d.m',kk); phi2=load(ss);    
ss=sprintf('/Users/wujingwen/research/fluid_mechanics/SAV+multi-component incompressible fluid system/msav+nc/example93D/a/dataout/datac3%d.m',kk); phi3=load(ss);    
    
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

%because of ny=nz=nx,we use two iteration.in genergal,we use three
%iteration for i,j,k respectviely.
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
pause(0.06);
end 

