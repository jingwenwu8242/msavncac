clear all; 

xright=1.0; yright=1.0; xleft = 0; yleft = 0; nx=128*2; ny = 128*2;h= (xright-xleft)/nx;
x=linspace(xleft+0.5*h,xright-0.5*h,nx); y=linspace(yleft+0.5*h,yright-0.5*h,ny);
[xx,yy]=meshgrid(x,y);
time =4;

 
ss=sprintf('./stability001/outdata/datac.m'); c11 = load(ss);
ss=sprintf('./stability001/outdata/datac2.m'); c21 = load(ss);
ss=sprintf('./stability001/outdata/datac3.m'); c31 = load(ss);
ss=sprintf('./stability0005/outdata/datac.m'); c12 = load(ss);
ss=sprintf('./stability0005/outdata/datac2.m'); c22 = load(ss);
ss=sprintf('./stability0005/outdata/datac3.m'); c32 = load(ss);
ss=sprintf('./stability00025/outdata/datac.m'); c13 = load(ss);
ss=sprintf('./stability00025/outdata/datac2.m'); c23 = load(ss);
ss=sprintf('./stability00025/outdata/datac3.m'); c33 = load(ss);
ss=sprintf('./stability0001/outdata/datac.m'); c14 = load(ss);
ss=sprintf('./stability0001/outdata/datac2.m'); c24 = load(ss);
ss=sprintf('./stability0001/outdata/datac3.m'); c34 = load(ss);




for kk = 1:20
    
    
A1 = c11((kk-1)*nx+1:kk*nx,:);
A1 = A1';
B1 = c21((kk-1)*nx+1:kk*nx,:);
B1 =B1';
C1 = c31((kk-1)*nx+1:kk*nx,:);
C1 = C1';


A2 = c12((kk-1)*nx+1:kk*nx,:);
A2 = A2';
B2 = c22((kk-1)*nx+1:kk*nx,:);
B2 =B2';
C2 = c32((kk-1)*nx+1:kk*nx,:);
C2 = C2';

A3 = c13((kk-1)*nx+1:kk*nx,:);
A3 = A3';
B3 = c23((kk-1)*nx+1:kk*nx,:);
B3 =B3';
C3 = c33((kk-1)*nx+1:kk*nx,:);
C3 = C3';

A4 = c14((kk-1)*nx+1:kk*nx,:);
A4 = A4';
B4 = c24((kk-1)*nx+1:kk*nx,:);
B4 =B4';
C4 = c34((kk-1)*nx+1:kk*nx,:);
C4 = C4';



massA1 = 0; massB1 = 0;massC1 =0;
massA2 = 0; massB2 = 0;massC2 =0;
massA3 = 0; massB3 = 0;massC3 =0;
massA4 = 0; massB4 = 0;massC4 =0;


for ii = 1:ny
    for jj = 1:nx
        massA1 = massA1 + A1(ii,jj);
        massB1 = massB1 + B1(ii,jj);
        massC1 = massC1 +C1(ii,jj);
        
        massA2 = massA2 + A2(ii,jj);
        massB2 = massB2 + B2(ii,jj);
        massC2 = massC2 +C2(ii,jj);
        
        massA3 = massA3 + A3(ii,jj);
        massB3 = massB3 + B3(ii,jj);
        massC3 = massC3 +C3(ii,jj);
        
        massA4 = massA4 + A4(ii,jj);
        massB4 = massB4 + B4(ii,jj);
        massC4 = massC4 +C4(ii,jj);
    end
end
MMA1(kk) = massA1*h*h;
MMB1(kk) = massB1*h*h;
MMC1(kk) = massC1*h*h;

MMA2(kk) = massA2*h*h;
MMB2(kk) = massB2*h*h;
MMC2(kk) = massC2*h*h;

MMA3(kk) = massA3*h*h;
MMB3(kk) = massB3*h*h;
MMC3(kk) = massC3*h*h;

MMA4(kk) = massA4*h*h;
MMB4(kk) = massB4*h*h;
MMC4(kk) = massC4*h*h;





end


fig=figure(137);
hold on;
t = linspace(0,4,20);
plot(t,MMA1,'rx--','markersize',15,'linewidth',1);hold on;
plot(t,MMB1,'rx--','markersize',15,'linewidth',1);hold on;
plot(t,MMC1,'rx--','markersize',15,'linewidth',1);hold on;
plot(t,MMA2,'go-.','markersize',15,'linewidth',1);hold on;
plot(t,MMB2,'go-.','markersize',15,'linewidth',1);hold on;
plot(t,MMC2,'go-.','markersize',15,'linewidth',1);hold on;
plot(t,MMA3,'bh:','markersize',15,'linewidth',1);hold on;
plot(t,MMB3,'bh:','markersize',15,'linewidth',1);hold on;
plot(t,MMC3,'bh:','markersize',15,'linewidth',1);hold on;
plot(t,MMA4,'k<-','markersize',15,'linewidth',1);hold on;
plot(t,MMB4,'k<-','markersize',15,'linewidth',1);hold on;
plot(t,MMC4,'k<-','markersize',15,'linewidth',1);hold on;
%title('Modified and Original energy curves for S =2 and four different dt')
h=legend('The mass of c_{1}: \Delta t =0.01','The mass of c_{2}: \Delta t =0.01','The mass of c_{3}: \Delta t =0.01','The mass of c_{1}: \Delta t =0.005','The mass of c_{2}: \Delta t =0.005','The mass of c_{3}: \Delta t =0.005','The mass of c_{1}: \Delta t =0.0025','The mass of c_{2}: \Delta t =0.0025','The mass of c_{3}: \Delta t =0.0025','The mass of c_{1}: \Delta t =0.001','The mass of c_{2}: \Delta t =0.001','The mass of c_{3}: \Delta t =0.001');
xlabel('Time');
ylabel('Mass');
set(gca,'fontsize',25);
box on;

ss = sprintf('fig2mass.eps');
print(fig,'-depsc',ss);




