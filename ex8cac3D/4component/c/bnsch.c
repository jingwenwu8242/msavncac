/******************************************************
x:periodic ,y:periodic,z:Neunman

cac
   with variable density and variable viscosity

******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "bnsch.h"
#define NR_END 1
#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

int nx,ny,nz, n_level,  c_relax, count, Num;
float pi, xleft, xright, yleft, yright,zleft,zright,volume,h, dt, gam, Cahn, M, Ev, oEv, acount, int_r, or, r, nrr, EnEO, EnEM,
        ***worku, ***workv, ***workw,***workp, ***mc, SS,CC, ***ct, ***sc, ***smu, ***mu, ***H1, ***H2, ***H3,***H4 ,***alp, ***intc1, ***intc2,
        ***intc3,***intc4, ***c4, ***oc4, ***ic, ***ic2, ***ic3,***ic4;
char bufferc[2000]={0}, bufferc2[2000]={0}, bufferc3[2000]={0},bufferc4[2000]={0};
int main()
{
    extern int nx,ny,nz, n_level,  c_relax, count, Num;
    extern float pi, xleft, xright, yleft, yright,zleft,zright,volume,h, dt, gam, Cahn, M, Ev, oEv, acount, int_r, or, r, nrr, EnEO, EnEM,
        ***worku, ***workv, ***workw,***workp, ***mc, SS,CC, ***ct, ***sc, ***smu, ***mu, ***H1, ***H2, ***H3,***H4 ,***alp, ***intc1, ***intc2,
        ***intc3,***intc4 ,***c4, ***oc4, ***ic, ***ic2, ***ic3,***ic4;
    extern char  bufferc[2000], bufferc2[2000], bufferc3[2000],bufferc4[2000];
    int i,j,k, it, max_it, ns;
    float ***c, ***nc, ***oc, ***c2, ***nc2, ***oc2,***c3, ***nc3, ***oc3;
    
    FILE  *fc, *fc2, *fc3,*fc4, *myr, *myq;
    
    
    /******/
     
    nx = gnx;
    ny = gny;
    nz = gnz;
    n_level = (int)(log(nx)/log(2)-0.9);
    
    c_relax = 5;
    
    pi = 4.0*atan(1.0);
    
    
    xleft = 0.0, xright = 1.0*pi;
    yleft = 0.0, yright = 1.0*pi;
    zleft = 0.0, zright = 1.0*pi;

    volume = (xright-xleft)*(yright-yleft)*(zright-zleft);
    
    count = 0;
    
    /**************************/
    
    h = xright/(float)nx;

    dt = 0.01*h;

    gam = 4.0*h/(4.0*sqrt(2.0)*atanh(0.9));
    Cahn = pow(gam,2);
    M = 1;

    SS = 2;

    CC =0.0;
    
        
    max_it = 600;
    ns = max_it/30;


    Num = 4;  
  

    /**************************/
    
    printf("nx = %d, ny = %d\n,nz=%d\n", nx, ny,nz);
    printf("n_level = %d\n", n_level);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);
    printf("M = %f\n",M);
    
 
    intc1 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    intc2 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    intc3 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    intc4 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    


    ic = cube(0, nx+1, 0, ny+1, 0, nz+1);
    ic2 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    ic3 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    ic4 = cube(0, nx+1, 0, ny+1, 0, nz+1);

    oc = cube(0, nx+1, 0, ny+1, 0, nz+1);
    c = cube(0, nx+1, 0, ny+1, 0, nz+1);
    nc = cube(0, nx+1, 0, ny+1, 0, nz+1);

    oc2 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    c2 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    nc2 = cube(0, nx+1, 0, ny+1, 0, nz+1);

    oc3 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    c3 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    nc3 = cube(0, nx+1, 0, ny+1, 0, nz+1);

    c4 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    oc4 = cube(0, nx+1, 0, ny+1, 0, nz+1);

    H1 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    H2 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    H3 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    H4 = cube(0, nx+1, 0, ny+1, 0, nz+1);
    alp = cube(0, nx+1, 0, ny+1, 0, nz+1);
    
    worku = cube(0, nx, 1, ny, 1, nz);
    workv = cube(1, nx, 0, ny, 1, nz);
    workw = cube(1, nx, 1, ny, 0, nz);
    workp = cube(1, nx, 1, ny, 1, nz);
   
    ct = cube(1, nx, 1, ny, 1, nz);
    sc = cube(1, nx, 1, ny, 1, nz);
    smu = cube(1, nx, 1, ny, 1, nz);
    mu = cube(1, nx, 1, ny, 1, nz);
    
    zero_cube(mu, 1, nx, 1, ny,1,nz); 
 
    
    initialization(c, c2, c3,c4);


    cube_copy(oc, c, 1, nx, 1, ny, 1, nz);
    cube_copy(oc2, c2, 1, nx, 1, ny, 1, nz);
    cube_copy(oc3, c3, 1, nx, 1, ny, 1, nz);

    cube_copy(nc, c, 1, nx, 1, ny, 1, nz);
    cube_copy(nc2, c2, 1, nx, 1, ny, 1, nz);
    cube_copy(nc3, c3, 1, nx, 1, ny, 1, nz);
    
    print_data(c, c2,c3, c4,1);
    
    augmenc(c, nx, ny,nz);
    augmenc(c2, nx, ny,nz);
    augmenc(c3, nx, ny,nz);
    augmenc(c4, nx, ny,nz);
    
    
     EnEO = 0.0;
     
     ijkloop{
     
      EnEO = EnEO + h*h*h*( 1/Cahn*0.25*pow(c[i][j][k],2)*pow(c[i][j][k]-1.0,2) 
         + 1/Cahn*0.25*pow(c2[i][j][k],2)*pow(c2[i][j][k]-1.0,2)
         + 1/Cahn*0.25*pow(c3[i][j][k],2)*pow(c3[i][j][k]-1.0,2)  
         + 1/Cahn*0.25*pow(c4[i][j][k],2)*pow(c4[i][j][k]-1.0,2) 
         + 0.5*( pow(c[i+1][j][k]-c[i][j][k],2)/(h*h) + pow(c[i][j+1][k]-c[i][j][k],2)/(h*h)+ pow(c[i][j][k+1]-c[i][j][k],2)/(h*h) )
         + 0.5*( pow(c2[i+1][j][k]-c2[i][j][k],2)/(h*h) + pow(c2[i][j+1][k]-c2[i][j][k],2)/(h*h) + pow(c2[i][j][k+1]-c2[i][j][k],2)/(h*h))
          + 0.5*( pow(c3[i+1][j][k]-c2[i][j][k],2)/(h*h) + pow(c3[i][j+1][k]-c3[i][j][k],2)/(h*h) + pow(c3[i][j][k+1]-c3[i][j][k],2)/(h*h))
         + 0.5*( pow(c4[i+1][j][k]-c4[i][j][k],2)/(h*h) + pow(c4[i][j+1][k]-c4[i][j][k],2)/(h*h)+ pow(c4[i][j][k+1]-c4[i][j][k],2)/(h*h))    );
         
     }
       
    myr = fopen("/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/Mene.m","w");
    myq = fopen("/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/Oene.m","w");
    
    fprintf(myr,"%16.12f \n",EnEO);
    fprintf(myq,"%16.12f \n",EnEO);

 
      

  for( it = 1; it <= 1; it++){

       Ev = Ef(c, c2, c3,c4, nx, ny,nz);

       r = Ev;
       or = r;
       


    ijkloop{

      H1[i][j][k] = c[i][j][k]*(c[i][j][k]-0.5)*(c[i][j][k]-1.0)/Ev;
      H2[i][j][k] = c2[i][j][k]*(c2[i][j][k]-0.5)*(c2[i][j][k]-1.0)/Ev;
      H3[i][j][k] = c3[i][j][k]*(c3[i][j][k]-0.5)*(c3[i][j][k]-1.0)/Ev;
      H4[i][j][k] = c4[i][j][k]*(c4[i][j][k]-0.5)*(c4[i][j][k]-1.0)/Ev;

      

      }


    lagrangeH(H1,nx,ny,nz);
    lagrangeH(H2,nx,ny,nz);
    lagrangeH(H3,nx,ny,nz);
    lagrangeH(H4,nx,ny,nz);

    ijkloop{

        alp[i][j][k] = -(1.0/Num)*(H1[i][j][k] + H2[i][j][k] + H3[i][j][k]+ H4[i][j][k] );


    }






    
          cahn_ini(c, H1, nc);

    
          cahn_ini(c2, H2, nc2);

          cahn_ini(c3, H3, nc3);

   

     ijkloop{ c4[i][j][k] = 1.0 - nc[i][j][k] - nc2[i][j][k]-nc3[i][j][k]; }


   

     acount = 0.0;

    ijkloop{

        acount = acount + h*h*h*( H1[i][j][k]*r*(nc[i][j][k]-c[i][j][k])/dt
               + H2[i][j][k]*r*(nc2[i][j][k]-c2[i][j][k])/dt 
               + H3[i][j][k]*r*(nc3[i][j][k]-c3[i][j][k])/dt
               + H4[i][j][k]*r*( (1.0-nc[i][j][k]-nc2[i][j][k]-nc3[i][j][k])-(1.0-c[i][j][k]-c2[i][j][k]-c3[i][j][k]) )/dt );

     }


        nrr = r + dt*acount;


         r = nrr;
  
        cube_copy(c, nc, 1, nx, 1, ny,1,nz);
        cube_copy(c2, nc2, 1, nx, 1, ny,1,nz);
        cube_copy(c3, nc3, 1, nx, 1, ny,1,nz);
        

   }



     


    for (it=2; it<=max_it; it++) {

   
     int_r = 2.0*r - or;
     

    ijkloop{

       intc1[i][j][k] = 2.0*c[i][j][k] - oc[i][j][k]; 
       intc2[i][j][k] = 2.0*c2[i][j][k] - oc2[i][j][k];
       intc3[i][j][k] = 2.0*c3[i][j][k] - oc3[i][j][k];
       intc4[i][j][k] = 1.0 - intc1[i][j][k] - intc2[i][j][k]- intc3[i][j][k];

      }
      


       Ev = Ef(intc1, intc2, intc3,intc4, nx,ny,nz);
 

    ijkloop{

      H1[i][j][k] = intc1[i][j][k]*(intc1[i][j][k]-0.5)*(intc1[i][j][k]-1.0)/Ev;
      H2[i][j][k] = intc2[i][j][k]*(intc2[i][j][k]-0.5)*(intc2[i][j][k]-1.0)/Ev;
      H3[i][j][k] = intc3[i][j][k]*(intc3[i][j][k]-0.5)*(intc3[i][j][k]-1.0)/Ev;
      H4[i][j][k] = intc4[i][j][k]*(intc4[i][j][k]-0.5)*(intc4[i][j][k]-1.0)/Ev;

      

      }
      lagrangeH(H1,nx,ny,nz);
      lagrangeH(H2,nx,ny,nz);
      lagrangeH(H3,nx,ny,nz);
      lagrangeH(H4,nx,ny,nz);

      ijkloop{

          alp[i][j][k] = -(1.0/Num)*(H1[i][j][k] + H2[i][j][k] + H3[i][j][k]+ H4[i][j][k] );

      }
  

    
          cahn(c, oc, H1, nc);

    
          cahn(c2, oc2, H2, nc2);
    
          cahn(c3, oc3, H3, nc3);

   

     ijkloop{ c4[i][j][k] = 1.0 - nc[i][j][k] - nc2[i][j][k]- nc3[i][j][k]; }


   

     acount = 0.0;

    ijkloop{

        acount = acount + h*h*h*( H1[i][j][k]*int_r*(3.0*nc[i][j][k]-4.0*c[i][j][k]+oc[i][j][k])/(2.0*dt)
               + H2[i][j][k]*int_r*(3.0*nc2[i][j][k]-4.0*c2[i][j][k]+oc2[i][j][k])/(2.0*dt)
               + H3[i][j][k]*int_r*(3.0*nc3[i][j][k]-4.0*c3[i][j][k]+oc3[i][j][k])/(2.0*dt)
               + H4[i][j][k]*int_r*(3.0*(1.0-nc[i][j][k]-nc2[i][j][k]-nc3[i][j][k])-4.0*(1.0-c[i][j][k]-c2[i][j][k]-c3[i][j][k])+(1.0-oc[i][j][k]-oc2[i][j][k]-oc3[i][j][k]))/(2.0*dt) );

     }


        nrr = ( 4.0*r - or + 2.0*dt*acount )/3.0;
   
         or = r;
         r = nrr; 
         
         augmenc(nc, nx,ny,nz);
         augmenc(nc2, nx,ny,nz);
         augmenc(nc3, nx,ny,nz);
         augmenc(c4, nx,ny,nz);
         
     
     EnEO = 0.0;
     
     ijkloop{
     
      EnEO = EnEO + h*h*h*( 1/Cahn*0.25*pow(nc[i][j][k],2)*pow(nc[i][j][k]-1.0,2) 
        + 1/Cahn*0.25*pow(nc2[i][j][k],2)*pow(nc2[i][j][k]-1.0,2)
         + 1/Cahn*0.25*pow(nc3[i][j][k],2)*pow(nc3[i][j][k]-1.0,2) 
         + 1/Cahn*0.25*pow(c4[i][j][k],2)*pow(c4[i][j][k]-1.0,2) 
         + 0.5*( pow(nc[i+1][j][k]-nc[i][j][k],2)/(h*h) + pow(nc[i][j+1][k]-nc[i][j][k],2)/(h*h) + pow(nc[i][j][k+1]-nc[i][j][k],2)/(h*h))
         + 0.5*( pow(nc2[i+1][j][k]-nc2[i][j][k],2)/(h*h) + pow(nc2[i][j+1][k]-nc2[i][j][k],2)/(h*h) + pow(nc2[i][j][k+1]-nc2[i][j][k],2)/(h*h))
         + 0.5*( pow(nc3[i+1][j][k]-nc3[i][j][k],2)/(h*h) + pow(nc3[i][j+1][k]-nc3[i][j][k],2)/(h*h) + pow(nc3[i][j][k+1]-nc3[i][j][k],2)/(h*h))
         + 0.5*( pow(c4[i+1][j][k]-c4[i][j][k],2)/(h*h) + pow(c4[i][j+1][k]-c4[i][j][k],2)/(h*h) + pow(c4[i][j][k+1]-c4[i][j][k],2)/(h*h))  );
         
     }
     
     
     EnEM = 0.0;   

    ijkloop{

     ic[i][j][k] = 2.0*nc[i][j][k] - c[i][j][k];
     ic2[i][j][k] = 2.0*nc2[i][j][k] - c2[i][j][k];
     ic3[i][j][k] = 2.0*nc3[i][j][k] - c3[i][j][k];
     ic4[i][j][k] = 1.0-ic[i][j][k]-ic2[i][j][k]-ic3[i][j][k];

    }

        augmenc(ic, nx,ny,nz);
         augmenc(ic2, nx,ny,nz);
         augmenc(ic3, nx,ny,nz);
         augmenc(ic4, nx,ny,nz);
  
     
     ijkloop{
     
      EnEM = EnEM + h*h*h*( 
          0.25*( pow(nc[i+1][j][k]-nc[i][j][k],2)/(h*h) + pow(nc[i][j+1][k]-nc[i][j][k],2)/(h*h)+ pow(nc[i][j][k+1]-nc[i][j][k],2)/(h*h) )
         + 0.25*( pow(nc2[i+1][j][k]-nc2[i][j][k],2)/(h*h) + pow(nc2[i][j+1][k]-nc2[i][j][k],2)/(h*h)+ pow(nc2[i][j][k+1]-nc2[i][j][k],2)/(h*h) )
         + 0.25*( pow(nc3[i+1][j][k]-nc3[i][j][k],2)/(h*h) + pow(nc3[i][j+1][k]-nc3[i][j][k],2)/(h*h)+ pow(nc3[i][j][k+1]-nc3[i][j][k],2)/(h*h) )
         + 0.25*( pow(c4[i+1][j][k]-c4[i][j][k],2)/(h*h) + pow(c4[i][j+1][k]-c4[i][j][k],2)/(h*h)+ pow(c4[i][j][k+1]-c4[i][j][k],2)/(h*h) ) 
         
         + 0.25*( pow(ic[i+1][j][k]-ic[i][j][k],2)/(h*h) + pow(ic[i][j+1][k]-ic[i][j][k],2)/(h*h)+ pow(ic[i][j][k+1]-ic[i][j][k],2)/(h*h) )  
         + 0.25*( pow(ic2[i+1][j][k]-ic2[i][j][k],2)/(h*h) + pow(ic2[i][j+1][k]-ic2[i][j][k],2)/(h*h)+ pow(ic2[i][j][k+1]-ic2[i][j][k],2)/(h*h) )  
         + 0.25*( pow(ic3[i+1][j][k]-ic3[i][j][k],2)/(h*h) + pow(ic3[i][j+1][k]-ic3[i][j][k],2)/(h*h)+ pow(ic3[i][j][k+1]-ic3[i][j][k],2)/(h*h) ) 
         + 0.25*( pow(ic4[i+1][j][k]-ic4[i][j][k],2)/(h*h) + pow(ic4[i][j+1][k]-ic4[i][j][k],2)/(h*h)+ pow(ic4[i][j][k+1]-ic4[i][j][k],2)/(h*h) )
         + 0.5*SS/Cahn*( pow(nc[i][j][k]-c[i][j][k],2)  
         + pow(nc2[i][j][k]-c2[i][j][k],2)
         + pow(nc3[i][j][k]-c3[i][j][k],2) 
          + pow(1.0-nc[i][j][k]-nc2[i][j][k]-nc3[i][j][k] - (1.0-c[i][j][k]-c2[i][j][k]-c3[i][j][k]),2)   )  );
         
     }
     
     
     EnEM = EnEM + 0.5/Cahn*(3.0*r-or)-CC/Cahn;
      

        cube_copy(oc, c, 1, nx, 1, ny, 1, nz);
        cube_copy(oc2, c2, 1, nx, 1, ny, 1, nz);
         cube_copy(oc3, c3, 1, nx, 1, ny, 1, nz);
        


        cube_copy(c, nc,  1, nx, 1, ny, 1, nz);
        cube_copy(c2, nc2, 1, nx, 1, ny, 1, nz);
         cube_copy(c3, nc3, 1, nx, 1, ny, 1, nz);
    
        
        printf("it = %d\n", it);
        
        if (it % ns == 0) {
            count++;
            print_data(nc,nc2,nc3,c4,count);
            fprintf(myr,"%16.12f \n",EnEM);
            fprintf(myq,"%16.12f \n",EnEO);
            printf("print out counts %d\n", count);
        }
    }
    
    return 0;
}



void initialization(float ***c, float ***c2, float ***c3, float ***c4)
{
    extern float xright, yright,zright, h, gam, pi;
    
    int i, j,k;
    float x, y,z;
    
    ijkloop {
        x = ((float)i-0.5)*h;
        y = ((float)j-0.5)*h;
        z = ((float)k-0.5)*h;
        
   c[i][j][k] = 1.0/6.0 + 0.01*(2.0*rand()/(float)RAND_MAX - 1.0);
   
   c2[i][j][k] = 1.0/6.0 + 0.01*(2.0*rand()/(float)RAND_MAX - 1.0);
    c3[i][j][k] = 1.0/6.0 + 0.01*(2.0*rand()/(float)RAND_MAX - 1.0);

   c4[i][j][k] = 1.0 - c[i][j][k] - c2[i][j][k] - c3[i][j][k];
        
    }
    
     
}


void augmenc(float ***c, int nxt, int nyt,int nzt)
{
       int i, j, k;
   
    for (j=1; j<=nyt; j++)
        for (k=1; k<=nzt; k++) {
            
            c[0][j][k] = c[nxt][j][k];
            c[nxt+1][j][k] = c[1][j][k]; }
    
    for (i=0; i<=nx+1; i++)
        for (k=1; k<=nzt; k++) {
            
            c[i][0][k] = c[i][nyt][k];
            c[i][nyt+1][k] = c[i][1][k]; }
    
    for (i=0; i<=nxt+1; i++)
        for (j=0; j<=nyt+1; j++) {
            
            c[i][j][0] = c[i][j][1];
            c[i][j][nzt+1] = c[i][j][nzt]; }
            
    
    
}


/***************** 1order phase field *******************/

void cahn_ini(float ***c_old,  float ***H, float ***c_new)
{
    extern int nx, ny,nz;
    extern float ***ct, ***sc, ***smu, ***mu;
    
    int it_mg = 1, max_it_CH = 200;
    float resid = 1.0, tol = 1.0e-6;
    
    cube_copy(ct, c_old, 1, nx, 1, ny,1,nz);
    
    source_ini(c_old, H, sc, smu);
    
    while (it_mg <= max_it_CH && resid > tol) {
        
        vcycle_ini(c_new, mu, sc, smu, nx ,ny, nz,1);
        resid = error(ct, c_new, nx,ny,nz);
        cube_copy(ct, c_new, 1, nx, 1, ny,1,nz);
        
        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg-1);
}

void source_ini(float ***c_old, float ***H, float ***src_c, float ***src_mu)
{
    extern int nx,ny,nz;
    extern float dt, h, SS, Cahn, ***alp, r;
    
    int i,j,k;
   

    ijkloop {

        src_c[i][j][k] = c_old[i][j][k]/dt;
        
        src_mu[i][j][k] = -1/Cahn*(H[i][j][k] + alp[i][j][k])*r +SS/Cahn*c_old[i][j][k];
    }

       
}



void vcycle_ini(float ***uf_new, float ***wf_new, float ***su, float ***sw, int nxf, int nyf, int nzf,int ilevel)
{
    extern int n_level;
    
    relax_ini(uf_new, wf_new, su, sw, ilevel, nxf, nyf,nzf);
    
   if (ilevel < n_level) {
        
        int nxc, nyc, nzc;
        float ***uc_new, ***wc_new, ***duc, ***dwc, 
               ***uc_def, ***wc_def, ***uf_def, ***wf_def;
        
        nxc = nxf/2, nyc = nyf/2, nzc = nzf/2;
        
        uc_new = cube(1, nxc, 1, nyc, 1, nzc);
        wc_new = cube(1, nxc, 1, nyc, 1, nzc);
        duc = cube(1, nxc, 1, nyc, 1, nzc);
        dwc = cube(1, nxc, 1, nyc, 1, nzc);
        uc_def = cube(1, nxc, 1, nyc, 1, nzc);
        wc_def = cube(1, nxc, 1, nyc, 1, nzc);
        uf_def = cube(1, nxf, 1, nyf, 1, nzf);
        wf_def = cube(1, nxf, 1, nyf, 1, nzf);

        defect_ini(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf,nzf);

          zero_cube(uc_def, 1, nxc, 1, nyc, 1, nzc);
          zero_cube(wc_def, 1, nxc, 1, nyc, 1, nzc);
        
        vcycle_ini(uc_def,  wc_def, duc, dwc, nxc, nyc,nzc, ilevel+1);
        
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc, nzc);
        
        cube_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf, 1, nzf);
        
        relax_ini(uf_new,  wf_new, su, sw, ilevel, nxf, nyf,nzf);
        
        free_cube(uc_new, 1, nxc, 1, nyc, 1, nzc);
        free_cube(wc_new, 1, nxc, 1, nyc, 1, nzc);
        free_cube(duc, 1, nxc, 1, nyc, 1, nzc);
        free_cube(dwc, 1, nxc, 1, nyc, 1,nzc);
        free_cube(uc_def, 1, nxc, 1, nyc, 1, nzc);
        free_cube(wc_def, 1, nxc, 1, nyc, 1, nzc);
        free_cube(uf_def, 1, nxf, 1, nyf, 1, nzf);
        free_cube(wf_def, 1, nxf, 1, nyf, 1, nzf);
    }
   
}

void relax_ini(float ***c_new, float ***mu_new, float ***su, float ***sw, int ilevel, int nxt, int nyt,int nzt)
{
    extern int c_relax;
    extern float xright, dt, Cahn, SS,M;
    
    int i,j,k, iter;
    float ht2, a[4], f[2], det;
    
    ht2 = pow(xright/(float)nxt,2);
    
    for (iter=1; iter<=c_relax; iter++) {
        
        ijkloopt {

            a[0] = 1.0/dt;
            
            a[1] = -M;

            a[2] = 4.0/ht2 +SS/Cahn;
            if (k> 1)
                a[2] += 1/ ht2;
            if (k < nzt)
                a[2] += 1 / ht2;


         
            a[3] = 1.0;

          
            f[0] = su[i][j][k];
        
            
            f[1] = sw[i][j][k];
            if (i > 1)   f[1] += c_new[i-1][j][k]/ht2;
            else         f[1] += c_new[nxt][j][k]/ht2;
            
            if (i < nxt) f[1] += c_new[i+1][j][k]/ht2;
            else         f[1] += c_new[1][j][k]/ht2;
            
            if (j > 1)   f[1] += c_new[i][j-1][k]/ht2;
            else         f[1] += c_new[i][nyt][k]/ht2;
            
            if (j < nyt) f[1] += c_new[i][j+1][k]/ht2;
            else         f[1] += c_new[i][1][k]/ht2;

            if (k > 1)   f[1] += c_new[i][j][k-1]/ht2;
            
            if (k < nzt) f[1] += c_new[i][j][k+1]/ht2;
            
            det = a[0]*a[3] - a[1]*a[2]; 
            
            c_new[i][j][k] = (a[3]*f[0] - a[1]*f[1])/det;
            mu_new[i][j][k] = (-a[2]*f[0] + a[0]*f[1])/det;
        }
    }
    
}

void defect_ini(float ***duc, float ***dwc, float ***uf_new, float ***wf_new,
            float ***suf, float ***swf, int nxf, int nyf,int nzf)
{
    float ***ruf, ***rwf, ***ruc, ***rwc, ***rruf, ***rrwf;
    
    
    ruf = cube(1, nxf, 1, nyf, 1, nzf);
    rwf = cube(1, nxf, 1, nyf, 1, nzf);
    rruf = cube(1, nxf/2, 1, nyf/2, 1, nzf/2);
    rrwf = cube(1, nxf/2, 1, nyf/2, 1, nzf/2);
        
    nonL_ini(ruf, rwf, uf_new, wf_new, nxf, nyf,nzf);
    
    cube_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf, 1, nzf);
    
    restrict2(ruf, duc, rwf, dwc, nxf/2, nyf/2, nzf/2);
    
    free_cube(ruf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rwf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rruf, 1, nxf/2, 1, nyf/2, 1, nzf/2);
    free_cube(rrwf, 1, nxf/2, 1, nyf/2, 1, nzf/2);
}

void nonL_ini(float ***ru, float ***rw, float ***c_new, float ***mu_new, int nxt, int nyt,int nzt)
{
    extern float xright, dt, SS, Cahn,M;
    
    int i,j,k;
    float ss, ht2, ***lap_c, ***lap_mu;
    
    lap_c = cube(1, nxt, 1, nyt, 1, nzt);
    
     laplace_ch(c_new, lap_c, nxt, nyt,nzt);
     
     ht2 = pow(xright/(float)nxt,2);
    
    ijkloopt {

        ru[i][j][k] = c_new[i][j][k]/dt - M*mu_new[i][j][k];
        rw[i][j][k] = mu_new[i][j][k] -lap_c[i][j][k] +SS/Cahn*c_new[i][j][k];
    }
    
    free_cube(lap_c, 1, nxt, 1, nyt, 1, nzt);
    
}



/************************ 2 order phase field part **************************/


void cahn(float ***c_old,  float ***cc_old, float ***H, float ***c_new)
{
    extern int nx,ny,nz;
    extern float ***ct, ***sc, ***smu, ***mu;
    
    int it_mg = 1, max_it_CH = 200;
    float resid = 1.0, tol = 1.0e-6;
    
    cube_copy(ct, c_old, 1, nx, 1, ny,1,nz);
    
    source(c_old, cc_old, H, sc, smu);
    
    while (it_mg <= max_it_CH && resid > tol) {
        
        vcycle(c_new, mu, sc, smu, nx ,ny, nz,1);
        resid = error(ct, c_new, nx,ny,nz);
        cube_copy(ct, c_new, 1, nx, 1, ny, 1, nz);
        
        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg-1);
}

void source(float ***c_old, float ***cc_old, float ***H, float ***src_c, float ***src_mu)
{
    extern int nx,ny,nz;
    extern float dt, h, int_r, gam, SS,Cahn, ***alp;
    
    int i,j,k;
   
    augmenc(c_old, nx,ny,nz);

    ijkloop {

        src_c[i][j][k] = (4.0*c_old[i][j][k]- cc_old[i][j][k])/(2.0*dt);
        
        src_mu[i][j][k] = -1/Cahn*(H[i][j][k] + alp[i][j][k])*int_r +SS/Cahn*(2.0*c_old[i][j][k]-cc_old[i][j][k]);
    }

       
}


void laplace_ch(float ***a, float ***lap_a, int nxt, int nyt,int nzt)
{
    extern float xright, yright, zright;
    int i,j,k;
    float ht2, dadx_L, dadx_R, dady_B, dady_T, dadz_D, dadz_U;
    
    ht2 = pow(xright/(float)nxt,2);
    
    ijkloopt {
        
        if (i > 1)
            dadx_L = a[i][j][k] - a[i-1][j][k];
        else
            dadx_L = a[i][j][k] - a[nxt][j][k];
        
        if (i < nxt)
            dadx_R = a[i+1][j][k] - a[i][j][k];
        else
            dadx_R = a[1][j][k] - a[i][j][k];

        
        if (j > 1)
            dady_B = a[i][j][k] - a[i][j-1][k];
        else
            dady_B = a[i][j][k] - a[i][nyt][k];
        
        if (j < nyt)
            dady_T = a[i][j+1][k] - a[i][j][k];
        else
            dady_T = a[i][1][k] - a[i][j][k];


       if (k>1)
        dadz_D = a[i][j][k] - a[i][j][k-1];
        else   
		dadz_D = 0.0;   
             
        if (k<nzt)
        dadz_U = a[i][j][k+1] - a[i][j][k];
        else 
		dadz_U = 0.0; 
        
        lap_a[i][j][k] = (dadx_R - dadx_L + dady_T - dady_B+ dadz_U - dadz_D)/ht2;
    }
    
}

void vcycle(float ***uf_new, float ***wf_new, float ***su, float ***sw, int nxf, int nyf, int nzf,int ilevel)
{
    extern int n_level;
    
    relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf,nzf);
    
   if (ilevel < n_level) {
        
        int nxc, nyc, nzc;
        float ***uc_new, ***wc_new, ***duc, ***dwc, 
               ***uc_def, ***wc_def, ***uf_def, ***wf_def;
        
        nxc = nxf/2, nyc = nyf/2; nzc = nzf/2;
        
        uc_new = cube(1, nxc, 1, nyc, 1, nzc);
        wc_new = cube(1, nxc, 1, nyc, 1, nzc);
        duc = cube(1, nxc, 1, nyc, 1, nzc);
        dwc = cube(1, nxc, 1, nyc, 1, nzc);
        uc_def = cube(1, nxc, 1, nyc, 1, nzc);
        wc_def = cube(1, nxc, 1, nyc, 1, nzc);
        uf_def = cube(1, nxf, 1, nyf, 1, nzf);
        wf_def = cube(1, nxf, 1, nyf, 1, nzf);

        defect(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf,nzf);

          zero_cube(uc_def, 1, nxc, 1, nyc, 1, nzc);
          zero_cube(wc_def, 1, nxc, 1, nyc, 1, nzc);
        
        vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, nzc, ilevel+1);
        
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc, nzc);
        
        cube_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf, 1, nzf);
        
        relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf, nzf);
        
        free_cube(uc_new, 1, nxc, 1, nyc, 1, nzc);
        free_cube(wc_new, 1, nxc, 1, nyc, 1, nzc);
        free_cube(duc, 1, nxc, 1, nyc, 1, nzc);
        free_cube(dwc, 1, nxc, 1, nyc, 1,nzc);
        free_cube(uc_def, 1, nxc, 1, nyc, 1, nzc);
        free_cube(wc_def, 1, nxc, 1, nyc, 1, nzc);
        free_cube(uf_def, 1, nxf, 1, nyf, 1, nzf);
        free_cube(wf_def, 1, nxf, 1, nyf, 1, nzf);
    }
   
}

void relax(float ***c_new, float ***mu_new, float ***su, float ***sw, int ilevel, int nxt, int nyt,int nzt)
{
    extern int c_relax;
    extern float xright, dt, Cahn, SS,M;
    
    int i,j,k, iter;
    float ht2, a[4], f[2], det;
    
    ht2 = pow(xright/(float)nxt,2);
    
    for (iter=1; iter<=c_relax; iter++) {
        
        ijkloopt {

            a[0] = 3.0/(2.0*dt);
            
            a[1] = -M;

                   
            a[2] = 4.0/ht2 +SS/Cahn;

            if (k > 1)
                a[2] += 1 / ht2;
            if (k < nzt)
                a[2] += 1/ ht2;
         
            a[3] = 1.0;

          
            f[0] = su[i][j][k];
            
            
            f[1] = sw[i][j][k];
            if (i > 1)   f[1] += c_new[i-1][j][k]/ht2;
            else         f[1] += c_new[nxt][j][k]/ht2;
            
            if (i < nxt) f[1] += c_new[i+1][j][k]/ht2;
            else         f[1] += c_new[1][j][k]/ht2;
            
            if (j > 1)   f[1] += c_new[i][j-1][k]/ht2;
            else         f[1] += c_new[i][nyt][k]/ht2;
            
            if (j < nyt) f[1] += c_new[i][j+1][k]/ht2;
            else         f[1] += c_new[i][1][k]/ht2;

            if (k > 1)   f[1] += c_new[i][j][k-1]/ht2;
        
            
            if (k < nzt) f[1] += c_new[i][j][k+1]/ht2;
        


            
            det = a[0]*a[3] - a[1]*a[2];
            
            c_new[i][j][k] = (a[3]*f[0] - a[1]*f[1])/det;
            mu_new[i][j][k] = (-a[2]*f[0] + a[0]*f[1])/det;
        }
    }
    
}

void defect(float ***duc, float ***dwc, float ***uf_new, float ***wf_new,
            float ***suf, float ***swf, int nxf, int nyf,int nzf)
{
    float ***ruf, ***rwf, ***ruc, ***rwc, ***rruf, ***rrwf;
    
    ruf = cube(1, nxf, 1, nyf, 1, nzf);
    rwf = cube(1, nxf, 1, nyf, 1, nzf);
    rruf = cube(1, nxf/2, 1, nyf/2, 1, nzf/2);
    rrwf = cube(1, nxf/2, 1, nyf/2, 1, nzf/2);
    
    nonL(ruf, rwf, uf_new, wf_new, nxf, nyf,nzf);
    
    cube_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf, 1, nzf);
    
    restrict2(ruf, duc, rwf, dwc, nxf/2, nyf/2, nzf/2);
    
    free_cube(ruf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rwf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rruf, 1, nxf/2, 1, nyf/2, 1, nzf/2);
    free_cube(rrwf, 1, nxf/2, 1, nyf/2, 1, nzf/2);
}

void nonL(float ***ru, float ***rw, float ***c_new, float ***mu_new, int nxt, int nyt,int nzt)
{
    extern float xright, dt, Cahn, SS,M;
    
    int i,j,k;
    float ss, ht2, ***lap_c, ***lap_mu;
    
    lap_c = cube(1, nxt, 1, nyt, 1, nzt);
  
    
     laplace_ch(c_new, lap_c, nxt, nyt,nzt);
    

     ht2 = pow(xright/(float)nxt,2);
    
    ijkloopt {

        ru[i][j][k] = 3.0*c_new[i][j][k]/(2.0*dt) - M*mu_new[i][j][k];
        rw[i][j][k] = mu_new[i][j][k] - lap_c[i][j][k] + SS/Cahn*c_new[i][j][k];
    }
    
    free_cube(lap_c, 1, nxt, 1, nyt, 1, nzt);
  
}

void restrict2(float ***uf, float ***uc, float ***vf, float ***vc, int nxt, int nyt,int nzt)
{
    int i, j, k;
    
    ijkloopt {
         uc[i][j][k] = 0.125*(uf[2*i][2*j][2*k]+uf[2*i-1][2*j][2*k]
                            + uf[2*i][2*j-1][2*k]+uf[2*i-1][2*j-1][2*k]
                            + uf[2*i][2*j][2*k-1]+uf[2*i-1][2*j][2*k-1]
                            + uf[2*i][2*j-1][2*k-1]+uf[2*i-1][2*j-1][2*k-1]);

          vc[i][j][k] = 0.125*(vf[2*i][2*j][2*k]+vf[2*i-1][2*j][2*k]
                            + vf[2*i][2*j-1][2*k]+vf[2*i-1][2*j-1][2*k]
                            + vf[2*i][2*j][2*k-1]+vf[2*i-1][2*j][2*k-1]
                            + vf[2*i][2*j-1][2*k-1]+vf[2*i-1][2*j-1][2*k-1]);
    }
    
}


void restrict1(float ***vf, float ***vc, int nxt, int nyt, int nzt)
{
    int i,j,k;
    
    ijkloopt {
        vc[i][j][k] = 0.125*(vf[2*i][2*j][2*k]+vf[2*i-1][2*j][2*k]
                            + vf[2*i][2*j-1][2*k]+vf[2*i-1][2*j-1][2*k]
                            + vf[2*i][2*j][2*k-1]+vf[2*i-1][2*j][2*k-1]
                            + vf[2*i][2*j-1][2*k-1]+vf[2*i-1][2*j-1][2*k-1]);
    }
    
}

void prolong_ch(float ***uc, float ***uf, float ***vc, float ***vf, int nxc, int nyc, int nzc)
{
   int i, j, k;

   for (i=1; i<=nxc; i++) 
      for (j=1; j<=nyc; j++)
         for (k=1; k<=nzc; k++){

         uf[2*i][2*j][2*k] = uc[i][j][k];
         uf[2*i-1][2*j][2*k] = uc[i][j][k];
         uf[2*i][2*j-1][2*k] = uc[i][j][k]; 
         uf[2*i-1][2*j-1][2*k] = uc[i][j][k];

         uf[2*i][2*j][2*k-1] = uc[i][j][k];
         uf[2*i-1][2*j][2*k-1] = uc[i][j][k];
         uf[2*i][2*j-1][2*k-1] = uc[i][j][k]; 
         uf[2*i-1][2*j-1][2*k-1] = uc[i][j][k];

         vf[2*i][2*j][2*k] = vc[i][j][k];
         vf[2*i-1][2*j][2*k] = vc[i][j][k];
         vf[2*i][2*j-1][2*k] = vc[i][j][k]; 
         vf[2*i-1][2*j-1][2*k] = vc[i][j][k];

         vf[2*i][2*j][2*k-1] = vc[i][j][k];
         vf[2*i-1][2*j][2*k-1] = vc[i][j][k];
         vf[2*i][2*j-1][2*k-1] = vc[i][j][k]; 
         vf[2*i-1][2*j-1][2*k-1] = vc[i][j][k];


   }
}
float error(float ***c_old, float ***c_new, int nxt, int nyt,int nzt)
{
    float ***r, res;
    
    r = cube(1, nxt, 1, nyt, 1, nzt);
    
    cube_sub(r, c_new, c_old, 1, nxt, 1, nyt, 1, nzt);
    res = cube_max(r, 1, nxt, 1, nyt, 1, nzt);
    
    free_cube(r, 1, nxt, 1, nyt, 1, nzt);
    
    return res;
}



/****************************************************/

float Ef(float ***c, float ***c2, float ***c3, float ***c4,int nxt, int nyt, int nzt)
{
    extern float Cahn, h, eta, theta, ksi, alpha1, alpha2, beta, delta, psis,CC;
    float r, res;
    int i,j,k;

    
      r = 0.0;

   ijkloopt{ r = r + h*h*h*( 0.25*pow(c[i][j][k],2)*pow(c[i][j][k]-1.0,2) 
   + 0.25*pow(c2[i][j][k],2)*pow(c2[i][j][k]-1.0,2) 
         + 0.25*pow(c3[i][j][k],2)*pow(c3[i][j][k]-1.0,2) 
         + 0.25*pow(c4[i][j][k],2)*pow(c4[i][j][k]-1.0,2)  );  }

   res = r+CC;
      
    return res;
}


void lagrangeH(float ***H , int nxt,int nyt, int nzt)
{
    extern float h, volume;
    float r;
    int i,j,k;
    r =0.0;
    ijkloopt{
        r = r+h*h*h*H[i][j][k];
    }
    r= r/volume;

    ijkloopt{
        H[i][j][k] = H[i][j][k]-r;
    }
}



/***************** util ******************/
float ***cube(int xl, int xr, int yl, int yr, int zl, int zr)
{
  int i,j,k,nrow=xr-xl+1, ncol=yr-yl+1, ndep=zr-zl+1;
  float ***t;
  
  t = (float ***) malloc(((nrow+1)*sizeof(float**)));
  memset(t, 0, (nrow+1) * sizeof(float**));
  t += 1;
  t -= xl;

  t[xl]=(float **) malloc(((nrow*ncol+1)*sizeof(float *)));
  memset(t[xl], 0, (nrow*ncol+1) *sizeof(float *));
  t[xl] += 1;
  t[xl] -= yl;

  t[xl][yl] = (float *) malloc(((nrow*ncol*ndep+1)*sizeof(float))); 
  memset(t[xl][yl], 0, (nrow*ncol*ndep+1)*sizeof(float));
  t[xl][yl] += 1;
  t[xl][yl] -= zl;

  for(j=yl+1; j<=yr; j++)
     t[xl][j]=t[xl][j-1]+ndep;

  for(i=xl+1; i<=xr; i++)
    {
       t[i] = t[i-1]+ncol;
       t[i][yl] = t[i-1][yl]+ncol*ndep;
    for(j=yl+1; j<=yr; j++)
       t[i][j] = t[i][j-1]+ndep;
    }

  return t;
} 



void free_cube(float ***t, int xl, int xr, int yl, int yr, int zl, int zr)
{
   free((char *) (t[xl][yl]+zl-1));
   free((char *) (t[xl]+yl-1));
   free((char *) (t+xl-1));
}

void zero_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i,j, k;

   for (i=xl; i<=xr; i++)
      for (j=yl; j<=yr; j++)
         for (k=zl; k<=zr; k++){

        a[i][j][k] = 0.0;
	 
  }

   return;   
}

void cube_copy(float ***a, float ***b, int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++) 
        for (j=yl; j<=yr; j++) 
            for (k=zl; k<=zr; k++) {
            a[i][j][k] = b[i][j][k];
         
    }
    
    return;
}

void cube_copy2(float ***a, float ***b, float ***a2, float ***b2,
               int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++) 
        for (j=yl; j<=yr; j++) 
           for (k=zl; k<=zr; k++) {
            a[i][j][k] = b[i][j][k];
            a2[i][j][k] = b2[i][j][k];
            }
        
    
    
    return;
}

void cube_add(float ***a, float ***b, float ***c,
             int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++) 
        for (j=yl; j<=yr; j++) 
            for (k=zl; k<=zr; k++) {
            a[i][j][k] = b[i][j][k]+c[i][j][k];
           
        
    }
    
    return;
}

void cube_add2(float ***a, float ***b, float ***c,
              float ***a2, float ***b2, float ***c2,
              int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++) 
        for (j=yl; j<=yr; j++) 
             for (k=zl; k<=zr; k++) {
            a[i][j][k] = b[i][j][k]+c[i][j][k];
            a2[i][j][k] = b2[i][j][k]+c2[i][j][k];
           }
    
    
    return;
}

void cube_sub(float ***a, float ***b, float ***c, 
              int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i=xl; i<=xr; i++)
      for (j=yl; j<=yr; j++)
         for (k=zl; k<=zr; k++){

        a[i][j][k]=b[i][j][k]-c[i][j][k];
	 
        }

   return;
}


void cube_sub2(float ***a, float ***b, float ***c, 
               float ***a2, float ***b2, float ***c2, 
               int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;

   for (i=xl; i<=xr; i++)
      for (j=yl; j<=yr; j++)
         for (k=zl; k<=zr; k++){

        a[i][j][k]=b[i][j][k]-c[i][j][k];
        a2[i][j][k]=b2[i][j][k]-c2[i][j][k];	 

        }

   return;
}

float cube_max(float ***a, 
               int xl, int xr, int yl, int yr, int zl, int zr)
{
   int i, j, k;
   float x = 0.0;

   for (i=xl; i<=xr; i++)
      for (j=yl; j<=yr; j++)
         for (k=zl; k<=zr; k++){

         if (fabs(a[i][j][k]) > x)
            x = fabs(a[i][j][k]);
      }

   return x;
}


 void print_cube(FILE *fptr,
               float ***a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
 {
    int i, j, k;
 
    for(k = ndl; k <= ndh; k++) {
        for(j = ncl; j <= nch; j++) {
 		  for(i = nrl; i <= nrh; i++) {
          fprintf(fptr, " %f", a[i][j][k]);
    
       fprintf(fptr, "\n");
    }}}

    return;


 }




void print_data(float ***c1,float ***c2,float ***c3, float ***c4, int kk)
{
   extern int nx,ny,nz;

   extern char bufferc[2000],bufferc2[2000],bufferc3[2000],bufferc4[2000];
   FILE *fc,*fc2,*fc3,*fc4;
   int i,j,k;

    sprintf(bufferc,"/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/datac1%d.m",kk);
    sprintf(bufferc2,"/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/datac2%d.m",kk);
    sprintf(bufferc3,"/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/datac3%d.m",kk);
    sprintf(bufferc4,"/Users/wujingwen/research/fluid_mechanics/msavncac/msav+ncac/example9cac3D/4component/c/dataout/datac4%d.m",kk);


   fc = fopen(bufferc,"w");
    fc2 = fopen(bufferc2,"w");
    fc3 = fopen(bufferc3,"w");
     fc4 = fopen(bufferc4,"w");


   print_cube(fc,c1,1,nx,1,ny,1,nz);
   print_cube(fc2,c2,1,nx,1,ny,1,nz);
   print_cube(fc3,c3,1,nx,1,ny,1,nz); 
   print_cube(fc4,c4,1,nx,1,ny,1,nz); 

    fclose(fc);
    fclose(fc2);
    fclose(fc3);
    fclose(fc4);

   
   return;
   
}