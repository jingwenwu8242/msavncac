/************************************************

x:periodic ,y:Neunman
for compare with SAV liquid lens
cac
   with variable density and variable viscosity
   // phase seperation,liquid lens

************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdexcept>
#include <string.h>
#include <time.h>
#include "util.h"
#include "bnsch.h"
#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)
#define NR_END 1

int nx, ny, n_level,  c_relax, count, Num;
double pi, xleft, xright, yleft, yright,volume, h, dt, gam, Cahn,M,  Ev, oEv, acount, int_r, orc, r, nr, EnEO, EnEM,
        **worku, **workv, **workp, **mc, SS, **ct, **sc, **smu, **mu, **H1, **H2, **H3, **alp, **intc1, **intc2,
        **intc3, **c3, **oc3 ,**ic, **ic2, **ic3;
char bufferc[1000]={0}, bufferc2[1000]={0}, bufferc3[1000]={0};
int main()
{
    extern int nx, ny, n_level,  c_relax, count, Num;
    extern double pi, xleft, xright, yleft, yright,volume, h, dt, gam, Cahn, Ev,  oEv, acount, int_r, orc, r, nr, EnEO, EnEM,
                  **worku, **workv, **workp, **mc, SS, **ct, **sc, **smu, **mu, **H1, **H2, **H3, **alp, **intc1, **intc2,
                  **intc3, **c3, **oc3, **ic, **ic2, **ic3;

    extern char  bufferc[1000], bufferc2[1000], bufferc3[1000];
    int i,j, it, max_it, ns;
    double **c, **nc, **oc, **c2, **nc2, **oc2;
    
    FILE  *fc, *fc2, *fc3, *myr, *myq;
    
    
    /*****/
     
    nx = gnx;
    ny = gny;
    n_level = (int)(log(nx)/log(2)-0.9);
    
    c_relax = 5;
    
    pi = 4.0*atan(1.0);
    
    xleft = 0.0, xright =2.0;//xright=1.0 for accuracy
    yleft = 0.0, yright = 1.0;
    volume = (xright-xleft)*(yright-yleft);
    
    count = 0;
    
    /***********************/
    
    h = xright/(double)nx;
    dt =0.04*h;

    gam = 1.0*h/(4.0*sqrt(2.0)*atanh(0.9));
   

    Cahn = pow(gam,2);
    M =0.1;


    SS=2;
         
    max_it=1200;

    ns = max_it/40;


    Num = 3;  //number of components
  

    /***********************/
    
    printf("nx = %d, ny = %d\n", nx, ny);
    printf("n_level = %d\n", n_level);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);
    
 
    intc1 = dmatrix(0, nx+1, 0, ny+1);
    intc2 = dmatrix(0, nx+1, 0, ny+1);
    intc3 = dmatrix(0, nx+1, 0, ny+1);

    ic = dmatrix(0, nx+1, 0, ny+1);
    ic2 = dmatrix(0, nx+1, 0, ny+1);
    ic3 = dmatrix(0, nx+1, 0, ny+1);//for modified energy

    oc = dmatrix(0, nx+1, 0, ny+1);
    c = dmatrix(0, nx+1, 0, ny+1);
    nc = dmatrix(0, nx+1, 0, ny+1);

    oc2 = dmatrix(0, nx+1, 0, ny+1); 
    c2 = dmatrix(0, nx+1, 0, ny+1);
    nc2 = dmatrix(0, nx+1, 0, ny+1);

    c3 = dmatrix(0, nx+1, 0, ny+1);
    oc3 = dmatrix(0, nx+1, 0, ny+1);

    H1 = dmatrix(0, nx+1, 0, ny+1);
    H2 = dmatrix(0, nx+1, 0, ny+1);
    H3 = dmatrix(0, nx+1, 0, ny+1);
    alp = dmatrix(0, nx+1, 0, ny+1); 
    
    worku = dmatrix(0, nx, 1, ny);
    workv = dmatrix(1, nx, 0, ny);
    workp = dmatrix(1, nx, 1, ny);
   
    ct = dmatrix(1, nx, 1, ny);
    sc = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);
    mu = dmatrix(1, nx, 1, ny); 
    
    zero_matrix(mu, 1, nx, 1, ny); 
 
    sprintf(bufferc,"/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/eplison/data1/datac.m");
    sprintf(bufferc2,"/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/eplison/data1/datac2.m");
    sprintf(bufferc3,"/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/eplison/data1/datac3.m");
   
    fc = fopen(bufferc,"w");
    fc2 = fopen(bufferc2,"w");
    fc3 = fopen(bufferc3,"w");
  
    fclose(fc);
    fclose(fc2);
    fclose(fc3);
    
    initialization(c, c2, c3);


    mat_copy(oc, c, 1, nx, 1, ny);
    mat_copy(oc2, c2, 1, nx, 1, ny);

    mat_copy(nc, c, 1, nx, 1, ny);
    mat_copy(nc2, c2, 1, nx, 1, ny);
    
   save as row//
    
    augmenc(c, nx, ny);
    augmenc(c2, nx, ny);
    augmenc(c3, nx, ny);
    
    //calculate initial energy
     EnEO = 0.0;
     
     ijloop{
     
      EnEO = EnEO + h*h*(1/Cahn* 0.25*pow(nc[i][j],2)*pow(nc[i][j]-1.0,2) 
         + 1/Cahn*0.25*pow(nc2[i][j],2)*pow(nc2[i][j]-1.0,2) 
         + 1/Cahn*0.25*pow(c3[i][j],2)*pow(c3[i][j]-1.0,2) 
         + 0.5*( pow(nc[i+1][j]-nc[i][j],2)/(h*h) + pow(nc[i][j+1]-nc[i][j],2)/(h*h) )
         + 0.5*( pow(nc2[i+1][j]-nc2[i][j],2)/(h*h) + pow(nc2[i][j+1]-nc2[i][j],2)/(h*h) )
         + 0.5*( pow(c3[i+1][j]-c3[i][j],2)/(h*h) + pow(c3[i][j+1]-c3[i][j],2)/(h*h) )  );
         
     }
       
    myr = fopen("/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/eplison/data1/Mene.m","w");
    myq = fopen("/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/eplison/data1/Oene.m","w");
    
    fprintf(myr,"%16.12f \n",EnEO);
    fprintf(myq,"%16.12f \n",EnEO);

 
   // Backward Euler for the first step

  for( it = 1; it <= 1; it++){

       Ev = Ef(c, c2, c3, nx, ny);

       r = Ev;
       orc = r;
      /*  printf("Ev=%f\n",Ev);
       if (std::isnan(Ev))
        {
        throw std::invalid_argument("wrong");
        } */
       


   ijloop{

      H1[i][j] = c[i][j]*(c[i][j]-0.5)*(c[i][j]-1.0)/Ev;
      H2[i][j] = c2[i][j]*(c2[i][j]-0.5)*(c2[i][j]-1.0)/Ev;
      H3[i][j] = c3[i][j]*(c3[i][j]-0.5)*(c3[i][j]-1.0)/Ev;

     /* 
      if (std::isnan(alp[i][j]))
        {
        throw std::invalid_argument("wrong");
        } */
        }

    lagrangeH(H1,nx,ny);
    lagrangeH(H2,nx,ny);
    lagrangeH(H3,nx,ny);

      ijloop{
          alp[i][j] = -(1.0/Num)*(H1[i][j] + H2[i][j] + H3[i][j] );

      }


    //calculate the phase-field: c1
          cahn_ini(c, H1, nc);

    //calculate the phase-field: c2
          cahn_ini(c2, H2, nc2);

   //calculate c3:

     ijloop{ c3[i][j] = 1.0 - nc[i][j] - nc2[i][j]; }


   //update the auxiliary variable 

     acount = 0.0;

    ijloop{

        acount = acount + h*h*( H1[i][j]*r*(nc[i][j]-c[i][j])/dt
               + H2[i][j]*r*(nc2[i][j]-c2[i][j])/dt
               + H3[i][j]*r*( (1.0-nc[i][j]-nc2[i][j])-(1.0-c[i][j]-c2[i][j]) )/dt );

     }


        nr = r + dt*acount;
        orc =r;

         r = nr;

                 
  
        mat_copy(c, nc, 1, nx, 1, ny);
        mat_copy(c2, nc2, 1, nx, 1, ny);
        

   }


     // BDF2 from the second step


    for (it=2; it<=max_it; it++) {

   
     int_r = 2.0*r - orc;
     

    ijloop{

       intc1[i][j] = 2.0*c[i][j] - oc[i][j]; 
       intc2[i][j] = 2.0*c2[i][j] - oc2[i][j];
       intc3[i][j] = 1.0 - intc1[i][j] - intc2[i][j];

      }
      


       Ev = Ef(intc1, intc2, intc3, nx, ny);
 

    ijloop{

      H1[i][j] = intc1[i][j]*(intc1[i][j]-0.5)*(intc1[i][j]-1.0)/Ev;
      H2[i][j] = intc2[i][j]*(intc2[i][j]-0.5)*(intc2[i][j]-1.0)/Ev;
      H3[i][j] = intc3[i][j]*(intc3[i][j]-0.5)*(intc3[i][j]-1.0)/Ev;
     }
     lagrangeH(H1,nx,ny);
      lagrangeH(H2,nx,ny);
      lagrangeH(H3,nx,ny);

      ijloop{
          alp[i][j] = -(1.0/Num)*(H1[i][j] + H2[i][j] + H3[i][j] );

      }

      
  

    //calculate the phase-field: c1
          cahn(c, oc, H1, nc);

    //calculate the phase-field: c2
          cahn(c2, oc2, H2, nc2);

   //calculate c3:

     ijloop{ c3[i][j] = 1.0 - nc[i][j] - nc2[i][j]; }


   //update the auxiliary variable 

     acount = 0.0;

    ijloop{

        acount = acount + h*h*( H1[i][j]*int_r*(3.0*nc[i][j]-4.0*c[i][j]+oc[i][j])/(2.0*dt)
               + H2[i][j]*int_r*(3.0*nc2[i][j]-4.0*c2[i][j]+oc2[i][j])/(2.0*dt)
               + H3[i][j]*int_r*(3.0*(1.0-nc[i][j]-nc2[i][j])-4.0*(1.0-c[i][j]-c2[i][j])+(1.0-oc[i][j]-oc2[i][j]))/(2.0*dt) );

     }


        nr = ( 4.0*r - orc + 2.0*dt*acount )/3.0;
   
         orc = r;
         r = nr; 
         
         augmenc(nc, nx, ny);
         augmenc(nc2, nx, ny);
         augmenc(c3, nx, ny);
         
     //calcualte original energy
     EnEO = 0.0;
     
     ijloop{
     
      EnEO = EnEO + h*h*(1/Cahn* 0.25*pow(nc[i][j],2)*pow(nc[i][j]-1.0,2) + 1/Cahn*0.25*pow(nc2[i][j],2)*pow(nc2[i][j]-1.0,2) 
         + 1/Cahn*0.25*pow(c3[i][j],2)*pow(c3[i][j]-1.0,2) 
         + 0.5*( pow(nc[i+1][j]-nc[i][j],2)/(h*h) + pow(nc[i][j+1]-nc[i][j],2)/(h*h) )
         + 0.5*( pow(nc2[i+1][j]-nc2[i][j],2)/(h*h) + pow(nc2[i][j+1]-nc2[i][j],2)/(h*h) )
         + 0.5*( pow(c3[i+1][j]-c3[i][j],2)/(h*h) + pow(c3[i][j+1]-c3[i][j],2)/(h*h) )  );
         
     }
     
     //calcualte modified energy
          EnEM = 0.0; 
     ijloop{

     ic[i][j] = 2.0*nc[i][j] - c[i][j];
     ic2[i][j] = 2.0*nc2[i][j] - c2[i][j];
     ic3[i][j] = 1.0-ic[i][j]-ic2[i][j];

    }

        augmenc(ic, nx, ny);
         augmenc(ic2, nx, ny);
         augmenc(ic3, nx, ny);
  
     
     ijloop{
     
      EnEM = EnEM + h*h*( 
          0.25*( pow(nc[i+1][j]-nc[i][j],2)/(h*h) + pow(nc[i][j+1]-nc[i][j],2)/(h*h) )
         + 0.25*( pow(nc2[i+1][j]-nc2[i][j],2)/(h*h) + pow(nc2[i][j+1]-nc2[i][j],2)/(h*h) )
         + 0.25*( pow(c3[i+1][j]-c3[i][j],2)/(h*h) + pow(c3[i][j+1]-c3[i][j],2)/(h*h) ) 
         
         + 0.25*( pow(ic[i+1][j]-ic[i][j],2)/(h*h) + pow(ic[i][j+1]-ic[i][j],2)/(h*h) )  
         + 0.25*( pow(ic2[i+1][j]-ic2[i][j],2)/(h*h) + pow(ic2[i][j+1]-ic2[i][j],2)/(h*h) )  
         + 0.25*( pow(ic3[i+1][j]-ic3[i][j],2)/(h*h) + pow(ic3[i][j+1]-ic3[i][j],2)/(h*h) ) 
         + 0.5*SS/Cahn*( pow(nc[i][j]-c[i][j],2)  + pow(nc2[i][j]-c2[i][j],2)  + pow(1.0-nc[i][j]-nc2[i][j] - (1.0-c[i][j]-c2[i][j]),2)   )  );
         
     }
     
     
     EnEM = EnEM + 0.5/Cahn*(3.0*r-orc);
      
  
        mat_copy(oc, c, 1, nx, 1, ny);
        mat_copy(oc2, c2, 1, nx, 1, ny);
        mat_copy(c, nc, 1, nx, 1, ny);
        mat_copy(c2, nc2, 1, nx, 1, ny);
        
        printf("it = %d\n", it);
        
        if (it % ns == 0) {
            count++;
            print_data(nc,nc2,c3);
            fprintf(myr,"%16.12f \n",EnEM);
            fprintf(myq,"%16.12f \n",EnEO);
            printf("print out counts %d\n", count);
        }
    }
    
    return 0;
}



void initialization(double **c, double **c2, double **c3)
{
    extern double xright, yright, h, gam, pi;
    
    int i, j;
    double x, y;
    
    ijloop {
        x = ((double)i-0.5)*h;
        y = ((double)j-0.5)*h;

        /* c[i][j] = 0.33+0.01*cos(3*pi*x)+0.04*cos(5*pi*x);
        c2[i][j] =  0.33+0.01*cos(2*pi*x)+0.02*cos(4*pi*x);
        c3[i][j] = 1.0 - c[i][j] - c2[i][j ]; */  



    
     /* c[i][j] = 0.33 + 0.01*(2.0*rand()/(float)RAND_MAX - 1.0);
   
    c2[i][j] = 0.33 + 0.01*(2.0*rand()/(float)RAND_MAX - 1.0);

   c3[i][j] = 1.0 - c[i][j] - c2[i][j];  */

     /*   c[i][j] = 0.5 + 0.5 * tanh(min(sqrt(pow(x - 1, 2) + pow(y - 0.5, 2)) - 0.22, y - 0.5) / gam);
        c2[i][j] = 0.5 + 0.5 * tanh(-max(0.22 - sqrt(pow(x - 1, 2) + pow(y - 0.5, 2)), y - 0.5) / gam);
        c3[i][j] = 1 - c[i][j] - c2[i][j]; */

 

/* 
    c[i][j]=0.5+0.5*tanh((0.25-sqrt(pow(x-1.26,2)+pow(y-0.5,2)))/gam);
    c2[i][j]=0.5+0.5*tanh((0.25-sqrt(pow(x-0.74,2)+pow(y-0.5,2)))/gam);
    c3[i][j]=1-c[i][j]-c2[i][j]; */


     /*     c[i][j] = c2[i][j] = c3[i][j] = 0;

        if (x >= 0.25 && x <= 0.5 && y >= 0.25 && y <= 0.58)
            c[i][j] = 1.0;
        if (x > 0.5 && x <= 0.75 && y >= 0.25 && y <= 0.58)
            c2[i][j] = 1.0;
        if (x >= 0.25 && x <= 0.75 && y > 0.58 && y <= 0.75)
            c3[i][j] = 1.0;  //in this case ,we should add CC in r*/


/*         c[i][j] =0.5+0.5*tanh(min(sqrt(pow(x-1,2)+pow(y-0.5,2))-0.3,y-0.5)/gam);
        c2[i][j]=0.5+0.5*tanh(-max(0.3-sqrt(pow(x-1,2)+pow(y-0.5,2)),y-0.5)/gam);
        c3[i][j]= 1-c[i][j]-c2[i][j]; */


       /* c2[i][j]=0.5+0.5*tanh((0.15-sqrt(pow(x-0.5,2)+pow(y-0.5,2)))/(2*sqrt(2)*gam));

        c[i][j]= max(0.5+0.5*tanh((y-0.5)/(2*sqrt(2)*gam))-c2[i][j],0);
        c3[i][j]= 1-c[i][j]-c2[i][j];  */

         c[i][j]=0.5+0.5*tanh((0.25-sqrt(pow(x-1.26,2)+pow(y-0.5,2)))/gam);
    c2[i][j]=0.5+0.5*tanh((0.25-sqrt(pow(x-0.74,2)+pow(y-0.5,2)))/gam);
    c3[i][j]=1-c[i][j]-c2[i][j];


    }
    
     
}


void augmenc(double **c, int nxt, int nyt)
{
    int i, j;
    
    for (j=1; j<=nyt; j++) { 
        c[0][j] = c[nxt][j];
        c[nxt+1][j] = c[1][j];
    }
    
    for (i=0; i<=nxt+1; i++) { 
        c[i][0] = c[i][1];
        c[i][nyt+1] = c[i][nyt];
    }
    
}


/*************** 1order phase field *****************/

void cahn_ini(double **c_old,  double **H, double **c_new)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu, **mu;
    
    int it_mg = 1, max_it_CH = 200;
    double resid = 1.0, tol = 1.0e-6;
    
    mat_copy(ct, c_old, 1, nx, 1, ny);
    
    source_ini(c_old, H, sc, smu);
    
    while (it_mg <= max_it_CH && resid > tol) {
        
        vcycle_ini(c_new, mu, sc, smu, nx ,ny, 1);
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny);
        
        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg-1);
}

void source_ini(double **c_old, double **H, double **src_c, double **src_mu)
{
    extern int nx, ny;
    extern double dt, h, SS, Cahn, **alp, r;
    
    int i, j;
   

    ijloop{

        src_c[i][j] = c_old[i][j]/dt;
        
        src_mu[i][j] = -1/Cahn*(H[i][j] + alp[i][j])*r +SS/Cahn*c_old[i][j];
    }

       
}



void vcycle_ini(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel)
{
    extern int n_level;
    
    relax_ini(uf_new, wf_new, su, sw, ilevel, nxf, nyf);
    
   if (ilevel < n_level) {
        
        int nxc, nyc;
        double **uc_new, **wc_new, **duc, **dwc, 
               **uc_def, **wc_def, **uf_def, **wf_def;
        
        nxc = nxf/2, nyc = nyf/2;
        
        uc_new = dmatrix(1, nxc, 1, nyc);
        wc_new = dmatrix(1, nxc, 1, nyc);
        duc = dmatrix(1, nxc, 1, nyc);
        dwc = dmatrix(1, nxc, 1, nyc);
        uc_def = dmatrix(1, nxc, 1, nyc);
        wc_def = dmatrix(1, nxc, 1, nyc);
        uf_def = dmatrix(1, nxf, 1, nyf);
        wf_def = dmatrix(1, nxf, 1, nyf); 

        defect_ini(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf);

          zero_matrix(uc_def, 1, nxc, 1, nyc);
          zero_matrix(wc_def, 1, nxc, 1, nyc);
        
        vcycle_ini(uc_def,  wc_def, duc, dwc, nxc, nyc, ilevel+1);
        
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc);
        
        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf);
        
        relax_ini(uf_new,  wf_new, su, sw, ilevel, nxf, nyf);
        
        free_dmatrix(uc_new, 1, nxc, 1, nyc);
        free_dmatrix(wc_new, 1, nxc, 1, nyc);
        free_dmatrix(duc, 1, nxc, 1, nyc);
        free_dmatrix(dwc, 1, nxc, 1, nyc);
        free_dmatrix(uc_def, 1, nxc, 1, nyc);
        free_dmatrix(wc_def, 1, nxc, 1, nyc);
        free_dmatrix(uf_def, 1, nxf, 1, nyf);
        free_dmatrix(wf_def, 1, nxf, 1, nyf);
    }
   
}

void relax_ini(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt)
{
    extern int c_relax;
    extern double xright, dt, Cahn, SS,M;
    
    int i, j, iter;
    double ht2, a[4], f[2], det;
    
    ht2 = pow(xright/(double)nxt,2);
    
    for (iter=1; iter<=c_relax; iter++) {
        
        ijloopt {

            a[0] = 1.0/dt;
            
            a[1] = -M;
                   
            a[2] = 2.0/ht2 +SS/Cahn;
            if (j > 1)
                a[2] += 1/ ht2;
            if (j < nyt)
                a[2] += 1 / ht2;

         
            a[3] = 1.0;

          
            f[0] = su[i][j];
            
            
            f[1] = sw[i][j];
            if (i > 1)   f[1] += c_new[i-1][j]/ht2;
            else         f[1] += c_new[nxt][j]/ht2;
            
            if (i < nxt) f[1] += c_new[i+1][j]/ht2;
            else         f[1] += c_new[1][j]/ht2;
            
            if (j > 1)   f[1] += c_new[i][j-1]/ht2;

            
            if (j < nyt) f[1] += c_new[i][j+1]/ht2;
    
            
            det = a[0]*a[3] - a[1]*a[2];
            
            c_new[i][j] = (a[3]*f[0] - a[1]*f[1])/det;
            mu_new[i][j] = (-a[2]*f[0] + a[0]*f[1])/det;
        }
    }
    
}

void defect_ini(double **duc, double **dwc, double **uf_new, double **wf_new,
            double **suf, double **swf, int nxf, int nyf)
{
    double **ruf, **rwf, **ruc, **rwc, **rruf, **rrwf;
    
    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf/2, 1, nyf/2);
    rrwf = dmatrix(1, nxf/2, 1, nyf/2);
        
    nonL_ini(ruf, rwf, uf_new, wf_new, nxf, nyf);
    
    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf);
    
    restrict2(ruf, duc, rwf, dwc, nxf/2, nyf/2);
    
    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf/2, 1, nyf/2);
    free_dmatrix(rrwf, 1, nxf/2, 1, nyf/2);
}

void nonL_ini(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt)
{
    extern double xright, dt, SS, Cahn,M;
    
    int i, j;
    double ss, ht2, **lap_c, **lap_mu;
    
    lap_c = dmatrix(1, nxt, 1, nyt);
    
     laplace_ch(c_new, lap_c, nxt, nyt);


     ht2 = pow(xright/(double)nxt,2);
    
    ijloopt {

        ru[i][j] = c_new[i][j]/dt - M*mu_new[i][j];
        rw[i][j] = mu_new[i][j] -lap_c[i][j] +SS/Cahn*c_new[i][j];
    }
    
    free_dmatrix(lap_c, 1, nxt, 1, nyt);
  
}



/********************* 2 order phase field part ***********************/


void cahn(double **c_old,  double **cc_old, double **H, double **c_new)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu, **mu;
    
    int it_mg = 1, max_it_CH = 200;
    double resid = 1.0, tol = 1.0e-6;
    
    mat_copy(ct, c_old, 1, nx, 1, ny);
    
    source(c_old, cc_old, H, sc, smu);
    
    while (it_mg <= max_it_CH && resid > tol) {
        
        vcycle(c_new, mu, sc, smu, nx ,ny, 1);
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny);
        
        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg-1);
}

void source(double **c_old, double **cc_old, double **H, double **src_c, double **src_mu)
{
    extern int nx, ny;
    extern double dt, h, int_r, gam, SS,Cahn, **alp;
    
    int i, j;
   
    augmenc(c_old, nx, ny);

    ijloop {

        src_c[i][j] = (4.0*c_old[i][j] - cc_old[i][j])/(2.0*dt);
        
        src_mu[i][j] = -1/Cahn*(H[i][j] + alp[i][j])*int_r +SS/Cahn*(2.0*c_old[i][j]-cc_old[i][j]);
    }

       
}


void laplace_ch(double **a, double **lap_a, int nxt, int nyt)
{
    extern double xright;
    
    int i, j;
    double ht2, dadx_L, dadx_R, dady_B, dady_T;
    
    ht2 = pow(xright/(double)nxt,2);
    
    ijloopt {
        
        if (i > 1)
            dadx_L = a[i][j] - a[i-1][j];
        else
            dadx_L = a[i][j] - a[nxt][j];
        
        if (i < nxt)
            dadx_R = a[i+1][j] - a[i][j];
        else
            dadx_R = a[1][j] - a[i][j];
        
        if (j > 1)
            dady_B = a[i][j] - a[i][j-1];
        else
            dady_B = 0.0;
        
        if (j < nyt)
            dady_T = a[i][j+1] - a[i][j];
        else
            dady_T = 0.0;
        
        lap_a[i][j] = (dadx_R - dadx_L + dady_T - dady_B)/ht2;
    }
    
}

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel)
{
    extern int n_level;
    
    relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf);
    
   if (ilevel < n_level) {
        
        int nxc, nyc;
        double **uc_new, **wc_new, **duc, **dwc, 
               **uc_def, **wc_def, **uf_def, **wf_def;
        
        nxc = nxf/2, nyc = nyf/2;
        
        uc_new = dmatrix(1, nxc, 1, nyc);
        wc_new = dmatrix(1, nxc, 1, nyc);
        duc = dmatrix(1, nxc, 1, nyc);
        dwc = dmatrix(1, nxc, 1, nyc);
        uc_def = dmatrix(1, nxc, 1, nyc);
        wc_def = dmatrix(1, nxc, 1, nyc);
        uf_def = dmatrix(1, nxf, 1, nyf);
        wf_def = dmatrix(1, nxf, 1, nyf); 

        defect(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf); 

          zero_matrix(uc_def, 1, nxc, 1, nyc);
          zero_matrix(wc_def, 1, nxc, 1, nyc);
        
        vcycle(uc_def,  wc_def, duc, dwc, nxc, nyc, ilevel+1);
        
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc);
        
        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf);
        
        relax(uf_new,  wf_new, su, sw, ilevel, nxf, nyf);
        
        free_dmatrix(uc_new, 1, nxc, 1, nyc);
        free_dmatrix(wc_new, 1, nxc, 1, nyc);
        free_dmatrix(duc, 1, nxc, 1, nyc);
        free_dmatrix(dwc, 1, nxc, 1, nyc);
        free_dmatrix(uc_def, 1, nxc, 1, nyc);
        free_dmatrix(wc_def, 1, nxc, 1, nyc);
        free_dmatrix(uf_def, 1, nxf, 1, nyf);
        free_dmatrix(wf_def, 1, nxf, 1, nyf);
    }
   
}

void relax(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt)
{
    extern int c_relax;
    extern double xright, dt, Cahn, SS,M;
    
    int i, j, iter;
    double ht2, a[4], f[2], det;
    
    ht2 = pow(xright/(double)nxt,2);
    
    for (iter=1; iter<=c_relax; iter++) {
        
        ijloopt {

            a[0] = 3.0/(2.0*dt);
            
            a[1] = -M;

                   
            a[2] = 2.0/ht2 +SS/Cahn;

            if (j > 1)
                a[2] += 1 / ht2;
            if (j < nyt)
                a[2] += 1/ ht2;
         
            a[3] = 1.0;

          
            f[0] = su[i][j];
        
            
            f[1] = sw[i][j];
            if (i > 1)   f[1] += c_new[i-1][j]/ht2;
            else         f[1] += c_new[nxt][j]/ht2;
            
            if (i < nxt) f[1] += c_new[i+1][j]/ht2;
            else         f[1] += c_new[1][j]/ht2;
            
            if (j > 1)   f[1] += c_new[i][j-1]/ht2;
           
            if (j < nyt) f[1] += c_new[i][j+1]/ht2;
            
            det = a[0]*a[3] - a[1]*a[2];
            
            c_new[i][j] = (a[3]*f[0] - a[1]*f[1])/det;
            mu_new[i][j] = (-a[2]*f[0] + a[0]*f[1])/det;
        }
    }
    
}

void defect(double **duc, double **dwc, double **uf_new, double **wf_new,
            double **suf, double **swf, int nxf, int nyf)
{
    double **ruf, **rwf, **ruc, **rwc, **rruf, **rrwf;
    
    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf/2, 1, nyf/2);
    rrwf = dmatrix(1, nxf/2, 1, nyf/2);
        
    nonL(ruf, rwf, uf_new, wf_new, nxf, nyf);
    
    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf);
    
    restrict2(ruf, duc, rwf, dwc, nxf/2, nyf/2);
    
    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf/2, 1, nyf/2);
    free_dmatrix(rrwf, 1, nxf/2, 1, nyf/2);
}

void nonL(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt)
{
    extern double xright, dt, Cahn, SS,M;
    
    int i, j;
    double ss, ht2, **lap_c, **lap_mu;
    
    lap_c = dmatrix(1, nxt, 1, nyt);

    
     laplace_ch(c_new, lap_c, nxt, nyt);

     ht2 = pow(xright/(double)nxt,2);
    
    ijloopt {

        ru[i][j] = 3.0*c_new[i][j]/(2.0*dt) - M*mu_new[i][j];
        rw[i][j] = mu_new[i][j] - lap_c[i][j] + SS/Cahn*c_new[i][j];
    }
    
    free_dmatrix(lap_c, 1, nxt, 1, nyt);
}

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt)
{
    int i, j;
    
    ijloopt {
        uc[i][j] = 0.25*(uf[2*i-1][2*j-1] + uf[2*i-1][2*j] + uf[2*i][2*j-1] + uf[2*i][2*j]);
        vc[i][j] = 0.25*(vf[2*i-1][2*j-1] + vf[2*i-1][2*j] + vf[2*i][2*j-1] + vf[2*i][2*j]);
    }
    
}


void restrict1(double **vf, double **vc, int nxt, int nyt)
{
    int i, j;
    
    ijloopt {
        vc[i][j] = 0.25*(vf[2*i-1][2*j-1] + vf[2*i-1][2*j] + vf[2*i][2*j-1] + vf[2*i][2*j]);
    }
    
}

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt)
{
    int i, j;
    
    ijloopt {
        uf[2*i-1][2*j-1] = uf[2*i-1][2*j] = uf[2*i][2*j-1] = uf[2*i][2*j] = uc[i][j];
        vf[2*i-1][2*j-1] = vf[2*i-1][2*j] = vf[2*i][2*j-1] = vf[2*i][2*j] = vc[i][j];
    }
    
}

double error(double **c_old, double **c_new, int nxt, int nyt)
{
    double **r, res;
    
    r = dmatrix(1, nxt, 1, nyt);
    
    mat_sub(r, c_new, c_old, 1, nxt, 1, nyt);
    res = mat_max(r, 1, nxt, 1, nyt);
    
    free_dmatrix(r, 1, nxt, 1, nyt);
    
    return res;
}


double error2(double **c_old,double **c_new,  int nxt,int nyt) 
{
	
    int i,j;
	double **rr,res2,x=0.0;

    ijloopt {
		x= x+pow(c_new[i][j]-c_old[i][j],2);
	}

	res2=sqrt(x/(1.0*nxt*nyt));

	return res2;
}


/**********************************************/

double Ef(double **c, double **c2, double **c3, int nxt, int nyt)
{
    extern double Cahn, h, eta, theta, ksi, alpha1, alpha2, beta, delta, psis,CC;
    double r, res;
    int i,j;

    
      r = 0.0;

   ijloopt{ r = r + h*h*( 0.25*pow(c[i][j],2)*pow(c[i][j]-1.0,2) + 0.25*pow(c2[i][j],2)*pow(c2[i][j]-1.0,2) 
         + 0.25*pow(c3[i][j],2)*pow(c3[i][j]-1.0,2)  );  }

   res = r;
   //printf("res=%f\n",res);

      
    return res;
}
void lagrangeH(double **H , int nxt,int nyt)
{
    extern double h, volume;
    double r;
    int i,j;
    r =0.0;
    ijloopt{
        r = r+h*h*(H[i][j]);
    }
    r= r/volume;

    ijloopt{
        H[i][j] = H[i][j]-r;
    }
}



/*************** util ****************/
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;
    long i, nrow=nrh-nrl+1+NR_END, ncol=nch-ncl+1+NR_END;
    
    m=(double **) malloc((nrow)*sizeof(double*));
    memset(m, 0, (nrow) * sizeof(double *));
    m+=NR_END;
    m-=nrl;
    
    m[nrl]=(double *) malloc((nrow*ncol)*sizeof(double));
    memset(m[nrl], 0, (nrow * ncol) * sizeof(double));
    m[nrl]+=NR_END;
    m[nrl]-=ncl;
    
    for (i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
    
    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl]+ncl-NR_END);
    free(m+nrl-NR_END);
    
    return;
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++) {
        for (j=yl; j<=yr; j++) {
            a[i][j] = 0.0;
        }
    }
    
    return;
}

void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++) {
        for (j=yl; j<=yr; j++) {
            a[i][j] = b[i][j];
        }
    }
    
    return;
}

void mat_copy2(double **a, double **b, double **a2, double **b2,
               int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++) {
        for (j=yl; j<=yr; j++) {
            a[i][j] = b[i][j];
            a2[i][j] = b2[i][j];
        }
    }
    
    return;
}

void mat_add(double **a, double **b, double **c,
             int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++) {
        for (j=yl; j<=yr; j++) {
            a[i][j] = b[i][j]+c[i][j];
        }
    }
    
    return;
}

void mat_add2(double **a, double **b, double **c,
              double **a2, double **b2, double **c2,
              int xl, int xr, int yl, int yr)
{
    int i, j;
    
    for (i=xl; i<=xr; i++) {
        for (j=yl; j<=yr; j++) {
            a[i][j] = b[i][j]+c[i][j];
            a2[i][j] = b2[i][j]+c2[i][j];
        }
    }
    
    return;
}

void mat_sub(double **a, double **b, double **c,
             int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    
    for (i=nrl; i<=nrh; i++) {
        for (j=ncl; j<=nch; j++) {
            a[i][j] = b[i][j]-c[i][j];
        }
    }
    
    return;
}

void mat_sub2(double **a, double **b, double **c,
              double **a2, double **b2, double **c2,
              int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    
    for (i=nrl; i<=nrh; i++) {
        for (j=ncl; j<=nch; j++) {
            a[i][j] = b[i][j]-c[i][j];
            a2[i][j] = b2[i][j]-c2[i][j];
        }
    }
    
    return;
}

double mat_max(double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    double x = 0.0;
    
    for (i=nrl; i<=nrh; i++) {
        for (j=ncl; j<=nch; j++) {
            if (fabs(a[i][j]) > x)
                x = fabs(a[i][j]);
        }
    }
    
    return x;
}


void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    
    for (i=nrl; i<=nrh; i++) {
        for (j=ncl; j<=nch; j++)
            fprintf(fptr, "   %16.14f", a[i][j]);
        
        fprintf(fptr, "\n");
    }
    
    return;
}


void print_data(double **c, double **c2, double **c3)
{
    extern char bufferc[1000], bufferc2[1000], bufferc3[1000];
    int i, j;
   FILE *fc, *fc2, *fc3;
 
    fc = fopen(bufferc,"a");
    fc2 = fopen(bufferc2,"a");
    fc3 = fopen(bufferc3,"a");

   iloop {
        jloop {
 		  
           fprintf(fc, "  %16.14f", c[i][j]);
           fprintf(fc2, "  %16.14f", c2[i][j]);
           fprintf(fc3, "  %16.14f", c3[i][j]);

	   }

                   fprintf(fc, "\n");
                   fprintf(fc2, "\n");
                   fprintf(fc3, "\n");
    }
   
    fclose(fc);
    fclose(fc2);
    fclose(fc3);

  return;
}



