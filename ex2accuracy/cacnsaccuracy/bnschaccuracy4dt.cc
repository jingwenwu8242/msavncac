/************************************************
   cac
   with variable density and variable viscosity
   phase seperation
   gravity=0,We,Re,rho=1
   boundary conditon:x:periodic;y:neumann

 well parameter
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

int nx, ny, n_level, c_relax, p_relax, count, Num;
double pi, xleft, xright, yleft, yright,volume, h, dt, gam, Cahn, Re,M, Fr,We, rho1, rho2, rho3,vis1, vis2, Ec, oEc, acount,
    int_r, orc, r, nr, int_q, oq, q, nq, EnEO, EnEM,gravity,
    **intu, **intv, **intmu1, **intmu2, **intmu3, **H1, **H2, **H3, **intc1, **intc2, **intc3, **alp, **adv_u, **adv_v, **adv_c1, **adv_c2, **adv_c3, **fx, **fy, **fx1, **fy1, **fx2, **fy2, **fx3, **fy3,
    **tu, **tv, **worku, **workv, **nworku, **nworkv, **workp, rho, **vis, **mc, SS, **ct, **sc, **smu, **mu1, **mu2, **mu3, **omu1, **omu2, **omu3, **c3, **oc3, **ic, **ic2, **ic3;
char bufferu[2000] = {0}, bufferv[2000] = {0}, bufferp[2000] = {0}, bufferc1[2000] = {0}, bufferc2[2000] = {0}, bufferc3[2000] = {0};
int main()
{

    int i, j, it, max_it, ns;
    double **ou, **u, **nu, **ov, **v, **nv, **p, **np, **c1, **nc1, **oc1, **c2, **nc2, **oc2;

    FILE *fu, *fv, *fp, *fc1, *fc2, *fc3, *myr, *myq,*auxiliary;

    /*****/

    nx = gnx;
    ny = gny;
    n_level = (int)(log(nx) / log(2) - 0.9); 

    c_relax = 5; 
    p_relax = 7;

    pi = 4.0 * atan(1.0);

    xleft = 0.0, xright = 1.0;
    yleft = 0.0, yright = 1.0;
    volume = (xright-xleft)*(yright-yleft);

    count = 0;

    /***********************/

    h = xright / (double)nx;
    dt = 0.01*pow(h,2)  ;
    dt = 4*dt;
    gam = 4.0*h/(4.0*sqrt(2.0)*atanh(0.9));        
    Cahn = pow(gam, 2); 
    Fr = 1;
    gravity=0;
    SS = 2;
    M = 0.001;
    max_it = 128/4; 
    ns = max_it ;
    rho1 =1;//
    rho2 =1.0;
    rho3=1.0; 
    rho = rho1;
    vis1 = 1.0; 
    vis2 = 1.0; 
    Re = 1;
    We = 1;
    Num = 3; 

    /***********************/

    printf("nx = %d, ny = %d\n", nx, ny);
    printf("n_level = %d\n", n_level);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);

    

    intc1 = dmatrix(0, nx + 1, 0, ny + 1);
    intc2 = dmatrix(0, nx + 1, 0, ny + 1);
    intc3 = dmatrix(0, nx + 1, 0, ny + 1);
    intu = dmatrix(-1, nx + 1, 0, ny + 1);
    intv = dmatrix(0, nx + 1, -1, ny + 1);
    intmu1 = dmatrix(0, nx + 1, 0, ny + 1);
    intmu2 = dmatrix(0, nx + 1, 0, ny + 1);
    intmu3 = dmatrix(0, nx + 1, 0, ny + 1);


    ic = dmatrix(0, nx+1, 0, ny+1);
    ic2 = dmatrix(0, nx+1, 0, ny+1);
    ic3 = dmatrix(0, nx+1, 0, ny+1);



    
    
    ou = dmatrix(-1, nx + 1, 0, ny + 1); 
    u = dmatrix(-1, nx + 1, 0, ny + 1);  
    nu = dmatrix(-1, nx + 1, 0, ny + 1); 

    ov = dmatrix(0, nx + 1, -1, ny + 1); 
    v = dmatrix(0, nx + 1, -1, ny + 1);  
    nv = dmatrix(0, nx + 1, -1, ny + 1); 

    p = dmatrix(0, nx + 1, 0, ny + 1);  
    np = dmatrix(0, nx + 1, 0, ny + 1); 

    adv_u = dmatrix(0, nx, 1, ny);
    adv_v = dmatrix(1, nx, 0, ny); 

    adv_c1 = dmatrix(1, nx, 1, ny); 
    adv_c2 = dmatrix(1, nx, 1, ny); 
    adv_c3 = dmatrix(1, nx, 1, ny); 

    
    fx = dmatrix(0, nx, 1, ny);
    fy = dmatrix(1, nx, 0, ny);
    fx1 = dmatrix(0, nx, 1, ny);
    fy1 = dmatrix(1, nx, 0, ny);
    fx2 = dmatrix(0, nx, 1, ny);
    fy2 = dmatrix(1, nx, 0, ny);
    fx3 = dmatrix(0, nx, 1, ny);
    fy3 = dmatrix(1, nx, 0, ny);

    
    tu = dmatrix(-1, nx + 1, 0, ny + 1);
    tv = dmatrix(0, nx + 1, -1, ny + 1);

    
    worku = dmatrix(0, nx, 1, ny);
    workv = dmatrix(1, nx, 0, ny);
    workp = dmatrix(1, nx, 1, ny); 
    nworku = dmatrix(0, nx, 1, ny);
    nworkv = dmatrix(1, nx, 0, ny);
    

    oc1 = dmatrix(0, nx + 1, 0, ny + 1); 
    c1 = dmatrix(0, nx + 1, 0, ny + 1);  
    nc1 = dmatrix(0, nx + 1, 0, ny + 1); 

    oc2 = dmatrix(0, nx + 1, 0, ny + 1); 
    c2 = dmatrix(0, nx + 1, 0, ny + 1);  
    nc2 = dmatrix(0, nx + 1, 0, ny + 1); 

    c3 = dmatrix(0, nx + 1, 0, ny + 1);  
    oc3 = dmatrix(0, nx + 1, 0, ny + 1); 

    H1 = dmatrix(0, nx + 1, 0, ny + 1);  
    H2 = dmatrix(0, nx + 1, 0, ny + 1);  
    H3 = dmatrix(0, nx + 1, 0, ny + 1);  
    alp = dmatrix(0, nx + 1, 0, ny + 1); 

    ct = dmatrix(1, nx, 1, ny); 
    sc = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny); 
    mu1 = dmatrix(0, nx + 1, 0, ny + 1);  
    mu2 = dmatrix(0, nx + 1, 0, ny + 1);  
    mu3 = dmatrix(0, nx + 1, 0, ny + 1);  
    omu1 = dmatrix(0, nx + 1, 0, ny + 1); 
    omu2 = dmatrix(0, nx + 1, 0, ny + 1); 
    omu3 = dmatrix(0, nx + 1, 0, ny + 1); 


    vis = dmatrix(0, nx + 1, 0, ny + 1); 



    zero_matrix(mu1, 1, nx, 1, ny);  
    zero_matrix(mu2, 1, nx, 1, ny);  
    zero_matrix(mu3, 1, nx, 1, ny);  
    zero_matrix(omu1, 1, nx, 1, ny); 
    zero_matrix(omu2, 1, nx, 1, ny); 
    zero_matrix(omu3, 1, nx, 1, ny); 

    sprintf(bufferu, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datau.m");

    sprintf(bufferv, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datav.m");

    sprintf(bufferp, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datap.m");

    sprintf(bufferc1, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datac1.m");
    sprintf(bufferc2, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datac2.m");
    sprintf(bufferc3, "/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/datac3.m");

    fu = fopen(bufferu, "w");
    fv = fopen(bufferv, "w");
    fp = fopen(bufferp, "w");
    fc1 = fopen(bufferc1, "w");
    fc2 = fopen(bufferc2, "w");
    fc3 = fopen(bufferc3, "w");
 
    fclose(fu);
    fclose(fv);
    fclose(fp);
    fclose(fc1);
    fclose(fc2);
    fclose(fc3);

    initialization(u, v, c1, c2, c3, p); 

    mat_copy(nu, u, 0, nx, 1, ny); 
    mat_copy(nv, v, 1, nx, 0, ny);
    mat_copy(ou, u, 0, nx, 1, ny);
    mat_copy(ov, v, 1, nx, 0, ny);

    mat_copy(oc1, c1, 1, nx, 1, ny); 
    mat_copy(oc2, c2, 1, nx, 1, ny); 

    mat_copy(nc1, c1, 1, nx, 1, ny); 
    mat_copy(nc2, c2, 1, nx, 1, ny);  


    Ec = Ef(c1, c2, c3, nx, ny); 
    r = Ec;
    orc = r;

    q = 1.0;
    oq = q; 
    auxiliary=fopen("/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/auq.m", "w");
    fprintf(auxiliary, "%16.12f \n", q);

    print_data(u, v, c1, c2, c3, p); 

    augmenc(c1, nx, ny);
    augmenc(c2, nx, ny);
    augmenc(c3, nx, ny);

  
    EnEO = 0.0;

    ijloop
    {

        EnEO = EnEO + h * h * (1/Cahn*0.25 * pow(c1[i][j], 2) * pow(c1[i][j] - 1.0, 2) + 1/Cahn*0.25 * pow(c2[i][j], 2) * pow(c2[i][j] - 1.0, 2) +1/Cahn* 0.25 * pow(c3[i][j], 2) * pow(c3[i][j] - 1.0, 2) 
        + 0.5  * (pow(c1[i + 1][j] - c1[i][j], 2) / (h * h) + pow(c1[i][j + 1] - c1[i][j], 2) / (h * h))
        + 0.5  * (pow(c2[i + 1][j] - c2[i][j], 2) / (h * h) + pow(c2[i][j + 1] - c2[i][j], 2) / (h * h)) 
        + 0.5  * (pow(c3[i + 1][j] - c3[i][j], 2) / (h * h) + pow(c3[i][j + 1] - c3[i][j], 2) / (h * h)));
    }
    i0jloop { EnEO = EnEO + We*0.5 * nu[i][j] * nu[i][j]; }

    ij0loop { EnEO = EnEO + We*0.5 * nv[i][j] * nv[i][j]; }
     EnEO = EnEO +1;

    myr = fopen("/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/Mene.m", "w");
    myq = fopen("/Users/wujingwen/research/fluid_mechanics/2msavncac/msav+ncac/cacnsaccuracy/data/accuracy4dt/Oene.m", "w");

    
    

    
  
  for (it = 1; it <= 1; it++)
    {

        cal_vis(oc1, vis);    
        augmenc(vis, nx, ny); //augement vis/
        


        full_step_ini(u, v, c1, c2,p,  mu1, mu2, nu, nv, nc1, nc2, np);

        

        acount = 0.0;

        ijloop
        {

            acount = acount + h * h * (H1[i][j] * r * (nc1[i][j] - c1[i][j]) / dt 
            + H2[i][j] * r * (nc2[i][j] - c2[i][j]) / dt 
            + H3[i][j] * r * ((1.0 - nc1[i][j] - nc2[i][j]) - (1.0 - c1[i][j] - c2[i][j])) / dt);
        }

        nr = r + dt * acount;

        r = nr;
        
        acount = 0.0;

        ijloop
        {

            acount = acount + h * h * (adv_c1[i][j] * mu1[i][j] + adv_c2[i][j] * mu2[i][j] + adv_c3[i][j] * mu3[i][j] 
            + (fx1[i][j] + fx2[i][j] + fx3[i][j]) * tu[i][j] 
            + (fy1[i][j] + fy2[i][j] + fy3[i][j]) * tv[i][j] 
            + We*adv_u[i][j] * tu[i][j] + We*adv_v[i][j] * tv[i][j]);
        }
        nq = q * acount * dt + q;
        q = nq;



        augmenc(nc1, nx, ny);
        augmenc(nc2, nx, ny);
        augmenc(c3, nx, ny); 

        
        EnEO = 0.0;

        ijloop
        {

            EnEO = EnEO + h * h * (1/Cahn*0.25 * pow(nc1[i][j], 2) * pow(nc1[i][j] - 1.0, 2) + 1/Cahn*0.25 * pow(nc2[i][j], 2) * pow(nc2[i][j] - 1.0, 2) + 1/Cahn*0.25 * pow(c3[i][j], 2) * pow(c3[i][j] - 1.0, 2) 
            + 0.5  * (pow(nc1[i + 1][j] - nc1[i][j], 2) / (h * h) + pow(nc1[i][j + 1] - nc1[i][j], 2) / (h * h)) 
            + 0.5  * (pow(nc2[i + 1][j] - nc2[i][j], 2) / (h * h) + pow(nc2[i][j + 1] - nc2[i][j], 2) / (h * h)) 
            + 0.5  * (pow(c3[i + 1][j] - c3[i][j], 2) / (h * h) + pow(c3[i][j + 1] - c3[i][j], 2) / (h * h)));
        }
        i0jloop { EnEO = EnEO + We*h * h * 0.5 * nu[i][j] * nu[i][j]; }

        ij0loop { EnEO = EnEO + We*h * h * 0.5 * nv[i][j] * nv[i][j]; }

         EnEO = EnEO +1;

        
        EnEM = 0.0;

        ijloop{

        ic[i][j] = 2.0*nc1[i][j] - oc1[i][j];
        ic2[i][j] = 2.0*nc2[i][j] - oc2[i][j];
        ic3[i][j] = 1.0-ic[i][j]-ic2[i][j];

    }
       augmenc(ic, nx, ny);
         augmenc(ic2, nx, ny);
         augmenc(ic3, nx, ny);

        ijloop
        {

            EnEM = EnEM + h*h*( 
          0.25*( pow(nc1[i+1][j]-nc1[i][j],2)/(h*h) + pow(nc1[i][j+1]-nc1[i][j],2)/(h*h) )
         + 0.25*( pow(nc2[i+1][j]-nc2[i][j],2)/(h*h) + pow(nc2[i][j+1]-nc2[i][j],2)/(h*h) )
         + 0.25*( pow(c3[i+1][j]-c3[i][j],2)/(h*h) + pow(c3[i][j+1]-c3[i][j],2)/(h*h) ) 

           + 0.25*( pow(ic[i+1][j]-ic[i][j],2)/(h*h) + pow(ic[i][j+1]-ic[i][j],2)/(h*h) )  
         + 0.25*( pow(ic2[i+1][j]-ic2[i][j],2)/(h*h) + pow(ic2[i][j+1]-ic2[i][j],2)/(h*h) )  
         + 0.25*( pow(ic3[i+1][j]-ic3[i][j],2)/(h*h) + pow(ic3[i][j+1]-ic3[i][j],2)/(h*h) ) 
        

         + 0.5*SS/Cahn*( pow(nc1[i][j]-oc1[i][j],2)  + pow(nc2[i][j]-oc2[i][j],2)  + pow(1.0-nc1[i][j]-nc2[i][j] - (1.0-oc1[i][j]-oc2[i][j]),2)   )  
         );
         
     }
       i0jloop { EnEM = EnEM + We*h * h * 0.25 * (pow(nu[i][j], 2) + pow((2 * nu[i][j] - u[i][j]), 2))+We*1/3*pow(dt,2)*pow(h,2)*pow(nworku[i][j],2); }
        ij0loop { EnEM = EnEM + We*h * h * 0.25 * (pow(nv[i][j], 2) + pow((2 * nv[i][j] - v[i][j]), 2))+We*1/3*pow(dt,2)*pow(h,2)*pow(nworkv[i][j],2); } 

         

        
        EnEM = EnEM + (3.0 * nr - r)*0.5/Cahn + (3.0* nq - q)*0.5;


        mat_copy(u, nu, 0, nx, 1, ny);
        mat_copy(v, nv, 1, nx, 0, ny); 
        mat_copy(c1, nc1, 1, nx, 1, ny); 
        mat_copy(c2, nc2, 1, nx, 1, ny); 
        mat_copy(p, np, 1, nx, 1, ny);   


        printf("it = %d\n", it);

        if (it % ns == 0){
        print_data(nu, nv, nc1, nc2, c3, np); 
        fprintf(myr, "%16.12f \n", EnEM);
        fprintf(myq, "%16.12f \n", EnEO);
        fprintf(auxiliary, "%16.12f \n", q);
        printf("print out counts %d\n", count);}
        
    }   
  

    

    for (it = 2; it <= max_it; it++)
    {

        cal_vis(oc1, vis);    
        augmenc(vis, nx, ny); //augement vis/

        full_step(u, v, ou, ov, c1, oc1, c2, oc2, p, omu1, mu1, omu2, mu2, nu, nv, nc1, nc2, np);

        

        acount = 0.0;

        ijloop
        {

            acount = acount + h * h * (H1[i][j] * int_r * (3.0 * nc1[i][j] - 4.0 * c1[i][j] + oc1[i][j]) / (2.0 * dt) 
            + H2[i][j] * int_r * (3.0 * nc2[i][j] - 4.0 * c2[i][j] + oc2[i][j]) / (2.0 * dt) 
            + H3[i][j] * int_r * (3.0 * (1.0 - nc1[i][j] - nc2[i][j]) - 4.0 * (1.0 - c1[i][j] - c2[i][j]) + (1.0 - oc1[i][j] - oc2[i][j])) / (2.0 * dt));
        }

        nr = (4.0 * r - orc + 2.0 * dt * acount) / 3.0;

        orc = r;
        r = nr;

        
        acount = 0.0;
   

        ijloop
        {

            acount = acount + h * h * (adv_c1[i][j] * mu1[i][j] + adv_c2[i][j] * mu2[i][j] + adv_c3[i][j] * (-mu1[i][j] - mu2[i][j]) 
            + (fx1[i][j] + fx2[i][j] + fx3[i][j]) * tu[i][j] + (fy1[i][j] + fy2[i][j] + fy3[i][j]) * tv[i][j] 
             +We*adv_u[i][j] * tu[i][j] + We*adv_v[i][j] * tv[i][j]);
        }
        nq = (int_q * acount * 2 * dt + 4 * q - oq) / 3.0;
        oq = q;
        q = nq;

   

        augmenc(nc1, nx, ny);
        augmenc(nc2, nx, ny);
        augmenc(c3, nx, ny); 

        
        EnEO = 0.0;

        ijloop
        {

            EnEO = EnEO + h * h * (0.25 * pow(nc1[i][j], 2) * pow(nc1[i][j] - 1.0, 2) + 0.25 * pow(nc2[i][j], 2) * pow(nc2[i][j] - 1.0, 2) + 0.25 * pow(c3[i][j], 2) * pow(c3[i][j] - 1.0, 2) 
            + 0.5 * Cahn * (pow(nc1[i + 1][j] - nc1[i][j], 2) / (h * h) + pow(nc1[i][j + 1] - nc1[i][j], 2) / (h * h)) 
            + 0.5 * Cahn * (pow(nc2[i + 1][j] - nc2[i][j], 2) / (h * h) + pow(nc2[i][j + 1] - nc2[i][j], 2) / (h * h)) 
            + 0.5 * Cahn * (pow(c3[i + 1][j] - c3[i][j], 2) / (h * h) + pow(c3[i][j + 1] - c3[i][j], 2) / (h * h)));
        }
        i0jloop { EnEO = EnEO + We*h * h * 0.5 * nu[i][j] * nu[i][j]; }

        ij0loop { EnEO = EnEO + We*h * h * 0.5 * nv[i][j] * nv[i][j]; }
        EnEO+=1;
        
        EnEM = 0.0;

        ijloop{

        ic[i][j] = 2.0*nc1[i][j] - oc1[i][j];
        ic2[i][j] = 2.0*nc2[i][j] - oc2[i][j];
        ic3[i][j] = 1.0-ic[i][j]-ic2[i][j];

    }
       augmenc(ic, nx, ny);
         augmenc(ic2, nx, ny);
         augmenc(ic3, nx, ny);

        ijloop
        {

            EnEM = EnEM + h*h*( 
          0.25*Cahn*( pow(nc1[i+1][j]-nc1[i][j],2)/(h*h) + pow(nc1[i][j+1]-nc1[i][j],2)/(h*h) )
         + 0.25*Cahn*( pow(nc2[i+1][j]-nc2[i][j],2)/(h*h) + pow(nc2[i][j+1]-nc2[i][j],2)/(h*h) )
         + 0.25*Cahn*( pow(c3[i+1][j]-c3[i][j],2)/(h*h) + pow(c3[i][j+1]-c3[i][j],2)/(h*h) ) 

           + 0.25*Cahn*( pow(ic[i+1][j]-ic[i][j],2)/(h*h) + pow(ic[i][j+1]-ic[i][j],2)/(h*h) )  
         + 0.25*Cahn*( pow(ic2[i+1][j]-ic2[i][j],2)/(h*h) + pow(ic2[i][j+1]-ic2[i][j],2)/(h*h) )  
         + 0.25*Cahn*( pow(ic3[i+1][j]-ic3[i][j],2)/(h*h) + pow(ic3[i][j+1]-ic3[i][j],2)/(h*h) ) 
        

         + 0.5*SS*( pow(nc1[i][j]-oc1[i][j],2)  + pow(nc2[i][j]-oc2[i][j],2)  + pow(1.0-nc1[i][j]-nc2[i][j] - (1.0-oc1[i][j]-oc2[i][j]),2)   )  
         );
         
     }
       i0jloop { EnEM = EnEM + We*h * h * 0.25 * (pow(nu[i][j], 2) + pow((2 * nu[i][j] - u[i][j]), 2))+We*1/3*pow(dt,2)*pow(h,2)*pow(nworku[i][j],2); }
        ij0loop { EnEM = EnEM + We*h * h * 0.25 * (pow(nv[i][j], 2) + pow((2 * nv[i][j] - v[i][j]), 2))+We*1/3*pow(dt,2)*pow(h,2)*pow(nworkv[i][j],2); } 

         

        
        EnEM = EnEM + (3.0 * nr - r)*0.5 + (3.0* nq - q)*0.5;


        mat_copy(ou, u, 0, nx, 1, ny);
        mat_copy(u, nu, 0, nx, 1, ny);
        mat_copy(ov, v, 1, nx, 0, ny);
        mat_copy(v, nv, 1, nx, 0, ny);

        mat_copy(omu1, mu1, 1, nx, 1, ny);
        mat_copy(omu2, mu2, 1, nx, 1, ny);
        mat_copy(omu3, mu3, 1, nx, 1, ny);

        mat_copy(oc1, c1, 1, nx, 1, ny); 
        mat_copy(oc2, c2, 1, nx, 1, ny); 
        mat_copy(c1, nc1, 1, nx, 1, ny); 
        mat_copy(c2, nc2, 1, nx, 1, ny); 
        mat_copy(p, np, 1, nx, 1, ny);

        printf("it = %d\n", it);

        if (it % ns == 0)
        {
            count++;
            print_data(nu, nv, nc1, nc2, c3, np); 
            fprintf(myr, "%16.12f \n", EnEM);
            fprintf(myq, "%16.12f \n", EnEO);
            fprintf(auxiliary, "%16.12f \n", q);

            printf("print out counts %d\n", count);
        }
    }
    return 0;
}

void initialization(double **u, double **v, double **c1, double **c2, double **c3, double **p)
{
    extern double xright, yright, h, gam, pi;

    int i, j;
    double x, y;
    
    ijloop
    {
        x = ((double)i - 0.5) * h;
        y = ((double)j - 0.5) * h;

         c1[i][j] = 0.33 + 0.01*cos(3*pi*x)+0.04*cos(5*pi*x);

        c2[i][j] = 0.33 + 0.01*cos(2*pi*x)+0.02*cos(4*pi*x);

        c3[i][j] = 1.0 - c1[i][j] - c2[i][j];
        p[i][j] = 0.0;
 
     /*    c1[i][j] =0.5+0.5*tanh(min(sqrt(pow(x-1,2)+pow(y-0.5,2))-0.22,y-0.5)/gam);
        c2[i][j]=0.5+0.5*tanh(-max(0.22-sqrt(pow(x-1,2)+pow(y-0.5,2)),y-0.5)/gam);
        c3[i][j]= 1-c1[i][j]-c2[i][j];
        p[i][j] = 0.0; */

        /*     c1[i][j] = c2[i][j] = c3[i][j] = 0;

        if (x >= 0.25 && x <= 0.5 && y >= 0.25 && y <= 0.58)
            c1[i][j] = 1.0;
        if (x > 0.5 && x <= 0.75 && y >= 0.25 && y <= 0.58)
            c2[i][j] = 1.0;
        if (x >= 0.25 && x <= 0.75 && y > 0.58 && y <= 0.75)
            c3[i][j] = 1.0;  //in this case ,we should add CC in r*/
/* 
            c1[i][j]=0.5+0.5*tanh((0.35-sqrt(pow(x-1,2)+pow(y-2.36,2)))/(0.5*gam));
            c2[i][j]=0.5+0.5*tanh((0.35-sqrt(pow(x-1,2)+pow(y-1.64,2)))/(0.5*gam));
            c3[i][j]= 1-c1[i][j]-c2[i][j];
            p[i][j]=0; */

    




    }
    zero_matrix(u, 0, nx, 1, ny);
    zero_matrix(v, 1, nx, 0, ny);
}


void cal_vis(double **c1, double **vis)
{
    extern double vis1, vis2;

    int i, j;

    ijloop
    {
        vis[i][j] = 1.0; //0.5*fabs(vis1*(1.0+c1[i][j])+vis2*(1.0-c1[i][j]))/vis1;
    }
}

/*************** 1order phase field *****************/
void full_step_ini(double **u, double **v, double **c1, double **c2, double **p, double **mu1, double **mu2, double **nu, double **nv, double **nc1, double **nc2, double **np)
{
    extern int nx, ny;
    extern double dt, **c3, **mu3, **adv_u, **adv_v, **adv_c1, **adv_c2, **adv_c3, **fx, **fy, **fx1, **fy1, **fx2, **fy2, **fx3, **fy3,
        **tu, **tv, **worku, **workv, **nworku, **nworkv, rho, q, int_r, r, Ec, **H1, **H2, **H3, **alp,rho1,rho2,rho3,gravity,We;

    int i, j;

    ijloop
    {
        c3[i][j] = 1.0 - c1[i][j] - c2[i][j];
    }

    advection_step(u, v, c1, c2, c3, adv_u, adv_v, adv_c1, adv_c2, adv_c3);
    

    
     sf_force(c1, mu1, fx1, fy1); 

    sf_force(c2, mu2, fx2, fy2); 

    sf_force(c3, mu3, fx3, fy3); 
    
    i0jloop { 
        fx[i][j] = fx1[i][j]/We + fx2[i][j]/We + fx3[i][j]/We; 
    
    }

    ij0loop { 
        fy[i][j] = fy1[i][j]/We + fy2[i][j]/We + fy3[i][j]/We;
    
     } 

    grad_p(p, worku, workv, nx, ny); 

    temp_uv_ini(c1,c2,c3,tu, tv, u, v, adv_u, adv_v, worku, workv, fx, fy); 
    augmentuv(tu, tv);

    Poisson_ini(tu, tv, p, np); 

    grad_p(np, nworku, nworkv, nx, ny); 

    i0jloop
    {
        nu[i][j] = tu[i][j] - dt * (nworku[i][j] - worku[i][j])/rho;
    }

    ij0loop
    {
        nv[i][j] = tv[i][j] - dt * (nworkv[i][j] - workv[i][j])/rho; 
    }  

    //get u_{i+1/2,j},v_{i,j+1/2}//

    ijloop
    {

        H1[i][j] = c1[i][j] * (c1[i][j] - 0.5) * (c1[i][j] - 1.0) / Ec;
        H2[i][j] = c2[i][j] * (c2[i][j] - 0.5) * (c2[i][j] - 1.0) / Ec;
        H3[i][j] = c3[i][j] * (c3[i][j] - 0.5) * (c3[i][j] - 1.0) / Ec;

    }
    lagrangeH(H1,nx,ny);
    lagrangeH(H2,nx,ny);
    lagrangeH(H3,nx,ny);

    ijloop{
          alp[i][j] = -(1.0/Num)*(H1[i][j] + H2[i][j] + H3[i][j] );

      }

    
    cahn_ini(c1, H1, adv_c1, nc1, mu1);

    
    cahn_ini(c2, H2, adv_c2, nc2, mu2);

    

    ijloop
    {
        c3[i][j] = 1.0 - nc1[i][j] - nc2[i][j];
        mu3[i][j] = -mu1[i][j] - mu2[i][j];
    }
}

void temp_uv_ini(double **c1,double **c2,double **c3, double **tu, double **tv, double **u, double **v, double **worku, double **workv, double **adv_u, double **adv_v, double **fx, double **fy)
{
    int i, j, it_mg = 1, max_it = 500;
    extern double h, dt, q,rho,rho1,rho2,rho3,gravity,Re;
    double residu = 1.0, residv = 1.0, tol = 1.0e-5, **soru, **sorv, **su, **sv;

    soru = dmatrix(0, nx, 1, ny);
    sorv = dmatrix(1, nx, 0, ny);

    su = dmatrix(0, nx, 1, ny);
    sv = dmatrix(1, nx, 0, ny);

    i0jloop
    {
        su[i][j] = u[i][j] / dt - q * adv_u[i][j] - worku[i][j]/rho - q * fx[i][j]/rho;
    } 

    ij0loop
    {
        sv[i][j] = v[i][j] / dt - q * adv_v[i][j] - workv[i][j]/rho - q * fy[i][j]/rho-gravity*(rho1*0.5*(c1[i][j+1]+c1[i][j])+rho2*0.5*(c2[i][j+1]+c2[i][j])+rho3*0.5*(c3[i][j+1]+c3[i][j])
    -rho)/rho;
    }

    while (it_mg <= max_it && residu >= tol && residv >= tol)
    {

        relax_uv_ini(tu, tv, su, sv, nx, ny);

        mat_sub(soru, soru, tu, 0, nx, 1, ny); 
        mat_sub(sorv, sorv, tv, 1, nx, 0, ny);
        residu = mat_max(soru, 0, nx, 1, ny);
        residv = mat_max(sorv, 1, nx, 0, ny);
        mat_copy(soru, tu, 0, nx, 1, ny);
        mat_copy(sorv, tv, 1, nx, 0, ny);

        it_mg++;
    }

    printf("Velocity1 iteration = %d ", it_mg - 1);

    free_dmatrix(soru, 0, nx, 1, ny);
    free_dmatrix(sorv, 1, nx, 0, ny);

    free_dmatrix(su, 0, nx, 1, ny);
    free_dmatrix(sv, 1, nx, 0, ny);
}
void relax_uv_ini(double **tu, double **tv, double **su, double **sv, int nx, int ny)
{
    extern double xright, dt, Re, h,rho;

    int i, j, iter;
    double h2, sorc, coef;

    h2 = pow(h, 2);

    for (iter = 1; iter <= 5; iter++)
    {

        augmentuv(tu, tv);

        
        i0jloop
        {

            sorc = su[i][j] + (tu[i + 1][j] + tu[i - 1][j] + tu[i][j + 1] + tu[i][j - 1]) / (Re * h2*rho); 
            coef = 1.0 / dt + 4.0 / (Re * h2*rho);

            tu[i][j] = sorc / coef;
        }

        ij0loop
        {

            sorc = sv[i][j] + (tv[i + 1][j] + tv[i - 1][j] + tv[i][j + 1] + tv[i][j - 1]) / (Re * h2*rho);
            coef = 1.0 / dt + 4.0 / (Re * h2*rho);

            tv[i][j] = sorc / coef;
        }
    }
}

void Poisson_ini(double **tu, double **tv, double **p, double **np)
{
    extern int nx, ny;
    extern double **workp;

    source_uv_ini(tu, tv, p, workp, nx, ny);

    MG_Poisson(np, workp);
}

void source_uv_ini(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt) 
{
    extern double dt, h,  lam,rho;

    int i, j;

    div_uv(tu, tv, divuv, nxt, nyt); 

    augmenc(p, nxt, nyt);

    ijloopt
    {
        divuv[i][j] = divuv[i][j] / dt + (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - 4.0 * p[i][j]) / (h * h*rho);
    }
}

void cahn_ini(double **c_old, double **P, double **adv_c, double **c_new, double **mu)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu;

    int it_mg = 1, max_it_CH = 200;
    double resid = 1.0, tol = 1.0e-6;

    mat_copy(ct, c_old, 1, nx, 1, ny); 

    source_ini(c_old, P, adv_c, sc, smu); 

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle_ini(c_new, mu, sc, smu, nx, ny, 1);
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny);

        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg - 1);
}

void source_ini(double **c_old, double **H, double **adv_c, double **src_c, double **src_mu) 
{
    extern int nx, ny;
    extern double dt, h, SS, Cahn, **alp, r, q;

    int i, j;

    ijloop
     {

       src_c[i][j] = c_old[i][j]/dt-q*adv_c[i][j];
        
        src_mu[i][j] = 1/Cahn*(H[i][j] + alp[i][j])*r -SS/Cahn*c_old[i][j];
    }  
}

void vcycle_ini(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel)
{
    extern int n_level;

    relax_ini(uf_new, wf_new, su, sw, ilevel, nxf, nyf); 

    if (ilevel < n_level)
    {

        int nxc, nyc;
        double **uc_new, **wc_new, **duc, **dwc,
            **uc_def, **wc_def, **uf_def, **wf_def;

        nxc = nxf / 2, nyc = nyf / 2;

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

        vcycle_ini(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1); 

        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc); 

        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf); 

        relax_ini(uf_new, wf_new, su, sw, ilevel, nxf, nyf); 

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

    ht2 = pow(xright / (double)nxt, 2);

    for (iter = 1; iter <= c_relax; iter++)
    {

        ijloopt
        {

             a[0] = 1.0/dt;
            
            a[1] = M;
                   
            a[2] = -2.0/ht2 -SS/Cahn;
            if (j > 1)
                a[2] += -1/ ht2;
            if (j < nyt)
                a[2] += -1 / ht2;


            a[3] = 1.0;

            f[0] = su[i][j];
        
           

            f[1] = sw[i][j];
            if (i > 1)   f[1] += -c_new[i-1][j]/ht2;
            else         f[1] += -c_new[nxt][j]/ht2;
            
            if (i < nxt) f[1] += -c_new[i+1][j]/ht2;
            else         f[1] += -c_new[1][j]/ht2;
            
            if (j > 1)   f[1] += -c_new[i][j-1]/ht2;

            
            if (j < nyt) f[1] += -c_new[i][j+1]/ht2; 
            

            det = a[0] * a[3] - a[1] * a[2];

            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}

void defect_ini(double **duc, double **dwc, double **uf_new, double **wf_new,
                double **suf, double **swf, int nxf, int nyf) 
{
    double **ruf, **rwf, **ruc, **rwc, **rruf, **rrwf;

    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf / 2, 1, nyf / 2);
    rrwf = dmatrix(1, nxf / 2, 1, nyf / 2);

    nonL_ini(ruf, rwf, uf_new, wf_new, nxf, nyf); 

    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf); 

    restrict2(ruf, duc, rwf, dwc, nxf / 2, nyf / 2); 

    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf / 2, 1, nyf / 2);
    free_dmatrix(rrwf, 1, nxf / 2, 1, nyf / 2);
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

        ru[i][j] = c_new[i][j]/dt + M*mu_new[i][j];
        rw[i][j] = mu_new[i][j] +lap_c[i][j] -SS/Cahn*c_new[i][j];
    }

    free_dmatrix(lap_c, 1, nxt, 1, nyt);
    
}

/********************* 2 order phase field part ***********************/
void full_step(double **u, double **v, double **ou, double **ov, double **c1, double **oc1, double **c2, double **oc2, double **p, double **omu1, double **mu1, double **omu2, double **mu2, double **nu, double **nv, double **nc1, double **nc2, double **np)
{
    extern int nx, ny;
    extern double dt, **c3, **oc3, **mu3, **omu3, **intu, **intv, **intc1, **intc2, **intc3, **intmu1, **intmu2, **intmu3,
        **adv_u, **adv_v, **adv_c1, **adv_c2, **adv_c3, **fx, **fy, **fx1, **fy1, **fx2, **fy2, **fx3, **fy3,
        **tu, **tv, **worku, **workv, **nworku, **nworkv, rho, int_q, q, oq, int_r, r, orc, Ec,
        **H1, **H2, **H3, **alp,rho1,rho2,rho3,gravity,We;

    int i, j;

    ijloop
    {
        c3[i][j] = 1.0 - c1[i][j] - c2[i][j];
        oc3[i][j] = 1.0 - oc1[i][j] - oc2[i][j];
    }

    int_q = 2.0 * q - oq;
    int_r = 2.0 * r - orc;
    
    i0jloop
    {
        intu[i][j] = 2.0 * u[i][j] - ou[i][j];
    }

    ij0loop
    {
        intv[i][j] = 2.0 * v[i][j] - ov[i][j];
    }

    ijloop
    {

        intc1[i][j] = 2.0 * c1[i][j] - oc1[i][j];
        intc2[i][j] = 2.0 * c2[i][j] - oc2[i][j];
        intc3[i][j] = 1.0 - intc1[i][j] - intc2[i][j];

        intmu1[i][j] = 2.0 * mu1[i][j] - omu1[i][j];
        intmu2[i][j] = 2.0 * mu2[i][j] - omu2[i][j];
        intmu3[i][j] = -intmu1[i][j] - intmu2[i][j];
    }
    

    advection_step(intu, intv, intc1, intc2, intc3, adv_u, adv_v, adv_c1, adv_c2, adv_c3);
    

    sf_force(intc1, intmu1, fx1, fy1); 

    sf_force(intc2, intmu2, fx2, fy2); 

    sf_force(intc3, intmu3, fx3, fy3); 
    

    i0jloop {
         fx[i][j] = fx1[i][j]/We + fx2[i][j]/We + fx3[i][j]/We; 
        
        }

    ij0loop { 
        fy[i][j] = fy1[i][j]/We + fy2[i][j]/We + fy3[i][j]/We; 
        
        } 

    grad_p(p, worku, workv, nx, ny); 

    temp_uv(intc1,intc2,intc3,tu, tv, u, v, ou, ov, adv_u, adv_v, worku, workv, fx, fy); 

    augmentuv(tu, tv);

    Poisson(tu, tv, p, np); 

    grad_p(np, nworku, nworkv, nx, ny); 

    i0jloop
    {
        nu[i][j] = tu[i][j] - 2.0 * dt * (nworku[i][j] - worku[i][j]) / (3.0*rho);
    }

    ij0loop
    {
        nv[i][j] = tv[i][j] - 2.0 * dt * (nworkv[i][j] - workv[i][j]) / (3.0*rho);
    }

    //get u_{i+1/2,j},v_{i,j+1/2}//

    Ec = Ef(intc1, intc2, intc3, nx, ny); 


    ijloop
    {

        H1[i][j] = intc1[i][j] * (intc1[i][j] - 0.5) * (intc1[i][j] - 1.0) / Ec;
        H2[i][j] = intc2[i][j] * (intc2[i][j] - 0.5) * (intc2[i][j] - 1.0) / Ec;
        H3[i][j] = intc3[i][j] * (intc3[i][j] - 0.5) * (intc3[i][j] - 1.0) / Ec;

       
    }
    lagrangeH(H1,nx,ny);
      lagrangeH(H2,nx,ny);
      lagrangeH(H3,nx,ny);

      ijloop{
          alp[i][j] = -(1.0/Num)*(H1[i][j] + H2[i][j] + H3[i][j] );

      }


    
    cahn(c1, oc1, H1, adv_c1, nc1, mu1);

    
    cahn(c2, oc2, H2, adv_c2, nc2, mu2);

    

    ijloop
    {
        c3[i][j] = 1.0 - nc1[i][j] - nc2[i][j];
        mu3[i][j] = -mu1[i][j] - mu2[i][j];
    }
}

void advection_step(double **u, double **v, double **c1, double **c2, double **c3, double **adv_u, double **adv_v, double **adv_c1, double **adv_c2, double **adv_c3)
{
    extern int nx, ny;
    extern double h, **ou, **ov;
    double a, b, d, **ux, **uy, **vx, **vy;

    int i, j, k;
    

    ux = dmatrix(-1, nx, 1, ny);
    uy = dmatrix(0, nx, 0, ny);
    vx = dmatrix(0, nx, 0, ny);
    vy = dmatrix(1, nx, -1, ny);

    augmenuv(u, v);
    augmenc(c1, nx, ny);
    augmenc(c2, nx, ny);
    augmenc(c3, nx, ny);

    

    i0jloop
    {
        if (u[i][j] > 0.0)
            adv_u[i][j] = u[i][j] * (u[i][j] - u[i - 1][j]) / h;
        else
            adv_u[i][j] = u[i][j] * (u[i + 1][j] - u[i][j]) / h;

        if (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j] > 0.0)
            adv_u[i][j] += 0.25 * (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j]) * (u[i][j] - u[i][j - 1]) / h;
        else
            adv_u[i][j] += 0.25 * (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j]) * (u[i][j + 1] - u[i][j]) / h;
    }

    ij0loop
    {
        if (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1] > 0.0)
            adv_v[i][j] = 0.25 * (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1]) * (v[i][j] - v[i - 1][j]) / h;
        else
            adv_v[i][j] = 0.25 * (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1]) * (v[i + 1][j] - v[i][j]) / h;

        if (v[i][j] > 0.0)
            adv_v[i][j] += v[i][j] * (v[i][j] - v[i][j - 1]) / h;
        else
            adv_v[i][j] += v[i][j] * (v[i][j + 1] - v[i][j]) / h;
    }

    ijloop
    { 
        adv_c1[i][j] = (u[i][j] * (c1[i + 1][j] + c1[i][j]) - u[i - 1][j] * (c1[i][j] + c1[i - 1][j])) / (2.0 * h) + (v[i][j] * (c1[i][j + 1] + c1[i][j]) - v[i][j - 1] * (c1[i][j] + c1[i][j - 1])) / (2.0 * h);

        adv_c2[i][j] = (u[i][j] * (c2[i + 1][j] + c2[i][j]) - u[i - 1][j] * (c2[i][j] + c2[i - 1][j])) / (2.0 * h) + (v[i][j] * (c2[i][j + 1] + c2[i][j]) - v[i][j - 1] * (c2[i][j] + c2[i][j - 1])) / (2.0 * h);

        adv_c3[i][j] = (u[i][j] * (c3[i + 1][j] + c3[i][j]) - u[i - 1][j] * (c3[i][j] + c3[i - 1][j])) / (2.0 * h) + (v[i][j] * (c3[i][j + 1] + c3[i][j]) - v[i][j - 1] * (c3[i][j] + c3[i][j - 1])) / (2.0 * h);
    }

    free_dmatrix(ux, -1, nx, 1, ny);
    free_dmatrix(uy, 0, nx, 0, ny);
    free_dmatrix(vx, 0, nx, 0, ny);
    free_dmatrix(vy, 1, nx, -1, ny);
}
void augmenuv(double **u, double **v)
{
    extern int nx, ny;

    int i, j;
    double aa;

    for (j = 1; j <= ny; j++)
    {
        aa = 0.5 * (u[0][j] + u[nx][j]);
        u[0][j] = u[nx][j] = aa;
        u[-1][j] = u[nx - 1][j];
        u[nx + 1][j] = u[1][j];
    }

    for (i = -1; i <= nx + 1; i++)
    {
        u[i][0] = -0.0 - u[i][1];
        u[i][ny + 1] = 0.0 - u[i][ny];
    }

    for (j = 0; j <= ny; j++)
    {
        v[0][j] = v[nx][j];
        v[nx + 1][j] = v[1][j];
    }

    for (i = 0; i <= nx + 1; i++)
    {
        v[i][0] = v[i][ny] = 0.0;
        v[i][-1] = -v[i][1];
        v[i][ny + 1] = -v[i][ny - 1];
    }
}

void augmentuv(double **u, double **v)
{
    extern int nx, ny;

    int i, j;
    double aa;

    for (j = 1; j <= ny; j++)
    {
        aa = 0.5 * (u[0][j] + u[nx][j]);
        u[0][j] = u[nx][j] = aa;
        u[-1][j] = u[nx - 1][j];
        u[nx + 1][j] = u[1][j];
    }

    for (i = -1; i <= nx + 1; i++)
    {
        u[i][0] = 0.0 - u[i][1];
        u[i][ny + 1] = 0.0 - u[i][ny];
    }

    for (j = 0; j <= ny; j++)
    {
        v[0][j] = v[nx][j];
        v[nx + 1][j] = v[1][j];
    }

    for (i = 0; i <= nx + 1; i++)
    {
        v[i][0] = v[i][ny] = 0.0;
        v[i][-1] = -v[i][1];
        v[i][ny + 1] = -v[i][ny - 1];
    }
}

void augmenc(double **c1, int nxt, int nyt) 
{
    int i, j;

    for (j = 1; j <= nyt; j++)
    {
        
        c1[0][j] = c1[nxt][j];
        c1[nxt + 1][j] = c1[1][j];
    }

    for (i = 0; i <= nxt + 1; i++)
    {
         c1[i][0] = c1[i][1];
        c1[i][nyt + 1] = c1[i][nyt];
        



    }
}
void sf_force(double **c1, double **mu, double **fx, double **fy) 
{
    extern int nx, ny;
    extern double h;

    int i, j;

    augmenc(c1, nx, ny);
    augmenc(mu, nx, ny);

    i0jloop
    {
        fx[i][j] = 0.5 * (c1[i + 1][j] + c1[i][j]) * (mu[i + 1][j] - mu[i][j]) / (h);
    }

    ij0loop
    {
        fy[i][j] = 0.5 * (c1[i][j + 1] + c1[i][j]) * (mu[i][j + 1] - mu[i][j]) / (h);
    }
}

void temp_uv(double **c1,double **c2,double **c3,double **tu, double **tv, double **u, double **v, double **ou, double **ov, double **adv_u, double **adv_v, double **worku, double **workv, double **fx, double **fy)
{
    int i, j, it_mg = 1, max_it = 500;
    extern double h, dt, int_q,rho,rho1,rho2,rho3,gravity;
    double residu = 1.0, residv = 1.0, tol = 1.0e-5, **soru, **sorv, **su, **sv;

    soru = dmatrix(0, nx, 1, ny);
    sorv = dmatrix(1, nx, 0, ny);

    su = dmatrix(0, nx, 1, ny);
    sv = dmatrix(1, nx, 0, ny);

    i0jloop
    {
        su[i][j] = (4.0 * u[i][j] - ou[i][j]) / (2.0 * dt) - int_q * adv_u[i][j] - worku[i][j] /rho- int_q * fx[i][j]/rho;
    } 

    ij0loop
    {
        sv[i][j] = (4.0 * v[i][j] - ov[i][j]) / (2.0 * dt) - int_q * adv_v[i][j] - workv[i][j]/rho - int_q * fy[i][j]/rho-gravity*(rho1*0.5*(c1[i][j+1]+c1[i][j])+rho2*0.5*(c2[i][j+1]+c2[i][j])+rho3*0.5*(c3[i][j+1]+c3[i][j])
    -rho)/rho;
    }

    while (it_mg <= max_it && residu >= tol && residv >= tol)
    {

        relax_uv(tu, tv, su, sv, nx, ny);

        mat_sub(soru, soru, tu, 0, nx, 1, ny); 
        mat_sub(sorv, sorv, tv, 1, nx, 0, ny);
        residu = mat_max(soru, 0, nx, 1, ny);
        residv = mat_max(sorv, 1, nx, 0, ny);
        mat_copy(soru, tu, 0, nx, 1, ny);
        mat_copy(sorv, tv, 1, nx, 0, ny);

        it_mg++;
    }

    printf("Velocity1 iteration = %d ", it_mg - 1);

    free_dmatrix(soru, 0, nx, 1, ny);
    free_dmatrix(sorv, 1, nx, 0, ny);

    free_dmatrix(su, 0, nx, 1, ny);
    free_dmatrix(sv, 1, nx, 0, ny);
}

void relax_uv(double **tu, double **tv, double **su, double **sv, int nx, int ny)
{
    extern double xright, dt, Re, h,rho;

    int i, j, iter;
    double h2, sorc, coef;

    h2 = pow(h, 2);

    for (iter = 1; iter <= 5; iter++)
    {

        augmentuv(tu, tv);

        
        i0jloop
        {

            sorc = su[i][j] + (tu[i + 1][j] + tu[i - 1][j] + tu[i][j + 1] + tu[i][j - 1]) / (Re * h2*rho); 
            coef = 3.0 / (2.0 * dt) + 4.0 / (Re * h2*rho);

            tu[i][j] = sorc / coef;
        }

        ij0loop
        {

            sorc = sv[i][j] + (tv[i + 1][j] + tv[i - 1][j] + tv[i][j + 1] + tv[i][j - 1]) / (Re * h2*rho);
            coef = 3.0 / (2.0 * dt) + 4.0 / (Re * h2*rho);

            tv[i][j] = sorc / coef;
        }
    }
}

void Poisson(double **tu, double **tv, double **p, double **np)
{
    extern int nx, ny;
    extern double **workp;

    source_uv(tu, tv, p, workp, nx, ny);

    MG_Poisson(np, workp);
}

void source_uv(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt) 
{
    extern double dt, h,  lam,rho;

    int i, j;

    div_uv(tu, tv, divuv, nxt, nyt); 

    augmenc(p, nxt, nyt);

    ijloopt
    {
        divuv[i][j] = 3.0  * divuv[i][j] / (2.0 * dt) + (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - 4.0 * p[i][j]) / (h * h*rho);
    }
}
void div_uv(double **tu, double **tv, double **divuv, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht;

    ht = xright / (double)nxt;

    ijloopt
    {
        divuv[i][j] = (tu[i][j] - tu[i - 1][j] + tv[i][j] - tv[i][j - 1]) / ht;
    }
}

void MG_Poisson(double **p, double **f)
{
    extern int nx, ny;
    extern double rho;

    int it_mg = 1, max_it = 100;
    double resid = 1.0, resid2 = 10.0, tol = 1.0e-5, **sor;

    sor = dmatrix(1, nx, 1, ny);

    mat_copy(sor, p, 1, nx, 1, ny); 

    while (it_mg <= max_it && resid >= tol)
    {

        vcycle_uv(p, f, rho, nx, ny, 1);

        pressure_update(p); 

        mat_sub(sor, sor, p, 1, nx, 1, ny);
        resid = mat_max(sor, 1, nx, 1, ny);
        mat_copy(sor, p, 1, nx, 1, ny);

        if (resid > resid2)
            it_mg = max_it;
        else
            resid2 = resid;

        it_mg++;
    }
    printf("Mac pressure iteration = %d   residual = %16.14f\n", it_mg - 1, resid);

    free_dmatrix(sor, 1, nx, 1, ny);

    return;
}

void vcycle_uv(double **uf, double **ff, double wf, int nxf, int nyf, int ilevel)
{
    extern int n_level;

    relax_p(uf, ff, wf, ilevel, nxf, nyf); 

    if (ilevel < n_level)
    {

        int nxc, nyc;
        double **rf, **fc, **uc;

        nxc = nxf / 2, nyc = nyf / 2; 

        rf = dmatrix(1, nxf, 1, nyf);         
        fc = dmatrix(1, nxc, 1, nyc);         
        uc = dmatrix(1, nxc, 1, nyc);         
    

        residual_den(rf, uf, ff, wf, nxf, nyf); 

        restrict1(rf, fc, nxc, nyc); 

     

        zero_matrix(uc, 1, nxc, 1, nyc);

        vcycle_uv(uc, fc, wf, nxc, nyc, ilevel + 1); 

        prolong(uc, rf, nxc, nyc); 

        mat_add(uf, uf, rf, 1, nxf, 1, nyf); 

        relax_p(uf, ff, wf, ilevel, nxf, nyf); 

        free_dmatrix(rf, 1, nxf, 1, nyf);
        free_dmatrix(fc, 1, nxc, 1, nyc);
        free_dmatrix(uc, 1, nxc, 1, nyc);
    }
}
void relax_p(double **p, double **f, double w, int ilevel, int nxt, int nyt)
{
    extern int ny, p_relax;
    extern double xright, Fr;

    int i, j, iter;
    double ht, ht2, a[4], sorc, coef;

    ht = xright / (double)nxt;
    ht2 = pow(ht, 2);

    for (iter = 1; iter <= p_relax; iter++)
    {

        ijloopt
        {
            a[0] = 1.0 / w;
            a[1] = 1.0 / w;
            a[2] = 1.0 / w;
            a[3] = 1.0 / w;

            sorc = f[i][j];
            coef = -(a[0] + a[1] ) / ht2;
            
            if (i == 1)
            {
                sorc -= (a[0] * p[i + 1][j] + a[1] * p[nxt][j]) / ht2;
            }
            else if (i == nxt)
            {
                sorc -= (a[0] * p[1][j] + a[1] * p[i - 1][j]) / ht2;
            }
            else
            {
                sorc -= (a[0] * p[i + 1][j] + a[1] * p[i - 1][j]) / ht2;
            }

            /
            if (j == 1)
            {
                sorc -= a[2]*p[i][j+1]/ht2;
                coef -= a[2]/ht2;

            }
            else if (j == nyt)
            {
                sorc -= a[3]*p[i][j-1]/ht2;
                coef -= a[3]/ht2;

            }
            else
            {
                sorc -= (a[2] * p[i][j + 1] + a[3] * p[i][j - 1]) / ht2;
                coef -= (a[2] + a[3]) / ht2;
            }

            p[i][j] = sorc / coef;
        }
    }
}

void residual_den(double **r, double **u, double **f, double den, int nxt, int nyt)
{
    int i, j;
    double **dpdx, **dpdy;

    dpdx = dmatrix(0, nxt, 1, nyt);
    dpdy = dmatrix(1, nxt, 0, nyt);

    grad_p(u, dpdx, dpdy, nxt, nyt); 

    i0jloopt
    {
        dpdx[i][j] = dpdx[i][j] / den;
    }

    ij0loopt
    {
        dpdy[i][j] = dpdy[i][j] / den;
    } //(dpdx,dpdy)=1/rho*grad_p//

    div_uv(dpdx, dpdy, r, nxt, nyt);  
    mat_sub(r, f, r, 1, nxt, 1, nyt); 

    free_dmatrix(dpdx, 0, nxt, 1, nyt);
    free_dmatrix(dpdy, 1, nxt, 0, nyt);
}

void grad_p(double **p, double **dpdx, double **dpdy, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht;

    ht = xright / (double)nxt;

    i0jloopt
    {
        if (i == 0)
            dpdx[0][j] = (p[1][j] - p[nxt][j]) / ht;
        else if (i == nxt)
            dpdx[nxt][j] = (p[1][j] - p[nxt][j]) / ht;
        else
            dpdx[i][j] = (p[i + 1][j] - p[i][j]) / ht;
    }

    ij0loopt
    {
        if (j == 0)
        {
            dpdy[i][0] = 0;
        }
        else if (j == nyt)
        {
            dpdy[i][nyt] = 0;
        }

        else
            dpdy[i][j] = (p[i][j + 1] - p[i][j]) / ht;
    }
}

void prolong(double **u_coarse, double **u_fine, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        u_fine[2 * i - 1][2 * j - 1] = u_fine[2 * i - 1][2 * j] =
            u_fine[2 * i][2 * j - 1] = u_fine[2 * i][2 * j] = u_coarse[i][j];
    }
}

void cahn(double **c_old, double **cc_old, double **P, double **adv_c1, double **c_new, double **mu)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu;

    int it_mg = 1, max_it_CH = 200; 
    double resid = 1.0, tol = 1.0e-6;

    mat_copy(ct, c_old, 1, nx, 1, ny); 

    source(c_old, cc_old, P, adv_c1, sc, smu); 

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle(c_new, mu, sc, smu, nx, ny, 1); 
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny); 

        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg - 1);
}

void source(double **c_old, double **cc_old, double **H, double **adv_c, double **src_c, double **src_mu) 
{
    extern int nx, ny;
    extern double dt, h, int_r, gam, SS, **alp, int_q;

    int i, j;

    augmenc(c_old, nx, ny);
    augmenc(cc_old, nx, ny);

    ijloop
    {

        src_c[i][j] = (4.0*c_old[i][j] - cc_old[i][j])/(2.0*dt)-int_q * adv_c[i][j];
        
        src_mu[i][j] = 1/Cahn*(H[i][j] + alp[i][j])*int_r -SS/Cahn*(2.0*c_old[i][j]-cc_old[i][j]);
        

    }  
}

void laplace_ch(double **a, double **lap_a, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht2, dadx_L, dadx_R, dady_B, dady_T;

    ht2 = pow(xright / (double)nxt, 2);

    ijloopt
    {

        if (i > 1)
            dadx_L = a[i][j] - a[i - 1][j];
        else
            dadx_L = a[i][j] - a[nxt][j];

        if (i < nxt)
            dadx_R = a[i + 1][j] - a[i][j];
        else
            dadx_R = a[1][j] - a[i][j];



        if (j > 1)
            dady_B = a[i][j] - a[i][j - 1];
        else
            dady_B = 0.0;

        if (j < nyt)
            dady_T = a[i][j + 1] - a[i][j];
        else
            dady_T = 0.0;

        lap_a[i][j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2;
    }
}

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel)
{
    extern int n_level;

    relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf); 

    if (ilevel < n_level)
    {

        int nxc, nyc;
        double **uc_new, **wc_new, **duc, **dwc,
            **uc_def, **wc_def, **uf_def, **wf_def;

        nxc = nxf / 2, nyc = nyf / 2;

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

        vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1); 

        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc); 

        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf); 

        relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf); 

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

    ht2 = pow(xright / (double)nxt, 2);

    for (iter = 1; iter <= c_relax; iter++)
    {

        ijloopt
        {

             a[0] = 3.0/(2.0*dt);           
            a[1] = M;
        a[2] = -2.0/ht2 -SS/Cahn;

            if (j > 1)
                a[2] += -1 / ht2;
            if (j < nyt)
                a[2] += -1/ ht2;
         
            a[3] = 1.0;
            f[0] = su[i][j];
            f[1] = sw[i][j];
            if (i > 1)   f[1] += -c_new[i-1][j]/ht2;
            else         f[1] += -c_new[nxt][j]/ht2;
            
            if (i < nxt) f[1] += -c_new[i+1][j]/ht2;
            else         f[1] += -c_new[1][j]/ht2;
            
            if (j > 1)   f[1] += -c_new[i][j-1]/ht2;
           
            if (j < nyt) f[1] += -c_new[i][j+1]/ht2;
             

            det = a[0] * a[3] - a[1] * a[2];

            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}

void defect(double **duc, double **dwc, double **uf_new, double **wf_new,
            double **suf, double **swf, int nxf, int nyf) 
{
    double **ruf, **rwf, **ruc, **rwc, **rruf, **rrwf;

    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf / 2, 1, nyf / 2);
    rrwf = dmatrix(1, nxf / 2, 1, nyf / 2);

    nonL(ruf, rwf, uf_new, wf_new, nxf, nyf); 

    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf); 

    restrict2(ruf, duc, rwf, dwc, nxf / 2, nyf / 2); 

    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf / 2, 1, nyf / 2);
    free_dmatrix(rrwf, 1, nxf / 2, 1, nyf / 2);
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

        ru[i][j] = 3.0*c_new[i][j]/(2.0*dt) + M*mu_new[i][j];
        rw[i][j] = mu_new[i][j] + lap_c[i][j] - SS/Cahn*c_new[i][j];
    } 
    
   
    
    free_dmatrix(lap_c, 1, nxt, 1, nyt);

}

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        uc[i][j] = 0.25 * (uf[2 * i - 1][2 * j - 1] + uf[2 * i - 1][2 * j] + uf[2 * i][2 * j - 1] + uf[2 * i][2 * j]);
        vc[i][j] = 0.25 * (vf[2 * i - 1][2 * j - 1] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i][2 * j]);
    }
}

void restrict1(double **vf, double **vc, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        vc[i][j] = 0.25 * (vf[2 * i - 1][2 * j - 1] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i][2 * j]);
    }
}

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        uf[2 * i - 1][2 * j - 1] = uf[2 * i - 1][2 * j] = uf[2 * i][2 * j - 1] = uf[2 * i][2 * j] = uc[i][j];
        vf[2 * i - 1][2 * j - 1] = vf[2 * i - 1][2 * j] = vf[2 * i][2 * j - 1] = vf[2 * i][2 * j] = vc[i][j];
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

/**********************************************/

double Ef(double **c1, double **c2, double **c3, int nxt, int nyt)
{
    extern double Cahn, h, eta, theta, ksi, alpha1, alpha2, beta, delta, psis;
    double r, res;
    int i, j;

    r = 0.0;

    ijloopt { r = r + h * h * (0.25 * pow(c1[i][j], 2) * pow(c1[i][j] - 1.0, 2) + 0.25 * pow(c2[i][j], 2) * pow(c2[i][j] - 1.0, 2) + 0.25 * pow(c3[i][j], 2) * pow(c3[i][j] - 1.0, 2)); }

    res = r;
    

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
    long i, nrow = nrh - nrl + 1 + NR_END, ncol = nch - ncl + 1 + NR_END;
    

    m = (double **)malloc((nrow) * sizeof(double *));
    memset(m, 0, (nrow) * sizeof(double *));
    m += NR_END;
    m -= nrl;

    m[nrl] = (double *)malloc((nrow * ncol) * sizeof(double));
    memset(m[nrl], 0, (nrow * ncol) * sizeof(double));
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - NR_END);
    free(m + nrl - NR_END);

    return;
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = 0.0;
        }
    }

    return;
}

void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
        }
    }

    return;
}

void mat_copy2(double **a, double **b, double **a2, double **b2,
               int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
            a2[i][j] = b2[i][j];
        }
    }

    return;
}

void mat_add(double **a, double **b, double **c1,
             int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j] + c1[i][j];
        }
    }

    return;
}

void mat_add2(double **a, double **b, double **c1,
              double **a2, double **b2, double **c2,
              int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j] + c1[i][j];
            a2[i][j] = b2[i][j] + c2[i][j];
        }
    }

    return;
}

void mat_sub(double **a, double **b, double **c1,
             int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - c1[i][j];
        }
    }

    return;
}

void mat_sub2(double **a, double **b, double **c1,
              double **a2, double **b2, double **c2,
              int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - c1[i][j];
            a2[i][j] = b2[i][j] - c2[i][j];
        }
    }

    return;
}

double mat_max(double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    double x = 0.0;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            if (fabs(a[i][j]) > x)
                x = fabs(a[i][j]);
        }
    }

    return x;
}

void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
            fprintf(fptr, "   %16.14f", a[i][j]);

        fprintf(fptr, "\n");
    }

    return;
}

void print_data(double **u, double **v, double **c1, double **c2, double **c3, double **p)
{
    extern char bufferu[2000], bufferv[2000], bufferp[2000], bufferc1[2000], bufferc2[2000], bufferc3[2000];
    int i, j;
    FILE *fu, *fv, *fp, *fc1, *fc2, *fc3;
    fu = fopen(bufferu, "a");
    fv = fopen(bufferv, "a");
    fp = fopen(bufferp, "a");
    fc1 = fopen(bufferc1, "a");
    fc2 = fopen(bufferc2, "a");
    fc3 = fopen(bufferc3, "a");
    if (fu == NULL)
    {
        printf("open failed %s", bufferu);
        exit(1);
    }

    if (fv == NULL)
    {
        printf("open failed %s", bufferv);
        exit(1);
    }

    if (fp == NULL)
    {
        printf("open failed %s", bufferp);
        exit(1);
    }

    if (fc1 == NULL)
    {
        printf("open failed %s", bufferc1);
        exit(1);
    }
    if (fc2 == NULL)
    {
        printf("open failed %s", bufferc2);
        exit(1);
    }
    if (fc3 == NULL)
    {
        printf("open failed %s", bufferc3);
        exit(1);
    }

    iloop
    {
        jloop
        {

            fprintf(fu, "  %16.14f", 0.5 * (u[i][j] + u[i - 1][j]));
            fprintf(fv, "  %16.14f", 0.5 * (v[i][j] + v[i][j - 1]));
            fprintf(fp, "  %16.14f", p[i][j]);
            fprintf(fc1, "  %16.14f", c1[i][j]);
            fprintf(fc2, "  %16.14f", c2[i][j]);
            fprintf(fc3, "  %16.14f", c3[i][j]);
        }

        fprintf(fu, "\n");
        fprintf(fv, "\n");
        fprintf(fp, "\n");
        fprintf(fc1, "\n");
        fprintf(fc2, "\n");
        fprintf(fc3, "\n");
    }

    fclose(fu);
    fclose(fv);
    fclose(fp);
    fclose(fc1);
    fclose(fc2);
    fclose(fc3);

    return;
}

void pressure_update(double **a)
{
    extern int nx, ny;

    int i, j;
    double ave = 0.0;

    ijloop
    {
        ave = ave + a[i][j];
    }
    ave /= (nx + 0.0) * (ny + 0.0);

    ijloop
    {
        a[i][j] -= ave;
    }

    return;
}
