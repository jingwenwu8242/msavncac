#define gnx 128
#define gny 128
#define gnz 128

#define iloop for(i=1;i<=gnx;i++)
#define i0loop for(i=0;i<=gnx;i++)
#define jloop for(j=1;j<=gny;j++)
#define j0loop for(j=0;j<=gny;j++)
#define kloop for(k=1;k<=gnz;k++)
#define k0loop for(k=0;k<=gnz;k++)

#define ijloop iloop jloop
#define i0jloop i0loop jloop
#define ij0loop iloop j0loop

#define ijkloop iloop jloop kloop
#define i0jkloop i0loop jloop kloop
#define ij0kloop iloop j0loop kloop
#define ijk0loop iloop jloop k0loop

#define iloopt for(i=1;i<=nxt;i++)
#define i0loopt for(i=0;i<=nxt;i++)
#define jloopt for(j=1;j<=nyt;j++)
#define j0loopt for(j=0;j<=nyt;j++)
#define kloopt for(k=1;k<=nzt;k++)
#define k0loopt for(k=0;k<=nzt;k++)

#define ijloopt iloopt jloopt 
#define i0jloopt i0loopt jloopt 
#define ij0loopt iloopt j0loopt 

#define ijkloopt iloopt jloopt kloopt
#define i0jkloopt i0loopt jloopt kloopt
#define ij0kloopt iloopt j0loopt kloopt
#define ijk0loopt iloopt jloopt k0loopt


float ***cube(int xl, int xr, int yl, int yr, int zl, int zr);


void free_cube(float ***t, int xl, int xr, int yl, int yr, int zl, int zr);

void zero_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr);


 void print_cube(FILE *fptr,
              float ***a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void print_data(float ***c1, float ***c2,float ***c3,float ***c4,int kk);




void cube_copy(float ***a, float ***b, int xl, int xr, int yl, int yr, int zl, int zr);

void cube_copy2(float ***a, float ***b, float ***a2, float ***b2,
               int xl, int xr, int yl, int yr, int zl, int zr);

void cube_add(float ***a, float ***b, float ***c,
             int xl, int xr, int yl, int yr, int zl, int zr);

void cube_add2(float ***a, float ***b, float ***c,
              float ***a2, float ***b2, float ***c2,
              int xl, int xr, int yl, int yr, int zl, int zr);

void cube_sub(float ***a, float ***b, float ***c, 
              int xl, int xr, int yl, int yr, int zl, int zr);

void cube_sub2(float ***a, float ***b, float ***c, 
               float ***a2, float ***b2, float ***c2, 
               int xl, int xr, int yl, int yr, int zl, int zr);

float cube_max(float ***a, 
               int xl, int xr, int yl, int yr, int zl, int zr);

