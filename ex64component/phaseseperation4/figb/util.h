#define gnx 128
#define gny nx

#define iloop for(i=1;i<=gnx;i++)
#define i0loop for(i=0;i<=gnx;i++)
#define jloop for(j=1;j<=gny;j++)
#define j0loop for(j=0;j<=gny;j++)
#define ijloop iloop jloop
#define i0jloop i0loop jloop
#define ij0loop iloop j0loop
#define iloopt for(i=1;i<=nxt;i++)
#define i0loopt for(i=0;i<=nxt;i++)
#define jloopt for(j=1;j<=nyt;j++)
#define j0loopt for(j=0;j<=nyt;j++)
#define ijloopt iloopt jloopt
#define i0jloopt i0loopt jloopt
#define ij0loopt iloopt j0loopt

double **dmatrix(long nrl, long nrh, long ncl, long nch);

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void zero_matrix(double **a, int xl, int xr, int yl, int yr);

void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr);

void mat_copy2(double **a, double **b, double **a2, double **b2,
               int xl, int xr, int yl, int yr);

void mat_add(double **a, double **b, double **c,
             int xl, int xr, int yl, int yr);

void mat_add2(double **a, double **b, double **c,
              double **a2, double **b2, double **c2,
              int xl, int xr, int yl, int yr);

void mat_sub(double **a, double **b, double **c,
             int nrl, int nrh, int ncl, int nch);

void mat_sub2(double **a, double **b, double **c,
              double **a2, double **b2, double **c2,
              int nrl, int nrh, int ncl, int nch);

double mat_max(double **a, int nrl, int nrh, int ncl, int nch);

void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch);

void print_data(double **c, double **c2, double **c3,double **c4);



