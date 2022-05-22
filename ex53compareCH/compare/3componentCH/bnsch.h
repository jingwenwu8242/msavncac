void initialization(double **c, double **c2, double **c3);

void augmenc(double **c, int nxt, int nyt);

/********* 1 order phase-field **********/

void cahn_ini(double **c_old, double **P, double **c_new);

void source_ini(double **c_old, double **P, double **src_c, double **src_mu);

void vcycle_ini(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel);

void relax_ini(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt);

void defect_ini(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf, double **swf, int nxf, int nyf);

void nonL_ini(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt);



/************ 2 order phase-field *******/
void cahn(double **c_old, double **cc_old, double **P, double **c_new);

void source(double **c_old, double **cc_old, double **P, double **src_c, double **src_mu);

void laplace_ch(double **a, double **lap_a, int nxt, int nyt);

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel);

void relax(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt);

void defect(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf, double **swf, int nxf, int nyf);

void nonL(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt);

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt);

void restrict1(double **vf, double **vc, int nxt, int nyt);

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt);

double error(double **c_old, double **c_new, int nxt, int nyt);
double error2(double **c_old,double **c_new,  int nxt,int nyt); 

double Ef(double **c, double **c2, double **c3, int nxt, int nyt);

