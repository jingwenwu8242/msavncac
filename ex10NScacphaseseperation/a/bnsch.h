void initialization(double **u, double **v, double **c1, double **c2, double **c3, double **p);
void cal_den(double **c1, double **rho);
void cal_vis(double **c1, double **vis);

double Ef(double **c, double **c2, double **c3, int nxt, int nyt);

void lagrangeH(double **H , int nxt,int nyt);


/********* 1 order phase-field **********/
void full_step_ini(double **u, double **v,  double **c1,  double **c2, double **p,double **mu1,double **mu2,  double **nu, double **nv, double **nc1, double **nc2,  double **np);


void temp_uv_ini(double **c1,double **c2,double **c3, double **tu, double **tv, double **u, double **v,  double **worku, double **workv,double **adv_u,double **adv_v, double **fx, double **fy);


void relax_uv_ini(double **tu, double **tv, double **su, double **sv, int nx, int ny);

void Poisson_ini(double **tu, double **tv, double **p, double **np);

void source_uv_ini(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt);
void cahn_ini(double **c_old, double **P,double **adv_c, double **c_new,double **mu);

void source_ini(double **c_old, double **P, double **adv_c,double **src_c, double **src_mu);

void vcycle_ini(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel);

void relax_ini(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt);

void defect_ini(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf, double **swf, int nxf, int nyf);

void nonL_ini(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt);



/************ 2 order phase-field *******/
void full_step(double **u, double **v, double **ou, double **ov, double **c1, double **oc1, double **c2, double **oc2, double **p,double **omu1, double **mu1, double **omu2, double **mu2, double **nu, double **nv, double **nc1, double **nc2,  double **np);
void advection_step(double **u, double **v, double **c1, double **c2, double **c3,  double **adv_u, double **adv_v, double **adv_c1,double **adv_c2,double **adv_c3);
void sf_force(double **c, double **mu, double **fx, double **fy);
void temp_uv(double **c1,double **c2,double **c3,double **tu, double **tv, double **u, double **v, double **ou, double **ov, double **adv_u,double **adv_v,double **worku, double **workv, double **fx, double **fy);
void relax_uv(double **tu, double **tv, double **su, double **sv, int nx, int ny);
void Poisson(double **tu, double **tv, double **p, double **np);
void source_uv(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt);
void div_uv(double **tu, double **tv, double **divuv, int nxt, int nyt);
void MG_Poisson(double **p, double **f);
void augmenc(double **c, int nxt, int nyt);
void augmenuv(double **u, double **v);
void augmentuv(double **u, double **v);
void vcycle_uv(double **uf, double **ff, double wf, int nxf, int nyf, int ilevel);
void relax_p(double **p, double **f, double w, int ilevel, int nxt, int nyt);
void residual_den(double **r, double **u, double **f, double den, int nxt, int nyt);
void grad_p(double **p, double **dpdx, double **dpdy, int nxt, int nyt);
void prolong(double **u_coarse, double **u_fine, int nxt, int nyt);
void cahn(double **c_old, double **cc_old, double **P, double **adv_c1, double **c_new, double **mu);

void source(double **c_old, double **cc_old, double **P, double **adv_c1, double **src_c, double **src_mu) ;

void laplace_ch(double **a, double **lap_a, int nxt, int nyt);

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf, int ilevel);

void relax(double **c_new, double **mu_new, double **su, double **sw, int ilevel, int nxt, int nyt);

void defect(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf, double **swf, int nxf, int nyf);

void nonL(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt);

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt);

void restrict1(double **vf, double **vc, int nxt, int nyt);

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt);

double error(double **c_old, double **c_new, int nxt, int nyt);


