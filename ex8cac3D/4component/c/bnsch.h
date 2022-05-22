void initialization(float ***c, float ***c2, float ***c3,float ***c4);

void augmenc(float ***c, int nxt, int nyt,int nzt);

void lagrangeH(float ***H , int nxt,int nyt, int nzt);
float Ef(float ***c, float ***c2, float ***c3,float ***c4, int nxt, int nyt, int nzt);

/************* 1 order phase-field ***************/

void cahn_ini(float ***c_old, float ***P, float ***c_new);

void source_ini(float ***c_old, float ***P, float ***src_c, float ***src_mu);

void vcycle_ini(float ***uf_new, float ***wf_new, float ***su, float ***sw, int nxf, int nyf,int nzf, int ilevel);

void relax_ini(float ***c_new, float ***mu_new, float ***su, float ***sw, int ilevel, int nxt, int nyt, int nzt);

void defect_ini(float ***duc, float ***dwc, float ***uf_new, float ***wf_new, float ***suf, float ***swf, int nxf, int nyf, int nzt);

void nonL_ini(float ***ru, float ***rw, float ***c_new, float ***mu_new, int nxt, int nyt,int nzt);



/****************** 2 order phase-field **********/
void cahn(float ***c_old, float ***cc_old, float ***P, float ***c_new);

void source(float ***c_old, float ***cc_old, float ***P, float ***src_c, float ***src_mu);

void laplace_ch(float ***a, float ***lap_a, int nxt, int nyt, int nzt);

void vcycle(float ***uf_new, float ***wf_new, float ***su, float ***sw, int nxf, int nyf,int nzt, int ilevel);

void relax(float ***c_new, float ***mu_new, float ***su, float ***sw, int ilevel, int nxt, int nyt,int nzt);

void defect(float ***duc, float ***dwc, float ***uf_new, float ***wf_new, float ***suf, float ***swf, int nxf, int nyf, int nzt);

void nonL(float ***ru, float ***rw, float ***c_new, float ***mu_new, int nxt, int nyt, int nzt);

void restrict2(float ***uf, float ***uc, float ***vf, float ***vc, int nxt, int nyt, int nzt);

void restrict1(float ***vf, float ***vc, int nxt, int nyt,int nzt);

void prolong_ch(float ***uc, float ***uf, float ***vc, float ***vf, int nxt, int nyt, int nzt);

float error(float ***c_old, float ***c_new, int nxt, int nyt ,int nzt);



