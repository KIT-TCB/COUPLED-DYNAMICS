/*
extern double xe[4*NNDIM][3], ze[4*NNDIM];  // common extchr
extern int ne;                              // common extchr
extern int lmax[MAXTYP];                    // common lmax
extern int izp[NNDIM];                      // common izp
extern double nel;                          // common izp
extern double qzero[MAXTYP], uhubb[MAXTYP]; // common mcharge
extern double skhtab[MAXTYP][MAXTYP][MAXTAB][10],
         skstab[MAXTYP][MAXTYP][MAXTAB][10],
	 skself[MAXTYP][3], dr[MAXTYP][MAXTYP]; // common sktab
extern int dim[MAXTYP][MAXTYP];             // common sktab
*/

// broyden.c
  void broyden(int niter, double alpha, int jtop, double *vecin, double *vecout, dftb_broyden_t *arrays);
// eglcao.c
  //extern int eglcao(int nn, double x[NNDIM][3], double *eel, int *miter, double qmat[NNDIM], int phase);
// fermi.c
  void fermi(int ndim, double *ev, double *occ, double *efermi, int nelectrons);
// gammamat.c
  void gammamatrix(int nat, dvec *rat, double **gammamat, double uhubb[DFTB_MAXTYPES], int *izp);
  double gam12(double r, double uhub1, double uhub2);
  double gamsub(double a, double b, double r, double rrc);
// mulliken.c
  void mulliken(int nn, double *qmat, double *qmulli, double *qtot, int ndim,
                double dacc, double *occ, double **a, double **overl, int *ind,
                int lmax[DFTB_MAXTYPES], int *izp);
  //void mulliken(int nn, double qmat[NNDIM], double qmulli[MDIM], double *qtot, int ndim, double dacc,
    //   double occ[MDIM], double a[MDIM][MDIM], double overl[MDIM][MDIM], int ind[NNDIM+1]);
// orthogonalize.c
  long orthogonalize(double **THamil, double **OverlF, double **THamilOr, long n, dftb_orthogo_t arrays);
// output.c
  void outeigenvectors(double **a, double *ev, int *ind, int nn, dftb_phase1_t dftb1);
  void outspec(int nn, int ndim, int *ind, double *ev, double *occ,
               double efermi, double *qmat, double *qmulli, dftb_t *dftb, dftb_phase1_t dftb1);
// shift.c
  //extern void hamilshift(int nat, double qmat[NNDIM], double gammamat[NNDIM][NNDIM], double shift[NNDIM]);
// skpar.c
  int skspar(int i, int j, double r2, double dd[13],
             int lmax[DFTB_MAXTYPES], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3],
             int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES]);
  int skhpar(int i, int j, double r2, double dd[13],
             int lmax[DFTB_MAXTYPES], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3],
             int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES]);
  double cubicspline(double f0, double f1, double f2, double x0, double x1,
       double xh, double hl, double dr);
  double spline5th(double f0, double f1, double f2, double x0, double x1, double x2,
       double xh, double dr, int mxind);
// slkode.c
  void slkmatrices(int i, int j, double (*xat)[3], double ham[LDIM][LDIM], double over[LDIM][LDIM],
                 int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
                 int *izp, tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void slkode(double dum[3], int i, int j, double em[LDIM][LDIM], int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
               int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES],
                 tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
               tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);

// slktrafo.c
  void skss(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void sksp(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void sksd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void skpp(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void skpd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void skdd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void selfs(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void selfp(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
  void selfd(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double[13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);

//gradient.c
void gammaall1(double r, double ui, double uj, double *dcdr);

