#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<f2c.h>
//#include"/usr/local/intel/mkl/10.0.3.020/include/mkl_lapack.h"
//#include"/usr/local/src/CLAPACK-3.1.1/INCLUDE/clapack.h"
//#include"/usr/include/atlas_enum.h"
//#include"/usr/include/clapack.h"
#include"maxima.h"
#include"external-and-prototypes.h"

/*     SUBROUTINE EGLCAO */
/*     ================= */

/*     Copyright 1997 by Peter Blaudeck, Dirk Porezag, Michael Haugk, */
/*                       Joachim Elsner */
/*     Bo Song for CMD-QM/MM based on Fragments, April 2007 */

/* ********************************************************************* */

/*     PROGRAM CHARACTERISTICS */
/*     ----------------------- */

/* eglcao calculates energy and gradient for dylcao as shown by Seifert. */
/* The determination of the occupation numbers has been changed to be */
/* also valid for metallic systems. */

/* PARAMETERS: */
/* nn      i  number of atoms */
/* x       r  coordinates (n3) */
/* eel     r  electronic energy */
/* miter   i  number of scf-iterations performed */
/* qmat    r */
/* phase   i  1 or 2 */

/* ********************************************************************* */

// lapack routine(s)
static long dsygv(long itype, char jobz, char uplo, long n, double *a, long lda,
             double *b, long ldb, double *w, double *work, long lwork)
{
  extern void dsygv_(long *, char *, char *, long *, double *, long *, double *,
                long *, double *, double *, long *, long *);
  long info;
  dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
  return info;
}

static long orthogonalize(double THamil[MFRG][MFRG], double OverlF[MFRG][MFRG],
             double THamilOr[MFRG][MFRG], long n)
{
  extern void dpotrf_(char *, long *, double *, long *, long *);
  extern void dpotri_(char *, long *, double *, long *, long *);
  extern void dsyevr_(char *, char *, char *, long *, double *, long *,
                double *, double *, long *, long *, double *, long *, double *,
		double *, long *, double *, double *, long *, double *, long *, long *);
  long i, j, k, m, info, eval_no, lwork, liwork;
  double *tij, *sij, *evec, *work, *iwork, *eval, *issupz;
  double abstol=1.0e-8;
  char itype='U';
  char jobz='V';
  char range='A';
  char uplo='U';
  double void_double;
  long void_long;

  // Allocate arrays
  tij = (double *) malloc(n*n*sizeof(double));
  sij = (double *) malloc(n*n*sizeof(double));
  evec = (double *) malloc(n*n*sizeof(double));
  lwork = 26*n;
  work = (double *) malloc(26*n*26*n*sizeof(double));
  liwork = 10*n;
  iwork = (double *) malloc(10*n*10*n*sizeof(double));
  eval = (double *) malloc(n*sizeof(double));
  issupz = (double *) malloc(2*n*sizeof(double));

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      sij[i+j*n] = OverlF[i][j];
      tij[i+j*n] = THamil[i][j];
    }
  
  // Choleski-factorize the overlap matrix
  dpotrf_(&itype, &n, sij, &n, &info);
  if (info) return info;

  // Invert the overlap matrix
  dpotri_(&itype, &n, sij, &n, &info);
  if (info) return info;

  for (i=0; i<n; i++)
    for (j=0; j<i; j++)
      sij[i+j*n] = sij[j+i*n];

  // Diagonalize the inverse of overlap matrix
  // evec - eigenvectors, eval - eigenvalues
  dsyevr_(&jobz, &range, &uplo, &n, sij, &n,
          &void_double, &void_double, &void_long, &void_long, &abstol, &eval_no, eval,
	  evec, &n, issupz, work, &lwork, iwork, &liwork, &info);
  if (info) return info;

  // sij = evec * sqrt(diag) * evec-transp
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      sij[i+j*n] = 0.0;
      for (k=0; k<n; k++)
        sij[i+j*n] += evec[i+k*n] * evec[j+k*n] * sqrt(eval[k]);
    }

  // THamilOr = sij * tij * sij
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      THamilOr[i][j] = 0.0;
      for (k=0; k<n; k++)
        for (m=0; m<n; m++)
	  THamilOr[i][j] += sij[i+k*n] * tij[k+m*n] * sij[m+j*n];
    }

  free(tij);
  free(sij);
  free(evec);
  free(work);
  free(iwork);
  free(eval);
  free(issupz);
  return 0;
}

int eglcao(int nn, double x[NNDIM][3], double *eel, int *miter, double qmat[NNDIM], int phase)
{
  const double scftol = 1.e-9;
  int ind[NNDIM+1], indj, indk, indj1, indk1;
  double ecoul, efermi, eelold, eext;
  double au[LDIM][LDIM], bu[LDIM][LDIM];
  static double a[MDIM][MDIM], b[MDIM][MDIM], gammamat[NNDIM][NNDIM];
  double ev[MDIM], occ[MDIM], qmulli[MDIM], qmold[NNDIM];
  static double hamil[MDIM][MDIM], overl[MDIM][MDIM];
  double shift[NNDIM];
  // for fragment analysis
  static double Taf[MDIM][MDIM], THamil[MDIM][MDIM];
  double tij[MFRG][MFRG], sij[MFRG][MFRG], THamilOrtho[MFRG][MFRG];
  //
  double shiftE[NNDIM];
  char str[80]; //originally 'char' !!!
  int nfrag, nfgs[MFRG], nfgn[MFRG], ifo, iao, totfrag;
  int ifrg[MFRG][MFAO], naco[MFRG], inf[MFRG];
  static double OverlF[MDIM][MDIM];   // overlap matrix of fragment orbitals by BoS
  //
  int i, j, k, m, n, li, lj, iao2fo, indfrg;
  int nmaofo, maxiter, niter, ii, jj, kk, ll, jfo, jao;
  double racc, dacc, r2;
  double almix, qtot, aux[3*MDIM];
  // lapack
  double *a_trans=NULL, *b_trans=NULL;
//  long la_itype=1;
//  char la_jobz='V', la_uplo='L';
//  long la_lwork, ier, ndim;
  long ier, ndim;
  FILE *f76, *f77, *f78, *f79;

  // machine accuracy
  racc = 1.0;
  while ((1.0 + racc) > 1.0)
    racc /= 2.0;
  racc *= 2.0;
       
  // calculation of indices for matrices H and S 
  ind[0] = 0;
  for (j=0; j<nn; j++)
    ind[j+1] = ind[j] + lmax[izp[j]-1]* lmax[izp[j]-1];

  // actual dimension of matrix
  ndim = ind[nn];
  if (ndim > MDIM) {
    printf(" eglcao: ndim = %ld > %d\n", ndim, MDIM);
    exit(-1);
  }

  // setup of charge-independent part of H and S
  for (j=0; j<nn; j++)
    for (k=0; k<=j; k++) {
      slkmatrices(j, k, x, au, bu);
      for (n=0; n<ind[k+1]-ind[k]; n++)
        for (m=0; m<ind[j+1]-ind[j]; m++) {
          hamil[ind[j]+m][ind[k]+n] = au[m][n];
          hamil[ind[k]+n][ind[j]+m] = au[m][n];      
          overl[ind[j]+m][ind[k]+n] = bu[m][n];
          overl[ind[k]+n][ind[j]+m] = bu[m][n];
	}
    }

  // Hamiltonian
  //printf("Hamiltonian\n");
  //for (i=0; i<10; i++) {
  //  for (j=0; j<10; j++) printf("%9.5f", hamil[i][j]);
  //  printf("\n");
  //}

  // Overlap
  //printf("Overlap\n");
  //for (i=0; i<10; i++) {
  //  for (j=0; j<10; j++) printf("%9.5f", overl[i][j]);
  //  printf("\n");
  //}

  // external charges
  for (j=0; j<nn; j++) {
    shiftE[j] = 0.0;
    for (k=0; k<ne; k++) {
      r2 = 0.0;
      for (i=0; i<3; i++)
        r2 += (x[j][i] - xe[k][i])*(x[j][i] - xe[k][i]);
      shiftE[j] += ze[k] / sqrt(r2);
    }
  }

  // calculation of complex within the Fragment Orbital Approach
  // will be done outside the SCC cycle
  if (phase == 2) {
    // add charge-dependent terms (Hubbard and consec. extcharges)
    gammamatrix(nn, x, gammamat);
    hamilshift(nn, qmat, gammamat, shift);
    for (i=0; i<nn; i++)
      shift[i] -= shiftE[i];

    // update the Hamiltonian matrix
    for (i=0; i<nn; i++)
      for (li=0; li<lmax[izp[i]-1]*lmax[izp[i]-1]; li++)
        for (j=0; j<=i; j++)
          for (lj=0; lj<lmax[izp[j]-1]*lmax[izp[j]-1]; lj++) {
            hamil[ind[i]+li][ind[j]+lj] += 0.5*overl[ind[i]+li][ind[j]+lj]*(shift[i]+shift[j]);
	    hamil[ind[j]+lj][ind[i]+li] = hamil[ind[i]+li][ind[j]+lj];
	  }

    // transform the H-matrix from AOs to FOs (Bo Song 06/03/07)
    printf("Fragment Analysis\n");
    f76 = fopen("FRAG.DEF", "r");
    fscanf(f76, "%s\n", str);
    fscanf(f76, "%d", &iao2fo);
    fscanf(f76, "%d", &nfrag);
    totfrag = 0;
    for (i=0; i<nfrag; i++) {
      fscanf(f76, "%d %d %d", nfgs+i, nfgn+i, naco+i);
      for (j=0; j<naco[i]; j++) {
        fscanf(f76, "%d", &(ifrg[i][j]));
	totfrag++;
      }
    }
    printf("totfrag = %d\n", totfrag);
    fclose(f76);

    nmaofo = ndim;
    for (ifo=0; ifo<nmaofo; ifo++)
      for (iao=0; iao<nmaofo; iao++)
        Taf[ifo][iao] = 0.0;
    
    f77 = fopen("AO2FO.DEF", "r");
    indfrg = 0;
    for (ii=0; ii<nfrag; ii++) {
      inf[ii] = indfrg; // orbital index of frag_ii
      for (i=nfgs[ii]-1; i<nfgn[ii]; i++)
      for (li=0; li<lmax[izp[i]-1]*lmax[izp[i]-1]; li++) {
        indfrg++;
	for (j=nfgs[ii]-1; j<nfgn[ii]; j++)
	for (lj=0; lj<lmax[izp[j]-1]*lmax[izp[j]-1]; lj++)
	  fscanf(f77, "%lf", &(Taf[ind[i]+li][ind[j]+lj]));
      }
    }
    fclose(f77);

    kk = 0;
    for (i=0; i<nfrag; i++)
    for (ii=0; ii<naco[i]; ii++) {
      ifo = ifrg[i][ii] + inf[i] - 1;
      ll = 0;
      for (j=0; j<nfrag; j++)
      for (jj=0; jj<naco[j]; jj++) {
        jfo = ifrg[j][jj] + inf[j] - 1;
	if (ifo <= jfo) {
	  THamil[ifo][jfo] = 0.0;
	  OverlF[ifo][jfo] = 0.0;
	  for (iao=0; iao<nmaofo; iao++)
	  for (jao=0; jao<nmaofo; jao++) {
	    THamil[ifo][jfo] += Taf[ifo][iao] * hamil[iao][jao] * Taf[jfo][jao];
	    OverlF[ifo][jfo] += Taf[ifo][iao] * overl[iao][jao] * Taf[jfo][jao];
	  }
	  tij[kk][ll] = THamil[ifo][jfo];
	  sij[kk][ll] = OverlF[ifo][jfo];
	  THamil[jfo][ifo] = THamil[ifo][jfo];
	  OverlF[jfo][ifo] = OverlF[ifo][jfo];
	  tij[ll][kk] = tij[kk][ll];
	  sij[ll][kk] = sij[kk][ll];
	}
	ll++;
      }
      kk++;
    }

    // orthogonalize the Hamiltonian
    ier = -512;
    ier = orthogonalize(tij, sij, THamilOrtho, totfrag);

    f78 = fopen("FRAGanalysis.DAT", "w");
    fprintf(f78, "  FOs   FOs        H-elements         Overlap-M   \n");
    fprintf(f78, "==================================================\n");
    f79 = fopen("FRAGanalysis.ortho.DAT", "w");
    fprintf(f79, " FO FO      H-elements  \n");
    fprintf(f79, "========================\n");
    for (i=0; i<totfrag; i++)
      for (ii=0; ii<naco[i]; ii++) {
        ifo = ifrg[i][ii] + inf[i] - 1;
        for (j=0; j<nfrag; j++)
          for (jj=0; jj<naco[j]; jj++) {
	    jfo = ifrg[j][jj] + inf[j] - 1;
	    fprintf(f78, "%3d%3d%3d%3d %16.12f %s %16.12f\n",
	                i, ifrg[i][ii], j, ifrg[j][jj],
			THamil[ifo][jfo]*27.2116, "eV", OverlF[ifo][jfo]);
	  }
      }
    for (i=0; i<totfrag; i++)
      for (j=0; j<totfrag; j++)
	fprintf(f79, "%3d%3d %16.12f Ha %16.12f eV\n",
	            i, j, THamilOrtho[i][j], THamilOrtho[i][j]*27.2116);
    fclose(f78);
    fclose(f79);

    outspec(nn, ndim, ev, occ, 0.0, qmat);
  
  } else {
    // calculation of individual fragment -- SCC

    // setup for SCC cycle
    almix = 0.2;
    maxiter = 70;
    eelold = 1.e10;
    a_trans = (double *) malloc(ndim * ndim * sizeof (double));
    b_trans = (double *) malloc(ndim * ndim * sizeof (double));

    // SCC cycle starts here
    for (niter=0; niter<maxiter; *miter=++niter) {
      // save old charges
      for (i=0; i<nn; i++)
        qmold[i] = qmat[i];
      
      // charge-independent part of H and S
      for (j=0; j<nn; j++) {
        indj = ind[j];
	indj1 = ind[j+1];
	for (k=0; k<nn; k++) {
	  indk = ind[k]; 
	  indk1 = ind[k+1];
	  for (n=0; n<indk1-indk; n++)
	    for (m=0; m<indj1-indj; m++) {
	      a[indj+m][indk+n] = hamil[indj+m][indk+n];
	      b[indj+m][indk+n] = overl[indj+m][indk+n];
	    }
	}
      }

      // zero shift
      for (i=0; i<nn; i++)
        shift[i] = 0.0;

      // add charge-dependent terms (Hubbard and consec. extcharges)
      if (*miter == 0)
        gammamatrix(nn, x, gammamat);
      hamilshift(nn, qmat, gammamat, shift);
      //
      //printf("shift\n");
      //for (i=0; i<10; i++) printf("%2d %12.8f%12.8f\n", i+1, shift[i], shiftE[i]);

      for (i=0; i<nn; i++)
        shift[i] -= shiftE[i];

      for (i=0; i<nn; i++)
        for (li=0; li<lmax[izp[i]-1]*lmax[izp[i]-1]; li++)
          for (j=0; j<=i; j++)
            for (lj=0; lj<lmax[izp[j]-1]*lmax[izp[j]-1]; lj++) 
	      a[ind[i]+li][ind[j]+lj] += 0.5*overl[ind[i]+li][ind[j]+lj]*(shift[i]+shift[j]);
      
      // A before ewevge
      //printf("A before ewevge\n");
      //for (i=0; i<10; i++) {
      //  for (j=0; j<10; j++) printf("%9.5f", a[i][j]);
      //  printf("\n");
      //}
      
      // transpose the arrays a and b
      for (j=0; j<ndim; j++)
        for (i=0; i<ndim; i++) {
	  a_trans[j*ndim+i] = a[i][j];
	  b_trans[j*ndim+i] = b[i][j];
	}

      ier = -512;
      ier = dsygv(1, 'V', 'L', ndim, a_trans, ndim, b_trans, ndim, ev, aux, 3*ndim);
      for (j=0; j<ndim; j++)
        for (i=0; i<ndim; i++)
	  a[i][j] = a_trans[j*ndim+i];

      //printf("\n");
      //printf("dsygv: ier = %ld\n", ier);
      //printf("\n");

      // machine accuracy
      dacc = 4*racc;

      // calculate occupation (occ) and Fermi energy (efermi),
      fermi(ndim, ev, occ, &efermi);

      // sum of occupied eigenvalues
      *eel = 0.0;
      for (i=0; i<ndim && occ[i]>dacc; i++)
        *eel += occ[i]*ev[i];

      // determine Mulliken charges, charge of the whole system and the mulliken
      mulliken(nn, qmat, qmulli, &qtot, ndim, dacc, occ, a, overl, ind);

      // complete calculation of electronic energy
      // charge-dependent contribution
      // warning: this will only lead to the right result if convergence has been reached
      ecoul = eext = 0.0;
      for (i=0; i<nn; i++) {
        ecoul += shift[i] * (qmat[i] + qzero[izp[i]-1]);
	eext += shiftE[i] * (qzero[izp[i]-1] - qmat[i]);
      }
      *eel += -0.5*ecoul + 0.5*eext;
      // remark: eel containts shiftE already via ev,
      // shift also contains -shiftE, i.e. ecoul also
      // contains contributions from EXT

      // print energy
      printf("iter: %d, E= %f, %f, %f\n", niter, *eel, ecoul, eext);

      // check convergence
      if (fabs(*eel-eelold) < scftol)
        break;
      eelold = *eel;

      // Broyden mixing
      broyden(niter, almix, nn, qmold, qmat);
      for (i=0; i<nn; i++)
        qmat[i] = qmold[i];

    } // end SCC cycle
    free(a_trans);
    free(b_trans);

    outspec(nn, ndim, ev, occ, efermi, qmat);
    outeigenvectors(a, ev, ind, nn);
  }

  return 0;
}
