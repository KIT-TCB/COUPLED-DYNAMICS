#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

#define QM_CHARGE(x) (-dftb1.qmat[(x)] + dftb->qzero1[dftb1.izp[(x)]])

/*      PROGRAM DYLCAO */
/*     ================ */

/*     Copyright 1991 by Peter Blaudeck, Dirk Porezag */

/* ********************************************************************* */

/*     PROGRAM CHARACTERISTICS */
/*     ----------------------- */

/* DYLCAO calculates the dynamics of various systems */
/* within a two-centre SETB formalism */

/* ********************************************************************* */
/*

     PHASE 1 FOR CHARGE TRANSFER -- CALC OF FRAGMENTS

*/
/* ********************************************************************* */


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

/* NDIM = NORB !!! */

int run_dftb1(charge_transfer_t *ct, dftb_t *dftb, int ibase) // i - nucleobase to be calculated
// int eglcao(int nn, double x[NNDIM][3], double *eel, int *miter, double qmat[NNDIM], int phase)
{
  const double scftol = 1.e-7;
  /* OLD VALUE WAS: const double scftol = 1.e-9;
   * we do not need it most likely,
   * because we do not compute the forces.
   * the loose criterion seems sufficient to obtain
   * converged orbital coefficients and orbital energies,
   * which are the only quantities that we need
   */
  const double almix = 0.2;
  const int maxiter = 70;
  
  int indj, indk, indj1, indk1;
  double eel, ecoul, efermi, eelold, eext;
  //
  int i, j, k, m, n, li, lj, niter;
  //int nmaofo, ii, jj, kk, ll, jfo, jao;
  //double r2;
  double qtot, sum;
  // lapack
  long ier; // ndim;

  int nn, ne;
  dvec *x, bond;
  dftb_phase1_t dftb1;

  // printf("Phase 1, site %d\n", ibase+1);

  dftb1 = dftb->phase1[ibase];
  nn = dftb1.nn;
  ne = dftb1.ne;
  x = dftb1.x;

  /* LENGTHS MUST BE CONVERTED TO BOHR! (5.292d-11 m) */

/*
  // write out the coordinates
  for (n=0; n<nn; n++) {
    printf("ATOM %d %d %f %f %f\n", n+1, dftb1.izp[n]+1, x[n][XX]*0.529177249, x[n][YY]*0.529177249, x[n][ZZ]*0.529177249);
  }

  // write out the external charges
  for (n=0; n<ne; n++) {
    printf("XTCH %f %f %f %f\n", dftb1.xe[n][XX]*0.529177249, dftb1.xe[n][YY]*0.529177249, dftb1.xe[n][ZZ]*0.529177249, dftb1.ze[n]);
  }
*/

  // set the initial charges
  if(ct->do_scc[ibase] < 2)  //with this option qmat(t) would be qmat(t - delta_t), otherwise it should be qzero before scc
  for (i=0; i<nn; i++)
    dftb1.qmat[i] = dftb->qzero1[dftb1.izp[i]];
  
  // initial setup
  eel = 0.0;

  // printf("  icycle iter niter   e(total)\n");
  // printf("====================================\n");

  // setup of charge-independent part of H and S
  for (j=0; j<nn; j++)
    for (k=0; k<=j; k++) {
      slkmatrices(j, k, x, dftb1.au, dftb1.bu, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
      for (n=0; n<dftb1.ind[k+1]-dftb1.ind[k]; n++)
        for (m=0; m<dftb1.ind[j+1]-dftb1.ind[j]; m++) {
          dftb1.hamil[dftb1.ind[j]+m][dftb1.ind[k]+n] = dftb1.au[m][n];
          dftb1.hamil[dftb1.ind[k]+n][dftb1.ind[j]+m] = dftb1.au[m][n];      
          dftb1.overl[dftb1.ind[j]+m][dftb1.ind[k]+n] = dftb1.bu[m][n];
          dftb1.overl[dftb1.ind[k]+n][dftb1.ind[j]+m] = dftb1.bu[m][n];
	}
    }

  /*
  printf("Hamiltonian matrix:\n");
  for (j=0; j<dftb1.ndim; j++) {
    for (k=0; k<dftb1.ndim; k++) printf("%9.5f", dftb1.hamil[j][k]);
    printf("\n");
  }
  printf("Overlap matrix:\n");
  for (j=0; j<dftb1.ndim; j++) {
    for (k=0; k<dftb1.ndim; k++) printf("%9.5f", dftb1.overl[j][k]);
    printf("\n");
  }
  // */

  // external charges
  if (ct->qmmm < 3)
    for (j=0; j<nn; j++) {
      dftb1.shiftE[j] = diprod(x[j],ct->efield); 
      for (k=0; k<ne; k++) {
        dvec_sub(x[j], dftb1.xe[k], bond);
        dftb1.shiftE[j] += dftb1.ze[k] / (dnorm(bond) * ct->esp_scaling_factor);  // apply the ESP scaling factor!
      }
      //printf("  shiftE %2d %f\n", j+1, dftb1.shiftE[j]);
    }

  // setup for SCC cycle

  eelold = 1.e10;

  // SCC cycle starts here
  for (niter=0; niter<maxiter; niter++) {

    /* calculate PME if desired */
    if (ct->qmmm == 3) {
      do_pme_for_dftb_part2(ct, dftb, ibase);
      /* save the computed values to shiftE */
      for (j=0; j<nn; j++){
        dftb1.shiftE[j] = diprod(x[j],ct->efield); 
        dftb1.shiftE[j] += dftb1.pot[j] / ct->esp_scaling_factor;
      }
    }

    // save old charges
    for (i=0; i<nn; i++)
      dftb1.qmold[i] = dftb1.qmat[i];
    
    /*
    printf("qmat - start SCF:\n");
    for (i=0; i<nn; i++)
      printf("%3d %f\n", i+1, dftb1.qmat[i]);
    // */

    // charge-independent part of H and S
    for (j=0; j<nn; j++) {
      indj = dftb1.ind[j];
      indj1 = dftb1.ind[j+1];
      for (k=0; k<nn; k++) {
        indk = dftb1.ind[k]; 
        indk1 = dftb1.ind[k+1];
        for (n=0; n<indk1-indk; n++)
          for (m=0; m<indj1-indj; m++) {
            dftb1.a[indj+m][indk+n] = dftb1.hamil[indj+m][indk+n];
            dftb1.b[indj+m][indk+n] = dftb1.overl[indj+m][indk+n];
          }
      }
    }

    // add charge-dependent terms (Hubbard and consec. extcharges)
    if (niter == 0) {
      // printf("call gammamatrix\n");
      gammamatrix(nn, x, dftb1.gammamat, dftb->uhubb1, dftb1.izp);
      // printf("gamma matrix:\n");
      // for (i=0; i<nn; i++) {
      //   for (j=0; j<nn; j++) printf("%9.5f", dftb1.gammamat[i][j]);
      //   printf("\n");
      // }
    }

    // calculate atomic hamilton shift (= sum over gamma*charge)
    for (i=0; i<nn; i++) {
      dftb1.shift[i] = 0.0;
      for (j=0; j<nn; j++)
        dftb1.shift[i] += (dftb1.qmat[j] - dftb->qzero1[dftb1.izp[j]]) * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
    }
    // hamilshift(nn, dftb1.qmat, dftb1.gammamat, dftb1.shift);

    /*
    printf("Shift:\n");
    for (i=0; i<nn; i++) printf("%10.6f", dftb1.shift[i]);
    printf("\n");
    // */

    for (i=0; i<nn; i++)
      dftb1.shift[i] -= dftb1.shiftE[i];

    for (i=0; i<nn; i++)
      for (li=0; li < SQR(dftb->lmax[dftb1.izp[i]]); li++)
        for (j=0; j<=i; j++)
          for (lj=0; lj < SQR(dftb->lmax[dftb1.izp[j]]); lj++) 
            dftb1.a[dftb1.ind[i]+li][dftb1.ind[j]+lj] += 0.5 * dftb1.overl[dftb1.ind[i]+li][dftb1.ind[j]+lj] * (dftb1.shift[i] + dftb1.shift[j]);
    
    // transpose the arrays a and b
    for (j=0; j<dftb1.ndim; j++)
      for (i=0; i<dftb1.ndim; i++) {
        dftb1.a_trans[j*dftb1.ndim+i] = dftb1.a[i][j];
        dftb1.b_trans[j*dftb1.ndim+i] = dftb1.b[i][j];
      }

    // print out the array a
    /*
    printf("A before dsygv\n");
    for (i=0; i<dftb1.ndim; i++) {
      for (j=0; j<dftb1.ndim; j++) printf ("%9.5f", dftb1.a_trans[i*dftb1.ndim+j]);
      printf("\n");
    }
    // */

    ier = -512;
    ier = dsygv(1, 'V', 'L', dftb1.ndim, dftb1.a_trans, dftb1.ndim, dftb1.b_trans, dftb1.ndim, dftb1.ev, dftb1.aux, 3*dftb1.ndim);
    if ((int) ier) {
      printf("\nDSYGV: ier = %d\nEXITING!\n\n", (int) ier);
      exit(-1);
    }
    for (j=0; j<dftb1.ndim; j++)
      for (i=0; i<dftb1.ndim; i++)
        dftb1.a[i][j] = dftb1.a_trans[j*dftb1.ndim+i];

    /*
    printf("\n");
    printf("ier = %d\n", (int) ier);
    printf("\n");
    */

    // print out the array a
    /*
    printf("A after dsygv\n");
    for (i=0; i<dftb1.ndim; i++) {
      for (j=0; j<dftb1.ndim; j++) printf ("%9.5f", dftb1.a[i][j]);
      printf("\n");
    }
    // */

    // calculate occupation (occ) and Fermi energy (efermi),
    fermi(dftb1.ndim, dftb1.ev, dftb1.occ, &efermi, dftb1.nel);
    // for (i=0; i<dftb1.ndim; i++)
    //   printf("%d: %f %f\n", i+1, dftb1.ev[i], dftb1.occ[i]);

    // sum of occupied eigenvalues
    eel = 0.0;
    for (i=0; i<dftb1.ndim && dftb1.occ[i] > dftb->dacc; i++)
      eel += dftb1.occ[i] * dftb1.ev[i];

    // determine Mulliken charges, charge of the whole system and the mulliken
    mulliken(nn, dftb1.qmat, dftb1.qmulli, &qtot, dftb1.ndim, dftb->dacc, dftb1.occ, dftb1.a, dftb1.overl, dftb1.ind, dftb->lmax, dftb1.izp);
    /*
    printf("qmat - after mulliken:\n");
    for (i=0; i<nn; i++)
      printf("%3d %f\n", i+1, dftb1.qmat[i]);
    */

    // complete calculation of electronic energy
    // charge-dependent contribution
    // warning: this will only lead to the right result if convergence has been reached
    ecoul = eext = 0.0;
    for (i=0; i<nn; i++) {
      ecoul += dftb1.shift[i] * (dftb1.qmat[i] + dftb->qzero1[dftb1.izp[i]]);
	eext += dftb1.shiftE[i] * (dftb->qzero1[dftb1.izp[i]] - dftb1.qmat[i]);
      }
    eel += -0.5*ecoul + 0.5*eext;
    // remark: eel containts shiftE already via ev,
    // shift also contains -shiftE, i.e. ecoul also
    // contains contributions from EXT

    // print energy
    // printf("iter: %d, E= %14.9f\n", niter, eel);

    // check convergence
    if (ct->do_scc[ibase] == 0){ break;}
    if (fabs(eel-eelold) < scftol) {
      printf("site %d converged after %d iterations. \n", ibase, niter);
      break;
    }else if (niter == maxiter-1)
      printf("site %d did not converge after %d iterations. \n", ibase, niter);

    eelold = eel;
    //printf("\n%2d %14.8f %14.8f ", niter, eel, dftb1.ev[24]);
    //for (i=0; i<dftb1.ndim; i++)
    //  printf("%12.8f ", dftb1.a[i][24]);

    // Broyden mixing
    broyden(niter, almix, nn, dftb1.qmold, dftb1.qmat, dftb->broyden + ibase);
    for (i=0; i<nn; i++)
      dftb1.qmat[i] = dftb1.qmold[i];

    /*
    printf("qmat - after Broyden:\n");
    for (i=0; i<nn; i++)
      printf("%3d %f\n", i+1, dftb1.qmat[i]);
    */

  } /* end of the SCC cycle */
  //free(a_trans);
  //free(b_trans);

  /**************************************************
   * WE NEED THE FOLLOWING STORED IN MEMORY:        *
   *  1. AO2FO - MO coefficients of the fragment    *
   *           - stored as dftb->phase1[i].a        *
   *  2. CHR.DAT - atom charges on the fragment     *
   *             - stored as dftb->phase1[i].qmat   *
   * WE HAVE IT ALL!                                *
   **************************************************/

  // write out the eigenvalues

/*
  for (i=0; i<dftb1.ndim; i++)
    printf("%2d %f\n", i+1, dftb1.ev[i]);
*/
   //outspec(nn, dftb1.ndim, dftb1.ind, dftb1.ev, dftb1.occ, efermi, dftb1.qmat, dftb1.qmulli, dftb, dftb1);
   //outeigenvectors(dftb1.a, dftb1.ev, dftb1.ind, nn, dftb1);

  // printf("%5d %5d / %2d %14.6f\n", 1, 1, niter, eel);
  // printf("\n***** end of dftb *****\n");


  return 0;
}



/* a variant of EGLCAO -- iterate in two cycles:
 * 1. vary the charges until self-consistence
 * 2. do PME with the self-cosistent charges, and go back to the cycle 1
 * IT LOOKS LIKE THIS IS NO PARTICULARLY GOOD IDEA...
 */
int run_dftb1_doublecycle(charge_transfer_t *ct, dftb_t *dftb, int ibase) // i - nucleobase to be calculated
{
  const double scftol = 1.e-9;
  const double almix = 0.2;
  const int maxiter = 70;
  const int max_pme_iter = 70;  

  int indj, indk, indj1, indk1;
  double eel, ecoul, efermi, eelold, eext, pme_eel, pme_eel_old;
  //
  int i, j, k, m, n, li, lj, niter, pme_iter;
  //int nmaofo, ii, jj, kk, ll, jfo, jao;
  //double r2;
  double qtot;
  // lapack
  long ier; // ndim;

  int nn, ne;
  dvec *x, bond;
  dftb_phase1_t dftb1;

  // printf("Phase 1, site %d\n", ibase+1);

  dftb1 = dftb->phase1[ibase];
  nn = dftb1.nn;
  ne = dftb1.ne;
  x = dftb1.x;

  /* LENGTHS MUST BE CONVERTED TO BOHR! (5.292d-11 m) */

/*
  // write out the coordinates
  for (n=0; n<nn; n++) {
    printf("%d %d %f %f %f\n", n+1, dftb1.izp[n]+1, x[n][XX]*0.529177249, x[n][YY]*0.529177249, x[n][ZZ]*0.529177249);
  }

  // write out the external charges
  for (n=0; n<ne; n++) {
    printf("%f %f %f %f\n", dftb1.xe[n][XX]*0.529177249, dftb1.xe[n][YY]*0.529177249, dftb1.xe[n][ZZ]*0.529177249, dftb1.ze[n]);
  }
*/
  
  // set the initial charges
  for (i=0; i<nn; i++)
    dftb1.qmat[i] = dftb->qzero1[dftb1.izp[i]];

  // initial setup
  eel = 0.0;

  // printf("  icycle iter niter   e(total)\n");
  // printf("====================================\n");

  // setup of charge-independent part of H and S
  for (j=0; j<nn; j++)
    for (k=0; k<=j; k++) {
      slkmatrices(j, k, x, dftb1.au, dftb1.bu, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
      for (n=0; n<dftb1.ind[k+1]-dftb1.ind[k]; n++)
        for (m=0; m<dftb1.ind[j+1]-dftb1.ind[j]; m++) {
          dftb1.hamil[dftb1.ind[j]+m][dftb1.ind[k]+n] = dftb1.au[m][n];
          dftb1.hamil[dftb1.ind[k]+n][dftb1.ind[j]+m] = dftb1.au[m][n];      
          dftb1.overl[dftb1.ind[j]+m][dftb1.ind[k]+n] = dftb1.bu[m][n];
          dftb1.overl[dftb1.ind[k]+n][dftb1.ind[j]+m] = dftb1.bu[m][n];
	}
    }

  /*
  printf("Hamiltonian matrix:\n");
  for (j=0; j<10; j++) {
    for (k=0; k<10; k++) printf("%9.5f", dftb1.hamil[j][k]);
    printf("\n");
  }
  printf("Overlap matrix:\n");
  for (j=0; j<10; j++) {
    for (k=0; k<10; k++) printf("%9.5f", dftb1.overl[j][k]);
    printf("\n");
  }
  */

  pme_eel_old = 1.e10;

  for (pme_iter = 0; pme_iter < max_pme_iter; pme_iter++) {

    /* calculate PME */
    do_pme_for_dftb_part2(ct, dftb, ibase);
    /* save the computed values to shiftE */
    for (j=0; j<nn; j++)
      dftb1.shiftE[j] = dftb1.pot[j];
 
 
    // setup for SCC cycle
    eelold = 1.e10;
 
    // SCC cycle starts here
    for (niter=0; niter<maxiter; niter++) {
 
      // save old charges
      for (i=0; i<nn; i++)
        dftb1.qmold[i] = dftb1.qmat[i];
      
      /*
      printf("qmat - start SCF:\n");
      for (i=0; i<nn; i++)
        printf("%3d %f\n", i+1, dftb1.qmat[i]);
      */
 
      // charge-independent part of H and S
      for (j=0; j<nn; j++) {
        indj = dftb1.ind[j];
        indj1 = dftb1.ind[j+1];
        for (k=0; k<nn; k++) {
          indk = dftb1.ind[k]; 
          indk1 = dftb1.ind[k+1];
          for (n=0; n<indk1-indk; n++)
            for (m=0; m<indj1-indj; m++) {
              dftb1.a[indj+m][indk+n] = dftb1.hamil[indj+m][indk+n];
              dftb1.b[indj+m][indk+n] = dftb1.overl[indj+m][indk+n];
            }
        }
      }
 
      // add charge-dependent terms (Hubbard and consec. extcharges)
      if (niter == 0) {
        // printf("call gammamatrix\n");
        gammamatrix(nn, x, dftb1.gammamat, dftb->uhubb1, dftb1.izp);
        // printf("gamma matrix:\n");
        // for (i=0; i<10; i++) {
        //   for (j=0; j<10; j++) printf("%9.5f", dftb1.gammamat[i][j]);
        //   printf("\n");
        // }
      }
 
      // calculate atomic hamilton shift (= sum over gamma*charge)
      for (i=0; i<nn; i++) {
        dftb1.shift[i] = 0.0;
        for (j=0; j<nn; j++)
          dftb1.shift[i] += (dftb1.qmat[j] - dftb->qzero1[dftb1.izp[j]]) * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
      }
      // hamilshift(nn, dftb1.qmat, dftb1.gammamat, dftb1.shift);
  
      /*
      printf("Shift:\n");
      for (i=0; i<10; i++) printf("%10.6f", dftb1.shift[i]);
      printf("\n");
      */
  
      for (i=0; i<nn; i++)
        dftb1.shift[i] -= dftb1.shiftE[i];
  
      for (i=0; i<nn; i++)
        for (li=0; li < SQR(dftb->lmax[dftb1.izp[i]]); li++)
          for (j=0; j<=i; j++)
            for (lj=0; lj < SQR(dftb->lmax[dftb1.izp[j]]); lj++) 
              dftb1.a[dftb1.ind[i]+li][dftb1.ind[j]+lj] += 0.5 * dftb1.overl[dftb1.ind[i]+li][dftb1.ind[j]+lj] * (dftb1.shift[i] + dftb1.shift[j]);
      
      // transpose the arrays a and b
      for (j=0; j<dftb1.ndim; j++)
        for (i=0; i<dftb1.ndim; i++) {
          dftb1.a_trans[j*dftb1.ndim+i] = dftb1.a[i][j];
          dftb1.b_trans[j*dftb1.ndim+i] = dftb1.b[i][j];
        }
  
      // print out the array a
      /*
      printf("A before dsygv\n");
      for (i=0; i<10; i++) {
        for (j=0; j<10; j++) printf ("%9.5f", dftb1.a_trans[i*dftb1.ndim+j]);
        printf("\n");
      }
      */
  
      ier = -512;
      ier = dsygv(1, 'V', 'L', dftb1.ndim, dftb1.a_trans, dftb1.ndim, dftb1.b_trans, dftb1.ndim, dftb1.ev, dftb1.aux, 3*dftb1.ndim);
      if ((int) ier) {
        printf("\nDSYGV: ier = %d\nEXITING!\n\n", (int) ier);
        exit(-1);
      }
      for (j=0; j<dftb1.ndim; j++)
        for (i=0; i<dftb1.ndim; i++)
          dftb1.a[i][j] = dftb1.a_trans[j*dftb1.ndim+i];
  
      /*
      printf("\n");
      printf("ier = %d\n", (int) ier);
      printf("\n");
      */
  
      // print out the array a
      /*
      printf("A after dsygv\n");
      for (i=0; i<10; i++) {
        for (j=0; j<10; j++) printf ("%9.5f", dftb1.a[i][j]);
        printf("\n");
      }
      */
    
      // calculate occupation (occ) and Fermi energy (efermi),
      fermi(dftb1.ndim, dftb1.ev, dftb1.occ, &efermi, dftb1.nel);
      // for (i=0; i<dftb1.ndim; i++)
      //   printf("%d: %f %f\n", i+1, dftb1.ev[i], dftb1.occ[i]);
    
      // sum of occupied eigenvalues
      eel = 0.0;
      for (i=0; i<dftb1.ndim && dftb1.occ[i] > dftb->dacc; i++)
        eel += dftb1.occ[i] * dftb1.ev[i];
    
      // determine Mulliken charges, charge of the whole system and the mulliken
      mulliken(nn, dftb1.qmat, dftb1.qmulli, &qtot, dftb1.ndim, dftb->dacc, dftb1.occ, dftb1.a, dftb1.overl, dftb1.ind, dftb->lmax, dftb1.izp);
      /*
      printf("qmat - after mulliken:\n");
      for (i=0; i<nn; i++)
        printf("%3d %f\n", i+1, dftb1.qmat[i]);
      */
    
      // complete calculation of electronic energy
      // charge-dependent contribution
      // warning: this will only lead to the right result if convergence has been reached
      ecoul = eext = 0.0;
      for (i=0; i<nn; i++) {
        ecoul += dftb1.shift[i] * (dftb1.qmat[i] + dftb->qzero1[dftb1.izp[i]]);
          eext += dftb1.shiftE[i] * (dftb->qzero1[dftb1.izp[i]] - dftb1.qmat[i]);
        }
      eel += -0.5*ecoul + 0.5*eext;
      // remark: eel containts shiftE already via ev,
      // shift also contains -shiftE, i.e. ecoul also
      // contains contributions from EXT
    
      // print energy
      // printf("iter: %d, E= %14.9f\n", niter, eel);
    
      // check convergence
      if (fabs(eel-eelold) < scftol)
        break;
      eelold = eel;
    
      // Broyden mixing
      broyden(niter, almix, nn, dftb1.qmold, dftb1.qmat, dftb->broyden + ibase);
      for (i=0; i<nn; i++)
      dftb1.qmat[i] = dftb1.qmold[i];

      /*
      printf("qmat - after Broyden:\n");
      for (i=0; i<nn; i++)
        printf("%3d %f\n", i+1, dftb1.qmat[i]);
      */

    } /* end of the inner cycle (without PME) */
    printf("%f ", eel);

    pme_eel = eel;
    if (fabs(pme_eel - pme_eel_old) < scftol)
      break;
    pme_eel_old = pme_eel;
  } /* end of the outer cycle (iterating PME) */
  printf("pme_iter=%d ", pme_iter+1);

  /**************************************************
   * WE NEED THE FOLLOWING STORED IN MEMORY:        *
   *  1. AO2FO - MO coefficients of the fragment    *
   *           - stored as dftb->phase1[i].a        *
   *  2. CHR.DAT - atom charges on the fragment     *
   *             - stored as dftb->phase1[i].qmat   *
   * WE HAVE IT ALL!                                *
   **************************************************/

  // write out the eigenvalues

/*
  for (i=0; i<dftb1.ndim; i++)
    printf("%2d %f\n", i+1, dftb1.ev[i]);
*/
  // outspec(nn, dftb1.ndim, dftb1.ind, dftb1.ev, dftb1.occ, efermi, dftb1.qmat, dftb1.qmulli, dftb, dftb1);
  // outeigenvectors(dftb1.a, dftb1.ev, dftb1.ind, nn, dftb1);

  // printf("%5d %5d / %2d %14.6f\n", 1, 1, niter, eel);
  // printf("\n***** end of dftb *****\n");

  return 0;
}

int run_esp_only(charge_transfer_t *ct, dftb_t *dftb, int ibase)
{
  if (ct->qmmm == 3) {
    do_pme_for_dftb_part2(ct, dftb, ibase);
  }

  return 0;
}
