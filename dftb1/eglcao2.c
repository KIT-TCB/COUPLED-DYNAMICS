#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

int run_dftb2(charge_transfer_t *ct, dftb_t *dftb) {
/* This function performs DFTB calculations of the complex consisting of several fragments.
   It takes the results from run_dftb1() and constructs an orthogonal "coarse-grained" Hamiltonian in the FO basis*/

// PARAMETERS:
// ct    = (in) main data structure with information about the charge transfer calculation
// dftb  = (in) main data structure for DFTB calculations
////////////////////////////////////////////////////////////////////////

  int i, j, k,l, m, n, li, lj;
  // lapack
  long ier;

  int counter, counter1,counter2, ifo, iao, ii, jj, kk, ll, jfo, jao, do_scc, has_neighbor;

  int nn, ne;
  dvec *x, bond;
  dftb_phase2_t dftb2;

  dftb2 = dftb->phase2;

  nn = dftb2.nn;
  ne = dftb2.ne;
  x = dftb2.x;

  double dist, dist_best;

  /* LENGTHS MUST BE CONVERTED TO BOHR! (5.292d-11 m) */

/*
  // write out the coordinates
  for (n=0; n<nn; n++) {
    printf("%d %d %f %f %f\n", n+1, dftb2.izp[n]+1, x[n][XX]*0.529177249, x[n][YY]*0.529177249, x[n][ZZ]*0.529177249);
  }
  printf("\n");

  // write out the external charges
  for (n=0; n<ne; n++) {
    printf("%f %f %f %f\n", dftb2.xe[n][XX]*0.529177249, dftb2.xe[n][YY]*0.529177249, dftb2.xe[n][ZZ]*0.529177249, dftb2.ze[n]);
  }
*/

  /* set the charges - read them from the phase 1! */
  counter = 0;
  for (i=0; i<ct->sites; i++) {
    for (j=0; j<ct->site[i].atoms; j++) {
      dftb2.qmat[counter] = ct->site[i].do_scc > 0 ? dftb->phase1[i].qmat[j] : dftb->qzero1[dftb->phase1[i].izp[j]] ; //use only DFTB2 Hamiltonian if the sites were also calculated using DFTB2
      //dftb2.qmat[counter] = dftb->phase1[i].qmat[j];
      counter++;
    }
  }


/* 
  // setup of charge-independent part of H and S (old version)
  for (j=0; j<nn; j++)
    for (k=0; k<=j; k++) {
      slkmatrices(j, k, x, dftb2.au, dftb2.bu, dftb->lmax, dftb->dim2, dftb->dr2, dftb2.izp, dftb->skstab2, dftb->skhtab2, dftb->skself2);
      for (n=0; n<dftb2.ind[k+1]-dftb2.ind[k]; n++)
        for (m=0; m<dftb2.ind[j+1]-dftb2.ind[j]; m++) {
          dftb2.hamil[dftb2.ind[j]+m][dftb2.ind[k]+n] = dftb2.au[m][n];
          dftb2.hamil[dftb2.ind[k]+n][dftb2.ind[j]+m] = dftb2.au[m][n];      
          dftb2.overl[dftb2.ind[j]+m][dftb2.ind[k]+n] = dftb2.bu[m][n];
          dftb2.overl[dftb2.ind[k]+n][dftb2.ind[j]+m] = dftb2.bu[m][n];
	}
    }
/*
  for (j=0; j<dftb2.norb; j++)
    for (k=0; k<=j; k++) {
          dftb2.hamil[j][k] = 0.0;
          dftb2.hamil[k][j] = 0.0;
          dftb2.overl[j][k] = 0.0;
          dftb2.overl[k][j] = 0.0;
   }

*/
///*

  // get neighbor list
  counter1=-1;
  for (i=0; i<ct->sites; i++)
  for (ii=0; ii<ct->site[i].atoms; ii++){
    counter1++;
    counter2=-1;
    for (j=0; j<ct->sites; j++)
    for (jj=0; jj<ct->site[j].atoms; jj++){
      counter2++;
      if (counter2 > counter1) {break;}
      dvec_sub(dftb->phase1[i].x[ii],dftb->phase1[j].x[jj], bond); 
      dist = dnorm(bond);
      if (dist < 20.0){
        dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=1;
      }else{
        dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=0;
      }  
    }
  }
/*
  for (i=0; i<dftb2.nn; i++){
  for (j=0; j<dftb2.nn; j++){
    printf("%d ", dftb->nl[i][j]);
    }
  printf("\n");
  }
*/ 
  
  // setup of charge-independent part of H and S ("double zeta" version)
  // in this version matrix elements between atoms of the same fragment are 
  counter1=-1;
  for (i=0; i<ct->sites; i++)
  for (ii=0; ii<ct->site[i].atoms; ii++){
    counter1++;
    counter2=-1;
    for (j=0; j<ct->sites; j++) 
    for (jj=0; jj<ct->site[j].atoms; jj++){
      counter2++;
      if (counter2 > counter1) {break;} //matrix is symmetric. other triangular part can be derived
      if (i==j){
        slkmatrices(counter1, counter2, x, dftb2.au, dftb2.bu, dftb->lmax, dftb->dim1, dftb->dr1, dftb2.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1); //use confined orbitals for AOs on the same site
      }else if (dftb->nl[counter1][counter2]){
          slkmatrices(counter1, counter2, x, dftb2.au, dftb2.bu, dftb->lmax, dftb->dim2, dftb->dr2, dftb2.izp, dftb->skstab2, dftb->skhtab2, dftb->skself2);
      }else{
        for (n=0; n<LDIM; n++)
        for (m=0; m<LDIM; m++) {
          dftb2.au[m][n]=0.0;
          dftb2.bu[m][n]=0.0;
        }
      }
      
      for (n=0; n<dftb2.ind[counter2+1]-dftb2.ind[counter2]; n++)
        for (m=0; m<dftb2.ind[counter1+1]-dftb2.ind[counter1]; m++) {
          dftb2.hamil[dftb2.ind[counter1]+m][dftb2.ind[counter2]+n] = dftb2.au[m][n];
          dftb2.hamil[dftb2.ind[counter2]+n][dftb2.ind[counter1]+m] = dftb2.au[m][n];      
          dftb2.overl[dftb2.ind[counter1]+m][dftb2.ind[counter2]+n] = dftb2.bu[m][n];
          dftb2.overl[dftb2.ind[counter2]+n][dftb2.ind[counter1]+m] = dftb2.bu[m][n];
	}
    }
  }
//*/
/*
  printf("Hamiltonian matrix:\n");
  for (j=0; j<20; j++) {
    for (k=0; k<20; k++) printf("%9.5f", dftb2.hamil[j][k]);
    printf("\n");
  }
  printf("Overlap matrix:\n");
  for (j=0; j<20; j++) {
    for (k=0; k<20; k++) printf("%9.5f", dftb2.overl[j][k]);
    printf("\n");
  }
// */


  /* calculate atomic hamilton shift */
  printf("phase2 before pme %f\n", (double) clock()/CLOCKS_PER_SEC);
  if (ct->qmmm == 3)
    do_pme_for_dftb_phase2(ct, dftb);

  // external charges and field
  for (j=0; j<nn; j++) {
    if (ct->qmmm == 3){ /* use the PME value */
      dftb2.shiftE[j] = diprod(x[j],ct->efield);
      dftb2.shiftE[j] += dftb2.pot[j] / ct->esp_scaling_factor;
    } else { /* loop over all of the MM atoms otherwise */
      dftb2.shiftE[j] = diprod(x[j],ct->efield);
      for (k=0; k<ne; k++) {
        dvec_sub(x[j], dftb2.xe[k], bond);
        dftb2.shiftE[j] += dftb2.ze[k] / (dnorm(bond) * ct->esp_scaling_factor);
      }
    }
  }
  for (i=0; i<nn; i++)
    dftb2.shift[i] = -dftb2.shiftE[i];


  // add charge-dependent terms (Hubbard)
  printf("phase2 before gamma %f\n", (double) clock()/CLOCKS_PER_SEC);
  do_scc=0;
  for (i=0; i<ct->sites; i++)
  if(ct->site[i].do_scc)
    do_scc=1;
  if(do_scc){
    gammamatrix(nn, x, dftb2.gammamat, dftb->uhubb2, dftb2.izp);
    for (i=0; i<nn; i++) {
      for (j=0; j<nn; j++)
        dftb2.shift[i] += (dftb2.qmat[j] - dftb->qzero2[dftb2.izp[j]]) * (i>j ? dftb2.gammamat[i][j] : dftb2.gammamat[j][i]);
    }
  }
  printf("phase2 after gamma%f\n", (double) clock()/CLOCKS_PER_SEC);

  /*
  printf("gamma matrix: nn= %d\n",nn);
  for (i=0; i<10; i++) {
    for (j=0; j<10; j++) printf("%9.5f", dftb2.gammamat[i][j]);
    printf("\n");
  }
  // */

  for (i=0; i<nn; i++)
    for (li=0; li < SQR(dftb->lmax[dftb2.izp[i]]); li++)
      for (j=0; j<=i; j++)
        for (lj=0; lj < SQR(dftb->lmax[dftb2.izp[j]]); lj++) {
          dftb2.hamil[dftb2.ind[i]+li][dftb2.ind[j]+lj] += 0.5 * dftb2.overl[dftb2.ind[i]+li][dftb2.ind[j]+lj] * (dftb2.shift[i] + dftb2.shift[j]);
          dftb2.hamil[dftb2.ind[j]+lj][dftb2.ind[i]+li] = dftb2.hamil[dftb2.ind[i]+li][dftb2.ind[j]+lj];
        }





  // calculate matrix which transforms the H-matrix from AOs to FOs:
  for (ifo=0; ifo<dftb2.norb; ifo++)
    for (iao=0; iao<dftb2.norb; iao++)
      dftb2.Taf[ifo][iao] = 0.0;

  // copy the orbital coefficients from phase 1
  // pay attention to correct indexing! - a[AO][MO] - the columns of matrix a[][] are MOs
  for (ii=0; ii<ct->sites; ii++)
    for (i=0; i<dftb->phase1[ii].norb; i++)
      for (j=0; j<dftb->phase1[ii].norb; j++){
        dftb2.Taf[dftb2.inf[ii]+i][dftb2.inf[ii]+j] = dftb->phase1[ii].a[i][j];
        //if (ct->jobtype == cteTDA)
        //  dftb2.Taf[dftb2.inf[ii]+i][dftb2.inf[ii]+j] = (i==j) ? 1 : 0 ; // dont transform backbone MOs for all atom calc
  }



  // calculate the elements "H_ij = sum_{mu nu} c_{i mu} c_{j nu} H_{mu nu}"
  printf("phase2_fast start assemble at %f\n", (double) clock()/CLOCKS_PER_SEC);
  for (i = 0; i < ct->sites; i++){
  for (j = i; j < ct->sites; j++){ 
    has_neighbor=0;                           

    for (m = 0; m < ct->site[i].atoms; m++)                                        
    for (n = 0; n < ct->site[j].atoms; n++){                                       
      if (dftb->nl[dftb2.atind[i]+m][dftb2.atind[j]+n]){
        has_neighbor=1;
        break;
      }
    } 
                                                                             
    for (ii = 0; ii < ct->site[i].homos; ii++) {
      counter1=ct->indFO[i]+ii;
    for (jj = 0; jj < ct->site[j].homos; jj++) {
      counter2=ct->indFO[j]+jj;

      dftb2.tij[counter1][counter2]=0.0;
      dftb2.sij[counter1][counter2]=0.0;

      if ( has_neighbor || i==j){ // calc hamilton only for close sites, SLKOs are tabulated up to 20 bohr // 
        ifo = ct->site[i].homo[ii] + dftb2.inf[i] - 1;
        jfo = ct->site[j].homo[jj] + dftb2.inf[j] - 1;
        for (iao=dftb2.inf[i]; iao<dftb2.inf[i+1]; iao++)
        for (jao=dftb2.inf[j]; jao<dftb2.inf[j+1]; jao++) {
          dftb2.tij[counter1][counter2] += dftb2.Taf[iao][ifo] * dftb2.hamil[iao][jao] * dftb2.Taf[jao][jfo];
          dftb2.sij[counter1][counter2] += dftb2.Taf[iao][ifo] * dftb2.overl[iao][jao] * dftb2.Taf[jao][jfo];
        }
      }
    }
    }
  }
  }
  for (i = 0; i < ct->dim; i++) 
  for (j = i; j < ct->dim; j++) {
    dftb2.tij[j][i]=dftb2.tij[i][j];
    dftb2.sij[j][i]=dftb2.sij[i][j];
  }
  printf("phase2_fast end assemble at %f\n", (double) clock()/CLOCKS_PER_SEC);

/* 
    // write out the obtained CG Hamiltonian
    printf("CG Hamiltonian (eV):\n");
    for (i=0; i<ct->dim; i++) {
      for (j=0; j<ct->dim; j++) printf("%12.6f", dftb2.tij[i][j]*27.21);
      printf("\n");
    }
    // write out the obtained CG overlap
    printf("CG overlap:\n");
    for (i=0; i<ct->dim; i++) {
      for (j=0; j<ct->dim; j++) printf("%12.6f", dftb2.sij[i][j]);
      printf("\n");
    }
//*/

///*
//purify S and H matrix (no longer useful with "double-zeta" version)
  kk=0;
  for (i = 0; i < ct->sites; i++) {
    for (j = 0; j < ct->site[i].homos; j++){
      dftb2.tij[kk][kk]=dftb->phase1[i].ev[ct->site[i].homo[j]-1]; //reintroduced. otherwise equal fragments have different energy if one interacts with QM neighbor and the other one with MM neighbor
      //for (ii = j+1; ii < ct->site[i].homos ; ii++) { //this loop is wrong // i guess now it's right
      //  dftb2.tij[kk][kk+ii] = dftb2.tij[kk+ii][kk] =0;
      //  dftb2.sij[kk][kk+ii] = dftb2.sij[kk+ii][kk] =0;
      //}
      kk++;
    } 
  }
//*/  

/*
    // write out the obtained CG Hamiltonian
    printf("CG Hamiltonian (eV):\n");
    for (i=0; i<ct->dim; i++) {
      for (j=0; j<ct->dim; j++) printf("%12.6f", dftb2.tij[i][j]*27.21);
      printf("\n");
    }
    // write out the obtained CG overlap
    printf("CG overlap:\n");
    for (i=0; i<ct->dim; i++) {
      for (j=0; j<ct->dim; j++) printf("%12.6f", dftb2.sij[i][j]);
      printf("\n");
    }
// */


    // orthogonalize the Hamiltonian
    ier = -512;
    ier = orthogonalize(dftb2.tij, dftb2.sij, dftb2.THamilOrtho, (long) ct->dim, dftb->orthogo);
    if ((int)ier) {
      printf("Orthogonalizer returned %d, exiting!\n", (int) ier);
      exit(-1);
    }
/*
    // write out the ortho CG Hamiltonian
    printf("ortho CG Hamiltonian (eV):\n");
    for (i=0; i<ct->dim; i++) {
      for (j=0; j<ct->dim; j++) printf("%12.6f", dftb2.THamilOrtho[i][j]*27.21);
      printf("\n");
    }
//  */

/*
    FILE *f78;
    FILE *f79;
    f78 = fopen("FRAGanalysis.DAT", "w");
    fprintf(f78, "  FOs   FOs        H-elements         Overlap-M   \n");
    fprintf(f78, "==================================================\n");
    f79 = fopen("FRAGanalysis.ortho.DAT", "w");
    fprintf(f79, " FO FO      H-elements  \n");
    fprintf(f79, "========================\n");
    for (i=0; i<ct->sites; i++) 
    for (ii=0; ii<ct->site[i].homos; ii++) {
      ifo = ct->site[i].homo[ii] + dftb2.inf[i] - 1;
      for (j=0; j<ct->sites; j++) 
      for (jj=0; jj<ct->site[j].homos; jj++) {
        jfo = ct->site[j].homo[jj] + dftb2.inf[j] - 1;
        fprintf(f78, "%3d%3d%3d%3d %16.12f %s %16.12f\n",
                i, ct->site[i].homo[ii], j, ct->site[j].homo[jj],
		dftb2.THamil[ifo][jfo]*27.2116, "eV", dftb2.OverlF[ifo][jfo]);
      }
    }
    for (i=0; i<ct->dim; i++)
      for (j=0; j<ct->dim; j++)
	fprintf(f79, "%3d%3d %16.12f Ha %16.12f eV\n",
	            i, j, dftb2.THamilOrtho[i][j], dftb2.THamilOrtho[i][j]*27.2116);
    fclose(f78);
    fclose(f79);
//*/

  return 0;
}
