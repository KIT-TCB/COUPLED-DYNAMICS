#define GMX_LIB_MPI
#define GMX_MPI

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

/* ********************************************************************* */
/*

     PHASE 2 FOR CHARGE TRANSFER -- CALC OF THE COMPLEX

*/
/* ********************************************************************* */

/* NDIM = NORB !!! */

int run_dftb2(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
  int i, j, k, m, n, li, lj;
  // lapack
  long ier;

  int counter, ifo, iao, ii, jj, kk, ll, jfo, jao;

  int nn, ne;
  dvec *x, bond;
  dftb_phase2_t dftb2;

  dftb2 = dftb->phase2;

  nn = dftb2.nn;
  ne = dftb2.ne;
  x = dftb2.x;

  double dist, dist_best;

  MPI_Status ct_mpi_status;

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
      dftb2.qmat[counter] = dftb->phase1[i].qmat[j];
      counter++;
    }
  }

  // setup of charge-independent part of H and S
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
  printf("Hamiltonian matrix:\n");
  for (j=0; j<10; j++) {
    for (k=0; k<10; k++) printf("%9.5f", dftb2.hamil[j][k]);
    printf("\n");
  }
  printf("Overlap matrix:\n");
  for (j=0; j<10; j++) {
    for (k=0; k<10; k++) printf("%9.5f", dftb2.overl[j][k]);
    printf("\n");
  }
// */

  /* call QM/MM if desired */
  if (ct->qmmm == 3)
    do_pme_for_dftb_phase2(ct, dftb);

  // external charges
  for (j=0; j<nn; j++) {
    if (ct->qmmm == 3) /* use the PME value */
      dftb2.shiftE[j] = dftb2.pot[j];
    else { /* loop over all of the MM atoms otherwise */
      dftb2.shiftE[j] = 0.0;
      for (k=0; k<ne; k++) {
        dvec_sub(x[j], dftb2.xe[k], bond);
        dftb2.shiftE[j] += dftb2.ze[k] / dnorm(bond);
      }
    }
    // apply the ESP scaling factor!
    dftb2.shiftE[j] /= ct->esp_scaling_factor;
  }

  // add charge-dependent terms (Hubbard and consec. extcharges)
  gammamatrix(nn, x, dftb2.gammamat, dftb->uhubb2, dftb2.izp);
/*
  printf("gamma matrix:\n");
  for (i=0; i<10; i++) {
    for (j=0; j<10; j++) printf("%9.5f", dftb2.gammamat[i][j]);
    printf("\n");
  }
*/
  // calculate atomic hamilton shift (= sum over gamma*charge)
  for (i=0; i<nn; i++) {
    dftb2.shift[i] = 0.0;
    for (j=0; j<nn; j++)
      dftb2.shift[i] += (dftb2.qmat[j] - dftb->qzero2[dftb2.izp[j]]) * (i>j ? dftb2.gammamat[i][j] : dftb2.gammamat[j][i]);
  }

  for (i=0; i<nn; i++)
    dftb2.shift[i] -= dftb2.shiftE[i];

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
        for (j=0; j<dftb->phase1[ii].norb; j++)
          dftb2.Taf[dftb2.inf[ii]+i][dftb2.inf[ii]+j] = dftb->phase1[ii].a[i][j];

  // send AO-Hamiltonian, overlap and Transformation matrix to all ranks, then calculate matrix elements parallel on different ranks
  if (ct_mpi_rank == 0){
    for (i=1; i < ct_mpi_size; i++) {
      MPI_Send(dftb2.hamil[0], SQR(dftb2.norb), MPI_DOUBLE, i,100+i, ct_mpi_comm);
      MPI_Send(dftb2.overl[0], SQR(dftb2.norb), MPI_DOUBLE, i,200+i, ct_mpi_comm);
      MPI_Send(dftb2.Taf[0], SQR(dftb2.norb), MPI_DOUBLE, i,300+i, ct_mpi_comm);
    }
  } else {
    MPI_Recv(dftb2.hamil[0], SQR(dftb2.norb), MPI_DOUBLE, 0,100+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
    MPI_Recv(dftb2.overl[0], SQR(dftb2.norb), MPI_DOUBLE, 0,200+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
    MPI_Recv(dftb2.Taf[0], SQR(dftb2.norb), MPI_DOUBLE, 0,300+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
  }

   printf("phase2_fast start at %f\n", (double) clock()/CLOCKS_PER_SEC);
    // calculate the elements "H_ij = sum_{mu nu} c_{i mu} c_{j nu} H_{mu nu}"
    // ifrg = # of orbital being considered, i.e. ct->homo
    counter=0;
    for (i = 0; i < ct->sites; i++) {
    for (j = i; j < ct->sites; j++) {
      dist_best=20.0; // = 20 bohr
    // check if sites are close enough for non-vanishing coupling //
      if (i < j) {  
        for (m = 0; m < ct->site[i].atoms; m++)
        for (n = 0; n < ct->site[j].atoms; n++){
          dvec_sub(dftb->phase1[i].x[m],dftb->phase1[j].x[n], bond);
          dist = dnorm(bond);
          if (dist < dist_best){ //find shortest distance 
            dist_best=dist;
          }
        }
      }
    //printf("phase2_fast check at %f\n", (double) clock()/CLOCKS_PER_SEC);
    for (ii = 0; ii < ct->site[i].homos; ii++) {
      ifo = ct->site[i].homo[ii] + dftb2.inf[i] - 1;
      for (jj = 0; jj < ct->site[j].homos; jj++) {
        jfo = ct->site[j].homo[jj] + dftb2.inf[j] - 1;
        dftb2.THamil[ifo][jfo] = 0.0;
        dftb2.OverlF[ifo][jfo] = 0.0;
        if (counter % ct_mpi_size == ct_mpi_rank){
          if ((i < j && dist_best < 20.0) || ifo==jfo){ //SLKOs are tabulated up to 20 bohr // only on-site energies, no on-site coupling
            for (iao=0; iao<dftb2.norb; iao++)
            for (jao=0; jao<dftb2.norb; jao++) {
              dftb2.THamil[ifo][jfo] += dftb2.Taf[iao][ifo] * dftb2.hamil[iao][jao] * dftb2.Taf[jao][jfo];
              dftb2.OverlF[ifo][jfo] += dftb2.Taf[iao][ifo] * dftb2.overl[iao][jao] * dftb2.Taf[jao][jfo];
            }
          }
          if (ct_mpi_rank != 0) { 
            MPI_Send(&dftb2.THamil[ifo][jfo], 1, MPI_DOUBLE,0,10000+counter, ct_mpi_comm);
            //printf("matrix element no. %d sent from rank %d\n", counter, ct_mpi_rank);
            MPI_Send(&dftb2.OverlF[ifo][jfo], 1, MPI_DOUBLE,0,20000+counter, ct_mpi_comm);
          }
        }
        else if (ct_mpi_rank == 0){ //rank 0 recieves data if it was calculated on another rank
          MPI_Recv(&dftb2.THamil[ifo][jfo], 1, MPI_DOUBLE, counter % ct_mpi_size ,10000+counter, ct_mpi_comm, &ct_mpi_status);
          //printf("matrix element no. %d recieved from rank %d\n", counter,  counter % ct_mpi_size);
          MPI_Recv(&dftb2.OverlF[ifo][jfo], 1, MPI_DOUBLE, counter % ct_mpi_size ,20000+counter, ct_mpi_comm, &ct_mpi_status);
        }
        counter++;
      }
    }
    }
    }
    printf("phase2_fast assemble at %f\n", (double) clock()/CLOCKS_PER_SEC);
    kk = 0;
    for (i = 0; i < ct->sites; i++) {
    for (ii = 0; ii < ct->site[i].homos; ii++) {
      ifo = ct->site[i].homo[ii] + dftb2.inf[i] - 1;
      ll = 0;
      for (j = 0; j < ct->sites; j++) {
      for (jj = 0; jj < ct->site[j].homos; jj++) { 
        jfo = ct->site[j].homo[jj] + dftb2.inf[j] - 1;
        if (ifo <= jfo) {
          dftb2.tij[kk][ll] = dftb2.tij[ll][kk] = dftb2.THamil[ifo][jfo];
          dftb2.sij[kk][ll] = dftb2.sij[ll][kk] = dftb2.OverlF[ifo][jfo];
          dftb2.THamil[jfo][ifo] = dftb2.THamil[ifo][jfo];
          dftb2.OverlF[jfo][ifo] = dftb2.OverlF[ifo][jfo];
        }
        ll++;
      }
      }
      kk++;
    }
    }
    printf("phase2_fast stop at %f\n", (double) clock()/CLOCKS_PER_SEC);


//set off-diagonals for each site to zero (and diagonals to one). (should be zero, deviations arise from the usage of different basis sets for phase1 and phase2)
// MAY CAUSE INCORRECT ORTHOGONALIZATION
  kk=0;
  for (i = 0; i < ct->sites; i++) {
    ll=0; 
    for (j = 0; j < ct->site[i].homos; j++){
      for (ii = 1; ii < ct->site[i].homos - ll ; ii++) {
        dftb2.tij[kk][kk+ii]=0;
        dftb2.sij[kk][kk+ii]=0;
        dftb2.tij[kk+ii][kk] = dftb2.tij[kk][kk+ii];
        dftb2.sij[kk+ii][kk] = dftb2.sij[kk][kk+ii];
      }
      ll++;
      kk++;
    } 
  }
  for (i = 0; i<ct->dim; i++)
    dftb2.sij[i][i]=1.0;


 
    // write out the obtained CG Hamiltonian
    if (ct_mpi_rank == 0){
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
  // used for projection the compressed basis, otherwise there is overlap between orbitals on the same site
  for (j=0; j<nn; j++)
    for (k=0; k<=j; k++) {
      slkmatrices(j, k, x, dftb2.au, dftb2.bu, dftb->lmax, dftb->dim2, dftb->dr2, dftb2.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
      for (n=0; n<dftb2.ind[k+1]-dftb2.ind[k]; n++)
        for (m=0; m<dftb2.ind[j+1]-dftb2.ind[j]; m++) {
          dftb2.overl_c[dftb2.ind[j]+m][dftb2.ind[k]+n] = dftb2.bu[m][n];
          dftb2.overl_c[dftb2.ind[k]+n][dftb2.ind[j]+m] = dftb2.bu[m][n];
        }
    }
*/
    // orthogonalize the Hamiltonian
    ier = -512;
    ier = orthogonalize(dftb2.tij, dftb2.sij, dftb2.THamilOrtho, (long) ct->dim, dftb->orthogo);
    if ((int)ier) {
      printf("Orthogonalizer returned %d, exiting!\n", (int) ier);
      exit(-1);
    }

/*
    f78 = fopen("FRAGanalysis.DAT", "w");
    fprintf(f78, "  FOs   FOs        H-elements         Overlap-M   \n");
    fprintf(f78, "==================================================\n");
    f79 = fopen("FRAGanalysis.ortho.DAT", "w");
    fprintf(f79, " FO FO      H-elements  \n");
    fprintf(f79, "========================\n");
    for (i=0; i<ct->sites; i++) {
      ifo = ct->homo[i] + dftb2.inf[i] - 1;
      for (j=0; j<ct->sites; j++) {
        jfo = ct->homo[j] + dftb2.inf[j] - 1;
        fprintf(f78, "%3d%3d%3d%3d %16.12f %s %16.12f\n",
                i, ct->homo[i], j, ct->homo[j],
		dftb2.THamil[ifo][jfo]*27.2116, "eV", dftb2.OverlF[ifo][jfo]);
      }
    }
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
	fprintf(f79, "%3d%3d %16.12f Ha %16.12f eV\n",
	            i, j, dftb2.THamilOrtho[i][j], dftb2.THamilOrtho[i][j]*27.2116);
    fclose(f78);
    fclose(f79);
*/

  return 0;
}
