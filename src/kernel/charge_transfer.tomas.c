/*************************
 * Charge transfer in DNA
 * Tomas Kubar
 *************************/

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include "charge_transfer.h"

#ifdef GMX_MPI
#define PRINTF(...) if (ct_mpi_rank==0) printf(__VA_ARGS__)
#else
#define PRINTF(...) printf(__VA_ARGS__)
#endif

void ct_get_index(int isize[], atom_id *index[]) { /* atom_id is typedef int! */
  /* these were parameters to get_index() - we do not need them here */
  char *fnm = "charge-transfer.ndx";
  /* maybe unnecessary?
  int ngrps = 1;
  char *grpnames;
  */
  char    ***gnames;
  t_blocka *grps = NULL; 
  int j;
  
  snew(gnames,1);

  grps=init_index(fnm,gnames);
  
  if (grps->nr != 1) {
    fprintf(stderr,"\nThe index file %s contains a number of groups different from 1! (%d)\nExiting!\n\n", fnm, grps->nr);
    exit(-1);
  }

  /* ct_read_group(grps,*gnames,isize,index); */

  fprintf(stderr,"   -- group %s has %d atoms\n", **gnames, grps->index[1] - grps->index[0]);

  isize[0]=grps->index[1]-grps->index[0];
  snew(index[0],isize[0]);
  for(j=0; j<isize[0]; j++)
    index[0][j]=grps->a[grps->index[0]+j];

  return;
}

int ct_atom_in_group(int atom, int *list, int list_size) {
  int i;

  for (i=0; i<list_size; i++)
    if (atom == list[i])
      return 1;

  return 0;
}

int negf_init_arrays(ct_negf_lorentz_t *negf, double *rk_timestep, double *wf)
{
  extern long negf_const_mat_(double *, double *, long *, long *);
  extern long negf_jacobi_(double *, double *, double *, long *);
  extern long negf_construct_nu_();
  extern long negf_calc_r_(double *, double *, double *, long *, double *, double *, double *, long *);
  extern long negf_create_hi_(double *, double *, double *, double *, long *,
                              double *, double *, double *, double *, long *, double *);
  //extern long negf_create_gam_(double_complex *, double_complex *, double_complex *, double_complex *, 
  //                             double_complex *, double_complex *, double_complex *, double_complex *,
  //	                       double *, double *, double *, double *, double *, double_complex *, double_complex *,
  //			       double *, double *, double *, double *, double *, double_complex *, double_complex *,
  //			       long *, long *, long *, long *, long *, double *);
  extern long negf_create_gam_(double *, double *, double *, double *, double *,
			       double *, double *, double *, double *, double *,
			       long *, long *, double *);
  extern long negf_set_initial_values_(double_complex *, double *);
  extern long negf_allocate_fortran_(long *, long *, long *, long *);
  extern long negf_rk_init_();

  long info;

  /* intialize */
  //ct->negf_arrays->kl[0] = ct->negf_arrays->n_poles_l[0] + ct->negf_arrays->n_lorentz[0];
  //ct->negf_arrays->kr[0] = ct->negf_arrays->n_poles_r[0] + ct->negf_arrays->n_lorentz[0];
  //snew(ct->negf_arrays->hi_left_p, ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->hi_left_m, ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->hi_right_p, ct->negf_arrays->kr[0]);
  //snew(ct->negf_arrays->hi_right_m, ct->negf_arrays->kr[0]);

  //snew(ct->negf_arrays->gam_greater_l_m, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->gam_greater_l_p, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->gam_lesser_l_m,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->gam_lesser_l_p,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
  //snew(ct->negf_arrays->gam_greater_r_m, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
  //snew(ct->negf_arrays->gam_greater_r_p, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
  //snew(ct->negf_arrays->gam_lesser_r_m,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
  //snew(ct->negf_arrays->gam_lesser_r_p,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);

  //snew(ct->negf_arrays->pi_l, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  (ct->negf_arrays->kl[0] + 1));
  //snew(ct->negf_arrays->pi_r, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  (ct->negf_arrays->kr[0] + 1));

  //snew(ct->negf_arrays->omega_ll1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_ll2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kl[0] + 1 - ct->negf_arrays->n_lorentz[0]));
  //snew(ct->negf_arrays->omega_ll3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kl[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_lr1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_lr2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kl[0] + 1 - ct->negf_arrays->n_lorentz[0]));
  //snew(ct->negf_arrays->omega_lr3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kl[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_rl1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_rl2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kr[0] + 1 - ct->negf_arrays->n_lorentz[0]));
  //snew(ct->negf_arrays->omega_rl3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kr[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_rr1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
  //snew(ct->negf_arrays->omega_rr2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kr[0] + 1 - ct->negf_arrays->n_lorentz[0]));
  //snew(ct->negf_arrays->omega_rr3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kr[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));

  negf->length_rkvec[0] = negf->n[0] * negf->n[0] * (1 + negf->kl[0] + negf->kr[0] +
                                      4 * negf->n_lorentz[0] * (negf->n_lorentz[0] + negf->n_poles_l[0] + negf->n_poles_r[0]));
  snew(negf->rkvec,  negf->length_rkvec[0]);
  snew(negf->drkvec, negf->length_rkvec[0]);

  negf->beta[0] = 1. / (BOLTZMANN_HARTREE_KELVIN * HARTREE_TO_EV * negf->temp[0]); /* make everything in eV within NEGF !!! */

  snew(negf->eig_val_mat_l,  2 * negf->n_poles_l[0]);
  snew(negf->eig_vect_mat_l, SQR(2 * negf->n_poles_l[0]));
  snew(negf->mat_l,          SQR(2 * negf->n_poles_l[0]));
  snew(negf->nu_l,           negf->n_poles_l[0]);
  snew(negf->r_alpha_l,      2 * negf->n_poles_l[0]);
  snew(negf->r_l,            negf->n_poles_l[0]);
  snew(negf->eig_val_mat_r,  2 * negf->n_poles_r[0]);
  snew(negf->eig_vect_mat_r, SQR(2 * negf->n_poles_r[0]));
  snew(negf->mat_r,          SQR(2 * negf->n_poles_r[0]));
  snew(negf->nu_r,           negf->n_poles_r[0]);
  snew(negf->r_alpha_r,      2 * negf->n_poles_r[0]);
  snew(negf->r_r,            negf->n_poles_r[0]);

  info = negf_allocate_fortran_(negf->n, negf->n_poles_l, negf->n_poles_r, negf->length_rkvec);

  info = negf_const_mat_(negf->mat_l, negf->mat_r, negf->n_poles_l, negf->n_poles_r);
  info = negf_jacobi_(negf->mat_l, negf->eig_val_mat_l, negf->eig_vect_mat_l, negf->n_poles_l);
  info = negf_jacobi_(negf->mat_r, negf->eig_val_mat_r, negf->eig_vect_mat_r, negf->n_poles_r);
  info = negf_calc_r_(negf->r_alpha_l, negf->eig_val_mat_l, negf->eig_vect_mat_l, negf->n_poles_l,
                      negf->r_alpha_r, negf->eig_val_mat_r, negf->eig_vect_mat_r, negf->n_poles_r);
  info = negf_construct_nu_(negf->eig_val_mat_l, negf->nu_l, negf->r_l, negf->r_alpha_l, negf->n_poles_l,
                            negf->eig_val_mat_r, negf->nu_r, negf->r_r, negf->r_alpha_r, negf->n_poles_r);

  //info = negf_create_hi_(negf->hi_left_m,  negf->hi_left_p,  negf->eps_l, negf->w0_l, negf->e_f_left,  negf->nu_l, negf->kl, negf->n_poles_l,
  //                       negf->hi_right_m, negf->hi_right_p, negf->eps_r, negf->w0_r, negf->e_f_right, negf->nu_r, negf->kr, negf->n_poles_r, negf->beta);
  //info = negf_create_gam_(negf->gam_greater_l_m, negf->gam_greater_l_p, negf->gam_lesser_l_m, negf->gam_lesser_l_p,
  //                        negf->gam_greater_r_m, negf->gam_greater_r_p, negf->gam_lesser_r_m, negf->gam_lesser_r_p,
//		          negf->gam_l, negf->w0_l, negf->eps_l, negf->e_f_left,  negf->r_l, negf->hi_left_m,  negf->hi_left_p,
//		          negf->gam_r, negf->w0_r, negf->eps_r, negf->e_f_right, negf->r_r, negf->hi_right_m, negf->hi_right_p,
//		          negf->n, negf->kl, negf->kr, negf->n_poles_l, negf->n_poles_r, negf->beta);
  //info = negf_set_initial_values_(negf->rho, negf->pi_l, negf->pi_r,
  //                                negf->omega_ll1, negf->omega_ll2, negf->omega_ll3, negf->omega_lr1, negf->omega_lr2, negf->omega_lr3,
  //                                negf->omega_rl1, negf->omega_rl2, negf->omega_rl3, negf->omega_rr1, negf->omega_rr2, negf->omega_rr3,
//				  negf->rkvec,
//				  negf->n, negf->n_lorentz, negf->kl, negf->kr, negf->length_rkvec);
  info = negf_create_hi_(negf->eps_l, negf->w0_l, negf->e_f_left,  negf->nu_l, negf->n_poles_l,
                         negf->eps_r, negf->w0_r, negf->e_f_right, negf->nu_r, negf->n_poles_r, negf->beta);
  info = negf_create_gam_(negf->gam_l, negf->w0_l, negf->eps_l, negf->e_f_left,  negf->r_l,
		          negf->gam_r, negf->w0_r, negf->eps_r, negf->e_f_right, negf->r_r,
		          negf->n_poles_l, negf->n_poles_r, negf->beta);
  info = negf_set_initial_values_(negf->rkvec, wf);
  info = negf_rk_init_();

  return 0;
}

int negf_propagate(charge_transfer_t *ct)
{
  extern void negf_create_h_(double *);
  extern void negf_do_step_(double_complex *, double_complex *, double *);
  extern void negf_calculate_current_(double_complex *, double *, double *);

  negf_create_h_(ct->hamiltonian[0]);
  negf_do_step_(ct->negf_arrays->rkvec, ct->negf_arrays->drkvec, &(ct->rk_timestep));
  negf_calculate_current_(ct->negf_arrays->rkvec, ct->negf_arrays->current, ct->wf);

  return 0;
}


/***************************************
 * INITIALIZE THE CHARGE TRANSFER CODE *
 ***************************************/

#ifdef GMX_MPI
void init_charge_transfer(gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path, int ct_mpi_rank) {
#else
void init_charge_transfer(gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path) {
#endif
  FILE *f;
  char *line1=NULL, *line2=NULL;
  size_t len;
  int i, j, counter, *counter_array, counter_cplx, environment, modif_cplx, counter_modif_cplx,
      mm_list_size, *mm_list;

#ifdef GMX_MPI
  printf("Initializing charge transfer at rank %d\n", ct_mpi_rank);
#else
  printf("Initializing charge transfer\n");
#endif

  /* Read the file charge-transfer.dat - list of sites */
  f = fopen("charge-transfer.dat", "r");
  if (f == NULL) {
    PRINTF("File charge-transfer.dat not accessible, exiting!\n");
    exit(-1);
  }
  PRINTF("Reading file charge-transfer.dat\n");
  //snew(ct, 1);

  fscanf(f, "%s\n", slko_path);
  PRINTF("SLKO files will be sought in %s\n", slko_path);

  /* Job type */
  getline(&line1, &len, f);
  if (strstr(line1, "SCC") || strstr(line1, "scc") || strstr(line1, "Scc")) {
    ct->jobtype = cteSCCDYNAMIC;
    PRINTF("Fully coupled electron-ion dynamics\n");
  } else {
    if (strstr(line1, "ADI") || strstr(line1, "adi") || strstr(line1, "Adi")) {
      ct->jobtype = cteADIABATIC;
      PRINTF("Adiabatic Born-Oppenheimer (SCF) dynamics of the hole\n");
    } else {
      if (strstr(line1, "NON") || strstr(line1, "non") || strstr(line1, "Non")) {
        ct->jobtype = cteNONSCCDYNAMIC;
        PRINTF("Uncoupled dynamics of the hole - w/o the polarization of solvent\n");
      } else {
        if (strstr(line1, "PAR") || strstr(line1, "par") || strstr(line1, "Par")) {
          ct->jobtype = ctePARAMETERS;
          PRINTF("Calculation of charge-transfer parameters only");
        } else {
          if (strstr(line1, "ADN") || strstr(line1, "adn") || strstr(line1, "Adn")) {
            ct->jobtype = cteADNONSCC; // adiabatic non-self-consistent-charges
            PRINTF("Adiabatic Born-Oppenheimer (SCF) dynamics of the hole w/o the polarization of solvent\n");
          } else {
            if (strstr(line1, "NOM") || strstr(line1, "nom") || strstr(line1, "Nom")) {
              ct->jobtype = cteNOMOVEMENT;
              PRINTF("Stationary charge, calculation of all contributions to hamiltonian\n");
            } else {
              if (strstr(line1, "SFH") || strstr(line1, "sfh") || strstr(line1, "Sfh")) {
                ct->jobtype = cteSURFACEHOPPING;
                PRINTF("Surface hopping between adiabatic surfaces, diabatic limit\n");
              } else {
                if (strstr(line1, "FER") || strstr(line1, "fer") || strstr(line1, "Fer")) {
                  ct->jobtype = cteFERMI;
                  PRINTF("Dynamics with Fermi-distribution-based combination of adiabatic states\n");
                } else {
                  if (strstr(line1, "FAD") || strstr(line1, "fad") || strstr(line1, "Fad")) {
                    ct->jobtype = cteFERMIADIABATIC;
                    PRINTF("Dynamics of the adiabatic ground state obtained from the Fermi-distribution-based combination\n");
                  } else {
                    if (strstr(line1, "FSH") || strstr(line1, "fsh") || strstr(line1, "Fsh")) {
                      ct->jobtype = cteFERMISFHOPPING;
                      PRINTF("Surface hopping between adiabatic states obtained from the Fermi-distribution-based combination, diabatic limit\n");
                    } else {
                      if (strstr(line1, "TFS") || strstr(line1, "tfs") || strstr(line1, "Tfs")) {
                        ct->jobtype = cteTULLYFEWESTSWITCHES;
                        PRINTF("Tully's fewest switches surface hopping between adiabatic states from Fermi-dist. combination\n");
                      } else {
		        if (strstr(line1, "PER") || strstr(line1, "per") || strstr(line1, "Per") || strstr(line1, "PED") || strstr(line1, "ped") || strstr(line1, "Ped")) {
			  ct->jobtype = ctePERSICOSFHOPPING;
			  PRINTF("Persico's locally diabatic surface hopping between adiabatic states from Fermi-dist. combination\n");
                          if (strstr(line1, "PED") || strstr(line1, "ped") || strstr(line1, "Ped")) {
                            ct->decoherence = 1;
			    PRINTF(" - correction for quantum decoherence switched on!\n");
                          } else {
                            ct->decoherence = 0;
                          }
			} else {
			  if (strstr(line1, "NGL") || strstr(line1, "ngl") || strstr(line1, "Ngl")) {
			    ct->jobtype = cteNEGFLORENTZ;
			    PRINTF("Calculation of electric current with non-equlibrium Green's function approach + Lorentzian functions\n");
			    PRINTF(" - populations of molecules mapped onto MD charges (self-consistent calculation)\n");
			  } else {
			    if (strstr(line1, "NGN") || strstr(line1, "ngn") || strstr(line1, "Ngn")) {
			      ct->jobtype = cteNEGFLORENTZNONSCC;
			      PRINTF("Calculation of electric current with non-equlibrium Green's function approach + Lorentzian functions\n");
			      PRINTF(" - no mapping of charges to MD (non-self-consistent calculation)\n");
			    } else {
			      if (strstr(line1, "ESP") || strstr(line1, "esp") || strstr(line1, "Esp")) {
			        ct->jobtype = cteESP;
			        PRINTF("Calculation of electrostatic potential only\n");
			      } else {
                                PRINTF("Unknown job type for charge transfer, use one of:\n");
                                PRINTF("   SCCDYNAMIC, ADIABATIC, NON S C C, PARAMETERS, ADNON S C C, NOMOVEMENT, SFH (surface hopping), FER (dynamic with Fermi distribution),\n");
			        PRINTF("   FAD (gruond state from Fermi mixture), FSH (Fermi-based surface hopping), TFS (Tully's fewest switches),\n");
			        PRINTF("   PER (Persico's surface hopping), PED (Persico's surface hopping with decoherence correction),\n");
			        PRINTF("   NGL (non-equilibrium Green's function calc. of current), NGN (dtto, non-self-consistent calculation),\n");
			        PRINTF("   ESP (calculation of electric potential only).\n");
                                PRINTF("Exiting!\n");
                                exit(-1);
                              }
                            }
                          }
			}
                      }
                    }
                  }
                }
              }
	    }
	  }
	}
      }
    }
  }

  /* how often should the parameters be calculated? */
  if (ct->jobtype == ctePARAMETERS || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteESP) {
    fscanf(f, "%d\n", &(ct->interval));
    if (ct->interval < 1) {
      PRINTF("\nRead a value of %d, setting to 1\nCalculation of charge-transfer parameters (or ESP) only", ct->interval);
      ct->interval = 1;
    }
    PRINTF(" every %d steps.\n", ct->interval);
  }


  /* get the electronic temperature for the Fermi distribution */

  /* QMMM yes/no */
  getline(&line2, &len, f);
  if (strstr(line2, "QMMM") || strstr(line2, "qmmm") || strstr(line2, "QM/MM") || strstr(line2, "qm/mm")) {
    ct->qmmm = 1;
    ct->esp_scaling_factor = 1.;
    PRINTF("QM/MM calculation - the charge-transfer hamiltonian affected by the electric field\n");
  } else {
    if (strstr(line2, "QMG") || strstr(line2, "qmg")) {
      ct->qmmm = 2;
      ct->esp_scaling_factor = 1.;
      PRINTF("QM/MM calculation - the list of MM atoms will be read from file charge-transfer.ndx\n");
      /* read the file "charge-transfer.ndx" here! */
      ct_get_index(&mm_list_size, &mm_list);
      //for (i=0; i<mm_list_size; i++) PRINTF(" %d", mm_list[i]); PRINTF("\n");
    } else {
      if (strstr(line2, "QMS") || strstr(line2, "qms")) {
        ct->qmmm = 1;
        fscanf(f, "%lf\n", &(ct->esp_scaling_factor));
        PRINTF("QM/MM calculation - the electrostatic interaction with MM atoms will be attenuated\n");
        PRINTF("                      by a factor of %f\n", ct->esp_scaling_factor);
      } else {
        if (strstr(line2, "QME") || strstr(line2, "qme")) {
	  ct->qmmm = 3;
	  ct->esp_scaling_factor = 1.;
	  PRINTF("QM/MM calculation with particle--mesh Ewald summation\n");
	} else {
	  if (strstr(line2, "QES") || strstr(line2, "qes")) {
	    ct->qmmm = 3;
	    fscanf(f, "%lf\n", &(ct->esp_scaling_factor));
	    PRINTF("QM/MM calculation with particle--mesh Ewald summation\n");
            PRINTF(" - the electrostatic interaction with MM atoms will be attenuated by a factor of %f\n", ct->esp_scaling_factor);
	  } else {
	    ct->qmmm = 0;
	    PRINTF("\"in vacuo\" calculation - no QM/MM\n");
	  }
	}
      }
    }
  }

  /* self-interaction correction yes or no */
  if (ct->jobtype != ctePARAMETERS && ct->jobtype != cteNOMOVEMENT && ct->jobtype != cteNEGFLORENTZ && ct->jobtype != cteNEGFLORENTZNONSCC &&
      ct->jobtype != cteESP) {
    getline(&line2, &len, f);
    if (strstr(line2, "SIC") || strstr(line2, "sic")) {
      fscanf(f, "%lf\n", &(ct->sic));
      PRINTF("Naive self-interaction correction applied, second-order term scaled by factor %f\n", ct->sic);

    } else {
      ct->sic = 1.0;
      PRINTF("No self-interaction correction\n");
    }
    getline(&line2, &len, f);
    if (strstr(line2, "L_I") || strstr(line2, "l_i")) {
      PRINTF("Emulation of inner-sphere reorganization energy applied\n");
      PRINTF("DNA: LAMBDA_I = %f \n",LAMBDA_I);
      PRINTF("TRP: LAMBDA_I = %f \n",LAMBDA_I_TRYPTOPHAN);
      PRINTF("TYR: LAMBDA_I = %f \n",LAMBDA_I_TYR);
      ct->do_lambda_i = 1;
    } else {
      PRINTF("No inner-sphere reorganization energy\n");
      ct->do_lambda_i = 0;
    }
  } else {
    ct->sic = 1.0;
  }

  /* Number of sites */
  fscanf(f, "%d", &(ct->sites));
  PRINTF("We have %d sites\n", ct->sites);
  snew(ct->site, ct->sites);
  //snew(ct->is_adenine, ct->sites);
  snew(ct->is_dna, ct->sites);
  snew(ct->atoms, ct->sites);
  snew(ct->homo, ct->sites);
  snew(ct->last_atom, ct->sites);
  snew(ct->atom, ct->sites);
  snew(ct->atomtype, ct->sites);
  snew(ct->extcharges, ct->sites);
  snew(ct->extcharge, ct->sites);
  snew(ct->modif_extcharge, ct->sites);
  //ct->hamiltonian = (double **) malloc(ct->sites * sizeof(double *));
  //  ct->hamiltonian[0] = (double *) malloc(SQR(ct->sites) * sizeof(double));
  snew(ct->hamiltonian, ct->sites);
    snew(ct->hamiltonian[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->hamiltonian[j] = ct->hamiltonian[0] + j * ct->sites;
  snew(ct->hamiltonian_mod, ct->sites);
  snew(ct->hamiltonian_adiab, SQR(ct->sites));
  snew(ct->ev_adiab, ct->sites);
  snew(ct->evec_adiab, SQR(ct->sites));
  snew(ct->work_adiab, 3*ct->sites);
  snew(ct->occupation, ct->sites);
  //ct->hubbard = (double **) malloc(ct->sites * sizeof(double *));
  //  ct->hubbard[0] = (double *) malloc(SQR(ct->sites) * sizeof(double));
  snew(ct->hubbard, ct->sites);
    snew(ct->hubbard[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->hubbard[j] = ct->hubbard[0] + j * ct->sites;
  snew(ct->lambda_i, ct->sites);
  snew(ct->delta_q, ct->sites);

  /* Runge-Kutta related */
  ct->rk_neq = 2 * ct->sites;
  snew(ct->wf, ct->rk_neq);
  snew(ct->wf_exc, ct->sites);
  snew(ct->dwf, ct->rk_neq);
  snew(ct->rk_ymax, ct->rk_neq);
  ct->rk_tol = 1.e-8;
  snew(ct->rk_thres, ct->rk_neq);
  for (i=0; i<ct->rk_neq; i++)
    ct->rk_thres[i] = ct->rk_tol;
  ct->rk_lenwrk = 32 * ct->rk_neq;
  snew(ct->rk_work, ct->rk_lenwrk);

  /* NEGF data */
  if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC) {
    snew(ct->negf_arrays, 1);
    ct->negf_arrays->n[0] = ct->sites;
    fscanf(f, "%ld\n", ct->negf_arrays->n_lorentz);
    fscanf(f, "%lf %lf\n", ct->negf_arrays->e_f_left, ct->negf_arrays->e_f_right);
    fscanf(f, "%lf %ld %ld\n", ct->negf_arrays->temp, ct->negf_arrays->n_poles_l, ct->negf_arrays->n_poles_r);
    fscanf(f, "%lf %lf %lf\n", ct->negf_arrays->gam_l, ct->negf_arrays->eps_l, ct->negf_arrays->w0_l);
    fscanf(f, "%lf %lf %lf\n", ct->negf_arrays->gam_r, ct->negf_arrays->eps_r, ct->negf_arrays->w0_r);
    PRINTF("Parameters read for the non-eq. Green's functions:\n");
    PRINTF("  n = %ld, n_lorentz = %ld\n", ct->negf_arrays->n[0], ct->negf_arrays->n_lorentz[0]);
    PRINTF("  e_f_left = %lf, e_f_right = %lf\n", ct->negf_arrays->e_f_left[0], ct->negf_arrays->e_f_right[0]);
    PRINTF("  temp = %lf, n_poles_l = %ld, n_poles_r = %ld\n", ct->negf_arrays->temp[0], ct->negf_arrays->n_poles_l[0], ct->negf_arrays->n_poles_r[0]);
    PRINTF("  gam_l = %lf, eps_l = %lf, w0_l = %lf\n", ct->negf_arrays->gam_l[0], ct->negf_arrays->eps_l[0], ct->negf_arrays->w0_l[0]);
    PRINTF("  gam_r = %lf, eps_r = %lf, w0_r = %lf\n", ct->negf_arrays->gam_r[0], ct->negf_arrays->eps_r[0], ct->negf_arrays->w0_r[0]);
  }

  /* Individual sites (including the number of extcharges) */
  switch (ct->qmmm) {
    case 1:
    case 3:
      ct->extcharges_cplx = top_global->natoms;
      break;
    case 2:
      ct->extcharges_cplx = mm_list_size;
      break;
    default:
      ct->extcharges_cplx = 0;
  }
  ct->atoms_cplx = 0;
  for (i=0; i<ct->sites; i++) {
    fscanf(f, "%d\n", ct->site + i);
    ct->site[i]--;
    if (ct->site[i] < 0 ||
        ct->site[i] >= top_global->moltype->atoms.nres ||
        // ((*(top_global->moltype->atoms.resname[ct->site[i]]))[0] != 'D' &&
        // strcmp((*(top_global->moltype->atoms.resname[ct->site[i]])), "TRP") && 
	// strcmp((*(top_global->moltype->atoms.resname[ct->site[i]])), "TYR")) ) {
        ((*(top_global->moltype->atoms.resinfo[ct->site[i]].name))[0] != 'D' &&
        strcmp((*(top_global->moltype->atoms.resinfo[ct->site[i]].name)), "TRP") && 
	strcmp((*(top_global->moltype->atoms.resinfo[ct->site[i]].name)), "TYR")) ) {
      PRINTF("Illegal site number #%d, exiting!\n", i+1);
      exit(-1);
    }
    // if (!strncmp(*(top_global->moltype->atoms.resname[ct->site[i]]), "DA", 2)) {
    if (!strncmp(*(top_global->moltype->atoms.resinfo[ct->site[i]].name), "DA", 2)) {
      //ct->is_adenine[i] = 1;
      ct->is_dna[i] = 1;
      ct->atoms[i] = 15;
      ct->homo[i] = 25;
      ct->lambda_i[i] = LAMBDA_I;
      ct->hubbard[i][i] = ct->sic * HUBBARD_ADENINE;
      if (ct->do_lambda_i == 1) ct->hubbard[i][i] -= ct->lambda_i[i];
      snew(ct->delta_q[i], ct->atoms[i]);
      ct->delta_q[i][0] =  0.0;
      ct->delta_q[i][1] = -0.0594;
      ct->delta_q[i][2] =  0.2290;
      ct->delta_q[i][3] =  0.0234;
      ct->delta_q[i][4] = -0.0199;
      ct->delta_q[i][5] =  0.1942;
      ct->delta_q[i][6] = -0.0515;
      ct->delta_q[i][7] =  0.2829;
      ct->delta_q[i][8] =  0.0043;
      ct->delta_q[i][9] =  0.0278;
      ct->delta_q[i][10]=  0.0098;
      ct->delta_q[i][11]=  0.1140;
      ct->delta_q[i][12]=  0.0444;
      ct->delta_q[i][13]=  0.1107;
      ct->delta_q[i][14]=  0.0903;

    } else { /* guanine */
    // if (!strncmp(*(top_global->moltype->atoms.resname[ct->site[i]]), "DG", 2)) {
    if (!strncmp(*(top_global->moltype->atoms.resinfo[ct->site[i]].name), "DG", 2)) {
      //ct->is_adenine[i] = 0;
      ct->is_dna[i] = 1;
      ct->atoms[i] = 16;
      ct->homo[i] = 28;
      ct->lambda_i[i] = LAMBDA_I;
      ct->hubbard[i][i] = ct->sic * HUBBARD_GUANINE;
      if (ct->do_lambda_i == 1) ct->hubbard[i][i] -= ct->lambda_i[i];
      snew(ct->delta_q[i], ct->atoms[i]);
      ct->delta_q[i][0] =  0.0;
      ct->delta_q[i][1] = -0.1041;
      ct->delta_q[i][2] =  0.2903;
      ct->delta_q[i][3] =  0.0112;
      ct->delta_q[i][4] = -0.0385;
      ct->delta_q[i][5] =  0.1941;
      ct->delta_q[i][6] = -0.0163;
      ct->delta_q[i][7] =  0.2286;
      ct->delta_q[i][8] = -0.0411;
      ct->delta_q[i][9] =  0.0391;
      ct->delta_q[i][10]=  0.0167;
      ct->delta_q[i][11]=  0.0868;
      ct->delta_q[i][12]=  0.0241;
      ct->delta_q[i][13]=  0.0291;
      ct->delta_q[i][14]=  0.0743;
      ct->delta_q[i][15]=  0.2057;
    } else { /* only A and G implemented; if C and T are needed, a modification is required! */
      /* amino-acid side chains are coming here! */
    // if (!strcmp((*(top_global->moltype->atoms.resname[ct->site[i]])), "TRP")) {
    if (!strcmp((*(top_global->moltype->atoms.resinfo[ct->site[i]].name)), "TRP")) {
      ct->is_dna[i] = 0;
      ct->atoms[i] = 19;
      ct->homo[i] = 25;
      ct->lambda_i[i] = LAMBDA_I_TRYPTOPHAN;
      ct->hubbard[i][i] = ct->sic * HUBBARD_TRYPTOPHAN;
      if (ct->do_lambda_i == 1) ct->hubbard[i][i] -= ct->lambda_i[i];
      snew(ct->delta_q[i], ct->atoms[i]);
      ct->delta_q[i][0] =  0.0;		/* C-alpha*/
      ct->delta_q[i][1] = -0.089202;	/* C-beta + H-dummy(C-alpha) */
      ct->delta_q[i][2] =  0.075236;
      ct->delta_q[i][3] =  0.075236;
      ct->delta_q[i][4] =  0.219781;
      ct->delta_q[i][5] =  0.025729;
      ct->delta_q[i][6] =  0.022458;
      ct->delta_q[i][7] =  0.32752;
      ct->delta_q[i][8] =  0.005057;
      ct->delta_q[i][9] = -0.110575;
      ct->delta_q[i][10]=  0.154573;
      ct->delta_q[i][11]=  0.019495;
      ct->delta_q[i][12]= -0.035344;
      ct->delta_q[i][13]=  0.044093;
      ct->delta_q[i][14]=  0.159628;
      ct->delta_q[i][15]=  0.018352;
      ct->delta_q[i][16]=  0.078786;
      ct->delta_q[i][17]=  0.022967;
      ct->delta_q[i][18]= -0.01379;
    } else {
    // if (!strcmp((*(top_global->moltype->atoms.resname[ct->site[i]])), "TYR")) {
    if (!strcmp((*(top_global->moltype->atoms.resinfo[ct->site[i]].name)), "TYR")) {
      ct->is_dna[i] = 0;
      ct->atoms[i] = 16;
      ct->homo[i] = 21;
      ct->lambda_i[i] = LAMBDA_I_TYR;
      ct->hubbard[i][i] = ct->sic * HUBBARD_TYR;
      if (ct->do_lambda_i == 1) ct->hubbard[i][i] -= ct->lambda_i[i];
      snew(ct->delta_q[i], ct->atoms[i]);
      ct->delta_q[i][0] =  0.0;         /* C-alpha*/
      ct->delta_q[i][1] = -0.037617;	/* C-beta + H-dummy(C-alpha) */
      ct->delta_q[i][2] =  0.059699;	/* HB2 */
      ct->delta_q[i][3] =  0.059699;	/* HB3 */
      ct->delta_q[i][4] =  0.132143;	/* CG  */
      ct->delta_q[i][5] =  0.074927;	/* CD1 */
      ct->delta_q[i][6] =  0.034173;	/* HD1 */
      ct->delta_q[i][7] =  0.05063;	/* CE1 */
      ct->delta_q[i][8] =  0.033783;	/* HE1 */
      ct->delta_q[i][9] =  0.224994;	/* CZ */
      ct->delta_q[i][10]=  0.114956;	/* OH */
      ct->delta_q[i][11]=  0.0591;	/* HH */
      ct->delta_q[i][12]=  0.05063;	/* CE2 */
      ct->delta_q[i][13]=  0.033783;	/* HE2 */
      ct->delta_q[i][14]=  0.074927;	/* CD2 */
      ct->delta_q[i][15]=  0.034173;	/* HD2 */
    } else { /* up to now only TRP and TYR residues for charge transfer*/
      PRINTF("\n  Only ADE, GUA, TRP and TYR implemented yet.\n  Exiting!\n\n");
      exit(-1);
    }}}}
    /* number of extcharges */
    switch (ct->qmmm) {
      case 1:
      case 3:
        ct->extcharges[i] = top_global->natoms - ct->atoms[i] - 1;
        ct->extcharges_cplx -= ct->atoms[i] + 1;
        break;
      case 2:
        /* NUMBER OF EXTCHARGES HERE! */
        ct->extcharges[i] = mm_list_size;
        /* DO NOT SUBTRACT ANYTHING HERE, YET! */
        break;
      default:
        ct->extcharges[i] = 0;
    }
    ct->modif_extcharge[i][0] = ct->modif_extcharge[i][1] = -1;
    ct->atoms_cplx += ct->atoms[i];
    ct->last_atom[i] = ct->atoms_cplx;
    /* print */
    PRINTF("Site %d: Residue #%d (%s)\n", i+1, ct->site[i]+1, (*(top_global->moltype->atoms.resinfo[ct->site[i]].name)));
  }

  /* read occupations if we have NOMOVEMENT */
  /* RATHER, READ IN THE WAVE FUNCTION ALWAYS! except ctePARAMETERS... */
  if (ct->jobtype != ctePARAMETERS && ct->jobtype != cteESP) {
    ct->survival = 0.0;
    for (i=0; i<ct->sites; i++) {
      fscanf(f, "%lf %lf\n", ct->wf + i, ct->wf + ct->sites + i);
      ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i + ct->sites]);
      ct->survival += ct->occupation[i];
    }
    PRINTF("Read the wavefunction:\n");
    for (i=0; i<ct->sites; i++)
      PRINTF(" Re_wf[%d] = %7.4f, Im_wf[%d] = %7.4f\n", i+1, ct->wf[i], i+1, ct->wf[i + ct->sites]);
    PRINTF("Sum of occupations = %7.5f\n", ct->survival);
  }
  /* NEGF initialization including initial density matrix */
  if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC) {
    PRINTF("Initializing the NEGF calculation\n");
#ifdef GMX_MPI
    if (ct_mpi_rank == 0)
#endif
    negf_init_arrays(ct->negf_arrays, &(ct->rk_timestep), ct->wf);
  }

  if (ct->jobtype == cteSCCDYNAMIC || ct->jobtype == cteNONSCCDYNAMIC || ct->jobtype == cteSURFACEHOPPING) {
    getline(&line2, &len, f);
    if (strstr(line2, "NIP") || strstr(line2, "nip")) {
      PRINTF("Charge taken out of the several sites by way of negative imaginary potential,\n");
      fscanf(f, "%d %lf\n", &(ct->neg_imag_pot), &(ct->neg_imag_pot_rate_constant));
      PRINTF("   this will be done for %d sites, with 1/tau = %e au-1 = %f ps-1.\n",
             ct->neg_imag_pot, ct->neg_imag_pot_rate_constant, ct->neg_imag_pot_rate_constant * PS_TO_AU);
      snew(ct->site_neg_imag_pot, ct->neg_imag_pot);
      snew(ct->site_annihilated_occupation, ct->neg_imag_pot);
      for (i=0; i<ct->neg_imag_pot; i++) {
        fscanf(f, "%d\n", ct->site_neg_imag_pot + i);
        ct->site_neg_imag_pot[i]--;
        ct->site_annihilated_occupation[i] = 0.0;
        PRINTF("   site no. %d for NIP occupation removal: residue %d\n",
          i+1, ct->site[ct->site_neg_imag_pot[i]]+1);
      }
//      fscanf(f, "%lf\n", &(ct->survival_threshold));
//      PRINTF("   and quit if the survival drops under %f.\n", ct->survival_threshold);
    } else {
      PRINTF("No negat. imag. potential\n");
      ct->neg_imag_pot = 0;
    }
  }

  PRINTF("Finished reading file charge-transfer.dat\n");
  fclose(f);

  snew(ct->atom_cplx, ct->atoms_cplx);
  snew(ct->atomtype_cplx, ct->atoms_cplx);
  snew(ct->modif_extcharge_cplx, ct->sites * 2);
  snew(ct->shift_extcharge_cplx, ct->sites * 2);
  for (i=0; i<ct->sites*2; i++)
    ct->modif_extcharge_cplx[i] = -1;
  counter_cplx = 0;
  /* Assign the atom numbers that the sites are composed of */
  for (i=0; i<ct->sites; i++) {
    snew(ct->atom[i], ct->atoms[i]);
    snew(ct->atomtype[i], ct->atoms[i]);
    counter = 0;
    PRINTF("Site %d:\n", i+1);
    for (j=0; j<top_global->moltype->atoms.nr; j++) {
      //if (top_global->moltype->atoms.atom[j].resnr == (*(top_global->moltype->atoms.resname[ct->site[i]]))) { /* atom j is in residue that site i corresponds to */
      if (top_global->moltype->atoms.atom[j].resind == ct->site[i]) { /* atom j is in residue that site i corresponds to */ /* TEST: resind instead of resnr */
        if ((ct->is_dna[i] &&
            (!strcmp((*(top_global->moltype->atoms.atomname[j])), "C1q") ||      /* include C1q, assumed to come first! */
	     (!strchr((*(top_global->moltype->atoms.atomname[j])), 'q') &&       /* exclude all other atoms ending with 'q' - backbone */
	      !strchr((*(top_global->moltype->atoms.atomname[j])), 'T') &&       /* exclude end hydrogens */
	      !strchr((*(top_global->moltype->atoms.atomname[j])), 'P')))) ||    /* exclude phosphate groups */
            (!ct->is_dna[i] &&
            ((*(top_global->moltype->atoms.atomname[j]))[1] &&                   /* exclude N, H, C, O */
	     (strcmp(*(top_global->moltype->atoms.atomname[j]), "HA"))))) {      /* exclude HA - this makes the backbone */
	  ct->atom[i][counter] = j;
	  ct->atom_cplx[counter_cplx] = j;
	  switch ((*(top_global->moltype->atoms.atomname[j]))[0]) {
	    case 'C' : ct->atomtype[i][counter] = 0; break;
	    case 'H' : ct->atomtype[i][counter] = 1; break;
	    case 'N' : ct->atomtype[i][counter] = 2; break;
	    case 'O' : ct->atomtype[i][counter] = 3; break;
	    default : PRINTF("Unknown atom type for atom %d (%s), exiting!\n", j, (*(top_global->moltype->atoms.atomname[j]))); exit(-1);
	  }
	  if ((ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "C1q")) || /* the sugar atom C1q will be substituted by a link hydrogen! */
             (!ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "CA")))    /* the backbone atom CA will become a link hydrogen */
	    ct->atomtype[i][counter] = 1;
	  ct->atomtype_cplx[counter_cplx] = ct->atomtype[i][counter];
          PRINTF("%5d (%5s, type %d)\n", ct->atom[i][counter], (*(top_global->moltype->atoms.atomname[ct->atom[i][counter]])), ct->atomtype[i][counter]+1);
	  counter++;
	  counter_cplx++;
	  if (counter > ct->atoms[i]) {
	    PRINTF("Site %d found to have %d atoms, which is more than the expected number of %d, exiting!\n", i, counter, ct->atoms[i]);
	    exit(-1);
	  }
	}
      }
    }
    PRINTF("\n");
  }

  /*
  if (ct->qmmm == 0) {
    for (i=0; i<ct->sites; i++)
      ct->extcharges[i] = 0;
    ct->extcharges_cplx = 0;
  }
  */
  if (ct->qmmm > 0) {
    /* Deal with the external charges */
    for (i=0; i<ct->sites; i++)
      snew(ct->extcharge[i], ct->extcharges[i]);
    snew(ct->extcharge_cplx, ct->extcharges_cplx);
    snew(counter_array, ct->sites);
    counter_cplx = counter_modif_cplx = 0;
 
    for (j=0; j<top_global->natoms; j++) 
    if (ct->qmmm == 1 || ct->qmmm == 3 || ct_atom_in_group(j, mm_list, mm_list_size)) { /* either normal QM/MM or group-QM/MM and j is in the group */
      environment = 1;
      modif_cplx = -1;
      for (i=0; i<ct->sites; i++) {
        if ((j >= top_global->moltype->atoms.nr) ||                              /* atom j is not in molecule 1 -> it is water of ion */
            (top_global->moltype->atoms.atom[j].resind != ct->site[i]) ||        /* atom j is in another residue than i */ /* TEST: resind instead of resnr */
            (ct->is_dna[i] &&
              ((strchr((*(top_global->moltype->atoms.atomname[j])), 'q') ||      /* otherwise we have to look at the name of atom j */
               strchr((*(top_global->moltype->atoms.atomname[j])), 'T') ||
               strchr((*(top_global->moltype->atoms.atomname[j])), 'P')) &&      /* DNA/RNA backbone atoms except C1q and H1q */
               strcmp((*(top_global->moltype->atoms.atomname[j])), "C1q") &&
               strcmp((*(top_global->moltype->atoms.atomname[j])), "H1q"))) ||
            (!ct->is_dna[i] &&
              (!(*(top_global->moltype->atoms.atomname[j]))[1]))) {              /* peptide backbone atoms N, H, C, O */
          ct->extcharge[i][counter_array[i]] = j;
          if (j < top_global->moltype->atoms.nr && top_global->moltype->atoms.atom[j].resind == ct->site[i]) { /* TEST: resind instead of resnr */
                                                                  /* atom j is in molecule 1 (DNA or protein) AND is in the right site */
            if (ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "O4q")) { /* charge of O4q will be modified by +0.06080, to neutralize */
              ct->modif_extcharge[i][0] = counter_array[i];
	      modif_cplx = i;
	    }
            if (ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "C2q")) { /* dtto for C2q */
              ct->modif_extcharge[i][1] = counter_array[i];
	      modif_cplx = i;
            }
            if (!ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "N")) {  /* charge of backbone N will be modified, to electro-neutralize */
              ct->modif_extcharge[i][0] = counter_array[i];
	      modif_cplx = i;
	    }
            if (!ct->is_dna[i] && !strcmp((*(top_global->moltype->atoms.atomname[j])), "C")) {  /* dtto for backbone C */
              ct->modif_extcharge[i][1] = counter_array[i];
	      modif_cplx = i;
            }
	  }
          counter_array[i]++;
          /* debug
          if (i==0 && j<top_global->moltype->atoms.nr && top_global->moltype->atoms.atom[j].resnr==2)
            printf("Atomno %d, resid %d, name %s\n", j+1, top_global->moltype->atoms.atom[j].resnr+1, (*(top_global->moltype->atoms.atomname[j])));
          */
        } else {
          environment = 0;
        }
      }
      if (environment) {
        ct->extcharge_cplx[counter_cplx] = j;
        if (modif_cplx > -1) {
          ct->modif_extcharge_cplx[counter_modif_cplx] = counter_cplx;
          ct->shift_extcharge_cplx[counter_modif_cplx] = ct->is_dna[modif_cplx] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP;
          counter_modif_cplx++;
        }
        counter_cplx++;
      }
    }
    /* check - the number of ext. charges */
    PRINTF("Number of external charges:\n");
    for (i=0; i<ct->sites; i++) {
      PRINTF("Site %2d: original group %d, (possibly) restricted to %d\n", i+1, ct->extcharges[i], counter_array[i]);
      ct->extcharges[i] = counter_array[i];
      if (ct->modif_extcharge[i][0] > -1) PRINTF("          modified extcharge no. %5d: atom %5d - %s (residue %d)\n", ct->modif_extcharge[i][0], ct->extcharge[i][ct->modif_extcharge[i][0]]+1, 
        *(top_global->moltype->atoms.atomname[ct->extcharge[i][ct->modif_extcharge[i][0]]]), top_global->moltype->atoms.atom[ct->extcharge[i][ct->modif_extcharge[i][0]]].resind+1); /* TEST: resind instead of resnr */
      if (ct->modif_extcharge[i][1] > -1) PRINTF("          modified extcharge no. %5d: atom %5d - %s (residue %d)\n", ct->modif_extcharge[i][1], ct->extcharge[i][ct->modif_extcharge[i][1]]+1, 
        *(top_global->moltype->atoms.atomname[ct->extcharge[i][ct->modif_extcharge[i][1]]]), top_global->moltype->atoms.atom[ct->extcharge[i][ct->modif_extcharge[i][1]]].resind+1); /* TEST: resind instead of resnr */
    }
    PRINTF("Complex: original group %d, (possibly) restricted to %d\n", ct->extcharges_cplx, counter_cplx);
    ct->extcharges_cplx = counter_cplx;
    //if (counter_modif_cplx != 2*ct->sites) {
    PRINTF("         number of atoms cutting QM/MM boundary = %d (there are %d sites)\n", counter_modif_cplx, ct->sites);
    //  exit(-1);
    //}
    for (j=0; j<counter_modif_cplx; j++) {
      PRINTF("          modified extcharge no. %5d: atom %5d - %s (residue %d)\n", ct->modif_extcharge_cplx[j], ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]+1, 
        *(top_global->moltype->atoms.atomname[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]]), top_global->moltype->atoms.atom[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]].resind+1); /* TEST: resind instead of resnr */
    }
  } /* if (ct->qmmm > 0) */

  if (ct->jobtype == cteSURFACEHOPPING || ct->jobtype == cteFERMISFHOPPING) {
    for (i=0; i<ct->sites; i++)
      ct->wf_exc[i] = 0.0;
    ct->surface = 0;
    snew(ct->surf_overlap, ct->sites);
    snew(ct->surf_massey, ct->sites);
    snew(ct->surf_prob, ct->sites);
  }

  /* DO HERE PREPARATIONS FOR TFS! */
  if (ct->jobtype == cteTULLYFEWESTSWITCHES) {
    ct->surface = 0;
    snew(ct->tfs_popul, 2*ct->sites); /* complex array: Re(0), Re(1), ..., Re(n-1), Im(0), Im(1), ..., Im(n-1) */
    /* initial conditions - ground state occupied */
    ct->tfs_popul[0] = 1.;
    for (i=1; i<2*ct->sites; i++)
      ct->tfs_popul[i] = 0.;
    snew(ct->tfs_popul_der, 2*ct->sites); /* complex array */
    snew(ct->tfs_vector, ct->sites); /* tfs_vector[n]: n-th eigenvector of the CG Hamiltonian */
    snew(ct->tfs_vector[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_vector[j] = ct->tfs_vector[0] + j * ct->sites;
    snew(ct->tfs_vector_old, ct->sites); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
    snew(ct->tfs_vector_old[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_vector_old[j] = ct->tfs_vector_old[0] + j * ct->sites;
    snew(ct->tfs_overlap, ct->sites); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
    snew(ct->tfs_overlap[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_overlap[j] = ct->tfs_overlap[0] + j * ct->sites;
    snew(ct->surf_prob, ct->sites);
    ct->tfs_initialization_step = 1;
  }

  /* DO HERE PREPARATIONS FOR PERSICO! */
  if (ct->jobtype == ctePERSICOSFHOPPING) {
    ct->surface = 0;
    snew(ct->tfs_popul, 2*ct->sites); /* complex array: Re(0), Re(1), ..., Re(n-1), Im(0), Im(1), ..., Im(n-1) */
    snew(ct->tfs_popul_old, 2*ct->sites);
    /* initial conditions - ground state occupied */
    ct->tfs_popul[0] = 1.;
    for (i=1; i<2*ct->sites; i++)
      ct->tfs_popul[i] = 0.;
    snew(ct->tfs_vector, ct->sites); /* tfs_vector[n]: n-th eigenvector of the CG Hamiltonian */
    snew(ct->tfs_vector[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_vector[j] = ct->tfs_vector[0] + j * ct->sites;
    snew(ct->tfs_vector_old, ct->sites); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
    snew(ct->tfs_vector_old[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_vector_old[j] = ct->tfs_vector_old[0] + j * ct->sites;
    snew(ct->tfs_overlap, ct->sites); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
    snew(ct->tfs_overlap[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->tfs_overlap[j] = ct->tfs_overlap[0] + j * ct->sites;
    snew(ct->per_propag_operator, ct->sites); /* */
    snew(ct->per_propag_operator[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->per_propag_operator[j] = ct->per_propag_operator[0] + j * ct->sites;
    snew(ct->per_transformator, ct->sites); /* */
    snew(ct->per_transformator[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->per_transformator[j] = ct->per_transformator[0] + j * ct->sites;
    snew(ct->per_diab_hamiltonian, ct->sites); /* */
    snew(ct->per_diab_hamiltonian[0], SQR(ct->sites));
    for(j = 1; j < ct->sites; j++)
      ct->per_diab_hamiltonian[j] = ct->per_diab_hamiltonian[0] + j * ct->sites;
    snew(ct->surf_prob, ct->sites);
    snew(ct->ev_adiab_old, ct->sites);
    ct->tfs_initialization_step = 1;
    /* auxiliary arrays for the orthogonalizer */
    snew(ct->per_arrays, 1);
    snew(ct->per_arrays->in, ct->sites*ct->sites);
    snew(ct->per_arrays->evec, ct->sites*ct->sites);
    ct->per_arrays->lwork = 26*ct->sites;
    snew(ct->per_arrays->work, 26*ct->sites*26*ct->sites);
    ct->per_arrays->liwork = 10*ct->sites;
    snew(ct->per_arrays->iwork, 10*ct->sites*10*ct->sites);
    snew(ct->per_arrays->eval, ct->sites);
    snew(ct->per_arrays->issupz, 2*ct->sites);
  }

  return;
}

/****************************
 * INITIALIZE THE DFTB CODE *
 ****************************/

#ifdef GMX_MPI
void init_dftb(gmx_mtop_t *top_global, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path, int ct_mpi_rank)
#else
void init_dftb(gmx_mtop_t *top_global, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path)
#endif
{
  // const char *slko_path = "/home/tomas/DFTB/slko/";
  const char *type_symbols = "chno";
  const char *suffix1 = "-c.spl";
  const char *suffix2 = "-8-7-c.spl";

  char filename[128];
  int i, j, k, l, izpj, counter;
  double espin, qzeroh[3], uhubbh[3], mass;
  FILE *f;

  //snew(dftb, 1);

  /* assign the number of shells in atom types */
  dftb->lmax[0] = 2;
  dftb->lmax[1] = 1;
  dftb->lmax[2] = 2;
  dftb->lmax[3] = 2;

  /*
  snew(dftb->qzero1, dftb_maxtypes);
  snew(dftb->uhubb1, dftb_maxtypes);
  snew(dftb->qzero2, dftb_maxtypes);
  snew(dftb->uhubb2, dftb_maxtypes);
  */

  for (i=0; i<DFTB_MAXTYPES; i++) {
    dftb->qzero1[i] = 0.0;
    dftb->uhubb1[i] = 0.0;
    dftb->qzero2[i] = 0.0;
    dftb->uhubb2[i] = 0.0;
    for (j=0; j<DFTB_MAXTYPES; j++) {
      PRINTF("Atomtype pair %d-%d\n", i+1, j+1);
      /* read the tables for DFTB phase 1 - calculation of monomers */
      sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix1);
      f = fopen(filename, "r");
      if (f == NULL) {
        PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
	exit(-1);
      }
      //printf("fscanf(f, \"%%lf %%d\", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]))\n");
      fscanf(f, "%lf %d", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]));
      if (i == j) {
        //printf("&(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin\n");
        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin,
                  uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
        dftb->uhubb1[i] = uhubbh[0];
        for (k=0; k<3; k++)
          dftb->qzero1[i] += qzeroh[k];
      }
      snew(dftb->skhtab1[i][j], dftb->dim1[i][j]);
      snew(dftb->skstab1[i][j], dftb->dim1[i][j]);
      for (k=0; k<dftb->dim1[i][j]; k++) {
        //printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skhtab1[i][j][k][l]))\n");
        for (l=0; l<10; l++) fscanf(f, "%lf", &(dftb->skhtab1[i][j][k][l]));
        //printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skstab1[i][j][k][l]))\n");
        for (l=0; l<10; l++) fscanf(f, "%lf", &(dftb->skstab1[i][j][k][l]));
      }
      PRINTF("skfile for pair %d-%d, phase 1: %s\n", i+1, j+1, filename);
      fclose(f);

      /* read the tables for DFTB phase 2 - calculation of monomers */
      sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix2);
      f = fopen(filename, "r");
      if (f == NULL) {
        PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
	exit(-1);
      }
      fscanf(f, "%lf %d", &(dftb->dr2[i][j]), &(dftb->dim2[i][j]));
      if (i == j) {
        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &(dftb->skself2[i][0]), &(dftb->skself2[i][1]), &(dftb->skself2[i][2]), &espin,
                  uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
        dftb->uhubb2[i] = uhubbh[0];
        for (k=0; k<3; k++)
          dftb->qzero2[i] += qzeroh[k];
      }
      snew(dftb->skhtab2[i][j], dftb->dim2[i][j]);
      snew(dftb->skstab2[i][j], dftb->dim2[i][j]);
      for (k=0; k<dftb->dim2[i][j]; k++) {
        for (l=0; l<10; l++) fscanf(f, "%lf", &(dftb->skhtab2[i][j][k][l]));
        for (l=0; l<10; l++) fscanf(f, "%lf", &(dftb->skstab2[i][j][k][l]));
      }
      PRINTF("skfile for pair %d-%d, phase 2: %s\n", i+1, j+1, filename);
      fclose(f);
    }
  }

  /* deal with phase1 and certain parts of phase2 */
  snew(dftb->phase1, ct->sites);
  //snew(dftb->phase2, 1);
  dftb->phase2.nn = 0;
  dftb->phase2.ne = ct->extcharges_cplx;
  dftb->phase2.nel = 0;
  dftb->phase2.norb = 0;

  /* prepare the broyden structures */
  snew(dftb->broyden, ct->sites);
  for (i=0; i<ct->sites; i++) {
    snew(dftb->broyden[i].f, ct->atoms[i]);
    snew(dftb->broyden[i].ui, ct->atoms[i]);
    snew(dftb->broyden[i].vti, ct->atoms[i]);
    snew(dftb->broyden[i].t1, ct->atoms[i]);
    snew(dftb->broyden[i].dumvi, ct->atoms[i]);
    snew(dftb->broyden[i].df, ct->atoms[i]);
    snew(dftb->broyden[i].vector, ct->atoms[i]);
    snew(dftb->broyden[i].unit31, ct->atoms[i]);
    snew(dftb->broyden[i].unit32, ct->atoms[i]);
  }

  /* determine the numbers of electrons and atom types (in arrays) */
  dftb->phase2.nel = 0;
  for (i=0; i<ct->sites; i++) {
    if (ct->is_dna[i]) {
      // if (!strncmp(*(top_global->moltype->atoms.resname[ct->site[i]]), "DA", 2)) { /* adenine */
      if (!strncmp(*(top_global->moltype->atoms.resinfo[ct->site[i]].name), "DA", 2)) { /* adenine */
        dftb->phase1[i].nel = 50;
        dftb->phase1[i].norb = 45;
      } else {                                                                     /* guanine */
        dftb->phase1[i].nel = 56;
        dftb->phase1[i].norb = 49;
      }
    } else {
      // if (!strcmp(*(top_global->moltype->atoms.resname[ct->site[i]]), "TRP")) {    /* tryptophan */
      if (!strcmp(*(top_global->moltype->atoms.resinfo[ct->site[i]].name), "TRP")) {    /* tryptophan */
        dftb->phase1[i].nel = 50;
        dftb->phase1[i].norb = 49;
      } else {                                                                     /* tyrosine */
        dftb->phase1[i].nel = 42;
        dftb->phase1[i].norb = 40;
      }
    }
    dftb->phase1[i].nn = ct->atoms[i];
    dftb->phase2.nn += ct->atoms[i];
    dftb->phase2.nel += dftb->phase1[i].nel;
    dftb->phase2.norb += dftb->phase1[i].norb;

    /* double arrays */
    snew(dftb->phase1[i].x, ct->atoms[i]);
    snew(dftb->phase1[i].xe, ct->extcharges[i]);
    snew(dftb->phase1[i].mass, ct->atoms[i]);
    snew(dftb->phase1[i].ze, ct->extcharges[i]);
    snew(dftb->phase1[i].qmat, ct->atoms[i]);
    snew(dftb->phase1[i].qmold, ct->atoms[i]);
    snew(dftb->phase1[i].qmulli, dftb->phase1[i].norb);
    snew(dftb->phase1[i].ev, dftb->phase1[i].norb);
    snew(dftb->phase1[i].occ, dftb->phase1[i].norb);
    //dftb->phase1[i].a = (double **) malloc(dftb->phase1[i].norb * sizeof(double *));
    //  dftb->phase1[i].a[0] = (double *) malloc(dftb->phase1[i].norb * dftb->phase1[i].norb * sizeof(double));
    snew(dftb->phase1[i].a, dftb->phase1[i].norb);
      snew(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb));
      for(j = 1; j < dftb->phase1[i].norb; j++)
        dftb->phase1[i].a[j] = dftb->phase1[i].a[0] + j * dftb->phase1[i].norb;
    snew(dftb->phase1[i].a_old, dftb->phase1[i].norb);
      snew(dftb->phase1[i].a_old[0], SQR(dftb->phase1[i].norb));
      for (j=0; j<SQR(dftb->phase1[i].norb); j++)
        dftb->phase1[i].a_old[0][j] = 0.;
      for(j = 1; j < dftb->phase1[i].norb; j++)
        dftb->phase1[i].a_old[j] = dftb->phase1[i].a_old[0] + j * dftb->phase1[i].norb;
    //dftb->phase1[i].b = (double **) malloc(dftb->phase1[i].norb * sizeof(double *));
    //  dftb->phase1[i].b[0] = (double *) malloc(dftb->phase1[i].norb * dftb->phase1[i].norb * sizeof(double));
    snew(dftb->phase1[i].b, dftb->phase1[i].norb);
      snew(dftb->phase1[i].b[0], SQR(dftb->phase1[i].norb));
      for(j = 1; j < dftb->phase1[i].norb; j++)
        dftb->phase1[i].b[j] = dftb->phase1[i].b[0] + j * dftb->phase1[i].norb;
    snew(dftb->phase1[i].a_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
    snew(dftb->phase1[i].b_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
    //dftb->phase1[i].hamil = (double **) malloc(dftb->phase1[i].norb * sizeof(double *));
    //  dftb->phase1[i].hamil[0] = (double *) malloc(dftb->phase1[i].norb * dftb->phase1[i].norb * sizeof(double));
    snew(dftb->phase1[i].hamil, dftb->phase1[i].norb);
      snew(dftb->phase1[i].hamil[0], SQR(dftb->phase1[i].norb));
      for(j = 1; j < dftb->phase1[i].norb; j++)
        dftb->phase1[i].hamil[j] = dftb->phase1[i].hamil[0] + j * dftb->phase1[i].norb;
    //dftb->phase1[i].overl = (double **) malloc(dftb->phase1[i].norb * sizeof(double *));
    //  dftb->phase1[i].overl[0] = (double *) malloc(dftb->phase1[i].norb * dftb->phase1[i].norb * sizeof(double));
    snew(dftb->phase1[i].overl, dftb->phase1[i].norb);
      snew(dftb->phase1[i].overl[0], SQR(dftb->phase1[i].norb));
      for(j = 1; j < dftb->phase1[i].norb; j++)
        dftb->phase1[i].overl[j] = dftb->phase1[i].overl[0] + j * dftb->phase1[i].norb;
    //dftb->phase1[i].gammamat = (double **) malloc(ct->atoms[i] * sizeof(double *));
    //  dftb->phase1[i].gammamat[0] = (double *) malloc(ct->atoms[i] * ct->atoms[i] * sizeof(double));
    snew(dftb->phase1[i].gammamat, ct->atoms[i]);
      snew(dftb->phase1[i].gammamat[0], SQR(ct->atoms[i]));
      for(j = 1; j < ct->atoms[i]; j++)
        dftb->phase1[i].gammamat[j] = dftb->phase1[i].gammamat[0] + j * ct->atoms[i];
    snew(dftb->phase1[i].shift, ct->atoms[i]);
    snew(dftb->phase1[i].shiftE, ct->atoms[i]);
    snew(dftb->phase1[i].aux, 3 * dftb->phase1[i].norb);

    if (ct->qmmm == 3) {
      snew(dftb->phase1[i].neighbors_pme, ct->atoms[i]);
      snew(dftb->phase1[i].neighbor_pme, ct->atoms[i]);
    }

    /* int arrays */
    snew(dftb->phase1[i].izp, ct->atoms[i]);
    snew(dftb->phase1[i].ind, ct->atoms[i] + 1);
    dftb->phase1[i].ind[0] = 0;
    for (j=0; j<ct->atoms[i]; j++) {
      dftb->phase1[i].izp[j] = ct->atomtype[i][j];
      izpj = dftb->phase1[i].izp[j];
      dftb->phase1[i].ind[j+1] = dftb->phase1[i].ind[j] + dftb->lmax[izpj]* dftb->lmax[izpj];
    }
    dftb->phase1[i].ndim = dftb->phase1[i].norb;
    dftb->phase1[i].ne = ct->extcharges[i];

    /* atom masses */
    mass = 0.0;
    for (j=0; j<ct->atoms[i]; j++) {
      dftb->phase1[i].mass[j] = mdatoms->massT[ct->atom[i][j]];
      mass += dftb->phase1[i].mass[j];
    }
    dftb->phase1[i].inv_tot_mass = 1.0 / mass;
  }

  /* the remainder of phase2 */
  snew(dftb->phase2.x, dftb->phase2.nn);
  snew(dftb->phase2.xe, dftb->phase2.ne);
  snew(dftb->phase2.mass, dftb->phase2.nn);
  snew(dftb->phase2.ze, dftb->phase2.ne);
  snew(dftb->phase2.qmat, dftb->phase2.nn);
  snew(dftb->phase2.shift, dftb->phase2.nn);
  snew(dftb->phase2.shiftE, dftb->phase2.nn);
  snew(dftb->phase2.ind, dftb->phase2.nn + 1);
  snew(dftb->phase2.inf, dftb->phase2.nn);
  snew(dftb->phase2.izp, dftb->phase2.nn);
  //dftb->phase2.hamil = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.hamil[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.hamil, dftb->phase2.norb);
    snew(dftb->phase2.hamil[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.hamil[j] = dftb->phase2.hamil[0] + j * dftb->phase2.norb;
  //dftb->phase2.overl = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.overl[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.overl, dftb->phase2.norb);
    snew(dftb->phase2.overl[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.overl[j] = dftb->phase2.overl[0] + j * dftb->phase2.norb;
  //dftb->phase2.Taf = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.Taf[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.Taf, dftb->phase2.norb);
    snew(dftb->phase2.Taf[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.Taf[j] = dftb->phase2.Taf[0] + j * dftb->phase2.norb;
  //dftb->phase2.THamil = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.THamil[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.THamil, dftb->phase2.norb);
    snew(dftb->phase2.THamil[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamil[j] = dftb->phase2.THamil[0] + j * dftb->phase2.norb;
  //dftb->phase2.OverlF = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.OverlF[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.OverlF, dftb->phase2.norb);
    snew(dftb->phase2.OverlF[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.OverlF[j] = dftb->phase2.OverlF[0] + j * dftb->phase2.norb;
  //dftb->phase2.THamilOrtho = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.THamilOrtho[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.THamilOrtho, dftb->phase2.norb);
    snew(dftb->phase2.THamilOrtho[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamilOrtho[j] = dftb->phase2.THamilOrtho[0] + j * dftb->phase2.norb;
  //dftb->phase2.tij = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.tij[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.tij, dftb->phase2.norb);
    snew(dftb->phase2.tij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.tij[j] = dftb->phase2.tij[0] + j * dftb->phase2.norb;
  //dftb->phase2.sij = (double **) malloc(dftb->phase2.norb * sizeof(double *));
  //  dftb->phase2.sij[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));
  snew(dftb->phase2.sij, dftb->phase2.norb);
    snew(dftb->phase2.sij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.sij[j] = dftb->phase2.sij[0] + j * dftb->phase2.norb;
  //dftb->phase2.gammamat = (double **) malloc(ct->atoms_cplx * sizeof(double *));
  //  dftb->phase2.gammamat[0] = (double *) malloc(SQR(ct->atoms_cplx) * sizeof(double));
  snew(dftb->phase2.gammamat, ct->atoms_cplx);
    snew(dftb->phase2.gammamat[0], SQR(ct->atoms_cplx));
    for(j = 1; j < ct->atoms_cplx; j++)
      dftb->phase2.gammamat[j] = dftb->phase2.gammamat[0] + j * ct->atoms_cplx;
  //snew(dftb->phase2.shift, dftb->phase2.nn); this was done above already!
  //snew(dftb->phase2.shiftE, dftb->phase2.nn);
  snew(dftb->phase2.ev, dftb->phase2.norb);
  snew(dftb->phase2.occ, dftb->phase2.norb);
  snew(dftb->phase2.aux, 3 * dftb->phase2.norb);

  if (ct->qmmm == 3) {
    snew(dftb->phase2.neighbors_pme, ct->atoms_cplx);
    snew(dftb->phase2.neighbor_pme, ct->atoms_cplx);
  }

  dftb->phase2.inf[0] = 0; // where do the orbitals of fragment i start?
  for (i=1; i<ct->sites; i++)
    dftb->phase2.inf[i] = dftb->phase2.inf[i-1] + dftb->phase1[i-1].norb;

  counter = 0;
  dftb->phase2.ind[0] = 0;
  mass = 0.0;
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->atoms[i]; j++) {
      dftb->phase2.izp[counter] = dftb->phase1[i].izp[j];
      izpj = dftb->phase2.izp[counter];
      dftb->phase2.ind[counter+1] = dftb->phase2.ind[counter] + dftb->lmax[izpj]* dftb->lmax[izpj];
      dftb->phase2.mass[counter] = mdatoms->massT[ct->atom[i][j]]; /* mass */
      mass += dftb->phase2.mass[counter];
      counter++;
    }
  dftb->phase2.inv_tot_mass = 1.0 / mass;
  if (counter != dftb->phase2.nn) {
    PRINTF("The number of atoms does not match: counter = %d, dftb->phase2.nn = %d !\n", counter, dftb->phase2.nn);
    exit(-1);
  }

  /* auxiliary arrays for function orthogonalize() */
  snew(dftb->orthogo.tij, ct->sites*ct->sites);
  snew(dftb->orthogo.sij, ct->sites*ct->sites);
  snew(dftb->orthogo.evec, ct->sites*ct->sites);
  dftb->orthogo.lwork = 26*ct->sites;
  snew(dftb->orthogo.work, 26*ct->sites*26*ct->sites);
  dftb->orthogo.liwork = 10*ct->sites;
  snew(dftb->orthogo.iwork, 10*ct->sites*10*ct->sites);
  snew(dftb->orthogo.eval, ct->sites);
  snew(dftb->orthogo.issupz, 2*ct->sites);

  /* machine accuracy */
  dftb->racc = 1.0;
  while ((1.0 + dftb->racc) > 1.0)
    dftb->racc /= 2.0;
  dftb->racc *= 2.0;
  dftb->dacc = 4 * dftb->racc;

  return;
}

#ifdef GMX_MPI
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir, int ct_mpi_rank)
#else
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir)
#endif
{
  int i, j, status;
  dftb->ewaldcoeff_pme = calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
  dftb->rcoulomb_pme = ir->rcoulomb;
  dftb->nstlist_pme = ir->nstlist;
  dftb->lastlist_pme = -1000000; /* to have the list generated before the first PME calculation (in any case) */
  /* for phase 1 -- loop over all sites
   * POSSIBLE MODIFICATION FOR PARALLEL RUNS
   *   -- DO IT ONLY FOR THOSE SITES THAT ARE TO BE CALCULATED ON THIS RESPECTIVE CT_MPI_RANK!
   */
  for (i=0; i<ct->sites; i++) {
    snew(dftb->phase1[i].q_pme, ct->atoms[i] + ct->extcharges[i]);
    snew(dftb->phase1[i].x_pme, ct->atoms[i] + ct->extcharges[i]);
    snew(dftb->phase1[i].pot, ct->atoms[i] + ct->extcharges[i]);
    snew(dftb->phase1[i].pot2, ct->atoms[i]);
    snew(dftb->phase1[i].pot3, ct->atoms[i]);
    snew(dftb->phase1[i].pot4, ct->atoms[i]);
    snew(dftb->phase1[i].nrnb_pme, 1);
    for (j=0; j<ct->atoms[i]; j++)
      dftb->phase1[i].qmat[j] = dftb->qzero1[dftb->phase1[i].izp[j]];
    snew(dftb->phase1[i].pmedata, 1);
    status = gmx_pme_init_dftb(dftb->phase1[i].pmedata, ir, ct->atoms[i] + ct->extcharges[i]);
  }
  /* for phase 2 -- POSSIBLY ONLY ON CT_MPI_RANK==0
   * ! */
  snew(dftb->phase2.q_pme, ct->atoms_cplx + ct->extcharges_cplx);
  snew(dftb->phase2.x_pme, ct->atoms_cplx + ct->extcharges_cplx);
  snew(dftb->phase2.pot, ct->atoms_cplx + ct->extcharges_cplx);
  snew(dftb->phase2.pot2, ct->atoms_cplx);
  snew(dftb->phase2.pot3, ct->atoms_cplx);
  snew(dftb->phase2.pot4, ct->atoms_cplx);
  snew(dftb->phase2.nrnb_pme, 1);
  snew(dftb->phase2.pmedata, 1);
  status = gmx_pme_init_dftb(dftb->phase2.pmedata, ir, ct->atoms_cplx + ct->extcharges_cplx);
  return;
}

/****************************
 * INITIALIZE THE DIIS CODE *
 ****************************/

void ct_init_diis(charge_transfer_t *ct, ct_diis_t *diis)
{
  int i;

  diis->n_elem = ct->sites;
  //diis->prev_q_input = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
  //  diis->prev_q_input[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
  snew(diis->prev_q_input, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_input[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (i=1; i<DIIS_MAX_PREV_VECTORS; i++)
      diis->prev_q_input[i] = diis->prev_q_input[0] + i * diis->n_elem;
  //diis->prev_q_diff = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
  //  diis->prev_q_diff[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
  snew(diis->prev_q_diff, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_diff[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (i=1; i<DIIS_MAX_PREV_VECTORS; i++)
      diis->prev_q_diff[i] = diis->prev_q_diff[0] + i * diis->n_elem;
  //diis->aa = (double *) malloc(SQR(DIIS_MAX_PREV_VECTORS + 1));
  //diis->bb = (double *) malloc(DIIS_MAX_PREV_VECTORS + 1);
  //diis->q_inp_result = (double *) malloc(ct->sites * sizeof(double));
  //diis->q_diff = (double *) malloc(ct->sites * sizeof(double));
  snew(diis->aa, SQR(DIIS_MAX_PREV_VECTORS + 1));
  snew(diis->bb, DIIS_MAX_PREV_VECTORS + 1);
  snew(diis->q_inp_result, ct->sites);
  snew(diis->q_diff, ct->sites);
  snew(diis->fermi_coef, ct->sites);

  return;
}

void ct_init_broyden(charge_transfer_t *ct, dftb_broyden_t *broyd)
{
  snew(broyd->f,      ct->sites);
  snew(broyd->ui,     ct->sites);
  snew(broyd->vti,    ct->sites);
  snew(broyd->t1,     ct->sites);
  snew(broyd->dumvi,  ct->sites);
  snew(broyd->df,     ct->sites);
  snew(broyd->vector, ct->sites);
  snew(broyd->unit31, ct->sites);
  snew(broyd->unit32, ct->sites);
  
  return;
}

/******************************************
 * PREREQUISITIES FOR A QM/MM CALCULATION *
 *    TO BE PERFORMED IN EVERY MD STEP    *
 ******************************************/

void prepare_charge_transfer(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct)
{
  int i, j, k, counter;
  ivec shift, shiftmin;
  double bond_length, mindist, curdist, sum;
  dvec bond, box, image, com, coord, masscoord, r;
  dftb_phase1_t dftb1;
  dftb_phase2_t dftb2;

//#ifdef GMX_MPI
//        for (i=0; i<10; i++)
//        {
//            printf("   INSIDE  rankX atom%d %12.7f%12.7f%12.7f\n", i,x_ct[i][XX]*10, x_ct[i][YY]*10, x_ct[i][ZZ]*10);
//        }
//        for (i=0; i<10; i++)
//        {
//            printf("   INSIDE  rankX atom%d Q= %12.7f\n", i, mdatoms->chargeA[i]);
//        }
//#endif

  /* debug begin
  printf("Information about system\n");
  printf("Number of atoms: %d\n", mdatoms->nr);
  printf("mdatoms->massA = %p, mdatoms->massT = %p, mdatoms->chargeA = %p\n", mdatoms->massA, mdatoms->massT, mdatoms->chargeA);
  printf("Selected atoms - massT, charge:\n");
  for (i=0; i<10; i++)
    printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
  for (i=7270; i<7280; i++)
    printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
  for (i=8270; i<8280; i++)
    printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
     debug end */

  // printf("prepare_charge_transfer\n");

  /* read the box dimensions */
  for (j=0; j<DIM; j++)
    box[j] = state_box[j][j] * NM_TO_BOHR;
//  printf("BOX     %12.7f %12.7f %12.7f\n", box[XX], box[YY], box[ZZ]);
  if (ct->qmmm == 3) {
    copy_mat(state_box, dftb->box_pme);
//    printf("BOX_PME %12.7f %12.7f %12.7f\n", dftb->box_pme[XX][XX], dftb->box_pme[YY][YY], dftb->box_pme[ZZ][ZZ]);
  }

  /* read coordinates of the quantum system */
  /* conversion from nanometer to bohr */
  counter = 0;
  for (i=0; i<ct->sites; i++) {
    for (j=0; j<ct->atoms[i]; j++) {
      for (k=0; k<3; k++)
        dftb->phase1[i].x[j][k] = x_ct[ct->atom[i][j]][k] * NM_TO_BOHR;
    }
    /* the link hydrogen atom */
    dvec_sub(dftb->phase1[i].x[0], dftb->phase1[i].x[1], bond);
    bond_length = dnorm(bond);
    for (k=0; k<DIM; k++)
      dftb->phase1[i].x[0][k] = dftb->phase1[i].x[1][k] + bond[k] *
        (ct->is_dna[i] ? NH_BOND_LENGTH : CH_BOND_LENGTH) / bond_length;
    /* copy to the phase2 - coordinates of the complex */
    for (j=0; j<ct->atoms[i]; j++) {
      copy_dvec(dftb->phase1[i].x[j], dftb->phase2.x[counter]);
      counter++;
    }
  }

  /* test - write out the coordinates
  //for (i=0; i<ct->sites; i++) {
  i=0;
    printf("Site %d - %d atoms\n", i+1, dftb->phase1[i].nn);
    for (j=0; j<dftb->phase1[i].nn; j++)
      //printf("%d %12.7f%12.7f%12.7f\n", dftb->phase1[i].izp[j]+1, dftb->phase1[i].x[j][0]*0.52, dftb->phase1[i].x[j][1]*0.52, dftb->phase1[i].x[j][2]*0.52);
      printf("C %12.7f%12.7f%12.7f\n", dftb->phase1[i].x[j][0]*0.52, dftb->phase1[i].x[j][1]*0.52, dftb->phase1[i].x[j][2]*0.52);
  //}

  printf("Complex - %d atoms\n", dftb->phase2.nn);
  for (j=0; j<dftb->phase2.nn; j++) {
    switch (dftb->phase2.izp[j]) {
      case 0: c = 'C'; break;
      case 1: c = 'H'; break;
      case 2: c = 'N'; break;
      case 3: c = 'O'; break;
    }
    printf("%c %12.7f%12.7f%12.7f\n", c, dftb->phase2.x[j][0]*0.52, dftb->phase2.x[j][1]*0.52, dftb->phase2.x[j][2]*0.52);
  }
  */

  /* get the center of mass of every fragment as well as of the complex */
  for (i=0; i<ct->sites; i++) {
    clear_dvec(com);
    for (j=0; j<ct->atoms[i]; j++) {
      coord[XX] = dftb->phase1[i].x[j][XX];
      coord[YY] = dftb->phase1[i].x[j][YY];
      coord[ZZ] = dftb->phase1[i].x[j][ZZ];
      dsvmul(mdatoms->massT[ct->atom[i][j]], coord, masscoord);
      dvec_inc(com, masscoord);
    }
    // dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
    dsvmul(dftb->phase1[i].inv_tot_mass, com, dftb->phase1[i].com);
//    printf("COM base %d: %f %f %f\n", i+1, dftb->phase1[i].com[XX] * 0.52, dftb->phase1[i].com[YY] * 0.52, dftb->phase1[i].com[ZZ] * 0.52);
  }

  // construct the Hubbard matrix
  if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteSURFACEHOPPING || ct->jobtype == cteFERMI ||
      ct->jobtype == cteFERMIADIABATIC || ct->jobtype == cteFERMISFHOPPING || ct->jobtype == cteTULLYFEWESTSWITCHES || ct->jobtype == ctePERSICOSFHOPPING) {
    // diag: Hubbard params of nucleobases (constant - do not calculate)
    // offdiag: 1/r_ij for nucleobases
    for (i=0; i<ct->sites; i++)
      for (j=i+1; j<ct->sites; j++) {
        // approximation of nucleobases as single charged points
        dvec_sub(dftb->phase1[i].com, dftb->phase1[j].com, bond);
        ct->hubbard[i][j] = ct->hubbard[j][i] = ct->sic / dnorm(bond);
      }
  } else {
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->hubbard[i][j] = 0.e0;
  }

  /* read coordinates and magnitudes of the external charges */
  /* attention - consider the extcharges to be in the nearest periodic image! */
  if (ct->qmmm > 0) { /* CONTINUE HERE FOR PME !!! */
    for (i=0; i<ct->sites; i++) {
      for (j=0; j<ct->extcharges[i]; j++) {
        /* coordinates */
        for (k=0; k<DIM; k++)
          dftb->phase1[i].xe[j][k] = NM_TO_BOHR * x_ct[ct->extcharge[i][j]][k];
	if (ct->qmmm < 3) {
	  /* when not doing particle--mesh Ewald: */
          /* identify the nearest periodic image */
          /* attention - doesn't necessarily keep molecules whole
              - possibly induces artificial dipole !!! */
          /* correction - it seems that Gromacs keeps molecules whole... */
        mindist = 1.e10;
        clear_ivec(shiftmin);
        for (shift[XX]=-1; shift[XX]<=1; shift[XX]++)
          for (shift[YY]=-1; shift[YY]<=1; shift[YY]++)
            for (shift[ZZ]=-1; shift[ZZ]<=1; shift[ZZ]++) {
              	image[XX] = dftb->phase1[i].xe[j][XX] + shift[XX] * box[XX];
              	image[YY] = dftb->phase1[i].xe[j][YY] + shift[YY] * box[YY];
              	image[ZZ] = dftb->phase1[i].xe[j][ZZ] + shift[ZZ] * box[ZZ];
              	dvec_sub(image, dftb->phase1[i].com, bond);
              	curdist = dnorm(bond);
              	if (curdist < mindist) {
              	  mindist = curdist;
                        copy_ivec(shift, shiftmin);
              	}
            }
        dftb->phase1[i].xe[j][XX] += shiftmin[XX] * box[XX];
        dftb->phase1[i].xe[j][YY] += shiftmin[YY] * box[YY];
        dftb->phase1[i].xe[j][ZZ] += shiftmin[ZZ] * box[ZZ];
	}
        /* magnitude, then add 0.06080 to O4a and C2q */
        dftb->phase1[i].ze[j] = mdatoms->chargeA[ct->extcharge[i][j]];
      }
      if (ct->modif_extcharge[i][0] > -1)
        dftb->phase1[i].ze[ct->modif_extcharge[i][0]] += ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP;
      if (ct->modif_extcharge[i][1] > -1)
        dftb->phase1[i].ze[ct->modif_extcharge[i][1]] += ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP;
    }
  }

  /* test - write out the extcharges for nucleobase 1
  i=0;
    for (j=0; j<10; j++) // ct->extcharges[i]; j++)
      //printf("H %12.7f%12.7f%12.7f\n", dftb->phase1[i].xe[j][XX]*0.52, dftb->phase1[i].xe[j][YY]*0.52, dftb->phase1[i].xe[j][ZZ]*0.52);
      printf("atom%d %12.7f%12.7f%12.7f\n", j, x_ct[j][XX]*10, x_ct[j][YY]*10, x_ct[j][ZZ]*10);
 end test */

  /* center of mass of the complex */
  clear_dvec(com);
  for (j=0; j<ct->atoms_cplx; j++) {
    coord[XX] = dftb->phase2.x[j][XX];
    coord[YY] = dftb->phase2.x[j][YY];
    coord[ZZ] = dftb->phase2.x[j][ZZ];
    dsvmul(mdatoms->massT[ct->atom_cplx[j]], coord, masscoord);
    dvec_inc(com, masscoord);
  }
  // dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
  dsvmul(dftb->phase2.inv_tot_mass, com, dftb->phase2.com);
//  printf(" PHASE 2 COM  %12.7f %12.7f %12.7f\n", dftb->phase2.com[XX], dftb->phase2.com[YY], dftb->phase2.com[ZZ]);

  /* coordinates and magnitude of external charges for the complex */
  if (ct->qmmm > 0) {
//    printf("  phase2.xe = %p\n", dftb->phase2.xe);
//    printf("  image     = %p\n", image          );
//    printf("  bond      = %p\n", bond           );
    for (j=0; j<ct->extcharges_cplx; j++) {
      for (k=0; k<DIM; k++) {
        dftb->phase2.xe[j][k] = NM_TO_BOHR * x_ct[ct->extcharge_cplx[j]][k];
      }
      //if (j<10)
      //       printf("          %12.7f %12.7f %12.7f\n", dftb->phase2.xe[j][XX], dftb->phase2.xe[j][YY], dftb->phase2.xe[j][ZZ]);
      if (ct->qmmm < 3) {
      mindist = 1.e10;
      clear_ivec(shiftmin);
      for (shift[XX]=-1; shift[XX]<=1; shift[XX]++)
        for (shift[YY]=-1; shift[YY]<=1; shift[YY]++)
          for (shift[ZZ]=-1; shift[ZZ]<=1; shift[ZZ]++) {
            image[XX] = dftb->phase2.xe[j][XX] + shift[XX] * box[XX];
            image[YY] = dftb->phase2.xe[j][YY] + shift[YY] * box[YY];
            image[ZZ] = dftb->phase2.xe[j][ZZ] + shift[ZZ] * box[ZZ];
            dvec_sub(image, dftb->phase2.com, bond);
            // printf("          %12.7f %12.7f %12.7f\n", bond[XX], bond[YY], bond[ZZ]);
            curdist = dnorm(bond);
            if (curdist < mindist) {
              mindist = curdist;
              copy_ivec(shift, shiftmin);
            }
          }
      dftb->phase2.xe[j][XX] += shiftmin[XX] * box[XX];
      dftb->phase2.xe[j][YY] += shiftmin[YY] * box[YY];
      dftb->phase2.xe[j][ZZ] += shiftmin[ZZ] * box[ZZ];
      }
      /* magnitude, then add 0.06080 to O4a and C2q */
      dftb->phase2.ze[j] = mdatoms->chargeA[ct->extcharge_cplx[j]];
    }
    for (j=0; j<ct->sites*2; j++)
      if (ct->modif_extcharge_cplx[j] > -1)
        dftb->phase2.ze[ct->modif_extcharge_cplx[j]] += ct->shift_extcharge_cplx[j];
  }
  /* begin debug - check of sum of extcharges
  printf("Site   sum of extcharges\n");
  for (i=0; i<ct->sites; i++) {
    sum = 0.0;
    for (j=0; j<ct->extcharges[i]; j++)
      sum += dftb->phase1[i].ze[j];
    printf("%d %12.7f\n", i+1, sum);
  }
  sum = 0.0;
  for (j=0; j<ct->extcharges_cplx; j++)
    sum += dftb->phase2.ze[j];
  printf("Complex: sum of extcharges = %12.7f\n", sum);
     end debug */

  if (ct->qmmm == 3) {
    /* feed the cordinates into the PME arrays */
    for (i=0; i<ct->sites; i++) {
      dftb1 = dftb->phase1[i];
      for (j=0; j<dftb1.nn; j++) {
            dftb1.x_pme[j][0] = (real) dftb1.x[j][0] / NM_TO_BOHR;
            dftb1.x_pme[j][1] = (real) dftb1.x[j][1] / NM_TO_BOHR;
            dftb1.x_pme[j][2] = (real) dftb1.x[j][2] / NM_TO_BOHR;
            // dftb1.q_pme[j]    = 0.; // (real) (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
      }
      for (j=0; j<dftb1.ne; j++) {
            dftb1.x_pme[dftb1.nn + j][0] = (real) dftb1.xe[j][0] / NM_TO_BOHR;
            dftb1.x_pme[dftb1.nn + j][1] = (real) dftb1.xe[j][1] / NM_TO_BOHR;
            dftb1.x_pme[dftb1.nn + j][2] = (real) dftb1.xe[j][2] / NM_TO_BOHR;
            dftb1.q_pme[dftb1.nn + j]    = (real) dftb1.ze[j];
      }
    }
    dftb2 = dftb->phase2;
    for (j=0; j<dftb2.nn; j++) {
          dftb2.x_pme[j][0] = (real) dftb2.x[j][0] / NM_TO_BOHR;
          dftb2.x_pme[j][1] = (real) dftb2.x[j][1] / NM_TO_BOHR;
          dftb2.x_pme[j][2] = (real) dftb2.x[j][2] / NM_TO_BOHR;
          // dftb2.q_pme[j]    = 0.; // (real) (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]);
    }
    for (j=0; j<dftb2.ne; j++) {
          dftb2.x_pme[dftb2.nn + j][0] = (real) dftb2.xe[j][0] / NM_TO_BOHR;
          dftb2.x_pme[dftb2.nn + j][1] = (real) dftb2.xe[j][1] / NM_TO_BOHR;
          dftb2.x_pme[dftb2.nn + j][2] = (real) dftb2.xe[j][2] / NM_TO_BOHR;
          dftb2.q_pme[dftb2.nn + j]    = (real) dftb2.ze[j];
    }
  }

  /* calculate ESP (electro-static potentials) at the nucleobases */
  //printf("Site   ESP (V)\n");
  if (ct->qmmm == 1 || ct->qmmm == 2)
    for (i=0; i<ct->sites; i++) {
      sum = 0.0;
      for (j=0; j<ct->extcharges[i]; j++) {
        //dvec_sub(dftb->phase1[i].xe[j], dftb->phase1[i].com, r);
        dvec_sub(dftb->phase1[i].xe[j], dftb->phase1[i].x[0], r);
        sum += dftb->phase1[i].ze[j] / dnorm(r);
      }
      dftb->phase1[i].esp = sum / ct->esp_scaling_factor;
      //printf("%d %12.7f\n", i+1, sum * AU_OF_ESP_TO_VOLT);
    }

  return;
}

#ifdef GMX_MPI
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
  int i, j, k, nn, ne, status;
  dftb_phase1_t dftb1;
  dftb_phase2_t dftb2;
  real energy_pme, normbond;
  rvec bond, x_qm, x_mm;
  dvec dbond;
  t_pbc pbc;

  set_pbc(&pbc, epbcXYZ, dftb->box_pme);

  for (i=0; i<ct->sites; i++)
#ifdef GMX_MPI
  if (i % ct_mpi_size == ct_mpi_rank)
#endif
  {
    dftb1 = dftb->phase1[i];
    nn = dftb1.nn;
    ne = dftb1.ne;
    /* do PME for DFTB here!
     * do it with the atomic charges obtained in the previous MD step
     * it is possibly slightly inaccurate,
     * but probably OK as a first guess */

    for (j=0; j<nn; j++)
      dftb1.pot[j] = 0.0;

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase1[i].nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
        	  dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

    /* PME -- corrections */
    for (j=0; j<nn; j++) {
      /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
      dftb1.pot2[j] = 0.;
      for (k=0; k<nn; k++)
        if (j != k) {
          dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
          dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
	                  * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);  // this should be OK
        }
      /* the self-interaction of charge densities */
      dftb1.pot3[j] = - dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI); // this should be OK
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j=0; j<nn; j++) {
      dftb1.pot4[j] = 0.;
      x_qm[XX] = (real) dftb1.x[j][XX] / NM_TO_BOHR;
      x_qm[YY] = (real) dftb1.x[j][YY] / NM_TO_BOHR;
      x_qm[ZZ] = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
      /* loop over the neighbors */
      for (k=0; k<dftb1.neighbors_pme[j]; k++) {
        x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
        x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
        x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
        /* calculate the distance (considering PBC) */
        status = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
	normbond = norm(bond);
	if (normbond < dftb->rcoulomb_pme)
	  dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
      }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j=0; j<nn; j++) {
      /* convert kJ/mol/e to atomic units of ESP */
      dftb1.pot[j] *= KJMOL_TO_HARTREE;
      if (i==0 && j==0)
        printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
      dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
      //printf("%12.7f\n", dftb1.pot[j]);
    }

    /* save some an average value in the variable esp */
    dftb->phase1[i].esp = 0.;
    for (j=0; j<nn; j++)
      dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
    dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
    /* the following would be a very approximative alternative...
     * dftb->phase1[i].esp = dftb1.pot[0];
     */

    /* result is in: dftb.phase1[i].pot */
  }
#ifdef GMX_MPI
  if (ct_mpi_rank == 0) {
#endif
    /* do it here for the complex */
    /* POSSIBLE OPTIMIZATION -- USE THE DATA OBTAINED FOR THE FRAGMENTS?
     * BUT ATTENTION: NOT ALL OF THE DATA HAVE BEEN CALCULATED ON MPI_RANK==0 !!!
     */
    dftb2 = dftb->phase2;
    nn = dftb2.nn;
    ne = dftb2.ne;

    for (j=0; j<nn; j++)
      dftb2.pot[j] = 0.0;

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase2.nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
        	  dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

    /* PME -- corrections */
    for (j=0; j<nn; j++) {
      /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
      dftb2.pot2[j] = 0.;
      for (k=0; k<nn; k++)
        if (j != k) {
          dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
          dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
	                  * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
        }
      /* the self-interaction of charge densities */
      dftb2.pot3[j] = - dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j=0; j<nn; j++) {
      dftb2.pot4[j] = 0.;
      x_qm[XX] = (real) dftb2.x[j][XX] / NM_TO_BOHR;
      x_qm[YY] = (real) dftb2.x[j][YY] / NM_TO_BOHR;
      x_qm[ZZ] = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
      /* loop over the neighbors */
      for (k=0; k<dftb2.neighbors_pme[j]; k++) {
        x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
        x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
        x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
        /* calculate the distance (considering PBC) */
        status = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
	normbond = norm(bond);
	if (normbond < dftb->rcoulomb_pme)
	  dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
      }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j=0; j<nn; j++) {
      /* convert kJ/mol/e to atomic units of ESP */
      dftb2.pot[j] *= KJMOL_TO_HARTREE;
      //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
      dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
      //printf("%12.7f\n", dftb1.pot[j]);
    }
  /* result is in: dftb.phase2.pot */
#ifdef GMX_MPI
  }
#endif
  return;
}

#ifdef GMX_MPI
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
  int i, j, k, nn, ne, status;
  dftb_phase1_t dftb1;
  dftb_phase2_t dftb2;
  real energy_pme, normbond;
  rvec bond, x_qm, x_mm;
  dvec dbond;
  t_pbc pbc;

  set_pbc(&pbc, epbcXYZ, dftb->box_pme);

  for (i=0; i<ct->sites; i++)
#ifdef GMX_MPI
  if (i % ct_mpi_size == ct_mpi_rank)
#endif
  {
    dftb1 = dftb->phase1[i];
    nn = dftb1.nn;
    ne = dftb1.ne;
    /* do PME for DFTB here! -- "preparatory" part, before actual SCC-DFTB calculation */

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j=0; j<nn; j++) {
      dftb1.pot4[j] = 0.;
      x_qm[XX] = (real) dftb1.x[j][XX] / NM_TO_BOHR;
      x_qm[YY] = (real) dftb1.x[j][YY] / NM_TO_BOHR;
      x_qm[ZZ] = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
      /* loop over the neighbors */
      for (k=0; k<dftb1.neighbors_pme[j]; k++) {
        x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
        x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
        x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
        /* calculate the distance (considering PBC) */
        status = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
	normbond = norm(bond);
	if (normbond < dftb->rcoulomb_pme)
	  dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
      }
    }
  }
  return;
}

void do_pme_for_dftb_part2(charge_transfer_t *ct, dftb_t *dftb, int i)
{
  int j, k, nn, ne, status;
  dftb_phase1_t dftb1;
  real energy_pme, normbond;
  rvec bond, x_qm, x_mm;
  dvec dbond;
  t_pbc pbc;

  set_pbc(&pbc, epbcXYZ, dftb->box_pme);

  dftb1 = dftb->phase1[i];
  nn = dftb1.nn;
  ne = dftb1.ne;

  /* transfer the QM charges to PME,
   * and reset the ESP array
   */
  for (j=0; j<nn; j++) {
    dftb1.pot[j] = 0.0;
    dftb1.q_pme[j] = -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]];
  }

  /* PME calculation */
  init_nrnb(dftb->phase1[i].nrnb_pme);
  gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
      	  dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

  /* PME -- corrections */
  for (j=0; j<nn; j++) {
    /* exclude the QM-QM interactions */
    dftb1.pot2[j] = 0.;
    for (k=0; k<nn; k++)
      if (j != k) {
        dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
        dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
                        * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
      }
    /* exclude the self-interaction of the charge densities */
    dftb1.pot3[j] = - dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI);
  }

  /* dftb1.pot4 was pre-calculated */

  //printf("Ewald potential kJ/mol/e = external shift\n");
  for (j=0; j<nn; j++) {
    /* convert kJ/mol/e to atomic units of ESP */
    dftb1.pot[j] *= KJMOL_TO_HARTREE;
    //if (i==0 && j==0)
    //  printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
    dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
    //if (i==0 && j==0) printf("%12.7f ", dftb1.pot[j]);
  }

  /* save some an average value in the variable esp */
  dftb->phase1[i].esp = 0.;
  for (j=0; j<nn; j++)
    dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
  dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
  /* this was too approximative:
    dftb->phase1[i].esp = dftb1.pot[0];
  */

  return;
}

/* PME for DFTB phase 2 -- this will be done only once per MD step */
void do_pme_for_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb)
{
  int i, j, k, nn, ne, status;
  dftb_phase2_t dftb2;
  real energy_pme, normbond;
  rvec bond, x_qm, x_mm;
  dvec dbond;
  t_pbc pbc;

  set_pbc(&pbc, epbcXYZ, dftb->box_pme);

  dftb2 = dftb->phase2;
  nn = dftb2.nn;
  ne = dftb2.ne;

  for (j=0; j<nn; j++) {
    dftb2.q_pme[j] = - dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]];
    dftb2.pot[j] = 0.0;
  }

  /* PME -- long-range component (reciprocal space) */
  init_nrnb(dftb->phase2.nrnb_pme);
  gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
      	  dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

  /* PME -- corrections */
  for (j=0; j<nn; j++) {
    /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
    dftb2.pot2[j] = 0.;
    for (k=0; k<nn; k++)
      if (j != k) {
        dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
        dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
                        * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
      }
    /* the self-interaction of charge densities */
    dftb2.pot3[j] = - dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
  }

  /* PME -- short-range component (real space) using the previously created neighbor list */
  for (j=0; j<nn; j++) {
    dftb2.pot4[j] = 0.;
    x_qm[XX] = (real) dftb2.x[j][XX] / NM_TO_BOHR;
    x_qm[YY] = (real) dftb2.x[j][YY] / NM_TO_BOHR;
    x_qm[ZZ] = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
    /* loop over the neighbors */
    for (k=0; k<dftb2.neighbors_pme[j]; k++) {
      x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
      x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
      x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
      /* calculate the distance (considering PBC) */
      status = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
      normbond = norm(bond);
      if (normbond < dftb->rcoulomb_pme)
        dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
    }
  }

  //printf("Ewald potential kJ/mol/e = external shift\n");
  for (j=0; j<nn; j++) {
    /* convert kJ/mol/e to atomic units of ESP */
    dftb2.pot[j] *= KJMOL_TO_HARTREE;
    //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
    dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
    //printf("%12.7f\n", dftb1.pot[j]);
  }
  /* result is in: dftb.phase2.pot */
  return;
}

#ifdef GMX_MPI
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x, int ct_mpi_rank, int ct_mpi_size)
#else
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x)
#endif
{
  int i, j, k, nn, ne, counter, status;
  dftb_phase1_t dftb1;
  dftb_phase2_t dftb2;
  rvec bond;
  t_pbc pbc;
  real rcoulomb2;

  set_pbc(&pbc, epbcXYZ, dftb->box_pme);
  rcoulomb2 = dftb->rcoulomb_pme * dftb->rcoulomb_pme;

  //printf("  creating PME DFTB neighborlists...");
  for (i=0; i<ct->sites; i++)
#ifdef GMX_MPI
  if (i % ct_mpi_size == ct_mpi_rank)
#endif
  {
    dftb1 = dftb->phase1[i];
    nn = dftb1.nn;
    ne = dftb1.ne;
    /* create a neighborlist for PME calculation for QM/MM DFTB
     * in a naive O(N^2) way:
     * double loop over QM atoms and over MM atoms */

    /* do it for every QM atom */
    for (j=0; j<nn; j++) {
      counter = 0;
      /* check every MM atom */
      for (k=0; k<ne; k++) {
        status = pbc_dx_aiuc(&pbc, x[ct->atom[i][j]], x[ct->extcharge[i][k]], bond);
	if (norm2(bond) < rcoulomb2) {
	  dftb1.neighbor_pme[j][counter] = k;
	  counter++;
	  // printf("Fragment %d, atom %d: neighbor no. %d is extcharge %d\n", i+1, j+1, counter, k+1);
	  if (counter == MAX_PME_NEIGHBORS) {
	    fprintf(stderr, "Too many PME neighbors for atom %d of fragment %d\n  Exiting !!!\n\n", j+1, i+1);
	    exit(-1);
	  }
	}
      }
      dftb1.neighbors_pme[j] = counter;
      //fprintf(stderr, "NS for PME/DFTB: atom %d in fragment %d has %d MM neighbors\n", j+1, i+1, counter);
    }
  }

#ifdef GMX_MPI
  if (ct_mpi_rank == 0) {
#endif
  /* do it here for the complex */
  dftb2 = dftb->phase2;
  nn = dftb2.nn;
  ne = dftb2.ne;
  /* do it for every QM atom */
  for (j=0; j<nn; j++) {
    counter = 0;
    /* check every MM atom */
    for (k=0; k<ne; k++) {
      status = pbc_dx_aiuc(&pbc, x[ct->atom_cplx[j]], x[ct->extcharge_cplx[k]], bond);
      if (norm2(bond) < rcoulomb2) {
        dftb2.neighbor_pme[j][counter] = k;
        counter++;
        if (counter == MAX_PME_NEIGHBORS) {
          fprintf(stderr, "Too many PME neighbors for atom %d of complex\n  Exiting !!!\n\n", j+1);
          exit(-1);
        }
      }
    }
    dftb2.neighbors_pme[j] = counter;
    //fprintf(stderr, "NS for PME/DFTB: atom %d in complex has %d MM neighbors\n", j+1, counter);
  }
#ifdef GMX_MPI
  }
#endif
  return;
}

#ifdef GMX_MPI
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
  int i, return_value;

  // printf("do_dftb_phase1 start at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);
/*
  if (ct_mpi_rank > 0) {
    ct_loc = (charge_transfer_t *) malloc(sizeof (charge_transfer_t));
    ct_loc->sites = ct->sites;
  }
*/
  for (i=0; i<ct->sites; i++)
    if (i % ct_mpi_size == ct_mpi_rank) {
      // printf("Doing nucleobase %d at rank %d\n", i, ct_mpi_rank);
      //if (ct->qmmm == 3)
      //  run_dftb1_doublecycle(ct, dftb, i);
      //else
      run_dftb1(ct, dftb, i);
      if (ct_mpi_rank != 0) {
        //printf("Transferring data on nucleobase %d at rank %d to master\n", i, ct_mpi_rank);
        return_value = MPI_Send(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb), MPI_DOUBLE, 0, 100 + i * 4, ct_mpi_comm);
        //printf("Transfer of array a, nucleobase %d at rank %d to master finished\n", i, ct_mpi_rank);
        return_value = MPI_Send(dftb->phase1[i].qmat, dftb->phase1[i].nn, MPI_DOUBLE, 0, 101 + i * 4, ct_mpi_comm);
        //printf("Transfer of array qmat, nucleobase %d at rank %d to master finished\n", i, ct_mpi_rank);
        return_value = MPI_Send(dftb->phase1[i].ev, dftb->phase1[i].norb, MPI_DOUBLE, 0, 102 + i * 4, ct_mpi_comm);
        //printf("Transfer of array ev, nucleobase %d at rank %d to master finished\n", i, ct_mpi_rank);
        if (ct->qmmm == 3) {
          return_value = MPI_Send(&(dftb->phase1[i].esp), 1, MPI_DOUBLE, 0, 103 + i * 4, ct_mpi_comm);
          //printf("Transfer of value esp, nucleobase %d at rank %d to master finished\n", i, ct_mpi_rank);
        }
      }
    }

  // printf("do_dftb_phase1 end at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);

  return;
}

void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
  int i, j, return_value;

  if (ct->qmmm == 3) {
    for (i=0; i<ct->sites; i++) {
      if (i % ct_mpi_size == ct_mpi_rank) {
        /* set the charges on "QM" atoms */
        for (j=0; j<dftb->phase1[i].nn; j++) {
          dftb->phase1[i].qmat[j] = - q[ct->atom[i][j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        dftb->phase1[i].qmat[0] += - q[ct->atom[i][0]+1];
        if (ct->modif_extcharge[i][0] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        if (ct->modif_extcharge[i][1] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        /* calculate the remaining components of PME */
        do_pme_for_dftb_part2(ct, dftb, i);
        if (ct_mpi_rank != 0) {
          return_value = MPI_Send(&(dftb->phase1[i].esp), 1, MPI_DOUBLE, 0, 103 + i * 4, ct_mpi_comm);
        }
      }
    }
  }
  return;
}
#else
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb)
{
  int i;

  // printf("do_dftb_phase1 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

  for (i=0; i<ct->sites; i++) {
    //if (ct->qmmm == 3)
    //  run_dftb1_doublecycle(ct, dftb, i);
    //else
    run_dftb1(ct, dftb, i);
  }
  //printf("%12.7f %12.7f %12.7f %12.7f %12.7f", dftb->phase1[0].pot[0], dftb->phase1[0].pot2[0], dftb->phase1[0].pot3[0], dftb->phase1[0].pot4[0] / NM_TO_BOHR,
  //                                             - dftb->phase1[0].qmat[0] + dftb->qzero1[dftb->phase1[0].izp[0]]);

  // printf("do_dftb_phase1 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);

  return;
}
void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q)
{
  int i, j;

  if (ct->qmmm == 3) {
      for (i=0; i<ct->sites; i++) {
        /* set the charges on "QM" atoms */
        for (j=0; j<dftb->phase1[i].nn; j++) {
          dftb->phase1[i].qmat[j] = - q[ct->atom[i][j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        dftb->phase1[i].qmat[0] += - q[ct->atom[i][0]+1];
        if (ct->modif_extcharge[i][0] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        if (ct->modif_extcharge[i][1] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        /* calculate the remaining components of PME */
        do_pme_for_dftb_part2(ct, dftb, i);
      }
  }

  return;
}
#endif

void do_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb)
{

  // printf("do_dftb_phase2 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

  run_dftb2(ct, dftb);

  // printf("do_dftb_phase2 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);

  return;
}

void ct_assemble_hamiltonian(charge_transfer_t *ct, dftb_t *dftb)
{
  int i, j;
  double value;

  for (i=0; i<ct->sites; i++) {
    ct->hamiltonian[i][i] = - dftb->phase1[i].ev[ct->homo[i]-1];
    for (j=i+1; j<ct->sites; j++) {
      value = dftb->phase2.THamilOrtho[i][j];
      ct->hamiltonian[i][j]
        = ct->hamiltonian[j][i]
        = value > 0 ? value : -value;
    }
  }

  printf("TB Hamiltonian:\n");
  for (i=0; i<ct->sites; i++) {
    for (j=i; j<ct->sites; j++) printf("%10.6f", ct->hamiltonian[i][j] * HARTREE_TO_EV);
  }
  printf("\n");

  return;
}

/*
void ct_integrate_tdse(charge_transfer_t *ct, double twant)
{
  int i, j;
  double value;

  return;
}
*/

// lapack routines
static long dsyev(int int_n, double *a, double *w, double *work, int int_lwork)
{
  extern void dsyev_(char *, char *, long *, double *, long *, double *, double *, long *, long *);
  long info, n, lwork;
  char jobz='V', uplo='U';
  n = (long) int_n;
  lwork = (long) int_lwork;

  dsyev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, &info);

  return info;
}

static long dgesv(int int_n, double *a, double *b, long ipiv[DIIS_MAX_PREV_VECTORS + 1])
{
  extern void dgesv_(long *, long *, double *, long *, long *, double *, long *, long *);
  long n, nrhs, lda, ldb, info;
  n = lda = ldb = (long) int_n;
  nrhs = 1L;

  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  return info;
}

void broyden(int niter, double alpha, int jtop, double *vecin, double *vecout, dftb_broyden_t *arrays);

double fermi_coef_sum(double x, int n, double *a, double *coeffs)
{
  int i;
  double sum;

  sum = 0.0;
  for (i=0; i<n; i++)
    sum += coeffs[i] = 1.0 / ( exp( (a[i]-x) / FERMI_KT ) + 1.0 );

  return sum;
}

int do_adiabatic(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
//int do_adiabatic(charge_transfer_t *ct, dftb_broyden_t *broyd, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, *ham, almix, fermi_energy, fermi_upper, fermi_lower;

  ham = ct->hamiltonian_adiab;
  almix = SIMPLE_ALMIX;

  diis->n_prev_vector = 1;
  diis->indx = 0;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations?
    if (step >= MAXITER_DIIS) {
      // copy the solution to the wave function
      for (i=0; i<ct->sites; i++)
        ct->wf[i] = sqrt(ct->q_act[i]);
      return 1;
    }

    printf("step %5d:", step);
    //fprintf(f, "step %5d:", step);

    // save and print old charges
    printf(" old q =");
    //fprintf(f, " old q =");
    for (i=0; i<ct->sites; i++)
      //fprintf(f, "%12.9f", ct->q_act[i]);
      printf("%12.9f", ct->q_act[i]);

    //  FOR BROYDEN:
    //  ct->q_old[i] = ct->q_act[i];

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // copy previous charge and the difference to the arrays
    for (i=0; i<ct->sites; i++) {
      diis->prev_q_input[diis->indx][i] = ct->q_act[i];
      diis->prev_q_diff[diis->indx][i] = - ct->q_act[i]; // here, the new occupation will be added
    }
/*
    // instead of Fermi:
    for (i=0; i<ct->sites; i++)
      diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);
*/
    // FERMI DISTRIBUTION (electronic temperature)
    // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, diis->fermi_coef) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, " %8.5f", diis->fermi_coef[i]);
    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++)   // i runs over the sites
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->sites + i]);

    // calculate the energy and check convergence
    old_energy = energy;
    energy = 0.e0;
    /* old - calculate energy by orbitals
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        energy += ham[i] * ham[j] * ct->hamiltonian[i][j]
	        + 0.5 * ct->q_act[i] * ct->q_act[j] * (ct->sic ? (SIC_COEFFICIENT * ct->hubbard[i][j]) : (ct->hubbard[i][j]));
    */
    for (i=0; i<ct->sites; i++)
      energy += diis->fermi_coef[i] * ct->ev_adiab[i];
    printf(" energy = %12.9f\n", energy);
    //fprintf(f, " energy = %12.9f\n", energy);
    if (fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //fprintf(f, "Converged!\n");
      printf("Converged!\n");
      // copy the solution to the wave function
      for (i=0; i<ct->sites; i++)
        //ct->wf[i] = sqrt(ct->q_act[i]);
        ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
        //ct->wf[i] = ham[i];
      break;
    }

    // DIIS - Direct Inversion in the Iterative Subspace
    // we need in memory:
    // p        -- matrix (# of bases) x (# of vectors)
    // delta p  -- matrix (# of bases) x (# of vectors)
    // b        -- matrix (# of vectors + 1) x (# of vectors + 1)
    
    // produce a new vector "ct->q_old", either by simple mixing (in first iterations) or by DIIS (later)
    if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS) {
      // SIMPLE MIXING in a couple of first steps
      for (i=0; i<ct->sites; i++)
        ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
        //ct->q_act[i] = (1.0 - DIIS_INIT_MIXING_PARAM) * ct->q_act[i] + DIIS_INIT_MIXING_PARAM * diis->q_inp_result[i];
      diis->n_prev_vector++;
    } else {
      // real DIIS here
      // what do the vectors "prev_q_diff" look like?
      //printf("\n DIFF VECTORS\n");
      //for (i=0; i<DIIS_MAX_PREV_VECTORS; i++) {
      //  for (j=0; j<ct->sites; j++)
      //  printf("%14.10f", diis->prev_q_diff[i][j]);
      //printf("\n");
      //}
      // construct the matrix of overlaps
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        for (j=0; j<DIIS_MAX_PREV_VECTORS; j++) {
	  diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
	  for (k=0; k<ct->sites; k++) {
	    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
	  }
	  //printf("aa[%d][%d] = %20.15f\n", i, j, diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
	}
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
      diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
      //printf("\n MATRIX AA:\n");
      //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++) {
      //  for (j=0; j<=DIIS_MAX_PREV_VECTORS; j++)
      //  printf("%14.10f", diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
      //printf("\n");
      //}
      // construct the vector of right hand sides
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->bb[i] = 0.e0;
      diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;
      //printf(" MATRIX BB:\n");
      //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
      //  printf("%14.10f", diis->bb[i]);
      //printf("\n");
      // solve it with lapack
      //printf(" DGESV: %ld", dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv));
      dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);
      // the column of solution is now in diis->bb -- print it!
      //printf(" DIIS coeffs");
      //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
      //  printf("%9.5f", diis->bb[i]);
      for (j=0; j<ct->sites; j++) {
        ct->q_act[j] = 0.e0;
        for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
	  ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
      }
    }

    // move the counter in a loop
    diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
        
/*
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in broyd->df
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, broyd->df) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    //fermi_coef2 = 1 / (exp((ct->ev_adiab[1] - ct->ev_adiab[0]) / FERMI_KT) + 1);
    //printf(" fermi coef2 = %7.4f", fermi_coef2);
    fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, " %8.5f", broyd->df[i]);

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += broyd->df[j] * SQR(ham[j * ct->sites + i]);
    }

    // mix it - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, *broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];
*/

/*
    // SIMPLE MIXING
    //
    // calculate the charges
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = SQR(ham[i]);
    // calculate oldness - the content of old vector in the new one
    oldness = 0.e0;
    for (i=0; i<ct->sites; i++)
      oldness += ct->q_act[i] * ct->q_old[i];
    // the orthogonal vector is now: new_ortho = new - oldness * old
    // calculate the new vector as: old + almix(t) * new_ortho
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = (1.0 - almix) * ct->q_old[i] + almix * (ct->q_act[i] - oldness * ct->q_old[i]);
    // attenuate the mixing coefficient
    almix *= ALMIX_ATTENUATOR;
*/

  }

  return 0;
}

int do_wf_decomp_evec(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, *ham, dotproduct_squared;

  ham = ct->hamiltonian_adiab;
  diis->n_prev_vector = 1;
  diis->indx = 0;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = SQR(ct->wf[i]);

  // SCC cycle
  for (step = 0, energy = 1.e10; ; step++) {

    if (step >= MAXITER_DIIS) {
      /* DO HERE SOMETHING, BUT DO NOT MODIFY ct->wf !
      for (i=0; i<ct->sites; i++)
        ct->wf[i] = sqrt(ct->q_act[i]);
      */
      return 1;
    }

    //printf("step %5d:", step);

    //printf(" old q =");
    //for (i=0; i<ct->sites; i++)
    //  printf("%9.5f", ct->q_act[i]);

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // copy previous charge and the difference to the arrays
    for (i=0; i<ct->sites; i++) {
      diis->prev_q_input[diis->indx][i] = ct->q_act[i];
      diis->prev_q_diff[diis->indx][i] = - ct->q_act[i]; // here, the new occupation will be added
    }
    // FERMI DISTRIBUTION (electronic temperature)
    // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
    /* SWITCH IT OFF FOR DECOMPOSITION !
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, diis->fermi_coef) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, " %8.5f", diis->fermi_coef[i]);
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->sites + i]);
    */
    for (i=0; i<ct->sites; i++)
      diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);

    old_energy = energy;
    /* AGAIN - DO NOT MIX THE ENERGY !
    energy = 0.e0;
    for (i=0; i<ct->sites; i++)
      energy += diis->fermi_coef[i] * ct->ev_adiab[i];
    */
    energy = ct->ev_adiab[0];
    //printf(" energy = %12.9f\n", energy);
    if (fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //printf("Converged!\n");
      /* NOW, WE HAVE THE ADIABATIC STATES IN ham !
      for (i=0; i<ct->sites; i++)
        ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
      */
      break;
    }

    if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS) {
      // SIMPLE MIXING in a couple of first steps
      for (i=0; i<ct->sites; i++)
        ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
      diis->n_prev_vector++;
    } else {
      // REAL DIIS here
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        for (j=0; j<DIIS_MAX_PREV_VECTORS; j++) {
	  diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
	  for (k=0; k<ct->sites; k++) {
	    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
	  }
	}
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
      diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;

      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->bb[i] = 0.e0;
      diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;

      dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);

      for (j=0; j<ct->sites; j++) {
        ct->q_act[j] = 0.e0;
        for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
	  ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
      }
    }
    diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
  }

  // print energies of adiabatic states
  //printf("Adiab E: ");
  for (i=0; i<ct->sites; i++)
    //fprintf(f, "%10.7f", ct->ev_adiab[i]);
    fprintf(f, "%10.7f", ct->ev_adiab[i] - ct->ev_adiab[0]);

  // calculate the expansion coefficients
  fprintf(f, "   WF exp'd in adiab states: ");
  /*
  for (i=0; i<ct->sites; i++) {
    // real part of a_i
    dotproduct_re = 0.e0;
    for (j=0; j<ct->sites; j++) dotproduct_re += ham[i * ct->sites + j] * ct->wf[j];
    // imaginary part of a_i
    dotproduct_im = 0.e0;
    for (j=0; j<ct->sites; j++) dotproduct_im += ham[i * ct->sites + j] * ct->wf[j + ct->sites];
    fprintf(f, "%9.6f", SQR(dotproduct_re)+SQR(dotproduct_im));
  }
  */
  for (i=0; i<ct->sites; i++) {
    // calculate directly the square of a_i
    dotproduct_squared = 0.e0;
    for (j=0; j<ct->sites; j++)
      for (k=0; k<ct->sites; k++)
        dotproduct_squared += ham[i * ct->sites + j] * ham[i * ct->sites + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->sites] * ct->wf[k + ct->sites]);
    fprintf(f, "%9.6f", dotproduct_squared);
  }
  fprintf(f, "\n");

  return 0;
}

int do_surface_hopping(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
  /* this should be done by the modification of do_adiabatic
   * followed by the testing of scalar product of the new wavefunction - ground and 1st excited state - with the old one
   */
  /* ct->surface = 0 (ground state) or 1 (1st excited state) */
  int i, j, k, step;
  double old_energy, energy, *ham, s01, s10, nonadiab_coupling;

  ham = ct->hamiltonian_adiab;
  diis->n_prev_vector = 1;
  diis->indx = 0;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = SQR(ct->wf[i]);

  // SCC cycle
  for (step = 0, energy = 1.e10; ; step++) {

    if (step >= MAXITER_DIIS) {
      /* DIIS NOT CONVERGED! */
      //printf("DIIS not converged !!!      ");
      //fprintf(f, "DIIS not converged !!!      ");
      fprintf(f, "DIIS   0 ");
      fprintf(f, "             ");
      for (i=0; i<ct->sites; i++)
        fprintf(f, "        ");
      fprintf(f, "                           DIIS not converged so    hopping not tested... ");
      fprintf(f, " %d", ct->surface);
      for (i=0; i<ct->sites; i++)
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
      fprintf(f, "\n");

      return 1;
    }

    //printf("step %5d:", step);

    //printf(" old q =");
    //for (i=0; i<ct->sites; i++)
      //printf("%9.5f", ct->q_act[i]);

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    /*
    // orbital shift:
    // 1. preliminary diagonalization
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    // 2. correction to the elements of the Hamiltonian
    //    Delta H = C * delta-Epsilon * C^T
    for (i=0; i<ct->sites; i++)
      ct->work_adiab[i] = ham[i];
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[j];
      }
    }
    // done with orbital shift
    */

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // copy previous charge and the difference to the arrays
    for (i=0; i<ct->sites; i++) {
      diis->prev_q_input[diis->indx][i] = ct->q_act[i];
      diis->prev_q_diff[diis->indx][i] = SQR(ham[i]) - ct->q_act[i];
    }

    old_energy = energy;
    energy = ct->ev_adiab[0];
    //printf(" energy = %12.9f\n", energy);
    if (fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //printf("DIIS Converged in %3d steps ", step);
      //fprintf(f, "DIIS Converged in %3d steps ", step);
      fprintf(f, "DIIS %3d ", step);
      break;
    }

    if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS) {
      // SIMPLE MIXING in a couple of first steps
      for (i=0; i<ct->sites; i++)
        ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
      diis->n_prev_vector++;
    } else {
      // REAL DIIS here
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        for (j=0; j<DIIS_MAX_PREV_VECTORS; j++) {
	  diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
	  for (k=0; k<ct->sites; k++) {
	    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
	  }
	}
      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
      diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;

      for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
        diis->bb[i] = 0.e0;
      diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;

      dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);

      for (j=0; j<ct->sites; j++) {
        ct->q_act[j] = 0.e0;
        for (i=0; i<DIIS_MAX_PREV_VECTORS; i++)
	  ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
      }
    }
    diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
  }

  // print energies of adiabatic states
  //printf(" Eigenvalues:");
  fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
    //fprintf(f, "%10.7f", ct->ev_adiab[i]);
    fprintf(f, "%8.5f", ct->ev_adiab[i]);
    //printf("%8.5f", ct->ev_adiab[i] - ct->ev_adiab[0]);
  }
  //printf("    ");
  fprintf(f, "    ");

  // check the energy gap between the ground and the first excited state
  if (ct->ev_adiab[1] - ct->ev_adiab[0] < SFH_LIMIT) {
    /* attempt surface hopping
     * a variation of the AS3 method in Fabiano et al., Chem. Phys. 351 (2008) 111-116
     */
    s01 = 0.e0;
    for (i=0; i<ct->sites; i++)
      s01 += (ct->surface == 0 ? ct->wf[i] : ct->wf_exc[i]) * ham[ct->sites + i];
    s10 = 0.e0;
    for (i=0; i<ct->sites; i++)
      s10 += (ct->surface == 0 ? ct->wf_exc[i] : ct->wf[i]) * ham[i];
    /* calculate the (quasi-)non-adiabatic coupling */
    nonadiab_coupling = 1 - pow(1 - (fabs(s01) + fabs(s10))*0.5, SFH_EXPONENT);
    /* print them out */
    //printf("wf.ground = %6.3f, wf.excited = %6.3f ", overlap_ground, overlap_excited);
    fprintf(f, "s01 = %7.4f, s10 = %7.4f, coupling = %6.3f ", s01, s10, nonadiab_coupling);

    /* check if hopping should be performed */
    if (nonadiab_coupling > rand()/(double)RAND_MAX) {
      fprintf(f, "hopping between states");
      if (ct->surface == 0) { // hopping 0->1
        for (i=0; i<ct->sites; i++) {
          ct->wf[i] = ham[ct->sites + i];
          ct->wf_exc[i] = ham[i];
        }
        ct->surface = 1;
      } else {                // hopping 1->0
        for (i=0; i<ct->sites; i++) {
          ct->wf[i] = ham[i];
          ct->wf_exc[i] = ham[ct->sites + i];
        }
        ct->surface = 0;
      }
    } else {
      fprintf(f, "no hopping now...     ");
      for (i=0; i<ct->sites; i++) {
        if (ct->surface == 0) { // system is remaining in state 0
          ct->wf[i] = ham[i];
          ct->wf_exc[i] = ham[ct->sites + i];
        } else {                // system is remaining in state 1
          ct->wf[i] = ham[ct->sites + i];
          ct->wf_exc[i] = ham[i];
        }
      }
    }
  } else {
    //printf("Egap > threshold, no surf hopping ");
    fprintf(f, "                                                hopping not tested... ");
    for (i=0; i<ct->sites; i++) {
      if (ct->surface == 0) { // system is remaining in state 0
        ct->wf[i] = ham[i];
        ct->wf_exc[i] = ham[ct->sites + i];
      } else {                // system is remaining in state 1
        ct->wf[i] = ham[ct->sites + i];
        ct->wf_exc[i] = ham[i];
      }
    }
  }
  fprintf(f, " %d", ct->surface);
  for (i=0; i<ct->sites; i++) {
    //printf("%8.4f", SQR(ct->wf[i]));
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }
  //printf("\n");
  fprintf(f, "\n");

  return 0;
}

int do_adiab_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;

  ham = ct->hamiltonian_adiab;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!\n");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    fprintf(f, " fermi E = %9.6f, coefs", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, "%8.5f", fermi_coeff[i]);

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // print new charges
    fprintf(f, " mixedQ =");
    for (i=0; i<ct->sites; i++)
      fprintf(f, "%9.6f", ct->q_act[i]);

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    fprintf(f, " E1=%10.7f E2=%10.7f totE=%13.10f", energy1, energy2, energy);

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      fprintf(f, " Converged!\n");
      // copy the solution to the wave function
      for (i=0; i<ct->sites; i++)
        ct->wf[i] = sqrt(ct->q_act[i]);
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

    // print mixed charges
    fprintf(f, " newQ =");
    for (i=0; i<ct->sites; i++)
      fprintf(f, "%9.6f", ct->q_act[i]);
    fprintf(f, "\n");

  } /* end SCC cycle */

  return 0;
}

int do_adiab_fermi_onestate(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;

  ham = ct->hamiltonian_adiab;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  fprintf(f, " %2d ", step);
  //fprintf(f, " %2d iters, SCC-Q =", step);
  for (i=0; i<ct->sites; i++)
    fprintf(f, "%9.6f", ct->q_act[i]);

  // copy the lowest pure adiabatic state to the wave function
  for (i=0; i<ct->sites; i++)
    ct->wf[i] = ham[i];

  // calculate and print the energies of pure adiabatic states (indexed by k)
  //fprintf(f, " E(adiab) =");
  for (k=0; k<ct->sites; k++) {
    energy = 0.0;
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        energy += ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j]
        // Hubbard / gamma terms
               +  0.5 * ct->hubbard[i][j] * SQR(ham[k * ct->sites + i]) * SQR(ham[k * ct->sites + j]);
      }
    fprintf(f, "%10.7f", energy);
  }
  return 0;
}

int do_fermi_surface_hopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
  double random_number; // in the interval (0,1)
  double cumulated_probability; // sum of ct->surf_prob[0..k]

  ham = ct->hamiltonian_adiab;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //fprintf(f, "Broyd %2d ", step);
      fprintf(f, " %2d ", step);
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
  //fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.5f", ct->ev_adiab[i]);
  }
  //fprintf(f, "    ");

  /* surface hopping tried according to Groenhof et al.
   * Adv. Quant. Chem. 59, 181-212 (2010).
   */
  for (k=0; k<ct->sites; k++) { // k spans all adiabatic states, which are candidates for surface hopping
    ct->surf_overlap[k] = 0.;
    for (j=0; j<ct->sites; j++)
      ct->surf_overlap[k] += ct->wf[j] * ham[k * ct->sites + j];
    fprintf(f, "%10.6f", ct->surf_overlap[k]);
  }
  // calculate the Massey parameter for all new states
  for (k=0; k<ct->sites; k++) {
    ct->surf_massey[k] = k == ct->surface ? 0. : fabs((ct->ev_adiab[k] - ct->ev_adiab[ct->surface]) * ct->rk_timestep / ct->surf_overlap[k]);
    // fprintf(f, "%10.2f", ct->surf_massey[k]);
  }
  // calculate the surface hopping probability for all new states!
  for (k=0; k<ct->sites; k++) {
    ct->surf_prob[k] = k == ct->surface ? 1. : exp(-M_PI / 4 * ct->surf_massey[k]);
    fprintf(f, "%9.6f", ct->surf_prob[k]);
  }
  // draw a random number and check if the *cumulative* probability is larger
  random_number = rand()/(double)RAND_MAX;
  cumulated_probability = 0.;
  for (k=0; k<ct->sites; k++) {
    if (k == ct->surface)
      continue;
    cumulated_probability += ct->surf_prob[k];
    if (cumulated_probability > random_number) {
      ct->surface = k;
      break;
    }
  }
  // assign the new wave function to the eigenvector ct->surface
  for (i=0; i<ct->sites; i++)
    ct->wf[i] = ham[ct->surface * ct->sites + i];
  
  // print out current state
  fprintf(f, " %d", ct->surface);
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }
  fprintf(f, "\n");
  return 0;
}

int do_tully_fewest_switches(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2)
{
  int i, j, k, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
  double random_number; // in the interval (0,1)
  double cumulated_probability; // sum of ct->surf_prob[0..k]
  double popul_norm; // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
  double decay_time, current_surface_factor; // in the decoherence algorithm
  double dot_product;

  ham = ct->hamiltonian_adiab;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //fprintf(f, "Broyd %2d ", step);
      fprintf(f, " %2d ", step);
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
  //fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.5f", ct->ev_adiab[i]);
  }
  //fprintf(f, "    ");

  /* surface hopping tried according to Tully JCP 93, 1061 (1990)
   * "fewest switches" algorithm
   * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
   */

  if (ct->tfs_initialization_step) {
    /* let us start with tfs_vector_old == tfs_vector */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ham[i * ct->sites + j];
    ct->tfs_initialization_step = 0;
  } else {
    /* push tfs_vector to tfs_vector_old */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
  }

  /* enter new eigenvectors to tfs_vector */
  for (i=0; i<ct->sites; i++) {
    /* calculate the dot product with tfs_vector_old[i] */
    dot_product = 0.;
    for (k=0; k<ct->sites; k++)
      dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->sites + k];
    /* for positive dot_product, take the vector as it is
     * for negative dot_product, take (-1) * the vector
     */
    for (k=0; k<ct->sites; k++)
      ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->sites + k] : - ham[i * ct->sites + k];
  }

  /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++) {
      ct->tfs_overlap[i][j] = 0.;
      for (k=0; k<ct->sites; k++)
        ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
    }

  /* print out the overlap matrix */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++)
      fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
  fprintf(f2, "\n");

  /* integrate tfs_popul for one step
   * (the procedure will use tfs_popul_der)
   */
  do_rksuite_tfs(ct);

  /* print out new populations */
  for (k=0; k<ct->sites; k++)
    fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->sites + k]));

  /* we are on ct->surface at the moment
   * now, go over all of the other states,
   * and calculate switching probability */
  /*
  ct->surf_prob[ct->surface] = 0.;
  for (j=0; j<ct->sites; j++) if (j != ct->surface)
    ct->surf_prob[j] = -2 * (ct->tfs_popul[ct->surface] * ct->tfs_popul[j] + ct->tfs_popul[ct->sites + ct->surface] * ct->tfs_popul[ct->sites + j])
                          / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]))
			  * ct->tfs_overlap[j][ct->surface];
  */
  /* use the trick by Hammes-Schiffer and Tully, JCP 101, 4657 (1994) */
  for (j=0; j<ct->sites; j++) if (j != ct->surface)
    ct->surf_prob[j] *= -2. * ct->rk_timestep * ct->tfs_overlap[j][ct->surface]
                      / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]));


  /* write out the probabilities */
  for (k=0; k<ct->sites; k++)
    fprintf(f, "%10.6f", ct->surf_prob[k]);
 
  /* generate a random number */
  random_number = rand()/(double)RAND_MAX;
  cumulated_probability = 0.;

  /* and determine the surface to continue on */
  for (j=0; j<ct->sites; j++) if (j != ct->surface && ct->surf_prob[j] > 0.) {
    cumulated_probability += ct->surf_prob[j];
    if (cumulated_probability > random_number) {
      ct->surface = j;
      break;
    }
  }

  /* decoherence correction
   * proposed by Truhlar JCTC 1, 527 (2005)
   * recommended by Persico 126, 134114 (2007)
   * simplified to:
   * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
   * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
   */
  /* WITHOUT DECOHERENCE!
  popul_norm = 0.;
  for (j=0; j<ct->sites; j++) if (j != ct->surface) {
    decay_time = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
    ct->tfs_popul[j]             *= exp( - ct->rk_timestep / decay_time);
    ct->tfs_popul[ct->sites + j] *= exp( - ct->rk_timestep / decay_time);
    popul_norm += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[ct->sites + j]);
  }
  current_surface_factor = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface])));
  ct->tfs_popul[ct->surface]             *= current_surface_factor;
  ct->tfs_popul[ct->sites + ct->surface] *= current_surface_factor;
  */

  /* assign the new wave function to the eigenvector ct->surface */
  for (i=0; i<ct->sites; i++)
    ct->wf[i] = ct->tfs_vector[ct->surface][i];
  
  /* print out current state */
  fprintf(f, " %d", ct->surface);
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }
  fprintf(f, "\n");
  return 0;
}

int do_persico_diabatic_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2)
{
  int i, j, k, l, n, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
  double prob_factor;
  double random_number; // in the interval (0,1)
  double cumulated_probability; // sum of ct->surf_prob[0..k]
  double popul_norm; // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
  double decay_time, current_surface_factor; // in the decoherence algorithm
  double dot_product, tot_prob;
  /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
  int surf, surf_old; /* originally: istati, istati_old */
  double popk0, popk1, rekk, rekl, rnum, rden;
  double *acoer0, *acoei0, *acoer1, *acoei1;
  twodoubles ac0, ac1, vkk, vkl;
  float maxval;
  const float eps_switch = 1.e-9;

  ham = ct->hamiltonian_adiab;
  n = ct->sites;
  surf_old = surf = ct->surface;
  acoer0 = ct->tfs_popul_old;
  acoei0 = ct->tfs_popul_old + n;
  acoer1 = ct->tfs_popul;
  acoei1 = ct->tfs_popul + n;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  // copy as the first vector

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //fprintf(f, "Broyd %2d ", step);
      fprintf(f, " %2d ", step);
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
  //fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.5f", ct->ev_adiab[i]);
  }
  //fprintf(f, "    ");

  /* surface hopping tried according to Tully JCP 93, 1061 (1990)
   * "fewest switches" algorithm
   * diabatization according to Persico JCP 114, 10608 (2001)
   * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
   */

  /* exp_imag_matrix - has to be written yet !!! */

  /* THIS WILL REMAIN! */
  if (ct->tfs_initialization_step) {
    /* let us start with tfs_vector_old == tfs_vector */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ham[i * ct->sites + j];
    for (i=0; i<ct->sites; i++)
      ct->ev_adiab_old[i] = ct->ev_adiab[i];
  } else {
    /* push tfs_vector to tfs_vector_old */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
  }

  /* PUSH THE OLD POPULATIONS TO tfs_popul_old */
  for (i=0; i<ct->sites; i++) {
    ct->tfs_popul_old[i] = ct->tfs_popul[i];
    ct->tfs_popul_old[ct->sites + i] = ct->tfs_popul[ct->sites + i];
  }

  /* THIS WILL REMAIN, TOO
     MAYBE JUST DO NOT PERFORM THE CHECKING FOR SIGN CHANGE! */
  /* enter new eigenvectors to tfs_vector */
  for (i=0; i<ct->sites; i++) {
    /* calculate the dot product with tfs_vector_old[i] */
    dot_product = 0.;
    for (k=0; k<ct->sites; k++)
      dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->sites + k];
    /* for positive dot_product, take the vector as it is
     * for negative dot_product, take (-1) * the vector
     */
    for (k=0; k<ct->sites; k++)
      ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->sites + k] : - ham[i * ct->sites + k];
    //ct->tfs_vector[i][k] = ham[i * ct->sites + k];
  }

  /* THIS WILL REMAIN - tfs_overlap[i][j] corresponds to matrix T[i][j] ! */
  /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++) {
      ct->tfs_overlap[i][j] = 0.;
      for (k=0; k<ct->sites; k++)
        ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
    }

  /* print out the overlap matrix */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++)
      fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
  fprintf(f2, "\n");

  /* DO THE PERSICO STUFF HERE!
   * BEGINNING FROM THE SECOND STEP OF THE SIMULATION
   */
  if (ct->tfs_initialization_step) {
    ct->tfs_initialization_step = 0;
  } else {

    /* SET UP THE AVERAGE DIABATIC HAMILTONIAN Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt) */
    for (i=0; i<ct->sites; i++) {
      for (j=0; j<ct->sites; j++) {
        if (i==j) {
          ct->per_diab_hamiltonian[i][i] = ct->ev_adiab_old[i] / 2.;
        } else {
          ct->per_diab_hamiltonian[i][j] = 0.;
        }
        for (k=0; k<ct->sites; k++) 
          ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k] / 2.;
      }
    }
 
    /* CONSTRUCT THE PROPAGATION OPERATOR exp[-i*Z*dt] ! */
    exp_imag_matrix(ct->per_diab_hamiltonian, ct->per_propag_operator, ct->rk_timestep, ct->sites, ct->per_arrays);
 
    /* CONSTRUCT THE TRANSFORMATION MATRIX U = T^t * exp[-i*Z*dt] */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        ct->per_transformator[i][j][0] = 0.;
        ct->per_transformator[i][j][1] = 0.;
        for (k=0; k<ct->sites; k++) {
          ct->per_transformator[i][j][0] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][0];
          ct->per_transformator[i][j][1] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][1];
        }
      }
 
    /* OBTAIN THE NEW POPULATIONS */
    for (i=0; i<ct->sites; i++) {
      ct->tfs_popul[i] = ct->tfs_popul[ct->sites + i] = 0.;
      for (j=0; j<ct->sites; j++) {
        ct->tfs_popul[i]             += ct->per_transformator[i][j][0] * ct->tfs_popul_old[j]
                                      - ct->per_transformator[i][j][1] * ct->tfs_popul_old[ct->sites + j];
        ct->tfs_popul[ct->sites + i] += ct->per_transformator[i][j][0] * ct->tfs_popul_old[ct->sites + j]
                                      + ct->per_transformator[i][j][1] * ct->tfs_popul_old[j];
      }
    }
 
    /* print out new populations */
    for (k=0; k<ct->sites; k++)
      fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->sites + k]));
 
    /* we are on ct->surface at the moment
     * now, go over all of the other states,
     * and calculate switching probability */

    /* the following is inspired/taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */

    /* save the value U_kk */
    vkk[0] = ct->per_transformator[surf][surf][0];
    vkk[1] = ct->per_transformator[surf][surf][1];

    /* has a switch of states just occured? */
    if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch) {
      for (j=0; j<n; j++)
        ct->surf_prob[j] = 0.;
      maxval = 0.;
      for (j=0; j<n; j++)
        if (fabs((float) ct->tfs_overlap[surf_old][j]) > maxval) {
          surf = j;
          maxval = fabs((float) ct->tfs_overlap[surf_old][j]);
        }
    }

    /* rename / evaluate simple variables */
    ac0[0] = acoer0[surf_old]; 
    ac0[1] = acoei0[surf_old]; 
    ac1[0] = acoer1[surf]; 
    ac1[1] = acoei1[surf]; 
    popk0 = SQR(ac0[0]) + SQR(ac0[1]);
    popk1 = SQR(ac1[0]) + SQR(ac1[1]);
    /* rekk = U_kk * A_k(0) * A_k^*(dt) */
    rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] + ac0[1]*ac1[0]);
    rnum = (popk1 - popk0) / popk0;
    rden = popk1 - rekk;

    /* if population of occupied surface increases,
     * then no hopping!
     */
    if (popk1 > popk0) {
      for (j=0; j<n; j++)
        ct->surf_prob[j] = -1.;
    /* otherwise, calculate the probabilities */
    } else {
      for (j=0; j<n; j++) {
        /* no hopping from state 'surf' to itself, obviously */
        if (j == surf) {
          ct->surf_prob[j] = -1.;
        /* calculate the probability for all of the states except 'surf' */
        } else {
          vkl[0] = ct->per_transformator[surf][j][0];
          vkl[1] = ct->per_transformator[surf][j][1];
          rekl = vkl[0] * (acoer0[j]*ac1[0] + acoei0[j]*ac1[1]) + vkl[1] * (acoer0[j]*ac1[1] - acoei0[j]*ac1[0]);
          ct->surf_prob[j] = - rekl * rnum / rden;
        }
      }
    }

    tot_prob = 0.;
    for (j=0; j<n; j++)
      tot_prob += ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.;
    if (tot_prob > 1.)
      for (j=0; j<n; j++)
        if (ct->surf_prob[j] > 0.)
          ct->surf_prob[j] /= tot_prob;

    /* END OF PROBABILITIES! */
 
//  /* CONTINUE HERE WITH SWITCHING PROBABILITIES! */
//  ct->surf_prob[ct->surface] = -1.; /* to skip this surface in the test */
//  if (SQR(ct->per_transformator[ct->surface][ct->surface][0]) + SQR(ct->per_transformator[ct->surface][ct->surface][1]) > 0.999999) {
//    /* the population of current state does not really change
//     * numerical problems are to be expected
//     * therefore, skip the calculation of probabilities
//     * no surface hop will be performed in this step!
//     *
//     * IS THIS REALLY NECESSARY ???
//     *
//     */
//    for (j=0; j<ct->sites; j++) if (j != ct->surface)
//      ct->surf_prob[j] = -0.1;
//  } else {
//    /* otherwise, there is change in the population
//     * and it is safe to calculate surface-hopping probabilities! */
//    /* prob_factor is the second factor on the r.h.s. in Eq. 19 */
//    prob_factor = (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]) - SQR(ct->tfs_popul_old[ct->surface]) - SQR(ct->tfs_popul_old[ct->sites + ct->surface]))
//               / ( SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface])
//                  - ( (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->surface] - ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->sites + ct->surface])
//                        * ct->tfs_popul[ct->surface]
//                    + (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->sites + ct->surface] + ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->surface])
//                        * ct->tfs_popul[ct->sites + ct->surface] ) );
//    /* then, what follows is contribution of state psi_j to the increment of population of currently occupied state (surface)
//     * B_KL according to Eq. 19
//     */
//    for (j=0; j<ct->sites; j++)
//      b_kl[j] = - ((ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[j] - ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[ct->sites + j]) * ct->tfs_popul[ct->surface]
//                          + (ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[ct->sites + j] + ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[j]) * ct->tfs_popul[ct->sites + ct->surface])
//                         * prob_factor;

//    /* and surface hopping probabilities from Eq. 17 */
//    for (j=0; j<ct->sites; j++) if (j != ct->surface)
//      ct->surf_prob[j] = - ((ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[j] - ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[ct->sites + j]) * ct->tfs_popul[ct->surface]
//                          + (ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[ct->sites + j] + ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[j]) * ct->tfs_popul[ct->sites + ct->surface])
//                         * prob_factor;
//      /* ct->surf_prob[j] = -2 * (ct->tfs_popul[ct->surface] * ct->tfs_popul[j] + ct->tfs_popul[ct->sites + ct->surface] * ct->tfs_popul[ct->sites + j])
//                            / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]))
//	  		  * ct->tfs_overlap[j][ct->surface]; */
//  }
//
    /* write out the probabilities */
    for (k=0; k<ct->sites; k++)
      fprintf(f, "%10.6f", ct->surf_prob[k]);
   
    /* generate a random number */
    random_number = rand()/(double)RAND_MAX;
    cumulated_probability = 0.;
  
    /* and determine the surface to continue on */
    for (j=0; j<ct->sites; j++) if (j != ct->surface && ct->surf_prob[j] > 0.) {
      cumulated_probability += ct->surf_prob[j];
      if (cumulated_probability > random_number) {
        ct->surface = j;
        break;
      }
    }
  
    /* decoherence correction
     * proposed by Truhlar JCTC 1, 527 (2005)
     * recommended by Persico 126, 134114 (2007)
     * simplified to:
     * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
     * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
     */

    if (ct->decoherence) {
      popul_norm = 0.;
      for (j=0; j<n; j++)
        if (j != ct->surface) {
          decay_time = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
          ct->tfs_popul[j]             *= exp( - ct->rk_timestep / decay_time);
          ct->tfs_popul[ct->sites + j] *= exp( - ct->rk_timestep / decay_time);
          popul_norm += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
        }
      current_surface_factor = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
      ct->tfs_popul[ct->surface]             *= current_surface_factor;
      ct->tfs_popul[ct->sites + ct->surface] *= current_surface_factor;
    }

  } /* end if initialization step else */

  /* assign the new wave function to the eigenvector ct->surface */
  for (i=0; i<ct->sites; i++)
    ct->wf[i] = ct->tfs_vector[ct->surface][i];
  
  /* print out current state */
  fprintf(f, " %d", ct->surface);
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }
  fprintf(f, "\n");

  /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
  for (i=0; i<ct->sites; i++)
  ct->ev_adiab_old[i] = ct->ev_adiab[i];

  return 0;
}

int do_persico_diabatic_sfhopping_new(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3)
{
  int i, j, k, l, n, step;
  double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
  double prob_factor;
  double random_number; // in the interval (0,1)
  double cumulated_probability; // sum of ct->surf_prob[0..k]
  double popul_norm; // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
  double decay_time, current_surface_factor; // in the decoherence algorithm
  double dot_product, tot_prob;
  /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
  int surf, surf_old; /* originally: istati, istati_old */
  double popk0, popk1, rekk, rekl, rnum, rden;
  double *acoer0, *acoei0, *acoer1, *acoei1;
  twodoubles ac0, ac1, vkk, vkl;
  float maxval;
  const float eps_switch = 1.e-9;
  static double **pop_flow=NULL;
  int lwork;
  static double *hdx=NULL, *ex=NULL, *work=NULL, **rotr=NULL, **roti=NULL, **rotsr=NULL, **rotsi=NULL, **rotr_new=NULL, **roti_new=NULL;
  double timetau, ld_sk, ld_ck, ld_fac;
  int ld_istep;
  const int ld_nstep = 20;

  ham = ct->hamiltonian_adiab;
  n = ct->sites;
  surf_old = surf = ct->surface;
  acoer0 = ct->tfs_popul_old;
  acoei0 = ct->tfs_popul_old + n;
  acoer1 = ct->tfs_popul;
  acoei1 = ct->tfs_popul + n;

  if (pop_flow == NULL) {
    snew(pop_flow, n);
    snew(pop_flow[0], SQR(n));
    for (i=1; i<n; i++)
      pop_flow[i] = pop_flow[0] + i * n;
  }

  if (hdx == NULL) {
    snew(hdx, SQR(n));
  }

  if (ex == NULL) {
    snew(ex, n);
  }

  if (work == NULL) {
    lwork = (n > 3) ? SQR(n) : 5;
    snew(work, lwork);
  }

  if (rotr == NULL) {
    snew(rotr, n);
    snew(rotr[0], SQR(n));
    for (i=1; i<n; i++)
      rotr[i]  = rotr[0] + i * n;
    snew(roti, n);
    snew(roti[0], SQR(n));
    for (i=1; i<n; i++)
      roti[i]  = roti[0] + i * n;
    snew(rotsr, n);
    snew(rotsr[0], SQR(n));
    for (i=1; i<n; i++)
      rotsr[i] = rotsr[0] + i * n;
    snew(rotsi, n);
    snew(rotsi[0], SQR(n));
    for (i=1; i<n; i++)
      rotsi[i] = rotsi[0] + i * n;
    snew(rotr_new, n);
    snew(rotr_new[0], SQR(n));
    for (i=1; i<n; i++)
      rotr_new[i]  = rotr_new[0] + i * n;
    snew(roti_new, n);
    snew(roti_new[0], SQR(n));
    for (i=1; i<n; i++)
      roti_new[i]  = roti_new[0] + i * n;
  }

  // calculate the charges from the wave function
  // and copy as the first vector
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy1 = energy2 = 0.0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy1 += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }
    energy = energy1 + energy2;

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      //fprintf(f, "Broyd %2d ", step);
      fprintf(f, " %2d ", step);
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
  //fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.5f", ct->ev_adiab[i]);
  }
  //fprintf(f, "    ");

  /* surface hopping tried according to Tully JCP 93, 1061 (1990)
   * "fewest switches" algorithm
   * diabatization according to Persico JCP 114, 10608 (2001)
   * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
   */

  /* exp_imag_matrix - has to be written yet !!! */

  /* 1. set up the state vectors */

  /*   a. back up state vectors from the previous step, tfs_vector_old*/

  if (ct->tfs_initialization_step) {
    /* let us start with tfs_vector_old == tfs_vector */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ham[i * ct->sites + j];
    for (i=0; i<ct->sites; i++)
      ct->ev_adiab_old[i] = ct->ev_adiab[i];
  } else {
    /* push tfs_vector to tfs_vector_old */
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++)
        ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
  }

  /*   b. back up populations (complex coefficients) of states from the previous step, tfs_popul_old */

  for (i=0; i<ct->sites; i++) {
    ct->tfs_popul_old[i]             = ct->tfs_popul[i];
    ct->tfs_popul_old[ct->sites + i] = ct->tfs_popul[ct->sites + i];
  }

  /*  c. save the new state vectors */

  for (i=0; i<ct->sites; i++) {
    /* calculate the dot product with tfs_vector_old[i],
       to check if we shall invert the sign
     */
    dot_product = 0.;
    for (k=0; k<ct->sites; k++)
      dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->sites + k];
    /* for positive dot_product, take the vector as it is
     * for negative dot_product, take (-1) * the vector
     */
    for (k=0; k<ct->sites; k++)
      ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->sites + k] : - ham[i * ct->sites + k];
  }
  
  /* 2. calculate the overlap of state vectors
          tfs_overlap[i][j] is the matrix T[i][j] = <tfs_vector_old[i] | tfs_vector[j]>
   */

  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++) {
      ct->tfs_overlap[i][j] = 0.;
      for (k=0; k<ct->sites; k++)
        ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
    }

  /* print out the overlap matrix */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++)
      fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
  fprintf(f2, "\n");

  /* print out the state vectors */
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++)
      fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
  fprintf(f3, "\n");

  /* check if the matrix is unitary - in our case orthogonal (real matrix) */
  printf("Testing the orthogonality of T - this should be a unit matrix:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      dot_product = 0.;
      for (k=0; k<n; k++)
        dot_product += ct->tfs_overlap[i][k] * ct->tfs_overlap[j][k];
      printf(" %9.6f", dot_product);
    }
    printf("\n");
  }

  /* 3. do the propagation of the wave function / populations
          this is equivalent to the propag_ld function in NewtonX */
	  
  /*    a. do not do it in the first step of the simulation */

  if (ct->tfs_initialization_step) {
    ct->tfs_initialization_step = 0;
  } else {

  /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
  // /*    b. set up the average diabatic hamiltonian Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt) */
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        ct->per_diab_hamiltonian[i][j] = 0.;
        for (k=0; k<n; k++) 
          ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k];
      }
    }

  /*    c. divide the propagation into smaller elementary steps */
    timetau = ct->rk_timestep / (double) ld_nstep;
    for (ld_istep=0; ld_istep<ld_nstep; ld_istep++) {
      /* matrix Hdx                 = (Hdia(istep*timetau) + Hdia((istep+1)*tau)) / 2
         matrix Hdia(istep*timetau) = E(0) + (2*istep-1)/(2*ld_nstep) * (Hdia(dt) - E(0))
       */
      ld_fac = (2. * ld_istep + 1) / (2 * ld_nstep);
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
	  hdx[i + j * n] = timetau * (ld_fac * ct->per_diab_hamiltonian[i][j] + (i==j ? (1. - ld_fac) * ct->ev_adiab_old[i] : 0.));
      /* eigenproblem of exp[-i * Hdx * timetau] */
      dsyev(n, hdx, ex, work, lwork);
      /* construct the rotation matrix (real + imaginary parts) for the elementary step */
      for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
	  rotsr[i][j] = 0.;
	  rotsi[i][j] = 0.;
	}
      for (k=0; k<n; k++) {
        ld_sk = sin(ex[k]);
        ld_ck = cos(ex[k]);
	for (i=0; i<n; i++)
          for (j=0; j<n; j++) {
	    rotsr[i][j] += hdx[i + k * n] * hdx[j + k * n] * ld_ck;
	    rotsi[i][j] -= hdx[i + k * n] * hdx[j + k * n] * ld_sk;
	  }
      }
      /* accumulate the overall rotation matrix (real + imaginary) */
      if (ld_istep == 0) {
	for (i=0; i<n; i++)
          for (j=0; j<n; j++) {
	    rotr[i][j] = rotsr[i][j];
	    roti[i][j] = rotsi[i][j];
	  }
      } else {
        /* tw  = rotsr * rotr
	   tx  = rotsi * roti
	   ty  = rotsr * roti
	   tz  = rotsi * rotr
	   rotr = tw - tx
	   roti = ty + tz
	 */
	for (i=0; i<n; i++)
          for (j=0; j<n; j++) {
	    rotr_new[i][j] = 0.;
	    roti_new[i][j] = 0.;
	    for (k=0; k<n; k++) {
	      rotr_new[i][j] += rotsr[i][k] * rotr[k][j] - rotsi[i][k] * roti[k][j];
	      roti_new[i][j] += rotsr[i][k] * roti[k][j] + rotsi[i][k] * rotr[k][j];
	    }
	  }
	for (i=0; i<n; i++)
          for (j=0; j<n; j++) {
	    rotr[i][j] = rotr_new[i][j];
	    roti[i][j] = roti_new[i][j];
	  }
      }
      /* print out the rotation matrix */
      printf("Rotation matrix for coefficients in the adiabatic representation:\n");
      for (i=0; i<n; i++) {
        for (j=0; j<n; j++)
          printf("  %9.6f + %9.6fi", rotr[i][j], roti[i][j]);
        printf("\n");
      }
  /*    d. update the coefficients in the adiabatic representation */
      for (i=0; i<n; i++) {
        ct->tfs_popul[i]     = 0.;
	ct->tfs_popul[i + n] = 0.;
	for (j=0; j<n; j++)
	  for (k=0; k<n; k++) {
	    ct->tfs_popul[i]     += ct->tfs_overlap[j][i] *
	                    (ct->tfs_popul_old[k] * rotr[j][k] - ct->tfs_popul_old[k + n] * roti[j][k]);
	    ct->tfs_popul[i + n] += ct->tfs_overlap[j][i] *
	                    (ct->tfs_popul_old[k] * roti[j][k] + ct->tfs_popul_old[k + n] * rotr[j][k]);
	  }
      }
    }
 

    /* GO ON HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */



    /* print out new populations */
    for (k=0; k<ct->sites; k++)
      fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->sites + k]));
 
    /* we are on ct->surface at the moment
     * now, go over all of the other states,
     * and calculate switching probability */

    /* the following is inspired/taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */

    /* save the value U_kk */
    vkk[0] = ct->per_transformator[surf][surf][0];
    vkk[1] = ct->per_transformator[surf][surf][1];

    /* has a switch of states just occured? */
    if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch) {
      for (j=0; j<n; j++)
        ct->surf_prob[j] = 0.;
      maxval = 0.;
      for (j=0; j<n; j++)
        if (fabs((float) ct->tfs_overlap[surf_old][j]) > maxval) {
          surf = j;
          maxval = fabs((float) ct->tfs_overlap[surf_old][j]);
        }
    }

    /* rename / evaluate simple variables */
    ac0[0] = acoer0[surf_old]; 
    ac0[1] = acoei0[surf_old]; 
    ac1[0] = acoer1[surf]; 
    ac1[1] = acoei1[surf]; 
    popk0 = SQR(ac0[0]) + SQR(ac0[1]);
    popk1 = SQR(ac1[0]) + SQR(ac1[1]);
    /* rekk = U_kk * A_k(0) * A_k^*(dt) */
    rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] + ac0[1]*ac1[0]);
    rnum = (popk1 - popk0) / popk0;
    rden = popk1 - rekk;

    /* if population of occupied surface increases,
     * then no hopping!
     */
    if (popk1 > popk0) {
      for (j=0; j<n; j++)
        ct->surf_prob[j] = -1.;
    /* otherwise, calculate the probabilities */
    } else {
      for (j=0; j<n; j++) {
        /* no hopping from state 'surf' to itself, obviously */
        if (j == surf) {
          ct->surf_prob[j] = -1.;
        /* calculate the probability for all of the states except 'surf' */
        } else {
          vkl[0] = ct->per_transformator[surf][j][0];
          vkl[1] = ct->per_transformator[surf][j][1];
          rekl = vkl[0] * (acoer0[j]*ac1[0] + acoei0[j]*ac1[1]) + vkl[1] * (acoer0[j]*ac1[1] - acoei0[j]*ac1[0]);
          ct->surf_prob[j] = - rekl * rnum / rden;
        }
      }
    }
    for (k=0; k<n; k++) {
      vkk[0] = ct->per_transformator[k][k][0];
      vkk[1] = ct->per_transformator[k][k][1];
      ac0[0] = acoer0[k]; 
      ac0[1] = acoei0[k]; 
      ac1[0] = acoer1[k]; 
      ac1[1] = acoei1[k]; 
      popk0 = SQR(ac0[0]) + SQR(ac0[1]);
      popk1 = SQR(ac1[0]) + SQR(ac1[1]);
      /* rekk = U_kk * A_k(0) * A_k^*(dt) */
      rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] + ac0[1]*ac1[0]);
      rnum = (popk1 - popk0) / popk0;
      rden = popk1 - rekk;
      for (l=0; l<n; l++) {
        vkl[0] = ct->per_transformator[k][l][0];
        vkl[1] = ct->per_transformator[k][l][1];
        rekl = vkl[0] * (acoer0[l]*ac1[0] + acoei0[l]*ac1[1]) + vkl[1] * (acoer0[l]*ac1[1] - acoei0[l]*ac1[0]);
        pop_flow[k][l] = - rekl * rnum * popk0 / rden;
        printf(" %8.5f", pop_flow[k][l]);
      }
    }
    printf("\n");

    tot_prob = 0.;
    for (j=0; j<n; j++)
      tot_prob += ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.;
    if (tot_prob > 1.)
      for (j=0; j<n; j++)
        if (ct->surf_prob[j] > 0.)
          ct->surf_prob[j] /= tot_prob;

    /* END OF PROBABILITIES! */
 
    /* write out the probabilities */
    for (k=0; k<ct->sites; k++)
      fprintf(f, "%10.6f", ct->surf_prob[k]);
   
    /* generate a random number */
    random_number = rand()/(double)RAND_MAX;
    cumulated_probability = 0.;
  
    /* and determine the surface to continue on */
    for (j=0; j<ct->sites; j++) if (j != ct->surface && ct->surf_prob[j] > 0.) {
      cumulated_probability += ct->surf_prob[j];
      if (cumulated_probability > random_number) {
        ct->surface = j;
        break;
      }
    }
  
    /* decoherence correction
     * proposed by Truhlar JCTC 1, 527 (2005)
     * recommended by Persico 126, 134114 (2007)
     * simplified to:
     * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
     * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
     */

    if (ct->decoherence) {
      popul_norm = 0.;
      for (j=0; j<n; j++)
        if (j != ct->surface) {
          decay_time = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
          ct->tfs_popul[j]             *= exp( - ct->rk_timestep / decay_time);
          ct->tfs_popul[ct->sites + j] *= exp( - ct->rk_timestep / decay_time);
          popul_norm += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
        }
      current_surface_factor = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
      ct->tfs_popul[ct->surface]             *= current_surface_factor;
      ct->tfs_popul[ct->sites + ct->surface] *= current_surface_factor;
    }

  } /* end if initialization step else */

  /* assign the new wave function to the eigenvector ct->surface */
  for (i=0; i<ct->sites; i++)
    ct->wf[i] = ct->tfs_vector[ct->surface][i];
  
  /* print out current state */
  fprintf(f, " %d", ct->surface);
  for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }
  fprintf(f, "\n");

  /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
  for (i=0; i<ct->sites; i++)
  ct->ev_adiab_old[i] = ct->ev_adiab[i];

  return 0;
}

int do_wf_decomp_evec_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
  int i, j, k, step;
  double old_energy, energy, *ham, fermi_energy, fermi_upper, fermi_lower, dotproduct_real, dotproduct_imag;

  ham = ct->hamiltonian_adiab;

  // calculate the charges from the wave function
  for (i=0; i<ct->sites; i++)
    ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

  energy = 1.e10;

  // SCC cycle
  for (step = 0; ; step++) {

    // too many iterations? break and indicate failure
    if (step >= MAXITER_BROYDEN) {
      fprintf(f, "Broyden/Fermi failed to converge!\n");
      return 1;
    }

    // construct the self-consistent Hamiltonian
    for (i=0; i<ct->sites; i++) {
      ham[i + i * ct->sites] = ct->hamiltonian[i][i];
      for (j=0; j<ct->sites; j++) {
        ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
        if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j];
      }
    }

    // diagonalize this Hamiltonian
    dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
    
    // FERMI DISTRIBUTION + BROYDEN MIXING

    // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, fermi_coeff) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->sites + i]);
    }

    // calculate the energy and check convergence
    old_energy = energy;
    energy = 0.e0;
    // indices i and j run over the sites
    for (i=0; i<ct->sites; i++)
      for (j=0; j<ct->sites; j++) {
        // TB Hamiltonian
        // index k runs over the adiabatic states (eigenvectors)
        for (k=0; k<ct->sites; k++)
          energy += fermi_coeff[k] * ham[k * ct->sites + i] * ham[k * ct->sites + j] * ct->hamiltonian[i][j];
        // Hubbard / gamma terms
        energy += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
      }

    if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL) {
      // Broyden mixing has converged
      break;
    }

    // mix the charges - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];

  } /* end SCC cycle */

  // print the energies of "adiabatic states"
  for (i=0; i<ct->sites; i++)
    fprintf(f, "%10.7f", ct->ev_adiab[i]);

  // calculate the expansion coefficients
  fprintf(f, " WF expansion:");
  for (i=0; i<ct->sites; i++) {
    dotproduct_real = 0.e0;
    for (j=0; j<ct->sites; j++)
      dotproduct_real += ham[i * ct->sites + j] * ct->wf[j];
    dotproduct_imag = 0.e0;
    for (j=0; j<ct->sites; j++)
      dotproduct_imag += ham[i * ct->sites + j] * ct->wf[j + ct->sites];
    fprintf(f, "%9.6f", SQR(dotproduct_real)+SQR(dotproduct_imag));
  }
  /*
  for (i=0; i<ct->sites; i++) {
    dotproduct_real = dotproduct_imag = 0.e0;
    for (j=0; j<ct->sites; j++)
      for (k=0; k<ct->sites; k++)
        dotproduct_squared += ham[i * ct->sites + j] * ham[i * ct->sites + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->sites] * ct->wf[k + ct->sites]);
    fprintf(f, "%9.6f", dotproduct_squared);
  }
  */
  fprintf(f, "\n");


  return 0;
}

void check_and_invert_orbital_phase(dftb_phase1_t *dftb1, int sites, int *homo, FILE *f, int step)
{
  int i, k;
  double dotprod;

  if (f)
    fprintf(f, "%10d", step);
  for (i=0; i<sites; i++) {
    /* look at orbitals in site i */
    /* calculate the dot product of HOMO with HOMO in the previous step */
    dotprod = 0.;
    for (k=0; k<dftb1[i].norb; k++)
      dotprod += dftb1[i].a_old[k][homo[i]-1] * dftb1[i].a[k][homo[i]-1];
    /* print out the dot product */
    if (f)
      fprintf(f, " %7.4f", dotprod);
    /* check the dot product:
     * if negative: invert a[..][homo[i]-1]
     */
    if (dotprod < 0.)
      for (k=0; k<dftb1[i].norb; k++)
        dftb1[i].a[k][homo[i]-1] = - dftb1[i].a[k][homo[i]-1];
    /* update the "old" array */
    for (k=0; k<dftb1[i].norb; k++)
      dftb1[i].a_old[k][homo[i]-1] = dftb1[i].a[k][homo[i]-1];
  }
  if (f)
    fprintf(f, "\n");

  return;
}

long exp_imag_matrix(double **in, twodoubles **out, double dt, long n, ct_per_orthogo_t *arrays)
{
  extern void dsyevr_(char *, char *, char *, long *, double *, long *,
                double *, double *, long *, long *, double *, long *, double *,
		double *, long *, double *, double *, long *, double *, long *, long *);
  double abstol=1.0e-8;
  char itype='U';
  char jobz='V';
  char range='A';
  char uplo='U';
  long i, j, k, m, info, eval_no;
  double void_double;
  long void_long;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      arrays->in[i+j*n] = in[i][j];
    }

/*
  // print the arrays
  printf("Orthogonalize - array sij:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.sij[i+j*n]);
    printf("\n");
  }
  printf("Orthogonalize - array tij:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.tij[i+j*n]);
    printf("\n");
  }
*/

  // Diagonalize the matrix
  // evec - eigenvectors, eval - eigenvalues
  dsyevr_(&jobz, &range, &uplo, &n, arrays->in, &n,
          &void_double, &void_double, &void_long, &void_long, &abstol, &eval_no, arrays->eval,
	  arrays->evec, &n, arrays->issupz, arrays->work, &(arrays->lwork), arrays->iwork, &(arrays->liwork), &info);
  //if (info) return info;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      out[i][j][0] = 0.;
      out[i][j][1] = 0.;
      for (k=0; k<n; k++) {
        out[i][j][0] += arrays->evec[i + n * k] * arrays->evec[j + n * k] * cos(arrays->eval[k] * dt);
        out[i][j][1] -= arrays->evec[i + n * k] * arrays->evec[j + n * k] * sin(arrays->eval[k] * dt);
      }
    }

  return info;
}
