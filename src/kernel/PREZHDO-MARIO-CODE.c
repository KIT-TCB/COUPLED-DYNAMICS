int do_prezhdo_sfhopping(charge_transfer_t *ct, int ext_step, FILE *f, FILE *f2, FILE *f3)
{
  int i, j, k, l, m, n, step, info;
  double q_sum; //Needed for normalization after fixing broyden for charges not in intervall 0,1
  double old_energy, energy, energy1, energy2, *ham;
  double prob_factor;
  double random_number; // in the interval (0,1)
  double cumulated_probability; // sum of ct->surf_prob[0..k]
  double popul_norm; // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
  double decay_time, current_surface_factor; // in the decoherence algorithm
  double dot_product, tot_prob;
  double dot_prod_r, dot_prod_i;
  /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
  int surf, surf_old; /* originally: istati, istati_old */
  double popk0, popk1, rekk, rekl, rnum, rden;
  double *acoer0, *acoei0, *acoer1, *acoei1;
  twodoubles ac0, ac1, vkk, vkl;
  float maxval;
  const float eps_switch = 1.e-9;
  static double **pop_flow=NULL;
  static int lwork=0;
  static double *hdx=NULL, *ex=NULL, *work=NULL, **rotr=NULL, **roti=NULL, **rotsr=NULL, **rotsi=NULL, **rotr_new=NULL, **roti_new=NULL;
  double timetau, ld_sk, ld_ck, ld_fac;
  int ld_istep;
  const int ld_nstep = 20;
  double diff_flow, gross_flow;

  ham = ct->hamiltonian_adiab;
  n = ct->sites;
  surf_old = surf = ct->surface;
  acoer0 = ct->tfs_popul_old;
  acoei0 = ct->tfs_popul_old + n;
  acoer1 = ct->tfs_popul;
  acoei1 = ct->tfs_popul + n;

  if (pop_flow == NULL) {
    pop_flow = (double**) calloc (ct->sites,sizeof(double));
    //snew(pop_flow, n);
    pop_flow[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(pop_flow[0], SQR(n));
    for (i=1; i<n; i++)
      pop_flow[i] = pop_flow[0] + i * n;
  }

  if (hdx == NULL) {
    hdx = (double*) calloc (SQR(ct->sites),sizeof(double));	  
    //snew(hdx, SQR(n));
  }

  if (ex == NULL) {
    ex = (double*) calloc (ct->sites,sizeof(double));
    //snew(ex, n);
  }

  if (work == NULL) {
    lwork = (n >= 3) ? SQR(n) : 5;
    work = (double*) calloc (lwork,sizeof(double));
    //snew(work, lwork);
  }

  if (rotr == NULL) {
    rotr = (double**) calloc (ct->sites,sizeof(double));
    //snew(rotr, n);
    rotr[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(rotr[0], SQR(n));
    for (i=1; i<n; i++)
      rotr[i]  = rotr[0] + i * n;
    roti = (double**) calloc (ct->sites,sizeof(double));
    //snew(roti, n);
    roti[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(roti[0], SQR(n));
    for (i=1; i<n; i++)
      roti[i]  = roti[0] + i * n;
    rotsr = (double**) calloc (ct->sites,sizeof(double));
    //snew(rotsr, n);
    rotsr[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(rotsr[0], SQR(n));
    for (i=1; i<n; i++)
      rotsr[i] = rotsr[0] + i * n;
    rotsi = (double**) calloc (ct->sites,sizeof(double));
    //snew(rotsi, n);
    rotsi[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(rotsi[0], SQR(n));
    for (i=1; i<n; i++)
      rotsi[i] = rotsi[0] + i * n;
    rotr_new = (double**) calloc (ct->sites,sizeof(double));
    //snew(rotr_new, n);
    rotr_new[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(rotr_new[0], SQR(n));
    for (i=1; i<n; i++)
      rotr_new[i]  = rotr_new[0] + i * n;
    roti_new = (double**) calloc (ct->sites,sizeof(double));
    //snew(roti_new, n);
    roti_new[0] = (double*) calloc (SQR(ct->sites),sizeof(double));
    //snew(roti_new[0], SQR(n));
    for (i=1; i<n; i++)
      roti_new[i]  = roti_new[0] + i * n;
  }

  // construct the self-consistent Hamiltonian
  for (i=0; i<ct->sites; i++) {
    for (j=0; j<ct->sites; j++) {
      ham[i + j * ct->sites] = ct->hamiltonian[i][j];
    }
  }
  
  // diagonalize this Hamiltonian
  info = dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
  if (info) {printf("dsyev returned %d, exiting\n", info); exit(-1);} 
  
  // mix the eigenvectors, with the appropriate coefficients
  for (i=0; i<ct->sites; i++) { // i runs over the sites
    ct->q_act[i] = SQR(ham[i]);
  }

  energy = ct->ev_adiab[0];

  // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
  //fprintf(f, " Eigenvalues:");
  for (i=0; i<ct->sites; i++) {
   if( ext_step % ct->ntout == 0 ) fprintf(f, "%8.5f", ct->ev_adiab[i]);
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
  if( ext_step % ct->ntout == 0 ) {
  for (i=0; i<ct->sites; i++)
    for (j=0; j<ct->sites; j++)
     fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
  fprintf(f2, "\n");
  }

  /* print out the state vectors */
//  for (i=0; i<ct->sites; i++)
//    for (j=0; j<ct->sites; j++)
//      fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
//  fprintf(f3, "\n");

  /* 3. do the propagation of the wave function / populations
          this is equivalent to the propag_ld function in NewtonX */
	  
  /*    a. do not do it in the first step of the simulation */

  if (ct->tfs_initialization_step) {
    ct->tfs_initialization_step = 0;
  } else {

  /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
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
      ld_fac = (double) (2. * ld_istep + 1) / (double) (2 * ld_nstep);
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

    /* print out new populations */
    for (k=0; k<ct->sites; k++)
      if( ext_step % ct->ntout == 0 ) fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->sites + k]));
 
    /* we are on ct->surface at the moment
     * now, go over all of the other states,
     * and calculate switching probability */

    /* if population of occupied surface increases,
     * then no hopping!
     */
    ac0[0] = acoer0[surf]; 
    ac0[1] = acoei0[surf]; 
    ac1[0] = acoer1[surf]; 
    ac1[1] = acoei1[surf]; 
    popk0 = SQR(ac0[0]) + SQR(ac0[1]);
    popk1 = SQR(ac1[0]) + SQR(ac1[1]);
    if (popk1 > popk0) {
      for (j=0; j<n; j++)
        ct->surf_prob[j] = 0.;

    /* otherwise, calculate the probabilities */
    } else {
      /* calculate the gross flow */
      gross_flow = 0.;
      for (j=0; j<n; j++) {
        diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
        if (diff_flow > 0.)
          gross_flow += diff_flow;
      }
      /* go over the states */        
      for (j=0; j<n; j++) {
        /* no hopping from state 'surf' to itself, obviously */
        if (j == surf) {
          ct->surf_prob[j] = 0.;
        /* calculate the probability for all of the states except 'surf' */
        } else {
          diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
	  /* no hopping to a state with a decreased population */
	  if (diff_flow < 0.) {
	    ct->surf_prob[j] = 0.;
	  /* otherwise, calculate probability according to Prezhdo */
	  } else {
	    ct->surf_prob[j] = diff_flow * (popk0 - popk1) / (popk0 * gross_flow);
	  }
        }
      }
    }
    
  
    tot_prob = 0.;
    for (j=0; j<n; j++)
      tot_prob += (ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.);
    if (tot_prob > 1.)
      for (j=0; j<n; j++)
        if (ct->surf_prob[j] > 0.)
          ct->surf_prob[j] /= tot_prob;

    /* END OF PROBABILITIES! */
 
    /* write out the probabilities */
    for (k=0; k<ct->sites; k++)
      if( ext_step % ct->ntout == 0 ) fprintf(f, "%10.6f", ct->surf_prob[k]);
   
    /* generate a random number */
    random_number = 100.;
    cumulated_probability = 0.;
  
    /* and determine the surface to continue on */
    for (j=0; j<ct->sites; j++) if (j != ct->surface && ct->surf_prob[j] > 0.) {
      cumulated_probability += ct->surf_prob[j];
      if (cumulated_probability > random_number) {
	fprintf(f3, "%10d %d %10.6f\n", ext_step, j, (ct->ev_adiab[j] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
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
  if( ext_step % ct->ntout == 0 ) {
    fprintf(f, " %d", ct->surface);
  /*for (i=0; i<ct->sites; i++) {
    fprintf(f, "%8.4f", SQR(ct->wf[i]));
  }*/
  fprintf(f, "\n");
  }

  /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
  for (i=0; i<ct->sites; i++)
  ct->ev_adiab_old[i] = ct->ev_adiab[i];

  return 0;
}
