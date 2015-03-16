/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "md_support.h"
#include "md_logging.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
#include "ionize.h"
#include "disre.h"
#include "orires.h"
#include "pme.h"
#include "mdatoms.h"
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "domdec_network.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "txtdump.h"
#include "string2.h"
#include "pme_loadbal.h"
#include "bondf.h"
#include "membed.h"
#include "types/nlistheuristics.h"
#include "types/iteratedconstraints.h"
#include "nbnxn_cuda_data_mgmt.h"

/* TOMAS KUBAR */
#include "charge_transfer.h"
#include <time.h>

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

static void reset_all_counters(FILE *fplog, t_commrec *cr,
                               gmx_large_int_t step,
                               gmx_large_int_t *step_rel, t_inputrec *ir,
                               gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                               gmx_runtime_t *runtime,
                               nbnxn_cuda_ptr_t cu_nbv)
{
    char sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    md_print_warn(cr, fplog, "step %s: resetting all time and cycle counters\n",
                  gmx_step_str(step, sbuf));

    if (cu_nbv)
    {
        nbnxn_cuda_reset_timings(cu_nbv);
    }

    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    if (DOMAINDECOMP(cr))
    {
        reset_dd_statistics_counters(cr->dd);
    }
    init_nrnb(nrnb);
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel      = 0;
    wallcycle_start(wcycle, ewcRUN);
    runtime_start(runtime);
    print_date_and_time(fplog, cr->nodeid, "Restarted time", runtime);
}

double do_md(FILE *fplog, t_commrec *cr, int nfile, const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose, gmx_bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int stepout, t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
             gmx_edsam_t ed, t_forcerec *fr,
             int repl_ex_nst, int repl_ex_nex, int repl_ex_seed, gmx_membed_t membed,
             real cpt_period, real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime)
{
    gmx_mdoutf_t   *outf;
    gmx_large_int_t step, step_rel;
    double          run_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM, bRerunMD, bNotLastFrame = FALSE,
                    bFirstStep, bStateFromCP, bStateFromTPX, bInitStep, bLastStep,
                    bBornRadii, bStartingFromCpt;
    gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
    gmx_bool          do_ene, do_log, do_verbose, bRerunWarnNoV = TRUE,
                      bForceUpdate = FALSE, bCPT;
    int               mdof_flags;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
    int               i, m;
    t_trxstatus      *status;
    rvec              mu_tot;
    t_vcm            *vcm;
    t_state          *bufstate = NULL;
    matrix           *scale_tot, pcoupl_mu, M, ebox;
    gmx_nlheur_t      nlh;
    t_trxframe        rerun_fr;
    gmx_repl_ex_t     repl_ex = NULL;
    int               nchkpt  = 1;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin = NULL;
    df_history_t      df_history;
    t_state          *state    = NULL;
    rvec             *f_global = NULL;
    int               n_xtc    = -1;
    rvec             *x_xtc    = NULL;
    gmx_enerdata_t   *enerd;
    rvec             *f = NULL;
    gmx_global_stat_t gstat;
    gmx_update_t      upd   = NULL;
    t_graph          *graph = NULL;
    globsig_t         gs;
    gmx_rng_t         mcrng = NULL;
    gmx_bool          bFFscan;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind, *ekind_save;
    gmx_shellfc_t     shellfc;
    int               count, nconverged = 0;
    real              timestep = 0;
    double            tcount   = 0;
    gmx_bool          bIonize  = FALSE;
    gmx_bool          bTCR     = FALSE, bConverged = TRUE, bOK, bSumEkinhOld, bExchanged;
    gmx_bool          bAppend;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    gmx_bool          bVV, bIterativeCase, bFirstIterate, bTemp, bPres, bTrotter;
    gmx_bool          bUpdateDoLR;
    real              mu_aver = 0, dvdl;
    int               a0, a1, gnx = 0, ii;
    atom_id          *grpindex = NULL;
    char             *grpname;
    t_coupl_rec      *tcr     = NULL;
    rvec             *xcopy   = NULL, *vcopy = NULL, *cbuf = NULL;
    matrix            boxcopy = {{0}}, lastbox;
    tensor            tmpvir;
    real              fom, oldfom, veta_save, pcurr, scalevir, tracevir;
    real              vetanew = 0;
    int               lamnew  = 0;
    /* for FEP */
    int               nstfep;
    real              rate;
    double            cycles;
    real              saved_conserved_quantity = 0;
    real              last_ekin                = 0;
    int               iter_i;
    t_extmass         MassQ;
    int             **trotter_seq;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/
    gmx_iterate_t     iterate;
    gmx_large_int_t   multisim_nsteps = -1;                        /* number of steps to do  before first multisim
                                                                      simulation stops. If equal to zero, don't
                                                                      communicate any more between multisims.*/
    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t pme_loadbal = NULL;
    double               cycles_pmes;
    gmx_bool             bPMETuneTry = FALSE, bPMETuneRunning = FALSE;

#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

    /* TOMAS KUBAR */
    charge_transfer_t *ct=NULL;
    dftb_t *dftb=NULL;
    char slko_path[MAX_PATH_LENGTH];
    dftb_broyden_t *ct_broyden=NULL;
    ct_diis_t *ct_diis=NULL;
    FILE *f_tb_hamiltonian=NULL, *f_tb_hamiltonian_hub=NULL, *f_tb_hubbard=NULL, *f_tb_occupation=NULL, *f_ct_esp=NULL, *f_ct_shift=NULL, *f_ct_energy=NULL,
         *f_ct_adiabatic=NULL, *f_ct_exp_adiab=NULL, *f_ct_surfacehopping=NULL, *f_ct_nonadiab_coupling=NULL, *f_ct_state_vectors=NULL, *f_ct_orbital=NULL, *f_ct_current=NULL,
         *f_ct_project_wf=NULL,*f_ct_project_wf_ref=NULL, *f_ct_overlap_fo=NULL,*f_ct_overlap_fo_ref=NULL, *f_ct_ortho_ev=NULL, *f_ct_spec=NULL, *f_ct_spec_evec=NULL, *f_ct_fo_ref=NULL, 
         *f_ct_tda=NULL, *f_tb_wfprops=NULL, *f_ct_startwf=NULL, *f_ct_active=NULL; 

    double shiftE;
    int j,k,l,n;
#ifdef GMX_MPI
    MPI_Comm ct_mpi_comm;
    MPI_Status ct_mpi_status;
    int *ct_mpi_gathered, nprocessed;
    if (DOMAINDECOMP(cr)){
      printf("CT-code has to be tested with Domain decomposition\n");
      //printf("CT-code works only with particle decomposition (mdrun -pd)\n");
      exit(-1);
    }      
#endif
    int ct_mpi_rank=0, ct_mpi_size=0, return_value;
    double ct_value, ct_energy, ct_energy1, ct_energy2, fraction_being_annihilated;
    dvec bond, std, msd;
    rvec *x_ct;
    int  counter, counter2, iao, jao;
    t_atoms *ct_atoms, *ct_atoms2, ct_atoms3;
    double tda, original_delta_t;
    static struct timespec time_1, time_end;
    /* END TOMAS KUBAR */

    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
    bIonize  = (Flags & MD_IONIZE);
    bFFscan  = (Flags & MD_FFSCAN);
    bAppend  = (Flags & MD_APPENDFILES);
    if (Flags & MD_RESETCOUNTERSHALFWAY)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf, top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */
    bIterativeCase = ((IR_NPH_TROTTER(ir) || IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    gmx_iterate_init(&iterate, FALSE); /* The default value of iterate->bIterationActive is set to
                                          false in this step.  The correct value, true or false,
                                          is set at each step, as it depends on the frequency of temperature
                                          and pressure control.*/
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || IR_NPH_TROTTER(ir) || IR_NVT_TROTTER(ir)));

    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr, fplog, ir, top_global);

    nstglobalcomm   = check_nstglobalcomm(fplog, cr, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    if (bRerunMD || bFFscan)
    {
        ir->nstxtcout = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog, cr, ir, oenv, &t, &t0, state_global->lambda,
            &(state_global->fep_state), lam0,
            nrnb, top_global, &upd,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, state_global, Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f, top_global->natoms);
    }

    /* lambda Monte carlo random number generator  */
    if (ir->bExpanded)
    {
        mcrng = gmx_rng_init(ir->expandedvals->lmc_seed);
    }
    /* copy the state into df_history */
    copy_df_history(&df_history, &state_global->dfhist);

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind);
    /* needed for iteration of constraints */
    snew(ekind_save, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind_save);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, n_flexible_constraints(constr),
                                 (ir->bContinuation ||
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);

    if (DEFORM(*ir))
    {
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }

    {
        double io = compute_io(ir, top_global->natoms, groups, mdebin->ebin->nener, 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        top = dd_init_local_top(top_global);

        snew(state, 1);
        dd_init_local_state(cr->dd, state_global, state);

        if (DDMASTER(cr->dd) && ir->nstfout)
        {
            snew(f_global, state_global->natoms);
        }
    }
    else
    {
        if (PAR(cr))
        {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog, top_global, ir, cr);

            pd_cg_range(cr, &fr->cg0, &fr->hcg);
            pd_at_range(cr, &a0, &a1);
        }
        else
        {
            top = gmx_mtop_generate_local_top(top_global, ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        forcerec_set_excl_load(fr, top, cr);

        state    = partdec_init_local_state(cr, state_global);
        f_global = f;

        atoms2md(top_global, ir, 0, NULL, a0, a1-a0, mdatoms);

        if (vsite)
        {
            set_vsite_top(vsite, top, mdatoms, cr);
        }

        if (ir->ePBC != epbcNONE && !fr->bMolPBC)
        {
            graph = mk_graph(fplog, &(top->idef), 0, top_global->natoms, FALSE, FALSE);
        }

        if (shellfc)
        {
            make_local_shells(cr, mdatoms, shellfc);
        }

        init_bonded_thread_force_reduction(fr, &top->idef);

        if (ir->pull && PAR(cr))
        {
            dd_make_local_pull_groups(NULL, ir->pull, mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            vsite, shellfc, constr,
                            nrnb, wcycle, FALSE);

    }

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    if (opt2bSet("-cpi", nfile, fnm))
    {
        bStateFromCP = gmx_fexist_master(opt2fn_master("-cpi", nfile, fnm, cr), cr);
    }
    else
    {
        bStateFromCP = FALSE;
    }

    if (MASTER(cr))
    {
        if (bStateFromCP)
        {
            /* Update mdebin with energy history if appending to output files */
            if (Flags & MD_APPENDFILES)
            {
                restore_energyhistory_from_state(mdebin, &state_global->enerhist);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                done_energyhistory(&state_global->enerhist);
                init_energyhistory(&state_global->enerhist);
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(&state_global->enerhist, mdebin);
    }

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG))
    {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd, state);
    }

    if (state->flags & (1<<estMC_RNG))
    {
        set_mc_state(mcrng, state);
    }

    /* Initialize constraints */
    if (constr)
    {
        if (!DOMAINDECOMP(cr))
        {
            set_constraints(constr, top, ir, mdatoms, cr);
        }
    }

    /* Check whether we have to GCT stuff */
    bTCR = ftp2bSet(efGCT, nfile, fnm);
    if (bTCR)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex, gnx);
        for (i = 0; (i < gnx); i++)
        {
            grpindex[i] = i;
        }
    }

    if (repl_ex_nst > 0)
    {
        /* We need to be sure replica exchange can only occur
         * when the energies are current */
        check_nst_param(fplog, cr, "nstcalcenergy", ir->nstcalcenergy,
                        "repl_ex_nst", &repl_ex_nst);
        /* This check needs to happen before inter-simulation
         * signals are initialized, too */
    }
    if (repl_ex_nst > 0 && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, cr->ms, state_global, ir,
                                        repl_ex_nst, repl_ex_nex, repl_ex_seed);
    }

    /* PME tuning is only supported with GPUs or PME nodes and not with rerun */
    if ((Flags & MD_TUNEPME) &&
        EEL_PME(fr->eeltype) &&
        ( (fr->cutoff_scheme == ecutsVERLET && fr->nbv->bUseGPU) || !(cr->duty & DUTY_PME)) &&
        !bRerunMD)
    {
        pme_loadbal_init(&pme_loadbal, ir, state->box, fr->ic, fr->pmedata);
        cycles_pmes = 0;
        if (cr->duty & DUTY_PME)
        {
            /* Start tuning right away, as we can't measure the load */
            bPMETuneRunning = TRUE;
        }
        else
        {
            /* Separate PME nodes, we can measure the PP/PME load balance */
            bPMETuneTry = TRUE;
        }
    }

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for (i = mdatoms->start; i < mdatoms->start+mdatoms->homenr; i++)
            {
                for (m = 0; m < DIM; m++)
                {
                    if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                    {
                        state->v[i][m] = 0;
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms, state, f,
                               graph, cr, nrnb, fr, top, shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, NULL,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, graph, cr, state->box);
        }
    }

    debug_gmx();

    /* set free energy calculation frequency as the minimum of nstdhdl, nstexpanded, and nstrepl_ex_nst*/
    nstfep = ir->fepvals->nstdhdl;
    if (ir->bExpanded && (nstfep > ir->expandedvals->nstexpanded))
    {
        nstfep = ir->expandedvals->nstexpanded;
    }
    if (repl_ex_nst > 0 && nstfep > repl_ex_nst)
    {
        nstfep = repl_ex_nst;
    }

    /* I'm assuming we need global communication the first time! MRS */
    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | ((ir->comm_mode != ecmNO) ? CGLO_STOPCM : 0)
                  | (bVV ? CGLO_PRESSURE : 0)
                  | (bVV ? CGLO_CONSTRAINT : 0)
                  | (bRerunMD ? CGLO_RERUNMD : 0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                    NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, NULL, FALSE, state->box,
                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld, cglo_flags);
    if (ir->eI == eiVVAK)
    {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                        NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, NULL, FALSE, state->box,
                        top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                        cglo_flags &~(CGLO_STOPCM | CGLO_PRESSURE));
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT))
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }
    if (ir->eI != eiVV)
    {
        enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                     and there is no previous step */
    }

    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterativeCase)
    {
        bufstate = init_bufstate(state);
    }
    if (bFFscan)
    {
        snew(xcopy, state->natoms);
        snew(vcopy, state->natoms);
        copy_rvecn(state->x, xcopy, 0, state->natoms);
        copy_rvecn(state->v, vcopy, 0, state->natoms);
        copy_mat(state->box, boxcopy);
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr, FALSE));
        }
        if (EI_STATE_VELOCITY(ir->eI))
        {
            fprintf(fplog, "Initial temperature: %g K\n", enerd->term[F_TEMP]);
        }
        if (bRerunMD)
        {
            fprintf(stderr, "starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name), opt2fn("-rerun", nfile, fnm));
            if (bVerbose)
            {
                fprintf(stderr, "Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
        {
            char tbuf[20];
            fprintf(stderr, "starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->nsteps >= 0)
            {
                sprintf(tbuf, "%8.1f", (ir->init_step+ir->nsteps)*ir->delta_t);
            }
            else
            {
                sprintf(tbuf, "%s", "infinite");
            }
            if (ir->init_step > 0)
            {
                fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps, sbuf), tbuf,
                        gmx_step_str(ir->init_step, sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr, "%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps, sbuf), tbuf);
            }
        }
        fprintf(fplog, "\n");
    }

    /* Set and write start time */
    runtime_start(runtime);
    print_date_and_time(fplog, cr->nodeid, "Started mdrun", runtime);
    wallcycle_start(wcycle, ewcRUN);
    if (fplog)
    {
        fprintf(fplog, "\n");
    }

    /* TOMAS KUBAR
     * initialize the charge transfer calculation
     */
#ifdef GMX_MPI
    ct_mpi_comm = MPI_COMM_WORLD;
    ct_mpi_size = cr->nnodes;
    ct_mpi_rank = cr->nodeid;
    if (ct_mpi_size>1 && ct_mpi_rank==0) snew(ct_mpi_gathered, ct_mpi_size);
    printf("MPI Initialization on node %d of %d done\n", ct_mpi_rank+1, ct_mpi_size);
#endif
// A.HECK //
    snew(ct_atoms,1);
#ifdef GMX_MPI
  if (ct_mpi_rank == 0) {
      *ct_atoms = gmx_mtop_global_atoms(top_global); //expands the topology. now arrays in struct t_atoms contain all atoms. in top_global they are grouped together in moltypes. only possible on rank0. relevant data has to be sent to other nodes.
      snew(ct_atoms->pdbinfo, ct_atoms->nr); //pdb files can be handy for special output where e.g. charge is written in bfactor which can be visualized with VMD 
    if(PROTEIN){ 
      printf("start protein preprocessing\n");
      //snew(ct_atoms2,1);
      //init_t_atoms(ct_atoms2, ct_atoms->nr, FALSE);
      //ct_atoms = protein_preprocessor(ct_atoms, state_global); //changes residues and atomnames for proteins
      printf("protein preprocessing done\n");
      for (j=0; j<ct_atoms->nr; j++){
        //printf("CHECK %s %s\n",*(ct_atoms->atomname[j]), *(ct_atoms->atomname[j]));
        //printf("CHECK %d %d\n",ct_atoms->atom[j].resind, ct_atoms->atom[j].resind );
       // printf("CHECK %s %s\n",*(ct_atoms->resinfo[ct_atoms->atom[j].resind].name), *(ct_atoms->resinfo[ct_atoms->atom[j].resind].name) );
      }/*
      FILE *f_ct_preprocessor=NULL;
      f_ct_preprocessor = fopen("preprocessor.gro", "w");
      fprintf(f_ct_preprocessor,"gro file with new residues\n %d\n", ct_atoms->nr);
      for (j=0; j<ct_atoms->nr; j++)
      fprintf(f_ct_preprocessor, "%5d%-5.5s%5.5s%5d %lf %lf %lf\n", ct_atoms->atom[j].resind%100000, "NEW", *(ct_atoms->atomname[j]), (j+1)%100000, state_global->x[j][XX], state_global->x[j][YY], state_global->x[j][ZZ]);
      fclose(f_ct_preprocessor);
*/
 //     write_sto_conf("preprocessor.gro", ".gro file to check residues", ct_atoms,  state_global->x, NULL, -1, NULL);
    }
    for (i=1; i < ct_mpi_size; i++){
      MPI_Send(&ct_atoms->nres, 1, MPI_INT, i, 200+i, ct_mpi_comm);
      MPI_Send(&ct_atoms->nr, 1, MPI_INT, i, 300+i, ct_mpi_comm);
      for (j=0; j < ct_atoms->nr; j++){
        MPI_Send(*ct_atoms->atomname[j],  6 , MPI_CHAR, i, (i)*1000000+j, ct_mpi_comm);
      }
      MPI_Send(ct_atoms->atom, ct_atoms->nr * sizeof(ct_atoms->atom[0]), MPI_CHAR, i, 400+i, ct_mpi_comm);
      MPI_Send(ct_atoms->resinfo, ct_atoms->nres * sizeof(ct_atoms->resinfo[0]), MPI_CHAR, i, 500+i, ct_mpi_comm);
    }
  } else {
      MPI_Recv(&ct_atoms->nr, 1, MPI_INT, 0, 300+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
   printf("nr %d\n",ct_atoms->nr);
      snew(ct_atoms->atom, ct_atoms->nr);
      snew(ct_atoms->atomname, ct_atoms->nr);
      for (j=0; j < ct_atoms->nr; j++){
        snew(ct_atoms->atomname[j], 1);
        snew(ct_atoms->atomname[j][0], 6);
      }
      MPI_Recv(&ct_atoms->nres, 1, MPI_INT, 0, 200+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
   //printf("nres %d\n",ct_atoms->nres);
      snew(ct_atoms->resinfo, ct_atoms->nres);

      for (j=0; j < ct_atoms->nr; j++){
        MPI_Recv(*ct_atoms->atomname[j],  6 , MPI_CHAR, 0, (ct_mpi_rank)*1000000+j, ct_mpi_comm, &ct_mpi_status);
   //printf("%5s\n",*(ct_atoms->atomname[j]));
      }
      MPI_Recv(ct_atoms->atom, ct_atoms->nr * sizeof(ct_atoms->atom[0]), MPI_CHAR, 0, 400+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
      MPI_Recv(ct_atoms->resinfo, ct_atoms->nres * sizeof(ct_atoms->resinfo[0]) , MPI_CHAR, 0, 500+ct_mpi_rank , ct_mpi_comm, &ct_mpi_status);
  }

#else
  *ct_atoms = gmx_mtop_global_atoms(top_global);
  if(PROTEIN){
    *ct_atoms = gmx_mtop_global_atoms(top_global); 
    ct_atoms = protein_preprocessor(ct_atoms,state_global); //changes residues for proteins
  }
#endif
// END A.HECK //
    snew(ct, 1);
    snew(dftb, 1);
    snew(x_ct, state_global->natoms);
    ct->rk_timestep = ir->delta_t * PS_TO_AU;
    printf("Time step for quantum mechanics: %f ps = %f a.u.\n", ct->rk_timestep / PS_TO_AU, ct->rk_timestep);
#ifdef GMX_MPI
    init_charge_transfer(ct_atoms, top_global, mdatoms, ct, slko_path, state, ct_mpi_rank);
    init_dftb(mdatoms, dftb, ct, slko_path, ct_mpi_rank);
#else
    init_charge_transfer(ct_atoms, top_global, mdatoms, ct, slko_path, state);
    init_dftb(mdatoms, dftb, ct, slko_path);
#endif

    if (ct->jobtype == cteBORNOPPENHEIMER)
      original_delta_t = ir->delta_t;
    if (ct->qmmm == 3) /* this means PME should be used */
#ifdef GMX_MPI
        init_dftb_pme(dftb, ct, ir, ct_mpi_rank);
#else
        init_dftb_pme(dftb, ct, ir);
#endif
#ifdef GMX_MPI
    if (ct_mpi_rank == 0) {
#endif

	/* optimize QM zone before starting MD loops */
        if(ct->opt_QMzone){
	  ct_collect_x(cr, state_global);
	  for (i=0; i<top_global->natoms; i++)
	    copy_rvec(state_global->x[i], x_ct[i]);
          search_starting_site(state->box, mdatoms, dftb, ct, x_ct, slko_path, top_global, state_global->x);
        }

        if (ct->jobtype != cteESP) {
          f_tb_hamiltonian = fopen("TB_HAMILTONIAN.xvg", "w");
          f_ct_project_wf = fopen("CT_PROJECT_WF.xvg", "w");
          f_ct_project_wf_ref = fopen("CT_PROJECT_WF_REF.xvg", "w");
          f_ct_ortho_ev = fopen("CT_ORTHO_EV.xvg", "w");
          f_ct_overlap_fo = fopen("CT_OVERLAP_FO.xvg", "w");
          f_ct_overlap_fo_ref = fopen("CT_OVERLAP_FO_REF.xvg", "w");
          if (ct->jobtype==ctePARAMETERS){
            f_ct_spec = fopen("CT_SPEC.xvg", "w");
            f_ct_spec_evec = fopen("CT_SPEC_EVEC.xvg", "w");
          }
          if (ct->jobtype==cteTDA)
            f_ct_tda = fopen("CT_TDA.xvg", "w");

        }
        if (ct->qmmm > 0) {
          f_ct_esp = fopen("CT_ESP.xvg", "w");
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype == cteSURFACEHOPPING || 
              ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC || ct->jobtype==ctePREZHDOSFHOPPING) {
          f_tb_occupation = fopen("TB_OCCUPATION.xvg", "w");
          f_tb_wfprops = fopen("TB_WFPROPS.xvg", "w");
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteSURFACEHOPPING || 
              ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING) {
          f_ct_energy = fopen("CT_ENERGY.xvg", "w");
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteSURFACEHOPPING ||
              ct->jobtype == cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING ||
              ct->jobtype==cteNEGFLORENTZ) {
          f_tb_hamiltonian_hub = fopen("TB_HAMILTONIAN_HUB.xvg", "w");
          f_tb_hubbard = fopen("TB_HUBBARD.xvg", "w");
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype == cteSURFACEHOPPING) {
          snew(ct_diis, 1);
          ct_init_diis(ct, ct_diis);
          snew(ct->q_act, ct->dim);
          snew(ct->q_old, ct->dim);
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || 
              ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNOMOVEMENT || ct->jobtype==ctePREZHDOSFHOPPING) { //we need here also NOMOVEMENT if no correct wave function is specified and one has to take lowest eigenvalue instead
          snew(ct_broyden, 1);
          ct_init_broyden(ct, ct_broyden);
          snew(ct->q_act, ct->dim);
          snew(ct->q_old, ct->dim);
        }
        if (ct->jobtype==cteFERMI) {
          f_ct_adiabatic = fopen("CT_BROYDEN.xvg", "w");
        }
        if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER) {
          f_ct_adiabatic = fopen("CT_ADIABATIC.xvg", "w");
        }
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteFERMI) {
          f_ct_exp_adiab = fopen("CT_EXPANS_IN_ADIAB_STATES.xvg", "w");
        }
        if (ct->jobtype==cteSURFACEHOPPING || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==ctePREZHDOSFHOPPING) {
          f_ct_surfacehopping = fopen("CT_SURFACE_HOPPING.xvg", "w");
        }
        if (ct->jobtype==cteTULLYFEWESTSWITCHES ||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==ctePREZHDOSFHOPPING) {
          f_ct_nonadiab_coupling = fopen("CT_NONADIAB_COUPLING.xvg", "w");
          f_ct_state_vectors = fopen("CT_STATE_VECTORS.xvg", "w");
        }
        if (ct->jobtype==ctePARAMETERS || ct->jobtype==cteNOMOVEMENT) {
          f_ct_orbital = fopen("CT_ORBITAL.xvg", "w");
        }
        if (ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC) {
          f_ct_current = fopen("CT_CURRENT.xvg", "w");
        }
        if (ct->pool_size > ct->sites){
          f_ct_active = fopen("CT_ACTIVE.xvg", "w");
        }
#ifdef GMX_MPI
    }
#endif
    /* END TOMAS KUBAR */

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
    chkpt_ret = fcCheckPointParallel( cr->nodeid,
                                      NULL, 0);
    if (chkpt_ret == 0)
    {
        gmx_fatal( 3, __FILE__, __LINE__, "Checkpoint error on step %d\n", 0 );
    }
#endif

    debug_gmx();
    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        rerun_fr.natoms = 0;
        if (MASTER(cr))
        {
            bNotLastFrame = read_first_frame(oenv, &status,
                                             opt2fn("-rerun", nfile, fnm),
                                             &rerun_fr, TRX_NEED_X | TRX_READ_V);
            if (rerun_fr.natoms != top_global->natoms)
            {
                gmx_fatal(FARGS,
                          "Number of atoms in trajectory (%d) does not match the "
                          "run input file (%d)\n",
                          rerun_fr.natoms, top_global->natoms);
            }
            if (ir->ePBC != epbcNONE)
            {
                if (!rerun_fr.bBox)
                {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %d time %f does not contain a box, while pbc is used", rerun_fr.step, rerun_fr.time);
                }
                if (max_cutoff2(ir->ePBC, rerun_fr.box) < sqr(fr->rlistlong))
                {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %d time %f has too small box dimensions", rerun_fr.step, rerun_fr.time);
                }
            }
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr, &rerun_fr, &bNotLastFrame);
        }

        if (ir->ePBC != epbcNONE)
        {
            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box, fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX    = !bStateFromCP;
    bInitStep        = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep        = FALSE;
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;

    init_global_signals(&gs, cr, ir, repl_ex_nst);

    step     = ir->init_step;
    step_rel = 0;

    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh, bGStatEveryStep, step);
    }

    if (MULTISIM(cr) && (repl_ex_nst <= 0 ))
    {
        /* check how many steps are left in other sims */
        multisim_nsteps = get_multisim_nsteps(cr, ir->nsteps);
    }


    /* and stop now if we should */
    bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps) ||
                 ((multisim_nsteps >= 0) && (step_rel >= multisim_nsteps )));
    while (!bLastStep || (bRerunMD && bNotLastFrame))
    {

        wallcycle_start(wcycle, ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        if (bRerunMD)
        {
            if (rerun_fr.bStep)
            {
                step     = rerun_fr.step;
                step_rel = step - ir->init_step;
		    }
		    if (rerun_fr.bTime)
		    {
			t = rerun_fr.time;
		    }
		    else
		    {
			t = step;
		    }
		}
		else
		{
		    bLastStep = (step_rel == ir->nsteps);
		    //t         = t0 + step*ir->delta_t; ALEX
                    //if (ct->jobtype == cteBORNOPPENHEIMER){
                    //  t += ir->delta_t;
		    //}else{
                      t         = t0 + step*ir->delta_t; 
                    //}
                    // END ALEX
		}

		if (ir->efep != efepNO || ir->bSimTemp)
		{
		    /* find and set the current lambdas.  If rerunning, we either read in a state, or a lambda value,
		       requiring different logic. */

		    set_current_lambdas(step, ir->fepvals, bRerunMD, &rerun_fr, state_global, state, lam0);
		    bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
		    bDoFEP       = (do_per_step(step, nstfep) && (ir->efep != efepNO));
		    bDoExpanded  = (do_per_step(step, ir->expandedvals->nstexpanded) && (ir->bExpanded) && (step > 0));
		}

		if (bSimAnn)
		{
		    update_annealing_target_temp(&(ir->opts), t);
		}

		if (bRerunMD)
		{
		    if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
		    {
			for (i = 0; i < state_global->natoms; i++)
			{
			    copy_rvec(rerun_fr.x[i], state_global->x[i]);
			}
			if (rerun_fr.bV)
			{
			    for (i = 0; i < state_global->natoms; i++)
			    {
				copy_rvec(rerun_fr.v[i], state_global->v[i]);
			    }
			}
			else
			{
			    for (i = 0; i < state_global->natoms; i++)
			    {
				clear_rvec(state_global->v[i]);
			    }
			    if (bRerunWarnNoV)
			    {
				fprintf(stderr, "\nWARNING: Some frames do not contain velocities.\n"
					"         Ekin, temperature and pressure are incorrect,\n"
					"         the virial will be incorrect when constraints are present.\n"
					"\n");
				bRerunWarnNoV = FALSE;
			    }
			}
		    }
		    copy_mat(rerun_fr.box, state_global->box);
		    copy_mat(state_global->box, state->box);

		    if (vsite && (Flags & MD_RERUN_VSITE))
		    {
			if (DOMAINDECOMP(cr))
			{
			    gmx_fatal(FARGS, "Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
			}
			if (graph)
			{
			    /* Following is necessary because the graph may get out of sync
			     * with the coordinates if we only have every N'th coordinate set
			     */
			    mk_mshift(fplog, graph, fr->ePBC, state->box, state->x);
			    shift_self(graph, state->box, state->x);
			}
			construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
					 top->idef.iparams, top->idef.il,
					 fr->ePBC, fr->bMolPBC, graph, cr, state->box);
			if (graph)
			{
			    unshift_self(graph, state->box, state->x);
			}
		    }
		}

		/* Stop Center of Mass motion */
		bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

		/* Copy back starting coordinates in case we're doing a forcefield scan */
		if (bFFscan)
		{
		    for (ii = 0; (ii < state->natoms); ii++)
		    {
			copy_rvec(xcopy[ii], state->x[ii]);
			copy_rvec(vcopy[ii], state->v[ii]);
		    }
		    copy_mat(boxcopy, state->box);
		}

		if (bRerunMD)
		{
		    /* for rerun MD always do Neighbour Searching */
		    bNS      = (bFirstStep || ir->nstlist != 0);
		    bNStList = bNS;
		}
		else
		{
		    /* Determine whether or not to do Neighbour Searching and LR */
		    bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

		    bNS = (bFirstStep || bExchanged || bNStList || bDoFEP ||
			   (ir->nstlist == -1 && nlh.nabnsb > 0));

		    if (bNS && ir->nstlist == -1)
		    {
			set_nlistheuristics(&nlh, bFirstStep || bExchanged || bDoFEP, step);
		    }
		}

		/* check whether we should stop because another simulation has
		   stopped. */
		if (MULTISIM(cr))
		{
		    if ( (multisim_nsteps >= 0) &&  (step_rel >= multisim_nsteps)  &&
			 (multisim_nsteps != ir->nsteps) )
		    {
			if (bNS)
			{
			    if (MASTER(cr))
			    {
				fprintf(stderr,
					"Stopping simulation %d because another one has finished\n",
					cr->ms->sim);
			    }
			    bLastStep         = TRUE;
			    gs.sig[eglsCHKPT] = 1;
			}
		    }
		}

		/* < 0 means stop at next step, > 0 means stop at next NS step */
		if ( (gs.set[eglsSTOPCOND] < 0 ) ||
		     ( (gs.set[eglsSTOPCOND] > 0 ) && ( bNS || ir->nstlist == 0)) )
		{
		    bLastStep = TRUE;
		}

		/* Determine whether or not to update the Born radii if doing GB */
		bBornRadii = bFirstStep;
		if (ir->implicit_solvent && (step % ir->nstgbradii == 0))
		{
		    bBornRadii = TRUE;
		}

		do_log     = do_per_step(step, ir->nstlog) || bFirstStep || bLastStep;
		do_verbose = bVerbose &&
		    (step % stepout == 0 || bFirstStep || bLastStep);

		if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
		{
		    if (bRerunMD)
		    {
			bMasterState = TRUE;
		    }
		    else
		    {
			bMasterState = FALSE;
			/* Correct the new box if it is too skewed */
			if (DYNAMIC_BOX(*ir))
			{
			    if (correct_box(fplog, step, state->box, graph))
			    {
				bMasterState = TRUE;
			    }
			}
			if (DOMAINDECOMP(cr) && bMasterState)
			{
			    dd_collect_state(cr->dd, state, state_global);
			}
		    }

		    if (DOMAINDECOMP(cr))
		    {
			/* Repartition the domain decomposition */
			wallcycle_start(wcycle, ewcDOMDEC);
			dd_partition_system(fplog, step, cr,
					    bMasterState, nstglobalcomm,
					    state_global, top_global, ir,
					    state, &f, mdatoms, top, fr,
					    vsite, shellfc, constr,
					    nrnb, wcycle,
					    do_verbose && !bPMETuneRunning);
			wallcycle_stop(wcycle, ewcDOMDEC);
			/* If using an iterative integrator, reallocate space to match the decomposition */
		    }
		}

		if (MASTER(cr) && do_log && !bFFscan)
		{
		    print_ebin_header(fplog, step, t, state->lambda[efptFEP]); /* can we improve the information printed here? */
		}

		if (ir->efep != efepNO)
		{
		    update_mdatoms(mdatoms, state->lambda[efptMASS]);
		}

		if ((bRerunMD && rerun_fr.bV) || bExchanged)
		{

		    /* We need the kinetic energy at minus the half step for determining
		     * the full step kinetic energy and possibly for T-coupling.*/
		    /* This may not be quite working correctly yet . . . . */
		    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
				    wcycle, enerd, NULL, NULL, NULL, NULL, mu_tot,
				    constr, NULL, FALSE, state->box,
				    top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
				    CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
		}
		clear_mat(force_vir);

		/* Ionize the atoms if necessary */
		if (bIonize)
		{
		    ionize(fplog, oenv, mdatoms, top_global, t, ir, state->x, state->v,
			   mdatoms->start, mdatoms->start+mdatoms->homenr, state->box, cr);
		}

		/* Update force field in ffscan program */
		if (bFFscan)
		{
		    if (update_forcefield(fplog,
					  nfile, fnm, fr,
					  mdatoms->nr, state->x, state->box))
		    {
			gmx_finalize_par();

			exit(0);
		    }
		}

		GMX_MPE_LOG(ev_timestep2);

		/* TOMAS KUBAR
		 * CHARGE TRANSFER
		 */
                clock_gettime(CLOCK_MONOTONIC, &time_1); //is here to measure the time of al QM stuff in one step

		/* do it all if we have a real simulation, or we want to calculate parameters in this step */
		/* now with NOMOVEMENT just the propagation is stopped.whith LIQM internal relaxation is still there every ct->interval step*/
		if ((ct->jobtype != ctePARAMETERS && ct->jobtype != cteNOMOVEMENT && ct->jobtype != cteESP && ct->jobtype != cteTDA) || step % ct->interval == 0)
		//if ((ct->jobtype != ctePARAMETERS && ct->jobtype != cteESP && ct->jobtype != cteTDA) || step % ct->interval == 0)
		{

		    /* RESTORE THE ORIGINAL CHARGES ON THE NUCLEOBASES - NO MAPPING OF THE HOLE
		     * ... IN ORDER TO AVOID DOUBLE COUNTING OF THE ELECTROSTATIC INTERACTION
		     * THIS WILL BE TAKEN INTO ACCOUNT WITHIN THE TIGHT-BINDING HAMILTONIAN!
		     */
		    for (i=0; i<ct->sites; i++)
			for (m=0; m<ct->site[i].atoms; m++)
			    mdatoms->chargeA[ct->site[i].atom[m]] = ct_atoms->atom[ct->site[i].atom[m]].q;    //A.HECK

		    /* collect the coordinates from all MPI threads */
		    ct_collect_x(cr, state_global);

		    /* collect the coordinates from all MPI threads */
		    if (!PAR(cr) || MASTER(cr))
		    {
			for (i=0; i<top_global->natoms; i++)
			    copy_rvec(state_global->x[i], x_ct[i]);
		    }
	#ifdef GMX_MPI
		    if (MASTER(cr))
		    {

                    do_pbc_mtop(NULL, ir->ePBC, state->box, top_global, x_ct); // make molecules whole -> no dipole , but now we no longer build pbc for every site // edit: for MD we will use site energies obtained in phase1 for CT-Hamiltonian. in phase 2 different sites have different environments (QM vs MM) and therefore enrgies (very bad!). //edit2: now we need again whole molecules to calculate correct COM distance for adaptive QM zones.
			for (i=1; i<ct_mpi_size; i++)
			{
			    MPI_Send(x_ct, 3 * top_global->natoms, GMX_MPI_REAL, i, 300+i, ct_mpi_comm);
	//                    printf("Data of size %d sent to rank %d\n", 3 * top_global->natoms, i);
			}
		    }
		    else
		    {
			return_value = MPI_Recv(x_ct, 3 * top_global->natoms, GMX_MPI_REAL, 0, 300 + ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
			if (return_value)
			    printf("Error receiving data on rank %d\n", ct_mpi_rank);
	//                else
	//                    printf("Data of size %d received on rank %d\n", 3 * top_global->natoms, ct_mpi_rank);
		    }
		    //for (i=0; i<10; i++)
		    //{
		    //    printf("rank%d atom%d %12.7f%12.7f%12.7f\n", ct_mpi_rank, i, x_ct[i][XX]*10, x_ct[i][YY]*10, x_ct[i][ZZ]*10);
		    //}
	#endif



   printf("start prepare at %f\n", (double) clock()/CLOCKS_PER_SEC);
	  prepare_charge_transfer(state->box, mdatoms, dftb, ct, x_ct);
   printf("stop prepare at %f\n", (double) clock()/CLOCKS_PER_SEC);

	#ifdef GMX_MPI
	   if (ct->qmmm == 3) {
	     if (step - dftb->lastlist_pme >= dftb->nstlist_pme) {
	       do_neighborlist_for_dftb(ct, dftb, state->x, ct_mpi_rank, ct_mpi_size);
	       dftb->lastlist_pme = step;
	     }
	     do_pme_for_dftb_part1(ct, dftb, ct_mpi_rank, ct_mpi_size);
   printf("stop pme-preparation at %f\n", (double) clock()/CLOCKS_PER_SEC);

	   }
	   if (ct->jobtype != cteESP) {
	     do_dftb_phase1(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
	   } else {
	     do_esp_only(ct, dftb, mdatoms->chargeA, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
	   }

           // broadcast the results if all calculations have finished. broadcaster is rank "i % ct_mpi_size" (there was phase1 performed) 
           MPI_Barrier(ct_mpi_comm);
           for (i=0; i<ct->sites; i++){
           //printf("Bcast at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);
           MPI_Bcast(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb),          MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);//needed to build the coarse grained hamilton matrix at main node
           MPI_Bcast(dftb->phase1[i].qmat, dftb->phase1[i].nn,                 MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);// "
           MPI_Bcast(dftb->phase1[i].ev, dftb->phase1[i].norb,                 MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);// "
           //printf("after Bcast at rank %d at %f\n",i % ct_mpi_size, (double) clock()/CLOCKS_PER_SEC);
           }
           MPI_Barrier(ct_mpi_comm); //wait until broadcasting has finished



	   if (ct_mpi_rank == 0 && ct->jobtype != cteESP){
             check_and_invert_orbital_phase(dftb->phase1, ct); //calc sign on master since we need later also the obtained overlap in project_wf_on_new_basis().
             do_dftb_phase2(ct, dftb);
           }
	#else
	   if (ct->qmmm == 3) {
	     if (step - dftb->lastlist_pme >= dftb->nstlist_pme) {
	       do_neighborlist_for_dftb(ct, dftb, state->x);
	       dftb->lastlist_pme = step;
	     }
	     do_pme_for_dftb_part1(ct, dftb);
	   }
	   if (ct->jobtype != cteESP) {
	     do_dftb_phase1(ct, dftb);
	     // check_and_invert_orbital_phase(dftb->phase1, ct->sites, ct->homo, f_ct_orbital, step); absorbed in sort_mobasis
	     /* copying of charges and PME for 2nd phase comes here !!! */
	     do_dftb_phase2(ct, dftb);
	   } else {
	     do_esp_only(ct, dftb, mdatoms->chargeA);
	   }
	#endif


///////////////////////////////////////////
// write interesting data in pbd file as bfactor in order to combine structre and property.
/*
for (i=0; i<ct->sites; i++)
for (j=0; j<ct->site[i].atoms; j++){
 ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac=0;
for (k=0; k<ct->site[i].homos; k++){
 //partial charges
 //ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac =(real) (dftb->phase1[i].qmat[j] - dftb->qzero1[dftb->phase1[i].izp[j]]);
 //delta charges
 ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac +=(real) ct->site[i].delta_q[k][j];
}
}
write_sto_conf("CT.pdb", "written by charge transfer code", ct_atoms, x_ct, NULL, 0, NULL);
//*/


//write_out_MOs(step, x_ct, ct_atoms, dftb, ct);
//////////////////////////////////////////

	#ifdef GMX_MPI
	   if (ct_mpi_rank == 0) {
	#endif

	   if (ct->jobtype != cteESP) {
	     ct_assemble_hamiltonian(ct, dftb); //now including averaging
/*	     for (i=0; i<ct->dim; i++) {
               ct->hamiltonian[i][i] = (ct->is_hole_transfer==1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
	       for (j=i+1; j<ct->dim; j++) {
                 //ct->hamiltonian[i][j] = ct->hamiltonian[j][i] = dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling; //  1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014)
                 ct->hamiltonian[i][j] = ct->hamiltonian[j][i] = (ct->is_hole_transfer==1) ? -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling:  dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling; //  1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014) // this is the right way to do it. without changing off-diagonal elements the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
               }
             }
*/



	     fprintf(f_tb_hamiltonian, "%10d", step);
	     for (i=0; i<ct->dim; i++) {
	       for (m=i; m<ct->dim; m++) fprintf(f_tb_hamiltonian, " %10.6f", ct->hamiltonian[i][m] * HARTREE_TO_EV);
	     }
	     fprintf(f_tb_hamiltonian, "\n");
	   }


	//overlap /////////////////////////////////////////////////////// 

	  // write overlap matrix for first site
	  i=0;
	  fprintf(f_ct_overlap_fo, "%10d ", step);
	  for (j = 0; j < ct->site[i].homos; j++)
	  for (k = j; k < ct->site[i].homos; k++)
      fprintf(f_ct_overlap_fo, " %10.6f ",ct->site[i].overlap[j][k]);
  fprintf(f_ct_overlap_fo, "  \n");
  i=0;
  fprintf(f_ct_overlap_fo_ref, "%10d ", step);
  for (j = 0; j < ct->site[i].homos; j++)
  for (k = j; k < ct->site[i].homos; k++)
      fprintf(f_ct_overlap_fo_ref, " %10.6f ",ct->site[i].overlap_ref[j][k]);
  fprintf(f_ct_overlap_fo_ref, "  \n");

//spectrum////////////////////
  if (ct->jobtype==ctePARAMETERS){ //&& ct->first_step == 1){
    get_spectrum(ct,dftb);
   //fprintf(f_ct_spec, "energy eigen values of FO Hamiltonian:\n");
   fprintf(f_ct_spec, "%10d", step);
   fprintf(f_ct_spec_evec, "%10d", step);

   for (i=0; i<ct->dim; i++)
     fprintf(f_ct_spec, " %10.6f ", ct->ev_spec[i] * HARTREE_TO_EV);
   fprintf(f_ct_spec, "\n");

   for (i=0; i<ct->dim; i++)
     for (j=0; j<ct->dim; j++)
     fprintf(f_ct_spec_evec, " %10.6f ", ct->evec_spec[j+i*ct->dim] ) ; //prints all evecs in one line
   fprintf(f_ct_spec_evec, "\n");
  }
/////////////////////////////////////////////////////////////


  // CONSTRUCT THE HUBBARD MATRIX //            now here because hubbard matrix depends on delta_q calculated by get_delta_q.

  if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteSURFACEHOPPING ||
       ct->jobtype == cteFERMI || ct->jobtype == cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype == cteFERMISFHOPPING || ct->jobtype == cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype == ctePERSICOSFHOPPING) {
    // diag: Hubbard params of nucleobases (constant - do not calculate)
/*
    // offdiag: 1/r_ij for nucleobases (original code)
    counter=-1;
    for (i=0; i<ct->sites; i++)
    for (k=0; k<ct->site[i].homos; k++){
      counter2=-1;
      counter++;
      for (j=0; j<ct->sites; j++) {
      dvec_sub(dftb->phase1[i].com, dftb->phase1[j].com, bond);
      for (l=0; l<ct->site[j].homos; l++){
        counter2++;
        if(j==i) continue;
        // approximation of nucleobases as single charged points
        ct->hubbard[counter][counter2] = ct->hubbard[counter2][counter] = ct->sic / dnorm(bond);
      }
      }
    } 
*/
    // offdiag: sum_AB[ (dq_iA * dq_jB)/r_AB  (improved code)
    counter=-1;
    for (i=0; i<ct->sites; i++)
    for (k=0; k<ct->site[i].homos; k++){
      counter2=-1;
      counter++;
      for (j=0; j<ct->sites; j++)
      for (l=0; l<ct->site[j].homos; l++){
          counter2++;
          if(j==i) continue;
          ct->hubbard[counter][counter2]=0;
          for (m=0; m < ct->site[i].atoms; m++)
          for (n=0; n < ct->site[j].atoms; n++){
            dvec_sub(dftb->phase1[i].x[m], dftb->phase1[j].x[n], bond);
            ct->hubbard[counter][counter2] += ct->sic * ct->site[i].delta_q[k][m] * ct->site[j].delta_q[l][n] / dnorm(bond);
            
          }
          if (ct->do_epol == 1)
            ct->hubbard[counter][counter2]= (1.0/EPSILON_OP)*ct->hubbard[counter][counter2];
      }
    } 
/*
//test: calculate diag elements on the fly    THIS IS WRONG (also overwrites fixed hubbard diagonals eg. lambda_i)
    counter=0;
    for (i=0; i<ct->sites; i++) 
    for (j=0; j<ct->site[i].homos; j++){
      ct->hubbard[counter][counter] = 0.0;
      for (k=0; k<ct->site[i].atoms; k++)
      for (l=k; l<ct->site[i].atoms; l++)
        ct->hubbard[counter][counter] += ct->site[i].delta_q[j][k]*ct->site[i].delta_q[j][l]* (k>l ? dftb->phase1[i].gammamat[k][l] : dftb->phase1[i].gammamat[l][k]);
      counter++;
    }
//end test */
  } else {
    for (i=0; i<ct->dim; i++)
      for (j=0; j<ct->dim; j++)
        ct->hubbard[i][j] = 0.e0;
  }

 /*
  printf("\n HUBBARD MATRIX \n");
  for (i=0; i< ct->dim; i++){
    for (j=0; j< ct->dim; j++)
      printf("%f ", ct->hubbard[i][j]);
    printf("\n");
  }
 // */

   if (ct->sic > 0.0 && (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteSURFACEHOPPING ||
         ct-> jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING ||
         ct->jobtype==cteNEGFLORENTZ) ) {
     fprintf(f_tb_hamiltonian_hub, "%10d", step);
     fprintf(f_tb_hubbard, "%10d", step);
     for (i=0; i<ct->dim; i++) {
       ct->hamiltonian_mod[i] = 0.0;
       for (m=0; m<ct->dim; m++)
         ct->hamiltonian_mod[i] += ct->hubbard[i][m] * ct->occupation[m];
       fprintf(f_tb_hamiltonian_hub, " %10.6f", ct->hamiltonian_mod[i] * HARTREE_TO_EV);
       for (m=i; m<ct->dim; m++)
         fprintf(f_tb_hubbard, " %10.6f", ct->hubbard[i][m]);
     }
     fprintf(f_tb_hamiltonian_hub, "\n");
     fprintf(f_tb_hubbard, "\n");
   }

   if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER) {
     fprintf(f_ct_adiabatic, "%10d ", step);
   }

   printf("start dynamic at %f\n", (double) clock()/CLOCKS_PER_SEC);
   if(step>=ct->n_avg_ham-1) //start propagation only after the running average is stable. steps starts from zero and we want to perform propagation in the first step if there is no averaging ct->n_avg=1
   switch (ct->jobtype) {
     case cteSCCDYNAMIC:
     case cteNONSCCDYNAMIC:
 
       /* evaluate the amount of annihilated charge due to neg. imag. potentials */
       for (i=0; i<ct->neg_imag_pot; i++)
         ct->site_annihilated_occupation[i] += 2 * ct->neg_imag_pot_rate_constant * PS_TO_AU * ct->occupation[ct->site_neg_imag_pot[i]] * ir->delta_t;
       /* additional option: if no wf is specified, take lowest eigenvector as starting function */
       ct->survival=0;
       for (i=0; i<ct->dim; i++) 
         ct->survival += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
       if (ct->survival < 0.99 && ct->first_step ){ //should be normalized, 0.01 tolerance
         printf("WARNING: no correct starting wave function was specified. Lowest eigenvector will be taken instead.\n");
         f_ct_startwf = fopen ("CT_STARTWF.xvg", "w" );
         ct->fermi_kt = BOLTZMANN_HARTREE_KELVIN * 300; // T=300K helps converging 
         if(do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df , f_ct_startwf) == 1){
           printf("ERROR: calculation of lowest eigenvalue did not converge.\n");
           exit(-1);
         }
         fclose(f_ct_startwf);
       }
       /* print out the wave function */
       if ((ir->nstxout > 0 && step % ir->nstxout == 0) || (ir->nstxtcout > 0 && step % ir->nstxtcout == 0)) {
         printf("Wave function before step %d:\n", step);
         for (i=0; i<ct->dim; i++)
           printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
       }
       /* SCC dynamics -- perform the integration with RKsuite
        * (the neg-imag-potential capability is included in do_rksuite() )
        */

       if(ct->dim > ct->sites){     //if there are degenerated orbitals 
         printf("start project at %f\n", (double) clock()/CLOCKS_PER_SEC);
         project_wf_on_new_basis(step, dftb, ct, f_ct_project_wf, f_ct_project_wf_ref );
         printf("stop project at %f\n", (double) clock()/CLOCKS_PER_SEC);
       }

       do_rksuite(ct);
       /* the coefficients of decomposition into adiabatic states */
       fprintf(f_ct_exp_adiab, "%10d ", step);
       //do_wf_decomp_evec(ct, ct_diis, f_ct_exp_adiab);
       do_wf_decomp_evec_fermi(ct, ct_broyden, ct_broyden->df, f_ct_exp_adiab);
       break;
     case cteADIABATIC:
     case cteADNONSCC:
       /* adiabatic dynamics -- find the wavefunction that minimizes the energy */
       //do_adiabatic(ct, ct_diis, f_ct_adiabatic);
       if (do_adiabatic(ct, ct_diis, f_ct_adiabatic)) {
         fprintf(stderr, "Adiabatic - DIIS not converged in step %d\n", step);
         //fprintf(f_ct_adiabatic, "Adiabatic - Broyden not converged\n");
       }
       break;
     case cteFERMI:
       if (do_adiab_fermi(ct, ct_broyden, ct_broyden->df, f_ct_adiabatic)) {
         fprintf(stderr, "Fermi distribution based Broyden not converged in step %d\n", step);
         fprintf(f_ct_exp_adiab, "%10d BROYDEN NOT CONVERGED\n", step);
       } else {
         fprintf(f_ct_exp_adiab, "%10d ", step);
         for (i=0; i<ct->dim; i++)
           fprintf(f_ct_exp_adiab, "%10.7f", ct->ev_adiab[i]);
         for (i=0; i<ct->dim; i++)
           fprintf(f_ct_exp_adiab, "%9.6f", ct_broyden->df[i]);
         fprintf(f_ct_exp_adiab, "\n");
       }
       break;
     case cteFERMIADIABATIC:
       if (do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df, f_ct_adiabatic)) {
         fprintf(stderr, "Fermi distribution based Broyden not converged in step %d\n", step);
       }
       break;
     case cteSURFACEHOPPING:
       /* evaluate the amount of annihilated charge due to neg. imag. potentials */
       for (i=0; i<ct->neg_imag_pot; i++) {
         fraction_being_annihilated = 2 * ct->neg_imag_pot_rate_constant * PS_TO_AU * ct->occupation[ct->site_neg_imag_pot[i]] * ir->delta_t;
         ct->site_annihilated_occupation[i] += fraction_being_annihilated;
         ct->survival -= fraction_being_annihilated;
       }
       /* (diabatic) surface hopping among adiabatic states */
       fprintf(f_ct_surfacehopping, "%10d ", step);
       do_surface_hopping(ct, ct_diis, f_ct_surfacehopping);
       /* in case there is annihilation, correct the wave function */
       if (ct->neg_imag_pot > 0)
         for (i=0; i<ct->dim; i++)
           ct->wf[i] *= sqrt(ct->survival);
       break;

     case cteBORNOPPENHEIMER:   //JJK
       do_born_oppenheimer(ct, ct_broyden, ct_broyden->df , f_ct_adiabatic); //scale timestep with overlap^2. i.e. take smaller steps if wf is changing strongly
       //ir->delta_t=do_born_oppenheimer(ct, ct_broyden, ct_broyden->df , f_ct_adiabatic); //scale timestep with overlap^2. i.e. take smaller steps if wf is changing strongly
//printf("delta_t %f\n", ir->delta_t);
       //if(ir->delta_t < 0){
       //  fprintf(stderr, "Fermi distribution based Broyden not converged in step %d\n", step);
       //} else {ir->delta_t *=original_delta_t;
//printf("original delta_t %f\n", original_delta_t);
       //}
       
       // maybe problematic since       t         = t0 + step*ir->delta_t;
       //ct->rk_timestep = ir->delta_t * PS_TO_AU;
     break;

     case cteFERMISFHOPPING:
       fprintf(f_ct_surfacehopping, "%10d ", step);
       do_fermi_surface_hopping(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping);
       break;
     case cteTULLYFEWESTSWITCHES:
       fprintf(f_ct_surfacehopping, "%10d ", step);
       fprintf(f_ct_nonadiab_coupling, "%10d ", step);
       do_tully_fewest_switches(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling);
       break;
     case cteTULLYLOC:
       printf(f_ct_surfacehopping, "%10d ", step);   //maybe replace step with time,mabye move into function                        
       fprintf(f_ct_nonadiab_coupling, "%10d ", step);
       do_tully_local(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling);
     case ctePERSICOSFHOPPING:
       fprintf(f_ct_surfacehopping, "%10d ", step);
       fprintf(f_ct_nonadiab_coupling, "%10d ", step);
       fprintf(f_ct_state_vectors, "%10d ", step);
       printf("%10d ", step);
       do_persico_diabatic_sfhopping_new(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling, f_ct_state_vectors); //new corrected version. adapted but not tested yet
       //do_persico_diabatic_sfhopping(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling);
       break;
     case cteNEGFLORENTZ:
     case cteNEGFLORENTZNONSCC:
       negf_propagate(ct);
       fprintf(f_ct_current, "%10d %12.5e %12.5e %12.5e\n", step, ct->negf_arrays->current[0], ct->negf_arrays->current[1], ct->negf_arrays->current[2]);
       break;
     case ctePARAMETERS:
       break;
     case cteNOMOVEMENT:
       if ( ct->first_step ){
         ct->survival=0;
         for (i=0; i<ct->dim; i++)
           ct->survival += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
         if (ct->survival < 0.99 ){ //should be normalized, 0.01 tolerance
           printf("WARNING: no correct starting wave function was specified. Lowest eigenvalue will be taken instead.\n");
           f_ct_startwf = fopen ("CT_STARTWF.xvg", "w");
           ct->fermi_kt = BOLTZMANN_HARTREE_KELVIN * 300; // T=300K helps converging 
           if(do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df , f_ct_startwf) == 1){
             printf("ERROR: calculation of lowest eigenvalue did not converge.\n");
             exit(-1);
           }
           fclose(f_ct_startwf);
         }
       }
       break;
     case cteESP:
       /* nothing to do */
       break;
     case cteTDA:
       tda=calc_tda(ct);
       fprintf(f_ct_tda, "%10d %15.12f \n", step, tda*HARTREE_TO_EV);
       break; 
     case ctePREZHDOSFHOPPING:
       printf("%10d ", step);
       fprintf(f_ct_surfacehopping, "%10d ", step);
       fprintf(f_ct_nonadiab_coupling, "%10d ", step);
       fprintf(f_ct_state_vectors, "%10d ", step);
       do_prezhdo_sfhopping(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling, f_ct_state_vectors);
       break;
   }
   printf("stop dynamic at %f\n", (double) clock()/CLOCKS_PER_SEC);

   if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER) {
     fprintf(f_ct_adiabatic, "\n");
   }

   if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteSURFACEHOPPING ||
         ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC || ct->jobtype==ctePREZHDOSFHOPPING) {
     /* calculate and print out the occupations */
     ct->survival = 0.0;
     fprintf(f_tb_occupation, "%10f", t*1000); //1000=converting ps to fs 
     for (i=0; i<ct->dim; i++) {
       ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
       fprintf(f_tb_occupation, "%12.8f", ct->occupation[i]);
       ct->survival += ct->occupation[i];
     }
     fprintf(f_tb_occupation, "%12.8f", ct->survival);
     if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC  || ct->jobtype==cteSURFACEHOPPING)
       for (i=0; i<ct->neg_imag_pot; i++)
         fprintf(f_tb_occupation, "%12.8f", ct->site_annihilated_occupation[i]);
     fprintf(f_tb_occupation, "\n");
   }

   /* print out the ESP */
   if (ct->qmmm > 0) {
     fprintf(f_ct_esp, "%10d", step);
     for (i=0; i<ct->sites; i++) {
       fprintf(f_ct_esp, " %9.5f", dftb->phase1[i].esp * AU_OF_ESP_TO_VOLT);
       /* For comparison (debug) purposes,
        * calculate the SCC-DFTB external shift
        * (be averaging over the atoms, mass weighted)
        */
       /* Maybe output this instead of ESP! */
//     shiftE = 0.0;
//     for (j=0; j<dftb->phase1[i].nn; j++)
//       shiftE += dftb->phase1[i].mass[j] * dftb->phase1[i].shiftE[j];
//     shiftE *= dftb->phase1[i].inv_tot_mass;
//     fprintf(f_ct_shift, " %9.5f", shiftE * AU_OF_ESP_TO_VOLT);
     }
     fprintf(f_ct_esp, "\n");
   }


   /* print out center, standard deviation (width) and mean square displacement (needed for mobility) of the charge */
   if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteSURFACEHOPPING ||
       ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES ||
       ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC || ct->jobtype==ctePREZHDOSFHOPPING) {
      // COC: <psi|R|psi> 				= sum_i occ_i(t)*R_i(t)  // R is operator that gives COM of site i (R_i)
      // MSD: <psi(t)|(R-R0)^2|psi(t)> 			= sum_i occ_i(t)*[R_i(t)-*R0)]^2    //R0 is position (COC) of the WF at t=0
      // STD: sqrt(<psi|R^2|psi> - (<psi|R|psi>)^2) 		= sqrt(  sum_i [occ_i(t)*R_i^2(t)- (occ_i(t)*R_i(t))^2]  )
     for (i=0; i<ct->dim; i++)
       ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
     counter=0;
     clear_dvec(ct->coc);
     for (i=0; i<ct->sites; i++)
     for (j=0; j<ct->site[i].homos; j++){
       ct->coc[XX]+=ct->occupation[counter]*dftb->phase1[i].com[XX];
       ct->coc[YY]+=ct->occupation[counter]*dftb->phase1[i].com[YY];
       ct->coc[ZZ]+=ct->occupation[counter]*dftb->phase1[i].com[ZZ];
       counter++;
     }
     if(ct->first_step)
       copy_dvec(ct->coc, ct->coc_start);
     
     counter=0;
     clear_dvec(std);
     clear_dvec(msd);
     for (i=0; i<ct->sites; i++)
     for (j=0; j<ct->site[i].homos; j++){
       std[XX]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[XX]) - SQR(ct->occupation[counter]*dftb->phase1[i].com[XX]);
       std[YY]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[YY]) - SQR(ct->occupation[counter]*dftb->phase1[i].com[YY]);
       std[ZZ]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[ZZ]) - SQR(ct->occupation[counter]*dftb->phase1[i].com[ZZ]);

       msd[XX]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[XX]-ct->coc_start[XX]);
       msd[YY]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[YY]-ct->coc_start[YY]);
       msd[ZZ]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[ZZ]-ct->coc_start[ZZ]);
       counter++;
     }
     for (i=0;i<3; i++ )
       std[i]=sqrt(std[i]);
   
     fprintf(f_tb_wfprops, "%d  %f %f %f    %f %f %f    %f %f %f\n",step, ct->coc[XX],ct->coc[YY],ct->coc[ZZ], msd[XX],msd[YY],msd[ZZ] , std[XX],std[YY],std[ZZ]); 
   }

#ifdef GMX_MPI
   } /* if (ct_mpi_rank == 0) */

   if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteSURFACEHOPPING || 
         ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING ||
         ct->jobtype==cteNEGFLORENTZ) {
     /* communicate the occupations to the slaves */
     if (ct_mpi_rank == 0) {
       for (i=1; i<ct_mpi_size; i++)
         MPI_Send(ct->occupation, ct->dim, MPI_DOUBLE, i, 200+i, ct_mpi_comm);
     } else {
       MPI_Recv(ct->occupation, ct->dim, MPI_DOUBLE, 0, 200+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
       ct->survival = 0.0;
       for (i=0; i<ct->dim; i++)
         ct->survival += ct->occupation[i];
     }
   }
#endif

  /* if we have neg. imag. potential and the survival is under the threshold, finish!
  if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC)
    if (ct->survival < ct->survival_threshold) {
#ifdef GMX_MPI
      if (ct_mpi_rank == 0)
#endif
      fprintf(stderr, "Survival of %f dropped under the threshold (%f) at step %d, terminating!\n", 
        ct->survival, ct->survival_threshold, step);
      bGotTermSignal = TRUE;
    }
  */

  } /* initial if */

  if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteSURFACEHOPPING || 
        ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ) {

    #ifdef GMX_MPI
    get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
    MPI_Barrier(ct_mpi_comm);
    for (i=0; i<ct->sites; i++){
      MPI_Bcast(ct->site[i].delta_q[0], ct->site[i].homos*ct->site[i].atoms, MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm); //needed to build the hubbard matrix at main node
      if ( ct->do_lambda_i==2 )
        MPI_Bcast(dftb->phase1[i].grad[0], 3*dftb->phase1[i].nn,  	       MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);//needed to add forces for explicit lambda_i //try to send rvec -> 3*double. have to  test
      if (ct->qmmm == 3) 
        MPI_Bcast(&(dftb->phase1[i].esp), 1,  MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
    }
    MPI_Barrier(ct_mpi_comm);
    #else
    get_MM_params(ct, dftb);
    #endif            

    /* map the hole on the atomic charges - first check the charges */
    counter=0;
    for (i=0; i<ct->sites; i++) {
      // printf("Site %d, charges top_global and mdatoms\n", i);
      for (m=0; m<ct->site[i].atoms; m++) {
        mdatoms->chargeA[ct->site[i].atom[m]] = ct_atoms->atom[ct->site[i].atom[m]].q;
//printf("q %f\n",mdatoms->chargeA[ct->site[i].atom[m]]);
      }
      for (j=0; j<ct->site[i].homos; j++){
        for (m=0; m<ct->site[i].atoms; m++){
          mdatoms->chargeA[ct->site[i].atom[m]] += ct->occupation[counter] * ct->site[i].delta_q[j][m];
        }
        counter++;
      }
      //for (m=0; m<ct->site[i].atoms; m++)
       // printf("rank%d %3d %8.5f %8.5f\n",ct_mpi_rank, ct->site[i].atom[m], ct_atoms.atom[ct->site[i].atom[m]].q, mdatoms->chargeA[ct->site[i].atom[m]]);
    }
  }


  /* print out active QM zone and set new QM zone if charge moved to the border */
  if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC){
    if(ct->pool_size > ct->sites){
      fprintf(f_ct_active, "%d  ", step);
      for (i=0; i<ct->sites; i++){
        fprintf(f_ct_active, " %d", ct->site[i].resnr);
      }
      fprintf(f_ct_active, "\n");

      while(adapt_QMzone(ct,x_ct, mdatoms, top_global, state->box, state_global->x)){} 
    }
  }

   clock_gettime(CLOCK_MONOTONIC, &time_end);
   if (ct_mpi_rank==0)
     print_time_difference("TOTAL QM TIME[ns]:", time_1, time_end);
  ct->first_step = 0; 
  /* END TOMAS KUBAR */

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
        bCPT = (((gs.set[eglsCHKPT] && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step && !bRerunMD);
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
        }

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        if (EI_VV(ir->eI) && (!bInitStep))
        {
            /* for vv, the first half of the integration actually corresponds
               to the previous step.  bCalcEner is only required to be evaluated on the 'next' step,
               but the virial needs to be calculated on both the current step and the 'next' step. Future
               reorganization may be able to get rid of one of the bCalcVir=TRUE steps. */

            bCalcEner = do_per_step(step-1, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && (do_per_step(step, ir->nstpcouple) || do_per_step(step-1, ir->nstpcouple)));
        }
        else
        {
            bCalcEner = do_per_step(step, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
        }

        /* Do we need global communication ? */
        bGStat = (bCalcVir || bCalcEner || bStopCM ||
                  do_per_step(step, nstglobalcomm) ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcVir  = TRUE;
            bCalcEner = TRUE;
            bGStat    = TRUE;
        }

        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                      (bGStat ? CGLO_GSTAT : 0)
                      );

        force_flags = (GMX_FORCE_STATECHANGED |
                       ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       GMX_FORCE_SEPLRF |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
                       );

        if (fr->bTwinRange)
        {
            if (do_per_step(step, ir->nstcalclr))
            {
                force_flags |= GMX_FORCE_DO_LR;
            }
        }

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count = relax_shell_flexcon(fplog, cr, bVerbose, bFFscan ? step+1 : step,
                                        ir, bNS, force_flags,
                                        bStopCM, top, top_global,
                                        constr, enerd, fcd,
                                        state, f, force_vir, mdatoms,
                                        nrnb, wcycle, graph, groups,
                                        shellfc, fr, bBornRadii, t, mu_tot,
                                        state->natoms, &bConverged, vsite,
                                        outf->fp_field);
            tcount += count;

            if (bConverged)
            {
                nconverged++;
            }
        }
        else
        {




            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ir, step, nrnb, wcycle, top, top_global, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, vsite, mu_tot, t, outf->fp_field, ed, bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);


/*********************************************************************************/
          //add QM forces to relax site geometry
/*********************************************************************************/
            if ( ct->do_lambda_i==2 || ct->do_lambda_i==3) 
            //if ( ct->do_lambda_i==2 ) 
              for(i=mdatoms->start;i<mdatoms->start+mdatoms->homenr;i++) //iterate over home atoms
              {
               for(j=0; j<ct->sites; j++)
               for(k=0; k<ct->site[j].atoms; k++)
               if(i==ct->site[j].atom[k]){ //atom is at home on this node
                 for(l=0; l<DIM; l++) {
//printf("atom i %d MM force %lf   QM force %lf \n", i, f[i][l], -(real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l]);
                     f[i][l]      -= (real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l]; 

             //      fshift[i][j] = (real) HARTREE_BOHR2MD * dftb1.grad[i][j]; TODO: what is fshift? this is used in QMMM gromacs code
                 }
               }
              }
       
            if ( ct->do_lambda_i==3 )
              for(i=mdatoms->start;i<mdatoms->start+mdatoms->homenr;i++){
               counter=0; 
               for(j=0; j<ct->sites; j++)
               for(k=0; k<ct->site[j].atoms; k++){
                 if(i==ct->site[j].atom[k]){ //atom is at home on this node
                   for(l=0; l<DIM; l++) {
//printf("atom i %d MM force %lf   QM force %lf \n", i, f[i][l], -(real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l]);
                     f[i][l]      -= (real) HARTREE_BOHR2MD * dftb->phase2.grad[counter][l];
                   }
                 }
                counter++;
               }
              }


 
/*********************************************************************************/

        }

        GMX_BARRIER(cr->mpi_comm_mygroup);

        if (bTCR)
        {
            mu_aver = calc_mu_aver(cr, state->x, mdatoms->chargeA,
                                   mu_tot, &top_global->mols, mdatoms, gnx, grpindex);
        }

        if (bTCR && bFirstStep)
        {
            tcr = init_coupling(fplog, nfile, fnm, cr, fr, mdatoms, &(top->idef));
            fprintf(fplog, "Done init_coupling\n");
            fflush(fplog);
        }

        if (bVV && !bStartingFromCpt && !bRerunMD)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            if (ir->eI == eiVV && bInitStep)
            {
                /* if using velocity verlet with full time step Ekin,
                 * take the first half step only to compute the
                 * virial for the first step. From there,
                 * revert back to the initial coordinates
                 * so that the input is actually the initial step.
                 */
                copy_rvecn(state->v, cbuf, 0, state->natoms); /* should make this better for parallelizing? */
            }
            else
            {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ1);
            }

            /* If we are using twin-range interactions where the long-range component
             * is only evaluated every nstcalclr>1 steps, we should do a special update
             * step to combine the long-range forces on these steps.
             * For nstcalclr=1 this is not done, since the forces would have been added
             * directly to the short-range forces already.
             */
            bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

            update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC,
                          f, bUpdateDoLR, fr->f_twin, fcd,
                          ekind, M, wcycle, upd, bInitStep, etrtVELOCITY1,
                          cr, nrnb, constr, &top->idef);

            if (bIterativeCase && do_per_step(step-1, ir->nstpcouple) && !bInitStep)
            {
                gmx_iterate_init(&iterate, TRUE);
            }
            /* for iterations, we save these vectors, as we will be self-consistently iterating
               the calculations */

            /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */

            /* save the state */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
            }

            bFirstIterate = TRUE;
            while (bFirstIterate || iterate.bIterationActive)
            {
                if (iterate.bIterationActive)
                {
                    copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
                    if (bFirstIterate && bTrotter)
                    {
                        /* The first time through, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */

                        veta_save = state->veta;
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ0);
                        vetanew     = state->veta;
                        state->veta = veta_save;
                    }
                }

                bOK = TRUE;
                if (!bRerunMD || rerun_fr.bV || bForceUpdate)     /* Why is rerun_fr.bV here?  Unclear. */
                {
                    dvdl = 0;

                    update_constraints(fplog, step, &dvdl, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, shake_vir, NULL,
                                       cr, nrnb, wcycle, upd, constr,
                                       bInitStep, TRUE, bCalcVir, vetanew);

                    if (!bOK && !bFFscan)
                    {
                        gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }

                }
                else if (graph)
                {
                    /* Need to unshift here if a do_force has been
                       called in the previous step */
                    unshift_self(graph, state->box, state->x);
                }

                /* if VV, compute the pressure and constraints */
                /* For VV2, we strictly only need this if using pressure
                 * control, but we really would like to have accurate pressures
                 * printed out.
                 * Think about ways around this in the future?
                 * For now, keep this choice in comments.
                 */
                /*bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir)); */
                /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));*/
                bPres = TRUE;
                bTemp = ((ir->eI == eiVV && (!bInitStep)) || (ir->eI == eiVVAK));
                if (bCalcEner && ir->eI == eiVVAK)  /*MRS:  7/9/2010 -- this still doesn't fix it?*/
                {
                    bSumEkinhOld = TRUE;
                }
                /* for vv, the first half of the integration actually corresponds to the previous step.
                   So we need information from the last step in the first half of the integration */
                if (bGStat || do_per_step(step-1, nstglobalcomm))
                {
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, state->box,
                                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                    cglo_flags
                                    | CGLO_ENERGY
                                    | (bTemp ? CGLO_TEMPERATURE : 0)
                                    | (bPres ? CGLO_PRESSURE : 0)
                                    | (bPres ? CGLO_CONSTRAINT : 0)
                                    | ((iterate.bIterationActive) ? CGLO_ITERATE : 0)
                                    | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                    | CGLO_SCALEEKIN
                                    );
                    /* explanation of above:
                       a) We compute Ekin at the full time step
                       if 1) we are using the AveVel Ekin, and it's not the
                       initial step, or 2) if we are using AveEkin, but need the full
                       time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                       b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                       EkinAveVel because it's needed for the pressure */
                }
                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (!bInitStep)
                {
                    if (bTrotter)
                    {
                        m_add(force_vir, shake_vir, total_vir); /* we need the un-dispersion corrected total vir here */
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ2);
                    }
                    else
                    {
                        if (bExchanged)
                        {

                            /* We need the kinetic energy at minus the half step for determining
                             * the full step kinetic energy and possibly for T-coupling.*/
                            /* This may not be quite working correctly yet . . . . */
                            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                            wcycle, enerd, NULL, NULL, NULL, NULL, mu_tot,
                                            constr, NULL, FALSE, state->box,
                                            top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
                        }
                    }
                }

                if (iterate.bIterationActive &&
                    done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                                   state->veta, &vetanew))
                {
                    break;
                }
                bFirstIterate = FALSE;
            }

            if (bTrotter && !bInitStep)
            {
                enerd->term[F_DVDL_BONDED] += dvdl;        /* only add after iterations */
                copy_mat(shake_vir, state->svir_prev);
                copy_mat(force_vir, state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI == eiVV)
                {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, NULL, (ir->eI == eiVV), FALSE, FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (bInitStep && ir->eI == eiVV)
            {
                copy_rvecn(cbuf, state->v, 0, state->natoms);
            }

            if (fr->bSepDVDL && fplog && do_log)
            {
                fprintf(fplog, sepdvdlformat, "Constraint", 0.0, dvdl);
            }
            enerd->term[F_DVDL_BONDED] += dvdl;

            GMX_MPE_LOG(ev_timestep1);
        }

        /* MRS -- now done iterating -- compute the conserved quantity */
        if (bVV)
        {
            saved_conserved_quantity = compute_conserved_from_auxiliary(ir, state, &MassQ);
            if (ir->eI == eiVV)
            {
                last_ekin = enerd->term[F_EKIN];
            }
            if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
            {
                saved_conserved_quantity -= enerd->term[F_DISPCORR];
            }
            /* sum up the foreign energy and dhdl terms for vv.  currently done every step so that dhdl is correct in the .edr */
            if (!bRerunMD)
            {
                sum_dhdl(enerd, state->lambda, ir->fepvals);
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        if (bDoExpanded)
        {
            /* perform extended ensemble sampling in lambda - we don't
               actually move to the new state before outputting
               statistics, but if performing simulated tempering, we
               do update the velocities and the tau_t. */

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, &df_history, step, mcrng, state->v, mdatoms);
        }
        /* ################## START TRAJECTORY OUTPUT ################# */

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        mdof_flags = 0;
        if (do_per_step(step, ir->nstxout))
        {
            mdof_flags |= MDOF_X;
        }
        if (do_per_step(step, ir->nstvout))
        {
            mdof_flags |= MDOF_V;
        }
        if (do_per_step(step, ir->nstfout))
        {
            mdof_flags |= MDOF_F;
        }
        if (do_per_step(step, ir->nstxtcout))
        {
            mdof_flags |= MDOF_XTC;
        }
        if (bCPT)
        {
            mdof_flags |= MDOF_CPT;
        }
        ;

#if defined(GMX_FAHCORE) || defined(GMX_WRITELASTSTEP)
        if (bLastStep)
        {
            /* Enforce writing positions and velocities at end of run */
            mdof_flags |= (MDOF_X | MDOF_V);
        }
#endif
#ifdef GMX_FAHCORE
        if (MASTER(cr))
        {
            fcReportProgress( ir->nsteps, step );
        }

        /* sync bCPT and fc record-keeping */
        if (bCPT && MASTER(cr))
        {
            fcRequestCheckPoint();
        }
#endif

        if (mdof_flags != 0)
        {
            wallcycle_start(wcycle, ewcTRAJ);
            if (bCPT)
            {
                if (state->flags & (1<<estLD_RNG))
                {
                    get_stochd_state(upd, state);
                }
                if (state->flags  & (1<<estMC_RNG))
                {
                    get_mc_state(mcrng, state);
                }
                if (MASTER(cr))
                {
                    if (bSumEkinhOld)
                    {
                        state_global->ekinstate.bUpToDate = FALSE;
                    }
                    else
                    {
                        update_ekinstate(&state_global->ekinstate, ekind);
                        state_global->ekinstate.bUpToDate = TRUE;
                    }
                    update_energyhistory(&state_global->enerhist, mdebin);
                    if (ir->efep != efepNO || ir->bSimTemp)
                    {
                        state_global->fep_state = state->fep_state; /* MRS: seems kludgy. The code should be
                                                                       structured so this isn't necessary.
                                                                       Note this reassignment is only necessary
                                                                       for single threads.*/
                        copy_df_history(&state_global->dfhist, &df_history);
                    }
                }
            }
            write_traj(fplog, cr, outf, mdof_flags, top_global,
                       step, t, state, state_global, f, f_global, &n_xtc, &x_xtc);
            if (bCPT)
            {
                nchkpt++;
                bCPT = FALSE;
            }
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj,
                 * because a checkpoint file will always be written
                 * at the last step.
                 */
                fprintf(stderr, "\nWriting final coordinates.\n");
                if (fr->bMolPBC)
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog, ir->ePBC, state->box, top_global, state_global->x);
                }
                write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                                    *top_global->name, top_global,
                                    state_global->x, state_global->v,
                                    ir->ePBC, state->box);
                debug_gmx();
            }
            wallcycle_stop(wcycle, ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);

        /* TOMAS KUBAR - then, this is not necessary...
        if (!(mdof_flags & (MDOF_X | MDOF_XTC)))
        {
            ct_collect_x(cr, state_global);
        }
         * but still we try to communicate coordinates:
         */
         /*
#ifdef GMX_MPI
        if (MASTER(cr))
        {
            for (i=1; i<ct_mpi_size; i++)
            {
                MPI_Send(state_global->x, 3 * top_global->natoms, GMX_MPI_REAL, i, 300+i, ct_mpi_comm);
                printf("Data of size %d sent to rank %d\n", 3 * top_global->natoms, i);
            }
        }
        else
        {
            return_value = MPI_Recv(state_global->x, 3 * top_global->natoms, GMX_MPI_REAL, 0, 300 + ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
            if (!return_value)
                printf("Data of size %d received on rank %d\n", 3 * top_global->natoms, ct_mpi_rank);
            else
                printf("Error receiving data on rank %d\n", ct_mpi_rank);
        }
        for (i=0; i<10; i++)
        {
            printf("rank%d atom%d %12.7f%12.7f%12.7f\n", ct_mpi_rank, i, state_global->x[i][XX]*10, state_global->x[i][YY]*10, state_global->x[i][ZZ]*10);
        }
#endif
*/
        /* END TOMAS KUBAR */

        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV)
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }
        /*  ################## END TRAJECTORY OUTPUT ################ */

        /* Determine the wallclock run time up till now */
        run_time = gmx_gettime() - (double)runtime->real;

        /* Check whether everything is still allright */
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREAD_MPI
            && MASTER(cr)
#endif
            )
        {
            /* this is just make gs.sig compatible with the hack
               of sending signals around by MPI_Reduce with together with
               other floats */
            if (gmx_get_stop_condition() == gmx_stop_cond_next_ns)
            {
                gs.sig[eglsSTOPCOND] = 1;
            }
            if (gmx_get_stop_condition() == gmx_stop_cond_next)
            {
                gs.sig[eglsSTOPCOND] = -1;
            }
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition = (int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.sig[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            run_time > max_hours*60.0*60.0*0.495)
        {
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        if (ir->nstlist == -1 && !bRerunMD)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir, fr, &top->cgs,
                                                     nlh.scale_tot, state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            run_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }

        /* at the start of step, randomize or scale the velocities (trotter done elsewhere) */
        if (EI_VV(ir->eI))
        {
            if (!bInitStep)
            {
                update_tcouple(fplog, step, ir, state, ekind, wcycle, upd, &MassQ, mdatoms);
            }
            if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
            {
                gmx_bool bIfRandomize;
                bIfRandomize = update_randomize_velocities(ir, step, mdatoms, state, upd, &top->idef, constr);
                /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
                if (constr && bIfRandomize)
                {
                    update_constraints(fplog, step, &dvdl, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir, NULL,
                                       cr, nrnb, wcycle, upd, constr,
                                       bInitStep, TRUE, bCalcVir, vetanew);
                }
            }
        }

        if (bIterativeCase && do_per_step(step, ir->nstpcouple))
        {
            gmx_iterate_init(&iterate, TRUE);
            /* for iterations, we save these vectors, as we will be redoing the calculations */
            copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
        }

        bFirstIterate = TRUE;
        while (bFirstIterate || iterate.bIterationActive)
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
            }

            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */

            /* #########   START SECOND UPDATE STEP ################# */
            GMX_MPE_LOG(ev_update_start);
            /* Box is changed in update() when we do pressure coupling,
             * but we should still use the old box for energy corrections and when
             * writing it to the energy file, so it matches the trajectory files for
             * the same timestep above. Make a copy in a separate array.
             */
            copy_mat(state->box, lastbox);

            bOK = TRUE;
            if (!(bRerunMD && !rerun_fr.bV && !bForceUpdate))
            {
                wallcycle_start(wcycle, ewcUPDATE);
                dvdl = 0;
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
                if (bTrotter)
                {
                    if (iterate.bIterationActive)
                    {
                        if (bFirstIterate)
                        {
                            scalevir = 1;
                        }
                        else
                        {
                            /* we use a new value of scalevir to converge the iterations faster */
                            scalevir = tracevir/trace(shake_vir);
                        }
                        msmul(shake_vir, scalevir, shake_vir);
                        m_add(force_vir, shake_vir, total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ3);
                    /* We can only do Berendsen coupling after we have summed
                     * the kinetic energy or virial. Since the happens
                     * in global_state after update, we should only do it at
                     * step % nstlist = 1 with bGStatEveryStep=FALSE.
                     */
                }
                else
                {
                    update_tcouple(fplog, step, ir, state, ekind, wcycle, upd, &MassQ, mdatoms);
                    update_pcouple(fplog, step, ir, state, pcoupl_mu, M, wcycle,
                                   upd, bInitStep);
                }

                if (bVV)
                {
                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    /* velocity half-step update */
                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, fcd,
                                  ekind, M, wcycle, upd, FALSE, etrtVELOCITY2,
                                  cr, nrnb, constr, &top->idef);
                }

                /* Above, initialize just copies ekinh into ekin,
                 * it doesn't copy position (for VV),
                 * and entire integrator for MD.
                 */

                if (ir->eI == eiVVAK)
                {
                    copy_rvecn(state->x, cbuf, 0, state->natoms);
                }
                bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                              bUpdateDoLR, fr->f_twin, fcd,
                              ekind, M, wcycle, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                wallcycle_stop(wcycle, ewcUPDATE);

                update_constraints(fplog, step, &dvdl, ir, ekind, mdatoms, state,
                                   fr->bMolPBC, graph, f,
                                   &top->idef, shake_vir, force_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   bInitStep, FALSE, bCalcVir, state->veta);

                if (ir->eI == eiVVAK)
                {
                    /* erase F_EKIN and F_TEMP here? */
                    /* just compute the kinetic energy at the half step to perform a trotter step */
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, lastbox,
                                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                    cglo_flags | CGLO_TEMPERATURE
                                    );
                    wallcycle_start(wcycle, ewcUPDATE);
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ4);
                    /* now we know the scaling, we can compute the positions again again */
                    copy_rvecn(cbuf, state->x, 0, state->natoms);

                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, fcd,
                                  ekind, M, wcycle, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                    wallcycle_stop(wcycle, ewcUPDATE);

                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure.
                     * For now, will call without actually constraining, constr=NULL*/
                    update_constraints(fplog, step, &dvdl, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir, force_vir,
                                       cr, nrnb, wcycle, upd, NULL,
                                       bInitStep, FALSE, bCalcVir,
                                       state->veta);
                }
                if (!bOK && !bFFscan)
                {
                    gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constrains");
                }

                if (fr->bSepDVDL && fplog && do_log)
                {
                    fprintf(fplog, sepdvdlformat, "Constraint dV/dl", 0.0, dvdl);
                }
                enerd->term[F_DVDL_BONDED] += dvdl;
            }
            else if (graph)
            {
                /* Need to unshift here */
                unshift_self(graph, state->box, state->x);
            }

            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_update_finish);

            if (vsite != NULL)
            {
                wallcycle_start(wcycle, ewcVSITECONSTR);
                if (graph != NULL)
                {
                    shift_self(graph, state->box, state->x);
                }
                construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
                                 top->idef.iparams, top->idef.il,
                                 fr->ePBC, fr->bMolPBC, graph, cr, state->box);

                if (graph != NULL)
                {
                    unshift_self(graph, state->box, state->x);
                }
                wallcycle_stop(wcycle, ewcVSITECONSTR);
            }

            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints  ############ */
            /* With Leap-Frog we can skip compute_globals at
             * non-communication steps, but we need to calculate
             * the kinetic energy one step before communication.
             */
            if (bGStat || (!EI_VV(ir->eI) && do_per_step(step+1, nstglobalcomm)))
            {
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    gs.sig[eglsNABNSB] = nlh.nabnsb;
                }
                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr,
                                bFirstIterate ? &gs : NULL,
                                (step_rel % gs.nstms == 0) &&
                                (multisim_nsteps < 0 || (step_rel < multisim_nsteps)),
                                lastbox,
                                top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                cglo_flags
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_ENERGY : 0)
                                | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
                                | (iterate.bIterationActive ? CGLO_ITERATE : 0)
                                | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                | CGLO_CONSTRAINT
                                );
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    nlh.nabnsb         = gs.set[eglsNABNSB];
                    gs.set[eglsNABNSB] = 0;
                }
            }
            /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
            /* #############  END CALC EKIN AND PRESSURE ################# */

            /* Note: this is OK, but there are some numerical precision issues with using the convergence of
               the virial that should probably be addressed eventually. state->veta has better properies,
               but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
               generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

            if (iterate.bIterationActive &&
                done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                               trace(shake_vir), &tracevir))
            {
                break;
            }
            bFirstIterate = FALSE;
        }

        /* only add constraint dvdl after constraints */
        enerd->term[F_DVDL_BONDED] += dvdl;
        if (!bVV || bRerunMD)
        {
            /* sum up the foreign energy and dhdl terms for md and sd. currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }
        update_box(fplog, step, ir, mdatoms, state, graph, f,
                   ir->nstlist == -1 ? &nlh.scale_tot : NULL, pcoupl_mu, nrnb, wcycle, upd, bInitStep, FALSE);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
        if (bFFscan && (shellfc == NULL || bConverged))
        {
            if (print_forcefield(fplog, enerd->term, mdatoms->homenr,
                                 f, NULL, xcopy,
                                 &(top_global->mols), mdatoms->massT, pres))
            {
                gmx_finalize_par();

                fprintf(stderr, "\n");
                exit(0);
            }
        }
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        if (bTCR)
        {
            /* Only do GCT when the relaxation of shells (minimization) has converged,
             * otherwise we might be coupling to bogus energies.
             * In parallel we must always do this, because the other sims might
             * update the FF.
             */

            /* Since this is called with the new coordinates state->x, I assume
             * we want the new box state->box too. / EL 20040121
             */
            do_coupling(fplog, oenv, nfile, fnm, tcr, t, step, enerd->term, fr,
                        ir, MASTER(cr),
                        mdatoms, &(top->idef), mu_aver,
                        top_global->mols.nr, cr,
                        state->box, total_vir, pres,
                        mu_tot, state->x, f, bConverged);
            debug_gmx();
        }

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && ir->eI == eiVV)
        {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

        if (bVV)
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
        }
        else
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + compute_conserved_from_auxiliary(ir, state, &MassQ);
        }
        /* Check for excessively large energies */
        if (bIonize)
        {
#ifdef GMX_DOUBLE
            real etot_max = 1e200;
#else
            real etot_max = 1e30;
#endif
            if (fabs(enerd->term[F_ETOT]) > etot_max)
            {
                fprintf(stderr, "Energy too large (%g), giving up\n",
                        enerd->term[F_ETOT]);
            }
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* TOMAS KUBAR
         *  energies
         */
#ifdef GMX_MPI
        if (ct_mpi_rank == 0) {
#endif
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteSURFACEHOPPING || 
              ct->jobtype==cteFERMI || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteFERMISFHOPPING || ct->jobtype==cteTULLYFEWESTSWITCHES||ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING) {
         /* calculate the QM energy
          * attention - the correct expression must be used for that!
          * THIS MAY BE WRONG FOR cteFERMI, WHERE WE MIX THE ADIABATIC STATES AND
          *     CONSTRUCT ARTIFICIAL ct->wf !
          */
         ct_energy = ct_energy1 = ct_energy2 = 0.0;
         for (i=0; i<ct->dim; i++)
           for (m=0; m<ct->dim; m++) {
             /* TB Hamiltonian */
             ct_energy1 += ct->hamiltonian[i][m] * (ct->wf[i] * ct->wf[m] + ct->wf[i + ct->dim] * ct->wf[m + ct->dim]);
             /* Hubbard / gamma terms */
             ct_energy2 += 0.5 * ct->hubbard[i][m] * (SQR(ct->wf[i]) + SQR(ct->wf[i + ct->dim])) * (SQR(ct->wf[m]) + SQR(ct->wf[m + ct->dim]));
           }
         ct_energy = ct_energy1 + ct_energy2;
         /* output both MM and QM energy */
         fprintf(f_ct_energy, "%10d %12.7f %12.7f %12.7f %12.7f %12.7f\n", step, ct_energy1, ct_energy2, ct_energy, enerd->term[F_ETOT] * KJMOL_TO_HARTREE, enerd->term[F_ETOT] * KJMOL_TO_HARTREE + ct_energy);
        }
        if (ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADNONSCC) {
         /* dtto, without Hubbard and MM energy */
         ct_energy = 0.0;
         for (i=0; i<ct->dim; i++)
           for (m=0; m<ct->dim; m++) {
             ct_energy += ct->hamiltonian[i][m] * (ct->wf[i] * ct->wf[m] + ct->wf[i + ct->dim] * ct->wf[m + ct->dim]);
           }
         /* output the QM energy */
         fprintf(f_ct_energy, "%10d %12.7f\n", step, ct_energy);
        }
#ifdef GMX_MPI
        }
#endif
        /* END TOMAS KUBAR */

        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep)
        {
            runtime_upd_proc(runtime);
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool do_dr, do_or;

            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog, ir->fepvals, ir->expandedvals, ir->bSimTemp ? ir->simtempvals : NULL,
                                          &df_history, state->fep_state, ir->nstlog, step);
            }
            if (!(bStartingFromCpt && (EI_VV(ir->eI))))
            {
                if (bCalcEner)
                {
                    upd_mdebin(mdebin, bDoDHDL, TRUE,
                               t, mdatoms->tmass, enerd, state,
                               ir->fepvals, ir->expandedvals, lastbox,
                               shake_vir, force_vir, total_vir, pres,
                               ekind, mu_tot, constr);
                }
                else
                {
                    upd_mdebin_step(mdebin);
                }

                do_dr  = do_per_step(step, ir->nstdisreout);
                do_or  = do_per_step(step, ir->nstorireout);

                print_ebin(outf->fp_ene, do_ene, do_dr, do_or, do_log ? fplog : NULL,
                           step, t,
                           eprNORMAL, bCompact, mdebin, fcd, groups, &(ir->opts));
            }
            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull, step, t);
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }
        if (bDoExpanded)
        {
            /* Have to do this part after outputting the logfile and the edr file */
            state->fep_state = lamnew;
            for (i = 0; i < efptNR; i++)
            {
                state_global->lambda[i] = ir->fepvals->all_lambda[i][lamnew];
            }
        }
        /* Remaining runtime */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal()) && !bPMETuneRunning)
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, runtime, step, ir, cr);
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step, repl_ex_nst))
        {
            bExchanged = replica_exchange(fplog, cr, repl_ex,
                                          state_global, enerd,
                                          state, step, t);

            if (bExchanged && DOMAINDECOMP(cr))
            {
                dd_partition_system(fplog, step, cr, TRUE, 1,
                                    state_global, top_global, ir,
                                    state, &f, mdatoms, top, fr,
                                    vsite, shellfc, constr,
                                    nrnb, wcycle, FALSE);
            }
        }

        bFirstStep       = FALSE;
        bInitStep        = FALSE;
        bStartingFromCpt = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & (1<<estPRES_PREV)) &&
            (bGStatEveryStep ||
             (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres, state->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if ( (membed != NULL) && (!bLastStep) )
        {
            rescale_membed(step_rel, membed, state_global->x);
        }

        if (bRerunMD)
        {
            if (MASTER(cr))
            {
                /* read next frame from input trajectory */
                bNotLastFrame = read_next_frame(oenv, status, &rerun_fr);
            }

            if (PAR(cr))
            {
                rerun_parallel_comm(cr, &rerun_fr, &bNotLastFrame);
            }
        }

        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (bPMETuneRunning || bPMETuneTry)
        {
            /* PME grid + cut-off optimization with GPUs or PME nodes */

            /* Count the total cycles over the last steps */
            cycles_pmes += cycles;

            /* We can only switch cut-off at NS steps */
            if (step % ir->nstlist == 0)
            {
                /* PME grid + cut-off optimization with GPUs or PME nodes */
                if (bPMETuneTry)
                {
                    if (DDMASTER(cr->dd))
                    {
                        /* PME node load is too high, start tuning */
                        bPMETuneRunning = (dd_pme_f_ratio(cr->dd) >= 1.05);
                    }
                    dd_bcast(cr->dd, sizeof(gmx_bool), &bPMETuneRunning);

                    if (bPMETuneRunning || step_rel > ir->nstlist*50)
                    {
                        bPMETuneTry     = FALSE;
                    }
                }
                if (bPMETuneRunning)
                {
                    /* init_step might not be a multiple of nstlist,
                     * but the first cycle is always skipped anyhow.
                     */
                    bPMETuneRunning =
                        pme_load_balance(pme_loadbal, cr,
                                         (bVerbose && MASTER(cr)) ? stderr : NULL,
                                         fplog,
                                         ir, state, cycles_pmes,
                                         fr->ic, fr->nbv, &fr->pmedata,
                                         step);

                    /* Update constants in forcerec/inputrec to keep them in sync with fr->ic */
                    fr->ewaldcoeff = fr->ic->ewaldcoeff;
                    fr->rlist      = fr->ic->rlist;
                    fr->rlistlong  = fr->ic->rlistlong;
                    fr->rcoulomb   = fr->ic->rcoulomb;
                    fr->rvdw       = fr->ic->rvdw;
                }
                cycles_pmes = 0;
            }
        }

        if (step_rel == wcycle_get_reset_counters(wcycle) ||
            gs.set[eglsRESETCOUNTERS] != 0)
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog, cr, step, &step_rel, ir, wcycle, nrnb, runtime,
                               fr->nbv != NULL && fr->nbv->bUseGPU ? fr->nbv->cu_nbv : NULL);
            wcycle_set_reset_counters(wcycle, -1);
            /* Correct max_hours for the elapsed time */
            max_hours                -= run_time/(60.0*60.0);
            bResetCountersHalfMaxH    = FALSE;
            gs.set[eglsRESETCOUNTERS] = 0;
        }
    }
    /* End of main MD loop */
    debug_gmx();

    /* Stop the time */
    runtime_end(runtime);

    if (bRerunMD && MASTER(cr))
    {
        close_trj(status);
    }

    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0 && !bRerunMD)
        {
            print_ebin(outf->fp_ene, FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, FALSE, mdebin, fcd, groups, &(ir->opts));
        }
    }

    done_mdoutf(outf);

    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog, "Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n", nlh.s1/nlh.nns, sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog, "Average number of atoms that crossed the half buffer length: %.1f\n\n", nlh.ab/nlh.nns);
    }

    if (pme_loadbal != NULL)
    {
        pme_loadbal_done(pme_loadbal, fplog);
    }

    if (shellfc && fplog)
    {
        fprintf(fplog, "Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog, "Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }

    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog, repl_ex);
    }

    runtime->nsteps_done = step_rel;

  /* TOMAS KUBAR */
#ifdef GMX_MPI
   if (ct_mpi_rank == 0) {
#endif
       if (f_ct_energy) fclose(f_ct_energy);
       if (f_ct_esp) fclose(f_ct_esp);
       if (f_ct_shift) fclose(f_ct_shift);
       if (f_tb_hamiltonian) fclose(f_tb_hamiltonian);
       if (f_tb_hamiltonian_hub) fclose(f_tb_hamiltonian_hub);
       if (f_tb_hubbard) fclose(f_tb_hubbard);
       if (f_tb_occupation) fclose(f_tb_occupation);
       if (f_ct_adiabatic) fclose(f_ct_adiabatic);
       if (f_ct_exp_adiab) fclose(f_ct_exp_adiab);
       if (f_ct_surfacehopping) fclose(f_ct_surfacehopping);
       if (f_ct_orbital) fclose(f_ct_orbital);
       if (f_ct_current) fclose(f_ct_current);
#ifdef GMX_MPI
   }
#endif
  /* END TOMAS KUBAR */

    return 0;
}
