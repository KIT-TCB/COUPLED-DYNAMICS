! negf method from Phys. Rev. B  80, 245311 (2009)
!**************************************************************
! Lorentz level-width functions 
! Ozaki spectral decomposition of the Fermi function


MODULE global
        implicit none 

!******** wire parameters: **************
        INTEGER         ::       n
        INTEGER         ::       n_lorentz
        
!       DOUBLE PRECISION, allocatable, dimension(:):: E         ! here also intersite coupling is time dependent
!       DOUBLE PRECISION, allocatable, dimension(:):: delta 

!******** lead parameters: *********
!       DOUBLE PRECISION                                :: E_F_left, E_F_right
!       INTEGER                                         :: N_poles_L, N_poles_R 
        
        INTEGER                                         :: kL, kR

!       DOUBLE PRECISION                                :: gam_L, eps_L, w0_L, gam_R, eps_R, w0_R
        
        DOUBLE COMPLEX, allocatable, dimension(:,:)     :: rho, H 
        
        DOUBLE COMPLEX, allocatable, dimension(:,:,:)   :: Pi_L    
        DOUBLE COMPLEX, allocatable, dimension(:,:,:)   :: Pi_R

        DOUBLE COMPLEX, allocatable, dimension(:,:,:,:) :: Omega_LL1, Omega_LL2, Omega_LL3
        DOUBLE COMPLEX, allocatable, dimension(:,:,:,:) :: Omega_LR1, Omega_LR2, Omega_LR3 

        DOUBLE COMPLEX, allocatable, dimension(:,:,:,:) :: Omega_RR1, Omega_RR2, Omega_RR3
        DOUBLE COMPLEX, allocatable, dimension(:,:,:,:) :: Omega_RL1, Omega_RL2, Omega_RL3

        DOUBLE COMPLEX, allocatable, dimension(:,:,:)   :: Gam_greater_L_m, Gam_greater_L_p, Gam_lesser_L_m, Gam_lesser_L_p 
        DOUBLE COMPLEX, allocatable, dimension(:,:,:)   :: Gam_greater_R_m, Gam_greater_R_p, Gam_lesser_R_m, Gam_lesser_R_p 
        
        DOUBLE COMPLEX, allocatable, dimension(:)       :: hi_left_m, hi_left_p, hi_right_m, hi_right_p

!       DOUBLE COMPLEX, allocatable, dimension(:)       :: rkvec, drkvec          
        INTEGER                                         :: length_rkvec
                
!******** matrices utilized to compute the poles of the Ozaki decomposition of the Fermi function                   
!     DOUBLE PRECISION, allocatable,dimension(:,:)      :: Mat_L, Mat_R, Eig_vect_mat_L, Eig_vect_mat_R
!     DOUBLE PRECISION, allocatable,dimension(:)        :: Eig_val_mat_L, Eig_val_mat_R, Nu_L, Nu_R, R_alpha_L, R_alpha_R, R_L, R_R

!******** others: ************************
!     DOUBLE PRECISION                                  :: t_0, t_end, t_step
      DOUBLE COMPLEX, parameter                         :: im       = dcmplx(0.d0,1.d0)  !imaginary unit
      DOUBLE COMPLEX, parameter                         :: hbar     = 0.658211928d0      !(eV*fs)
      DOUBLE PRECISION, parameter                       :: hbar_r   = 0.658211928d0      !(eV*fs)
      DOUBLE PRECISION, parameter                       :: ha_to_ev = 2.7211396132d1     !(convert Hartree to eV)
END MODULE global


MODULE rksuite_vars
      implicit none

      double precision                       :: TSTART !starting time (psi0)
      double precision                       :: TNOW   !output (time reached by rksuite)
      double precision                       :: TEND   !end time
      integer                                :: NEQ    !number of eq (*2 complex->real), defined in rk_init
      integer, parameter                     :: METHOD = 3 !(2 or 3)
      integer                                :: LENWRK !in this code set to 16*NEQ
      double precision, parameter            :: TOL = 1.d-6 !(rel error)
      character, parameter                   :: TASK = "U"
      logical, parameter                     :: ERRASS = .false.
      logical                                :: IsRKInit = .false.
      logical                                :: message = .false.
      integer                                :: OUTCH
      double precision                       :: MCHPES
      double precision                       :: DWARF
      double precision, allocatable          :: THRES(:)  !these 3 allocated in rk_init since number of sites parameter (-> NEQ)
      double precision, allocatable          :: WORK(:), YMAX(:)
      integer                                :: IFAIL
      logical                                :: expinit = .false.
END MODULE rksuite_vars

!!!!!!!!!!!!!!!!!!!!!!!!!
! wrapper routine -- one step of NEGF integration
!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE negf_do_step(rkvec, drkvec, t_step)
       use global
       use rksuite_vars
       implicit none
       double complex   :: rkvec(length_rkvec), drkvec(length_rkvec)
       double precision :: t_step, t_step_in_fs

!      call create_H
                     
       t_step_in_fs = t_step / 41.341373             ! AU of time to femtosecond
       TSTART = 0.d0
       TEND   = t_step_in_fs * 1.1d0
       CALL SETUP(NEQ,TSTART,rkvec,TEND,TOL,THRES,METHOD,TASK,ERRASS,TSTART,WORK,LENWRK,MESSAGE)
       call negf_solve_with_rk(t_step_in_fs, rkvec, drkvec)

!      call save_output(t_step, 90)

END SUBROUTINE negf_do_step

!!!!!!!!!!!!!!!!!!!!!!!
! preparations...
!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE negf_allocate_fortran(n_in, N_poles_L, N_poles_R, length_rkvec_in)
        use global
        implicit none
        integer :: n_in, N_poles_L, N_poles_R, length_rkvec_in
!---
        integer :: check_alloc, alloc_status

        n            = n_in
        length_rkvec = length_rkvec_in

        alloc_status = 0
        check_alloc  = 0

        allocate(H(n,n), stat= alloc_status)
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(rho(n,n), stat= alloc_status)            
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1        
        if(check_alloc .NE.0) write(*,*) 'Error in allocation of H/rho/G!'
        check_alloc = 0

        kL = N_poles_L + n_lorentz        
        kR = N_poles_R + n_lorentz

        allocate(hi_left_m(kL),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(hi_left_p(kL),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(hi_right_m(kR), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(hi_right_p(kR), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1  
        if(check_alloc .NE.0) write(*,*) 'Error in allocation of hi !'
        check_alloc = 0

        allocate(Gam_greater_L_m(n,n,kL), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_greater_L_p(n,n,kL), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_lesser_L_m(n,n,kL),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_lesser_L_p(n,n,kL),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_greater_R_m(n,n,kR), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_greater_R_p(n,n,kR), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_lesser_R_m(n,n,kR),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Gam_lesser_R_p(n,n,kR),  stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        if(check_alloc .NE.0) write(*,*) 'Error in allocation of Gam !'
        check_alloc = 0

        allocate(Pi_L(n,n,kL+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Pi_R(n,n,kR+1), stat= alloc_status)
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        if(check_alloc .NE.0) write(*,*) 'Error in allocation of Pi !'
        check_alloc = 0

        allocate(Omega_LL1(n,n,n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_LL2(n,n,n_lorentz,kL+1-n_lorentz), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_LL3(n,n,kL-n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_LR1(n,n,n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_LR2(n,n,n_lorentz,kL+1-n_lorentz), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_LR3(n,n,kL-n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1 
        allocate(Omega_RR1(n,n,n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_RR2(n,n,n_lorentz,kR+1-n_lorentz), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_RR3(n,n,kR-n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1 
        allocate(Omega_RL1(n,n,n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_RL2(n,n,n_lorentz,kR+1-n_lorentz), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1
        allocate(Omega_RL3(n,n,kR-n_lorentz,n_lorentz+1), stat= alloc_status)      
                                                        if (alloc_status .NE.0) check_alloc=check_alloc+1  
        if(check_alloc .NE.0) write(*,*) 'Error in allocation of Omega !'

END SUBROUTINE negf_allocate_fortran

SUBROUTINE negf_create_hi(eps_L, w0_L, E_F_left,  Nu_L, N_poles_L, &
                          eps_R, w0_R, E_F_right, Nu_R, N_poles_R, beta)

        use global
        implicit none
!       integer          :: kL, kR, N_poles_L, N_poles_R
        integer          :: N_poles_L, N_poles_R
        double precision :: beta
!       double complex   :: hi_left_m(kL), hi_left_p(kL), hi_right_m(kR), hi_right_p(kR)
        double precision :: eps_L, w0_L, E_F_left, Nu_L(N_poles_L), eps_R, w0_R, E_F_right, Nu_R(N_poles_R)
!---
        integer          :: i
                 
        hi_left_m(1) = eps_L - im*w0_L  
        hi_left_p(1) = eps_L + im*w0_L  

        hi_right_m(1) = eps_R - im*w0_R  
        hi_right_p(1) = eps_R + im*w0_R  


        DO i=2,kL               
             hi_left_m(i) = E_F_left - im / (Nu_L(i-1)*beta)
             hi_left_p(i) = E_F_left + im / (Nu_L(i-1)*beta)             
        END DO 

        DO i=2,kR
             hi_right_m(i) = E_F_right - im / (Nu_R(i-1)*beta) 
             hi_right_p(i) = E_F_right + im / (Nu_R(i-1)*beta)
        END DO
        
END SUBROUTINE negf_create_hi

SUBROUTINE negf_create_gam(gam_L, w0_L, eps_L, E_F_left,  R_L, &
                           gam_R, w0_R, eps_R, E_F_right, R_R, &
                           N_poles_L, N_poles_R, beta)

        use global
        implicit none
!       double complex   :: Gam_greater_L_m(n,n,kL), Gam_greater_L_p(n,n,kL), Gam_lesser_L_m(n,n,kL), Gam_lesser_L_p(n,n,kL)
!       double complex   :: Gam_greater_R_m(n,n,kR), Gam_greater_R_p(n,n,kR), Gam_lesser_R_m(n,n,kR), Gam_lesser_R_p(n,n,kR)
        double precision :: gam_L, w0_L, eps_L, E_F_left,  R_L(N_poles_L)
        double precision :: gam_R, w0_R, eps_R, E_F_right, R_R(N_poles_R)
!       double complex   :: hi_left_m(kL), hi_left_p(kL), hi_right_m(kR), hi_right_p(kR)
!       integer          :: n, kL, kR, N_poles_L, N_poles_R
        integer          :: N_poles_L, N_poles_R
        double precision :: beta
!---
        integer          :: i, j, k
        double complex   :: negf_fermi, negf_spectraldensity 
!       double complex   :: sigma_greater, sigma_lesser  

        DO i=1,n
           DO j=1,n

              IF( (i.EQ.1) .AND. (j.EQ.1)) THEN 

                   Gam_greater_L_m(i,j,1) = -dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_L*w0_L * negf_fermi(-(eps_L-im*w0_L)+E_F_left, beta)        
                   Gam_greater_L_p(i,j,1) = -dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_L*w0_L * negf_fermi(-(eps_L+im*w0_L)+E_F_left, beta)
                   Gam_lesser_L_m(i,j,1)  =  dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_L*w0_L * negf_fermi((eps_L-im*w0_L)-E_F_left, beta)
                   Gam_lesser_L_p(i,j,1)  =  dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_L*w0_L * negf_fermi((eps_L+im*w0_L)-E_F_left, beta)  

                 DO k=2,kL 
                    Gam_greater_L_m(i,j,k) = -dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_left_m(k), gam_L, w0_L, eps_L) * R_L(k-1)  
                    Gam_greater_L_p(i,j,k) =  dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_left_p(k), gam_L, w0_L, eps_L) * R_L(k-1) 
                    Gam_lesser_L_m(i,j,k)  = -dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_left_m(k), gam_L, w0_L, eps_L) * R_L(k-1) 
                    Gam_lesser_L_p(i,j,k)  =  dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_left_p(k), gam_L, w0_L, eps_L) * R_L(k-1)
                 END DO 
               
              ELSE 

                 DO k=1,kL
                    Gam_greater_L_m(i,j,k) = dcmplx(0.d0,0.d0)              
                    Gam_greater_L_p(i,j,k) = dcmplx(0.d0,0.d0)
                    Gam_lesser_L_m(i,j,k)  = dcmplx(0.d0,0.d0)
                    Gam_lesser_L_p(i,j,k)  = dcmplx(0.d0,0.d0)
                 END DO
              
              END IF 

           END DO
        END DO

        DO i=1,n
           DO j=1,n

              IF( (i.EQ.n) .AND. (j.EQ.n)) THEN 

                  Gam_greater_R_m(i,j,1) =  -dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_R*w0_R * negf_fermi(-(eps_R-im*w0_R)+E_F_right, beta)
                  Gam_greater_R_p(i,j,1) =  -dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_R*w0_R * negf_fermi(-(eps_R+im*w0_R)+E_F_right, beta)
                  Gam_lesser_R_m(i,j,1)  =   dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_R*w0_R * negf_fermi((eps_R-im*w0_R)-E_F_right, beta)  
                  Gam_lesser_R_p(i,j,1)  =   dcmplx(1.d0,0.d0)/dcmplx(2.d0,0.d0)*im*gam_R*w0_R * negf_fermi((eps_R+im*w0_R)-E_F_right, beta)  
       
                 DO k=2,kR
                    Gam_greater_R_m(i,j,k) = -dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_right_m(k), gam_R, w0_R, eps_R) * R_R(k-1)  
                    Gam_greater_R_p(i,j,k) =  dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_right_p(k), gam_R, w0_R, eps_R) * R_R(k-1) 
                    Gam_lesser_R_m(i,j,k)  = -dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_right_m(k), gam_R, w0_R, eps_R) * R_R(k-1) 
                    Gam_lesser_R_p(i,j,k)  =  dcmplx(1.d0,0.d0)/beta * negf_spectraldensity(hi_right_p(k), gam_R, w0_R, eps_R) * R_R(k-1) 
                 END DO 

              ELSE

                 DO k=1,kR 
                    Gam_greater_R_m(i,j,k) = dcmplx(0.d0,0.d0)               
                    Gam_greater_R_p(i,j,k) = dcmplx(0.d0,0.d0)
                    Gam_lesser_R_m(i,j,k)  = dcmplx(0.d0,0.d0)
                    Gam_lesser_R_p(i,j,k)  = dcmplx(0.d0,0.d0)
                 END DO

              END IF 

           END DO
        END DO

END SUBROUTINE negf_create_gam

SUBROUTINE negf_create_h(h_in)

        USE global
        implicit none
        double precision :: h_in(n,n)    ! Hamiltonian in atomic units - convert to eV
!---
!       double precision :: factor
        integer          :: i, j !, l, m

        do i = 1, n
          do j = 1, n
            H(i,j) = dcmplx(h_in(i,j) * ha_to_ev, 0.d0)
          end do
        end do
        
!       DO i=1,n
!         DO j=1,n
!            H(i,j)=0.d0
!         END DO 
!       END DO 

!       DO l=1,n
!          DO m=1,n

!            IF (l==m) THEN   
!               factor=E(l)
!            ELSE IF ( ((l+1)==m).AND.(l .NE. n) ) THEN
!               factor=delta(l)
!            ELSE IF ( (l==(m+1)).AND.(m .NE. n) ) THEN
!               factor=delta(m)
!            ELSE                                 
!               factor=0.d0
!            END IF 

!            H(l,m) = dcmplx(factor,0.d0)

!          END DO
!       END DO           
    
END SUBROUTINE negf_create_h

SUBROUTINE negf_set_initial_values(rkvec, wf)
!                                  rho, Pi_L, Pi_R, &
!                                  Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                                  Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                                  rkvec, &
!                                  n, n_lorentz, kL, kR, length_rkvec)

        use global
        implicit none
!       double complex   :: rho(n,n), Pi_L(n,n,kL+1), Pi_R(n,n,kR+1)
!       double complex   :: Omega_LL1(n,n,n_lorentz,n_lorentz+1), Omega_LL2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LL3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex   :: Omega_LR1(n,n,n_lorentz,n_lorentz+1), Omega_LR2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LR3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex   :: Omega_RR1(n,n,n_lorentz,n_lorentz+1), Omega_RR2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RR3(n,n,kR-n_lorentz,n_lorentz+1)
!       double complex   :: Omega_RL1(n,n,n_lorentz,n_lorentz+1), Omega_RL2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RL3(n,n,kR-n_lorentz,n_lorentz+1)
        double complex   :: rkvec(length_rkvec)
!       integer          :: n, n_lorentz, kL, kR, length_rkvec
        double precision :: wf(2*n)
!---
        integer:: i, j, k, p

        write (*,*) 'negf_set_initial_values. n =', n

        DO i=1,n 
          DO j=1,n

              IF (i==j) THEN  
                rho(i,j)=dcmplx(wf(i)**2 + wf(i+n)**2, 0.d0)
              ELSE
                rho(i,j)=dcmplx(0.d0,0.d0)
              END IF 
              write (*,*) 'initial rho:', i, j, rho(i,j)
            
              DO k=1,kL
                Pi_L(i,j,k)=dcmplx(0.d0,0.d0)
              END DO

              DO k=1,kR
                Pi_R(i,j,k)=dcmplx(0.d0,0.d0)
              END DO 
              
!>>>>>>>>>>>>>>>>>>>>

              DO k=1,n_lorentz
                 DO p=1,n_lorentz
                     Omega_LL1(i,j,k,p)=dcmplx(0.d0,0.d0)
                     Omega_LR1(i,j,k,p)=dcmplx(0.d0,0.d0)                     
                 END DO
              END DO 

              DO k=1,n_lorentz
                 DO p=1,kL-n_lorentz
                    Omega_LL2(i,j,k,p)=dcmplx(0.d0,0.d0)
                    Omega_LR2(i,j,k,p)=dcmplx(0.d0,0.d0)
                 END DO
              END DO

              DO k=1,kL-n_lorentz
                 DO p=1,n_lorentz
                    Omega_LL3(i,j,k,p)=dcmplx(0.d0,0.d0)
                    Omega_LR3(i,j,k,p)=dcmplx(0.d0,0.d0)
                 END DO
              END DO

              DO k=1,n_lorentz
                 DO p=1,n_lorentz
                    Omega_RR1(i,j,k,p)=dcmplx(0.d0,0.d0)
                    Omega_RL1(i,j,k,p)=dcmplx(0.d0,0.d0)     
                 END DO 
              END DO 

              DO k=1,n_lorentz
                 DO p=1,kR-n_lorentz
                    Omega_RR2(i,j,k,p)=dcmplx(0.d0,0.d0)
                    Omega_RL2(i,j,k,p)=dcmplx(0.d0,0.d0)
                 END DO
              END DO
                             
              DO k=1,kR-n_lorentz
                 DO p=1,n_lorentz
                    Omega_RR3(i,j,k,p)=dcmplx(0.d0,0.d0)
                    Omega_RL3(i,j,k,p)=dcmplx(0.d0,0.d0)
                 END DO
              END DO

          END DO 
        END DO 

        call negf_matrix_to_vector(rkvec)
!                             , length_rkvec, Pi_L, Pi_R, rho, &
!                             Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                             Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                             n, n_lorentz, kL, kR)
  
END SUBROUTINE negf_set_initial_values

!******************* eom ************************************************************************************
!                computes the equation of motion needed for the Runke-Kutta method
!                input: mat = rkmat from the Runge-Kutta method
!                output: dmat = all time derivatives in 1 array                  
!************************************************************************************************************

SUBROUTINE negf_eom(t, vec, dvec)

        use global
        implicit none
        double precision                        :: t        ! will not be used whatsover, but it is needed for RKsuite
        double complex, dimension(length_rkvec) :: vec, dvec
!       integer                                 :: length_rkvec, n, n_lorentz, kL, kR
!       double complex                          :: rho(n,n), Pi_L(n,n,kL+1), Pi_R(n,n,kR+1)
!       double complex                          :: Omega_LL1(n,n,n_lorentz,n_lorentz+1), Omega_LL2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LL3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_LR1(n,n,n_lorentz,n_lorentz+1), Omega_LR2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LR3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RR1(n,n,n_lorentz,n_lorentz+1), Omega_RR2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RR3(n,n,kR-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RL1(n,n,n_lorentz,n_lorentz+1), Omega_RL2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RL3(n,n,kR-n_lorentz,n_lorentz+1)
!---
        integer :: k,p
        
        call negf_vector_to_matrix(vec)
!                             , length_rkvec, Pi_L, Pi_R, rho, &
!                             Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                             Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                             n, n_lorentz, kL, kR)
      
        Pi_L(:,:,kL+1) = dcmplx(0.d0,0.d0)
        Pi_R(:,:,kR+1) = dcmplx(0.d0,0.d0)
           
        DO k=1,kL
            Pi_L(:,:,kL+1) = Pi_L(:,:,kL+1) + Pi_L(:,:,k)
        END DO 

        DO k=1,kR
            Pi_R(:,:,kR+1) = Pi_R(:,:,kR+1) + Pi_R(:,:,k)
        END DO                
                                
!>>>>>>>>>>>>>>>>>>>>

        DO k=1,n_lorentz
            Omega_LL1(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
            Omega_LR1(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
        END DO

        DO k=1,n_lorentz
            Omega_LL2(:,:,k,kL-n_lorentz+1)=dcmplx(0.d0,0.d0) 
            Omega_LR2(:,:,k,kL-n_lorentz+1)=dcmplx(0.d0,0.d0)   
        END DO 

        DO k=1,kL-n_lorentz
           Omega_LL3(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
           Omega_LR3(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0) 
        END DO 

        DO k=1,n_lorentz
           Omega_RR1(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
           Omega_RL1(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0) 
        END DO 

        DO k=1,n_lorentz
           Omega_RR2(:,:,k,kR-n_lorentz+1)=dcmplx(0.d0,0.d0) 
           Omega_RL2(:,:,k,kR-n_lorentz+1)=dcmplx(0.d0,0.d0)
        END DO 
 
        DO k=1,kR-n_lorentz
           Omega_RR3(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
           Omega_RL3(:,:,k,n_lorentz+1)=dcmplx(0.d0,0.d0)
        END DO 

!>>>>>>>>>>>>>>>>>>>>

        DO k=1,n_lorentz
           DO p=1,n_lorentz 
               Omega_LL1(:,:,k,n_lorentz+1) = Omega_LL1(:,:,k,n_lorentz+1) + Omega_LL1(:,:,k,p)
               Omega_LR1(:,:,k,n_lorentz+1) = Omega_LR1(:,:,k,n_lorentz+1) + Omega_LR1(:,:,k,p)               
           END DO
        END DO 

        DO k=1,n_lorentz
           DO p=1,kL-n_lorentz
              Omega_LL2(:,:,k,kL-n_lorentz+1) = Omega_LL2(:,:,k,kL-n_lorentz+1) + Omega_LL2(:,:,k,p)
              Omega_LR2(:,:,k,kL-n_lorentz+1) = Omega_LR2(:,:,k,kL-n_lorentz+1) + Omega_LR2(:,:,k,p)
           END DO 
        END DO 

        DO k=1,kL-n_lorentz
           DO p=1,n_lorentz
              Omega_LL3(:,:,k,n_lorentz+1) = Omega_LL3(:,:,k,n_lorentz+1) + Omega_LL3(:,:,k,p)
              Omega_LR3(:,:,k,n_lorentz+1) = Omega_LR3(:,:,k,n_lorentz+1) + Omega_LR3(:,:,k,p)
           END DO
        END DO

       DO k=1,n_lorentz
          DO p=1,n_lorentz
             Omega_RR1(:,:,k,n_lorentz+1) = Omega_RR1(:,:,k,n_lorentz+1) + Omega_RR1(:,:,k,p)
             Omega_RL1(:,:,k,n_lorentz+1) = Omega_RL1(:,:,k,n_lorentz+1) + Omega_RL1(:,:,k,p)
          END DO
       END DO

       DO k=1,n_lorentz
          DO p=1,kR-n_lorentz
             Omega_RR2(:,:,k,kR-n_lorentz+1) = Omega_RR2(:,:,k,kR-n_lorentz+1) + Omega_RR2(:,:,k,p)
             Omega_RL2(:,:,k,kR-n_lorentz+1) = Omega_RL2(:,:,k,kR-n_lorentz+1) + Omega_RL2(:,:,k,p)
          END DO 
       END DO    

       DO k=1,kR-n_lorentz
          DO p=1,n_lorentz
             Omega_RR3(:,:,k,n_lorentz+1) = Omega_RR3(:,:,k,n_lorentz+1) + Omega_RR3(:,:,k,p)
             Omega_RL3(:,:,k,n_lorentz+1) = Omega_RL3(:,:,k,n_lorentz+1) + Omega_RL3(:,:,k,p)
          END DO
       END DO 

!>>>>>>>>>>>>>>>>>>>>

        DO k=1,n_lorentz
             DO p=1,n_lorentz                      
                  Omega_LL1(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p) - Gam_lesser_L_m(:,:,p) ) , Pi_L(:,:,k) )   &
                                     + matmul( dconjg(transpose(Pi_L(:,:,p))) , ( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                                     - im * 1.d0/hbar * (hi_left_m(p) - hi_left_p(k)) * Omega_LL1(:,:,k,p)
             
                  Omega_LR1(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p) - Gam_lesser_R_m(:,:,p) ) , Pi_L(:,:,k) )  &
                                     + matmul( dconjg(transpose(Pi_R(:,:,p))) , ( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                                     - im * 1.d0/hbar * (hi_right_m(p) - hi_left_p(k)) * Omega_LR1(:,:,k,p) 
             END DO
        END DO

        DO k=1,n_lorentz
             DO p=1,kL-n_lorentz                      
                  Omega_LL2(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p+n_lorentz) - Gam_lesser_L_m(:,:,p+n_lorentz) ) , Pi_L(:,:,k) )   &
                                     + matmul( dconjg(transpose(Pi_L(:,:,p+n_lorentz))) , ( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                                     - im * 1.d0/hbar * (hi_left_m(p+n_lorentz) - hi_left_p(k)) * Omega_LL2(:,:,k,p)
             
                  Omega_LR2(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p+n_lorentz) - Gam_lesser_R_m(:,:,p+n_lorentz) ) , Pi_L(:,:,k) )  &
                                     + matmul( dconjg(transpose(Pi_R(:,:,p+n_lorentz))) , ( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                                     - im * 1.d0/hbar * (hi_right_m(p+n_lorentz) - hi_left_p(k)) * Omega_LR2(:,:,k,p)
             END DO
        END DO 

        DO k=1,kL-n_lorentz
             DO p=1,n_lorentz                      
                  Omega_LL3(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p) - Gam_lesser_L_m(:,:,p) ) , Pi_L(:,:,k+n_lorentz) )   &
                                     + matmul( dconjg(transpose(Pi_L(:,:,p))) , ( Gam_greater_L_p(:,:,k+n_lorentz) - Gam_lesser_L_p(:,:,k+n_lorentz) ) )  &
                                     - im * 1.d0/hbar * (hi_left_m(p) - hi_left_p(k+n_lorentz)) * Omega_LL3(:,:,k,p)
             
                  Omega_LR3(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p) - Gam_lesser_R_m(:,:,p) ) , Pi_L(:,:,k+n_lorentz) )  &
                                     + matmul( dconjg(transpose(Pi_R(:,:,p))) , ( Gam_greater_L_p(:,:,k+n_lorentz) - Gam_lesser_L_p(:,:,k+n_lorentz) ) )  &
                                     - im * 1.d0/hbar * (hi_right_m(p) - hi_left_p(k+n_lorentz)) * Omega_LR3(:,:,k,p)
             END DO
        END DO 

       DO k=1,n_lorentz
            DO p=1,n_lorentz                      
                 Omega_RR1(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p) - Gam_lesser_R_m(:,:,p) ) , Pi_R(:,:,k) )   &
                                    + matmul( dconjg(transpose(Pi_R(:,:,p))) , ( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                                    - im * 1.d0/hbar * (hi_right_m(p) - hi_right_p(k)) * Omega_RR1(:,:,k,p) 

                  
                 Omega_RL1(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p) - Gam_lesser_L_m(:,:,p) ) , Pi_R(:,:,k) )   &
                                    + matmul( dconjg(transpose(Pi_L(:,:,p))) , ( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                                    - im * 1.d0/hbar * (hi_left_m(p) - hi_right_p(k)) * Omega_RL1(:,:,k,p) 
            END DO
       END DO         

       DO k=1,n_lorentz
            DO p=1,kR-n_lorentz                      
                 Omega_RR2(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p+n_lorentz) - Gam_lesser_R_m(:,:,p+n_lorentz) ) , Pi_R(:,:,k) )   &
                                    + matmul( dconjg(transpose(Pi_R(:,:,p+n_lorentz))) , ( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                                    - im * 1.d0/hbar * (hi_right_m(p+n_lorentz) - hi_right_p(k)) * Omega_RR2(:,:,k,p)

                  
                 Omega_RL2(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p+n_lorentz) - Gam_lesser_L_m(:,:,p+n_lorentz) ) , Pi_R(:,:,k) )   &
                                    + matmul( dconjg(transpose(Pi_L(:,:,p+n_lorentz))) , ( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                                    - im * 1.d0/hbar * (hi_left_m(p+n_lorentz) - hi_right_p(k)) * Omega_RL2(:,:,k,p) 
            END DO
       END DO
        
       DO k=1,kR-n_lorentz
            DO p=1,n_lorentz                      
                 Omega_RR3(:,:,k,p) = matmul( ( Gam_greater_R_m(:,:,p) - Gam_lesser_R_m(:,:,p) ) , Pi_R(:,:,k+n_lorentz) )   &
                                    + matmul( dconjg(transpose(Pi_R(:,:,p))) , ( Gam_greater_R_p(:,:,k+n_lorentz) - Gam_lesser_R_p(:,:,k+n_lorentz) ) )  &
                                    - im * 1.d0/hbar * (hi_right_m(p) - hi_right_p(k+n_lorentz)) * Omega_RR3(:,:,k,p)

                  
                 Omega_RL3(:,:,k,p) = matmul( ( Gam_greater_L_m(:,:,p) - Gam_lesser_L_m(:,:,p) ) , Pi_R(:,:,k+n_lorentz) )   &
                                    + matmul( dconjg(transpose(Pi_L(:,:,p))) , ( Gam_greater_R_p(:,:,k+n_lorentz) - Gam_lesser_R_p(:,:,k+n_lorentz) ) )  &
                                    - im * 1.d0/hbar * (hi_left_m(p) - hi_right_p(k+n_lorentz)) * Omega_RL3(:,:,k,p)
            END DO
       END DO 

!>>>>>>>>>>>>>>>>>>>>

        DO k=1,n_lorentz
              Pi_L(:,:,k) = - im * Gam_lesser_L_p(:,:,k) - im * matmul( rho(:,:),( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                            - im * 1.d0/hbar * matmul( H(:,:) , Pi_L(:,:,k) ) + im * 1.d0/hbar * hi_left_p(k) * Pi_L(:,:,k) &  
                            - im * 1.d0/(hbar) * Omega_LL1(:,:,k,n_lorentz+1) - im * 1.d0/(hbar*hbar) * Omega_LL2(:,:,k,kL-n_lorentz+1) &
                            - im * 1.d0/(hbar) * Omega_LR1(:,:,k,n_lorentz+1) - im * 1.d0/(hbar*hbar) * Omega_LR2(:,:,k,kL-n_lorentz+1) 
        END DO
 
        DO k=n_lorentz+1,kL
              Pi_L(:,:,k) = - im * Gam_lesser_L_p(:,:,k) - im * matmul( rho(:,:),( Gam_greater_L_p(:,:,k) - Gam_lesser_L_p(:,:,k) ) )  &
                            - im * 1.d0/hbar * matmul( H(:,:) , Pi_L(:,:,k) ) + im * 1.d0/hbar * hi_left_p(k) * Pi_L(:,:,k) &
                            - im * 1.d0/(hbar) * Omega_LL3(:,:,k-n_lorentz,n_lorentz+1) & 
                            - im * 1.d0/(hbar) * Omega_LR3(:,:,k-n_lorentz,n_lorentz+1)
        END DO 
        
        DO k=1,n_lorentz
              Pi_R(:,:,k) = - im * Gam_lesser_R_p(:,:,k) - im * matmul( rho(:,:),( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                            - im * 1.d0/hbar * matmul( H(:,:) , Pi_R(:,:,k) ) + im * 1.d0/hbar * hi_right_p(k) * Pi_R(:,:,k) &
                            - im * 1.d0/(hbar) * Omega_RR1(:,:,k,n_lorentz+1) - im * 1.d0/(hbar*hbar) * Omega_RR2(:,:,k,kR-n_lorentz+1) &
                            - im * 1.d0/(hbar) * Omega_RL1(:,:,k,n_lorentz+1) - im * 1.d0/(hbar*hbar) * Omega_RL2(:,:,k,kR-n_lorentz+1) 
        END DO        

        DO k=n_lorentz+1,kR
              Pi_R(:,:,k) = - im * Gam_lesser_R_p(:,:,k) - im * matmul( rho(:,:),( Gam_greater_R_p(:,:,k) - Gam_lesser_R_p(:,:,k) ) )  &
                            - im * 1.d0/hbar * matmul( H(:,:) , Pi_R(:,:,k) ) + im * 1.d0/hbar * hi_right_p(k) * Pi_R(:,:,k) &
                            - im * 1.d0/(hbar) * Omega_RR3(:,:,k-n_lorentz,n_lorentz+1) &
                            - im * 1.d0/(hbar) * Omega_RL3(:,:,k-n_lorentz,n_lorentz+1)
        END DO

!>>>>>>>>>>>>>>>>>>>>

        rho(:,:) = - im * 1.d0/hbar * ( matmul(H(:,:),rho(:,:)) - matmul(rho(:,:),H(:,:)) )&
                   + 1.d0/(hbar) * ( Pi_L(:,:,kL+1) + dconjg(transpose(Pi_L(:,:,kL+1))) )&
                   + 1.d0/(hbar) * ( Pi_R(:,:,kR+1) + dconjg(transpose(Pi_R(:,:,kR+1))) )        

        call negf_matrix_to_vector(dvec)
!                             , length_rkvec, Pi_L, Pi_R, rho, &
!                             Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                             Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                             n, n_lorentz, kL, kR)
END SUBROUTINE negf_eom


! **************************** rk_init ****************************************
!                         initializes the Runge-Kutta integration with rksuite,
!                         necessary variables WORK, THRES and YMAX are allocated here,
!                         initializing subroutines ENVIRN and SETUP are called here
! *****************************************************************************

SUBROUTINE negf_rk_init
        use global
        USE rksuite_vars
        implicit none
!       double complex   :: rkvec(length_rkvec)
!---
        integer         :: alloc_status

!       NEQ = 9*nn*nn*max(N_Lorentzians_L+N_nu_k_Left+1,N_Lorentzians_R+N_nu_k_right+1)*2        
        NEQ = length_rkvec*2
        LENWRK = 16*NEQ
!       write(*,*) 'NEQ:',NEQ

!       *** allocation of missing arrays for rksuite, "missing" since they could not be defined in the MODULE rksuite_vars, because n (and nn) is not a parameter but acquired during the program
!       check_alloc = 0
        allocate(WORK(LENWRK), stat = alloc_status)
!       if (alloc_status .NE. 0) check_alloc = check_alloc+1
        if (alloc_status .NE. 0) write(*,*) 'Error in allocation WORK of RK Variables (negf_rk_init) '

!       check_alloc=0
        allocate(THRES(NEQ), stat = alloc_status)
!       if (alloc_status .NE. 0) check_alloc = check_alloc+1
        if (alloc_status .NE. 0) write(*,*) 'Error in allocation THRES of RK Variables (negf_rk_init) '
        
!       check_alloc=0
        allocate(YMAX(NEQ), stat = alloc_status)
!       if (alloc_status .NE. 0) check_alloc = check_alloc+1
        if (alloc_status .NE. 0) write(*,*) 'Error in allocation YMAX of RK Variables (negf_rk_init) '
!       if (check_alloc  .EQ. 0) write(*,*) 'rk allocation ok'

!       TSTART = 0.d0
!       TEND   = t_step / 41.341373 * 1.1   ! a.u. of time to femtosecond, plus buffer
        THRES  = TOL
        CALL ENVIRN(OUTCH,MCHPES,DWARF)

!       Do the setup rather before each time step?
!       CALL SETUP(NEQ,TSTART,rkvec,TEND,TOL,THRES,METHOD,TASK,ERRASS,TSTART,WORK,LENWRK,MESSAGE)

end subroutine negf_rk_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interface to RKsuite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE negf_solve_with_rk(TWANT, rkvec, drkvec)
        use global
        USE rksuite_vars
        implicit none
        double precision :: TWANT
        double complex   :: rkvec(length_rkvec), drkvec(length_rkvec)
        external negf_eom
!
        call UT(negf_eom,TWANT,TNOW,rkvec,drkvec,YMAX,WORK,ifail)
!
        if (ifail > 4) write(*,*) "rk error: ", ifail
!
!       *** on soft errors:
!       *** repeat until TWANT is reached (TNOW=TWANT)
        DO WHILE ( (IFAIL.eq.2 .or. IFAIL.eq.3 .or. IFAIL.eq.4).and. TNOW .ne. TWANT )
                IFAIL = -1
                call UT(negf_eom,TWANT,TNOW,rkvec,drkvec,YMAX,WORK,ifail)
                if (ifail > 4) write(*,*) "rk error: ", ifail
        END DO
        
END SUBROUTINE negf_solve_with_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!
! obtain the result
!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE negf_calculate_current(rkvec, curr, wf)
        USE global
        implicit none
        double complex               :: rkvec(length_rkvec)
        double precision             :: curr(3), wf(n)
        integer                      :: k

!       function
        double complex               :: negf_cdtrace

        call negf_vector_to_matrix(rkvec)
!                             , length_rkvec, Pi_L, Pi_R, rho, &
!                             Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                             Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                             n, n_lorentz, kL, kR)

        Pi_L(:,:,kL+1) = dcmplx(0.d0,0.d0)
        Pi_R(:,:,kR+1) = dcmplx(0.d0,0.d0)        
            
        DO k=1, kL
             Pi_L(:,:,kL+1) = Pi_L(:,:,kL+1) + Pi_L(:,:,k)
        END DO
         
        DO k=1, kR
             Pi_R(:,:,kR+1) = Pi_R(:,:,kR+1) + Pi_R(:,:,k)
        END DO                    
         
        ! curr_L
        curr(1) = 2. / hbar_r * real(negf_cdtrace(Pi_L(:,:,kL+1),n))
        ! curr_R
        curr(2) = 2. / hbar_r * real(negf_cdtrace(Pi_R(:,:,kR+1),n))
        ! total current
        curr(3) = 0.5d0 * (curr(1) - curr(2))
       
        curr(1) = 2.43413479624d-4 * curr(1)            !Current in Amps i.e. multiply with ( e [C] / hbar [fs] )
        curr(2) = 2.43413479624d-4 * curr(2)
        curr(3) = 2.43413479624d-4 * curr(3)

!       write (*,*) 'current', curr(1), curr(2), curr(3)

!       store the instantaneous density matrix
        do k=1, n
             if (RealPart(rho(k,k)) .gt. 0.) then
                 wf(k) = dsqrt(RealPart(rho(k,k)))
             else
                 wf(k) = 0.
             end if
        end do

!       write (*,*) 'density matrix', rho(1,1), rho(2,2), rho(3,3), rho(4,4)

END SUBROUTINE negf_calculate_current

!!!!!!!!!!!!!!!!!!!
! obsoleted...
!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE save_output(timestep, file)
!        USE global
!        implicit none
!        double precision, intent(in) :: timestep
!        double precision             :: curr_L, curr_R, curr
!        integer                      :: i, j, k, l, m, file
!        integer                      :: alloc_st, dimn,counter
!
!        !double precision, allocatable, dimension(:):: current_int
!        double complex :: negf_cdtrace
!
!        dimn = (n-1)*n/2
!
!        !allocate(current_int(dimn), stat= alloc_st)
!        !if (alloc_st.NE.0) write(*,*) "Error allocation"
!
!        call negf_vector_to_matrix(rkvec)
!!                             , length_rkvec, Pi_L, Pi_R, rho, &
!!                             Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!!                             Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!!                             n, n_lorentz, kL, kR)
!
!        Pi_L(:,:,kL+1) = dcmplx(0.d0,0.d0)
!        Pi_R(:,:,kR+1) = dcmplx(0.d0,0.d0)        
!            
!        DO k=1, kL
!             Pi_L(:,:,kL+1) = Pi_L(:,:,kL+1) + Pi_L(:,:,k)
!        END DO
!         
!        DO k=1, kR
!             Pi_R(:,:,kR+1) = Pi_R(:,:,kR+1) + Pi_R(:,:,k)
!        END DO                    
!         
!        
!        !curr_L=0.d0
!        curr_L= 2*(1/hbar)*real(negf_cdtrace(Pi_L(:,:,kL+1),n))
!        !curr_R=0.d0
!        curr_R= 2*(1/hbar)*real(negf_cdtrace(Pi_R(:,:,kR+1),n))
!        curr=0.5d0*(curr_L-curr_R)
!       
!        !DO i=1,dimn
!        !     current_int(i) = 0.d0
!        !END DO 
!
!        !IF (n .GT. 1) THEN
!        !  counter=0
!        !  DO i=1,n-1
!        !     DO j=i+1,n
!        !        counter=counter+1
!        !        current_int(counter) = im * 2.43413479624d-4 * abs(delta(i,j)) * (rho(j,i) - rho(i,j))
!        !     END DO
!        !  END DO  
!        !END IF
!
!        curr   = 2.43413479624d-4 * curr             !Current in Amps i.e. multiply with ( e [C] / hbar [fs] )
!        curr_L = 2.43413479624d-4 * curr_L
!        curr_R = 2.43413479624d-4 * curr_R       
!
!        write(file,'(f9.0,t16,15es16.8)') timestep, (real(rho(i,i)),i=1,n), curr, curr_L, curr_R  !, (current_int(j),j=1,dimn)
!
!END SUBROUTINE save_output


SUBROUTINE negf_matrix_to_vector(dvector)
!                           , length_rkvec, Pi_L, Pi_R, rho, &
!                           Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                           Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                           n, n_lorentz, kL, kR)
        use global
        implicit none
        double complex, dimension(length_rkvec) :: dvector
!       integer                                 :: length_rkvec, n, n_lorentz, kL, kR
!       double complex                          :: rho(n,n), Pi_L(n,n,kL+1), Pi_R(n,n,kR+1)
!       double complex                          :: Omega_LL1(n,n,n_lorentz,n_lorentz+1), Omega_LL2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LL3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_LR1(n,n,n_lorentz,n_lorentz+1), Omega_LR2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LR3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RR1(n,n,n_lorentz,n_lorentz+1), Omega_RR2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RR3(n,n,kR-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RL1(n,n,n_lorentz,n_lorentz+1), Omega_RL2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RL3(n,n,kR-n_lorentz,n_lorentz+1)
!---
        integer:: j,k,l,m,p

        m=1

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_LL1(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO 

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,kL-n_lorentz 
                        dvector(m) = Omega_LL2(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,kL-n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_LL3(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_LR1(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO 

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,kL-n_lorentz 
                        dvector(m) = Omega_LR2(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,kL-n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_LR3(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_RR1(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO 

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,kR-n_lorentz 
                        dvector(m) = Omega_RR2(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,kR-n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_RR3(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_RL1(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO 

       DO l=1,n
          DO j=1,n
                DO k=1,n_lorentz
                     DO p=1,kR-n_lorentz 
                        dvector(m) = Omega_RL2(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
       END DO

       DO l=1,n
          DO j=1,n
                DO k=1,kR-n_lorentz
                     DO p=1,n_lorentz 
                        dvector(m) = Omega_RL3(l,j,k,p)
                        m=m+1
                     END DO
               END DO
          END DO 
        END DO
 
        DO l=1,n
           DO j=1,n
                DO k=1,kL
                        dvector(m) = Pi_L(l,j,k)
                        m=m+1
                END DO
            END DO
        END DO                

        DO l=1,n
           DO j=1,n
                DO k=1,kR
                        dvector(m) = Pi_R(l,j,k)
                        m=m+1
                END DO
            END DO
        END DO
               
        DO l=1,n
           DO j=1,n
              dvector(m) = rho(l,j)
              m=m+1
           END DO
        END DO

return
END SUBROUTINE negf_matrix_to_vector

SUBROUTINE negf_vector_to_matrix(vector)
!                           , length_rkvec, Pi_L, Pi_R, rho, &
!                           Omega_LL1, Omega_LL2, Omega_LL3, Omega_LR1, Omega_LR2, Omega_LR3, &
!                           Omega_RL1, Omega_RL2, Omega_RL3, Omega_RR1, Omega_RR2, Omega_RR3, &
!                           n, n_lorentz, kL, kR)
        use global
        implicit none
        double complex, dimension(length_rkvec) :: vector
!       integer                                 :: length_rkvec, n, n_lorentz, kL, kR
!       double complex                          :: rho(n,n), Pi_L(n,n,kL+1), Pi_R(n,n,kR+1)
!       double complex                          :: Omega_LL1(n,n,n_lorentz,n_lorentz+1), Omega_LL2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LL3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_LR1(n,n,n_lorentz,n_lorentz+1), Omega_LR2(n,n,n_lorentz,kL+1-n_lorentz), Omega_LR3(n,n,kL-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RR1(n,n,n_lorentz,n_lorentz+1), Omega_RR2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RR3(n,n,kR-n_lorentz,n_lorentz+1)
!       double complex                          :: Omega_RL1(n,n,n_lorentz,n_lorentz+1), Omega_RL2(n,n,n_lorentz,kR+1-n_lorentz), Omega_RL3(n,n,kR-n_lorentz,n_lorentz+1)
!---
        integer:: j,k,l,p,m
        
        m=1

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,n_lorentz
                       Omega_LL1(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,kL-n_lorentz
                       Omega_LL2(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,kL-n_lorentz
                    DO p=1,n_lorentz
                       Omega_LL3(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,n_lorentz
                       Omega_LR1(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,kL-n_lorentz
                       Omega_LR2(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,kL-n_lorentz
                    DO p=1,n_lorentz
                       Omega_LR3(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,n_lorentz
                       Omega_RR1(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,kR-n_lorentz
                       Omega_RR2(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,kR-n_lorentz
                    DO p=1,n_lorentz
                       Omega_RR3(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,n_lorentz
                       Omega_RL1(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,n_lorentz
                    DO p=1,kR-n_lorentz
                       Omega_RL2(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
                DO k=1,kR-n_lorentz
                    DO p=1,n_lorentz
                       Omega_RL3(l,j,k,p) = vector(m)
                       m=m+1
                    END DO 
                END DO 
           END DO
        END DO

        DO l=1,n
           DO j=1,n
                DO k=1,kL
                        Pi_L(l,j,k) = vector(m)
                        m=m+1
                END DO
           END DO
        END DO


        DO l=1,n
           DO j=1,n
                DO k=1,kR
                        Pi_R(l,j,k) = vector(m)
                        m=m+1
                END DO
           END DO
        END DO 

        DO l=1,n
           DO j=1,n
              rho(l,j)=vector(m)
              m=m+1
           END DO 
        END DO 

return
END SUBROUTINE negf_vector_to_matrix


!********************************subroutines for constructing the poles of the Ozaki sum decomposition of the Fermi function********************************


SUBROUTINE negf_const_mat(Mat_L, Mat_R, N_poles_L, N_poles_R)
! This subroutine constructs the matrix used to compute the poles of the Ozaki sum decomposition of the Fermi function 
     implicit none
     double precision  :: Mat_L(2*N_poles_L,2*N_poles_L), Mat_R(2*N_poles_R,2*N_poles_R)
     integer           :: N_poles_L, N_poles_R 
! ---
     double precision  :: x
     integer           :: i, j
   
     do i=1,2*N_poles_L
       x=i
       do j=1,2*N_poles_L           
           if(i .EQ. j) then
           Mat_L(i,j)=0d0
           end if 
           if(i .EQ. j-1) then
           Mat_L(i,j)=1/(2*sqrt(4*x**2-1)) 
           end if 
           Mat_L(j,i)=Mat_L(i,j)
       end do
     end do

     do i=1,2*N_poles_R
       x=i
       do j=1,2*N_poles_R           
           if(i .EQ. j) then
           Mat_R(i,j)=0d0
           end if 
           if(i .EQ. j-1) then
           Mat_R(i,j)=1/(2*sqrt(4*x**2-1)) 
           end if 
           Mat_R(j,i)=Mat_R(i,j)
       end do
     end do

END SUBROUTINE negf_const_mat


SUBROUTINE negf_construct_nu(Eig_val_mat_L, Nu_L, R_L, R_alpha_L, N_poles_L, Eig_val_mat_R, Nu_R, R_R, R_alpha_R, N_poles_R)
        implicit none
        double precision :: Eig_val_mat_L(2*N_poles_L), Nu_L(N_poles_L), R_L(N_poles_L), R_alpha_L(2*N_poles_L)
        double precision :: Eig_val_mat_R(2*N_poles_R), Nu_R(N_poles_R), R_R(N_poles_R), R_alpha_R(2*N_poles_R)
        integer          :: N_poles_L, N_poles_R
! ---
        integer          :: i,j

     j=1
     DO i=1,2*N_poles_L
       IF(Eig_val_mat_L(i) .GT. 0) THEN
         Nu_L(j)=Eig_val_mat_L(i)
         R_L(j)=R_alpha_L(i)
         j=j+1
       END IF
     END DO

     j=1
     DO i=1,2*N_poles_R
       IF(Eig_val_mat_R(i) .GT. 0) THEN
         Nu_R(j)=Eig_val_mat_R(i)
         R_R(j)=R_alpha_R(i)
         j=j+1
       END IF
     END DO
 
END SUBROUTINE negf_construct_nu


subroutine negf_calc_r(R_alpha_L, Eig_val_mat_L, Eig_vect_mat_L, N_poles_L, R_alpha_R, Eig_val_mat_R, Eig_vect_mat_R, N_poles_R)
     implicit none

     double precision :: R_alpha_L(2*N_poles_L), Eig_val_mat_L(2*N_poles_L), Eig_vect_mat_L(2*N_poles_L,2*N_poles_L)
     double precision :: R_alpha_R(2*N_poles_R), Eig_val_mat_R(2*N_poles_R), Eig_vect_mat_R(2*N_poles_R,2*N_poles_R)
     integer          :: N_poles_L, N_poles_R
! ---
     integer          :: i
    
   do i=1,2*N_poles_L
     R_alpha_L(i) = Eig_vect_mat_L(1,i)**2 / (4 * Eig_val_mat_L(i)**2)
   end do
  
   do i=1,2*N_poles_R
     R_alpha_R(i) = Eig_vect_mat_R(1,i)**2 / (4 * Eig_val_mat_R(i)**2)
   end do  
  
end subroutine negf_calc_r


SUBROUTINE negf_jacobi(A,D,V,M)
!This subroutine calculates the eigenvalues and eigenvectors of a real symmetric square matrix A(N,N). 
!D(N) returns the eigenvalues of matrix A. V(N,N) contains the eigenvectors of A by columns. 

integer                      ::  ialloc, check_alloc
integer                      ::  M,N,NROT

double precision             ::  A(1:2*M,1:2*M),V(1:2*M,1:2*M),D(1:2*M)

double precision, pointer    ::  B(:), Z(:)
double precision             ::  c,g,h,s,sm,t,tau,theta,tresh
integer                      ::  i,j,iq,ip

 check_alloc=0
 N=2*M

 allocate(B(1:N),stat=ialloc)
 if (ialloc .NE.0) check_alloc=check_alloc+1
 allocate(Z(1:N),stat=ialloc)
 if (ialloc .NE.0) check_alloc=check_alloc+1
 if(check_alloc .NE.0) write(*,*) 'Error allocation'

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0 
    end do
      V(ip,ip)=1.d0
  end do  
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0    
  end do
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      end do
    end do
    if(sm==0.d0) return  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    end if
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
                .and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
                  A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
          h=D(iq)-D(ip)
          if(DABS(h)+g.eq.DABS(h)) then
            t=A(ip,iq)/h
          else
            theta=0.5d0*h/A(ip,iq)  
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
            if(theta.lt.0.d0) t=-t
          end if
          c=1.d0/DSQRT(1.d0+t**2)
          s=t*c
          tau=s/(1.d0+c)
          h=t*A(ip,iq)
          Z(ip)=Z(ip)-h
          Z(iq)=Z(iq)+h
          D(ip)=D(ip)-h
          D(iq)=D(iq)+h
          A(ip,iq)=0.d0
          do j=1, ip-1
            g=A(j,ip)
            h=A(j,iq)
            A(j,ip)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          end do
          do j=ip+1, iq-1
            g=A(ip,j)
            h=A(j,iq)
            A(ip,j)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          end do
          do j=iq+1, N
            g=A(ip,j)
            h=A(iq,j)
            A(ip,j)=g-s*(h+g*tau)
            A(iq,j)=h+s*(g-h*tau)
          end do
          do j=1, N
            g=V(j,ip)
            h=V(j,iq)
            V(j,ip)=g-s*(h+g*tau)
            V(j,iq)=h+s*(g-h*tau)
          end do
          NROT=NROT+1
        end if !if ((i.gt.4)...
      end do !main iq loop
    end do !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    end do
  end do !main i loop
!  pause ' 50 iterations !'
  return

END SUBROUTINE negf_jacobi

!**********************************************************************************************************************************************************

DOUBLE COMPLEX FUNCTION negf_fermi(energy, beta)
        implicit none
        double complex   :: energy
        double precision :: beta

        negf_fermi=dcmplx(1.d0,0.d0) / ( dcmplx(1.d0,0.d0) + cdexp( dcmplx(beta, 0.d0) * energy ) ) 

END FUNCTION negf_fermi


DOUBLE COMPLEX FUNCTION negf_spectraldensity(w, gam, w_0, eps)
        implicit none
        double complex   :: w
        double precision :: gam, w_0, eps 
        
        negf_spectraldensity = (dcmplx(w_0,0d0))**2 * dcmplx(gam,0d0) / (  (w-dcmplx(eps,0d0))**2 + (dcmplx(w_0,0d0))**2  )
         
END FUNCTION negf_spectraldensity


DOUBLE COMPLEX FUNCTION negf_cdtrace(A,B)
  implicit none
  integer, intent(IN)          :: B
  double complex, intent(IN)   :: A(B,B)
  integer                      :: i

  negf_cdtrace = dcmplx(0.d0,0.d0)
  do i=1,B
     negf_cdtrace = negf_cdtrace + A(i,i)
  end do
  return
END FUNCTION negf_cdtrace


!SUBROUTINE display_parameters(file)
!        USE global
!        implicit none
!        integer::i,file
!        write(file,'(t1,a,t20,i6,t30)') '%#NUMBER OF SITES:',n,'%#without spin'
!        write(file,'(t1,a,t20,f6.1)') '%#E_F_left:', E_F_left
!        write(file,'(t1,a,t20,f6.1)') '%#E_F_right:', E_F_right
!        write(file,'(t1,a,t20,10f6.1)') '%#E(i):', (E(i),i=1,n)
!        !write(file,'(t1,a,t20,10f6.3)') '%#delta(i)', (delta(i),i=1,n-1)        
!        write(file,'(t1,a,t20,2f6.2)') '%#Temp', Temp
!        write(file,'(t1,a,t20,2i6)')  '%#N_nu_k L/R:', N_poles_L, N_poles_R
!        write(file,'(t1,a,t20,f6.1)') '%#t_0:', t_0
!        write(file,'(t1,a,t20,f8.1)') '%#t_end:', t_end
!        write(file,'(t1,a,t20,f6.0)') '%#t_step:', t_step
!        write(file,'(t1,a,t20,f6.1)') '%#gam_L:', gam_L
!        write(file,'(t1,a,t20,f6.1)') '%#w0_L:',  w0_L 
!        write(file,'(t1,a,t20,f6.1)') '%#eps_L:', eps_L
!        write(file,'(t1,a,t20,f6.1)') '%#gam_R:', gam_R
!        write(file,'(t1,a,t20,f6.1)') '%#w0_R:',  w0_R
!        write(file,'(t1,a,t20,f6.1)') '%#eps_R:', eps_R
!        write(file,*)
!END SUBROUTINE display_parameters

