!*********************************************************************
!> \brief Initial solution
!> @author 
!> Mattia de' Michieli Vitturi
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!> \date 04/02/2017
!*********************************************************************

MODULE init

  USE parameters, ONLY : method_of_moments_flag
  
  USE parameters, ONLY : n_vars , n_cry , n_gas , n_mom

  USE parameters, ONLY : idx_p1 , idx_p2 , idx_u1 , idx_u2 , idx_T ,            &
       idx_xd_first , idx_xd_last , idx_alfa_first , idx_alfa_last ,            &
       idx_beta_first , idx_beta_last, idx_components_first,                    &
       idx_components_last

  USE parameters, ONLY : idx_mix_mass_eqn , idx_vol1_eqn , idx_mix_mom_eqn ,    &
       idx_rel_vel_eqn , idx_mix_engy_eqn , idx_dis_gas_eqn_first ,             &
       idx_dis_gas_eqn_last , idx_ex_gas_eqn_first , idx_ex_gas_eqn_last ,      &
       idx_cry_eqn_first , idx_cry_eqn_last, idx_components_eqn_first,          &
       idx_components_eqn_last
  
  IMPLICIT none

  REAL*8 :: alfa1_in         !< Inlet first phase volume fraction
  REAL*8 :: p1_in            !< Inlet first phase pressure 
  REAL*8 :: p2_in            !< Inlet second phase pressure 
  REAL*8 :: delta_p_in       !< Inlet pressure difference ( p2-p1 )
  REAL*8 :: u1_in            !< Inlet first phase velocity
  REAL*8 :: u2_in            !< Inlet second phase velocity
  REAL*8 :: T_in             !< Inlet temperature
  REAL*8, ALLOCATABLE :: beta_in(:)   !< Inlet crystal volume fraction
  REAL*8, ALLOCATABLE :: xd_md_in(:)  !< Inlet dissolved gas mass fraction
  REAL*8 :: frag_eff_in      !< Inlet fragmentation efficiency

  REAL*8 :: p_out            !< Outlet pressure

CONTAINS

  !******************************************************************************
  !> \brief Steady problem initialization
  !
  !> This subroutine initialize the state for a steady solution problem, given 
  !> the inlet velocity u_0.
  !> \param[in]   u_0       inlet velocity
  !> \param[out]  qp        physical variables at the inlet
  !> \date 18/09/2012
  !******************************************************************************

  SUBROUTINE init_steady(u_0,qp)

    ! external procedures
    USE constitutive, ONLY : f_beta_eq , f_xdis_eq , eval_densities
    USE constitutive, ONLY : f_alfa3, f_alfa

    ! external variables
    USE constitutive, ONLY : p_1 , p_2 , T , x_d_md, x_g, alfa_g_2
    USE constitutive, ONLY : x_ex_dis_in , x_d_md_eq, beta , beta_eq 
    USE constitutive, ONLY : rho_2 , rho_md, rho_g, rho_1, rho_c
    USE constitutive, ONLY : bar_p_m, gamma_m, cv_m
    USE constitutive, ONLY : bar_p_c, gamma_c, cv_c
    USE constitutive, ONLY : L0_cry_in , L_nucleus_in , mom_cry , cry_shape_factor_in
    USE constitutive, ONLY : beta0
    USE parameters, ONLY : n_components
    USE constitutive, ONLY : cry_init_solid_solution, rhoB_components
    USE melts_fit_module, ONLY : wt_tot_0,  wt_components_init, wt_components_fit
 
    IMPLICIT none

    REAL*8, INTENT(IN) :: u_0
    REAL*8, INTENT(INOUT) :: qp(n_vars)

    REAL*8 :: r_rho_g(1:n_gas)   !> exsolved gas density
    REAL*8 :: r_rho_2   !> total exsolved gas density
    REAL*8 :: r_rho_1   !> total exsolved melt density
    REAL*8 :: r_rho_c(1:n_cry)   !> crystals density

    REAL*8 :: r_rho_md  !> dis.gas+melt density

    REAL*8 :: r_u_1     !> melt+crystal phase velocity
    REAL*8 :: r_u_2     !> exsolved gas velocity

    REAL*8 :: r_frag_eff

    REAL*8 :: alfa2_in
    REAL*8 :: alfa_g_in(1:n_gas)

    COMPLEX*16 :: x_g_old(1:n_gas)

    REAL*8 :: xd_md_tot

    REAL*8 :: xtot_in

    INTEGER :: i,j,k
    INTEGER :: idx,iter,max_iter
    REAL*8 :: error_iter

    p2_in = p1_in + delta_p_in

    ! evaluate the initial crystal volume fraction

    p_1 = DCMPLX( p1_in , 0.D0 ) 
    p_2 = DCMPLX( p2_in , 0.D0 ) 
    T = DCMPLX( T_in , 0.D0 )

    ! evaluate the initial dissolved gas mass fraction

    r_rho_md = REAL( ( p_1 + bar_p_m ) / ( T * cv_m * ( gamma_m - 1.D0) ) )
    rho_md = DCMPLX(r_rho_md,0.0)
    
    r_rho_c(1:n_cry) = REAL( ( p_1 + bar_p_c(1:n_cry) ) / ( T * cv_c(1:n_cry) * &
         ( gamma_c(1:n_cry) - DCMPLX(1.D0,0.D0) ) ) )    
    
    r_rho_1 = r_rho_md
    
    rho_c(1:n_cry) = DCMPLX( r_rho_c(1:n_cry) , 0.D0 )
    rho_1 = DCMPLX( r_rho_1 , 0.D0 )
    
    DO i=1,n_gas
       alfa_g_2(i) = DCMPLX(1.0 / n_gas,0.0);
       x_g(i) = DCMPLX(1e-2,0.0)
    END DO
    
    iter = 1
    max_iter =1000
    error_iter = 1.0
    
    DO WHILE( error_iter .GT. 1e-15 .AND. (iter .LT. max_iter) )
       
       x_g_old = x_g
       
       CALL f_xdis_eq
       
       DO i=1,n_gas
          
          xd_md_in(i) = MIN( REAL(x_d_md_eq(i)) , x_ex_dis_in(i) )
          
       END DO
       
       ! required by eval_densities
       x_d_md(1:n_gas) = DCMPLX( xd_md_in(1:n_gas) , 0.D0 )
       
       IF( method_of_moments_flag ) THEN
       
          beta_in(1:n_cry) = beta0(1:n_cry)

       ELSE

          CALL f_beta_eq

          beta_in(1:n_cry) = REAL(beta_eq(1:n_cry))
      
       END IF       
       
       r_u_1 = u_0
       r_u_2 = u_0 + 1.D-10
       
       ! evaluate the volume fractions of the two phases
       
       beta = DCMPLX( beta_in(1:n_cry) , 0.D0 )
       
       CALL eval_densities
       
       r_rho_1 = REAL( rho_1 )
       r_rho_2 = REAL( rho_2 )
       r_rho_g = REAL( rho_g )
       
       r_rho_md = REAL( rho_md )
       
       xd_md_tot = SUM( xd_md_in(1:n_gas) )
       
       xtot_in = SUM( x_ex_dis_in(1:n_gas) ) 
       
       IF ( n_gas .EQ. 1 ) THEN
          
          CALL f_alfa( x_ex_dis_in(1:n_gas) , xd_md_in(1:n_gas) ,               &
               beta_in(1:n_cry) , r_rho_md , r_rho_2 , alfa_g_in(1:n_gas) )
          
       ELSE
          
          CALL f_alfa3( p2_in, x_ex_dis_in(1:n_gas) , beta_in(1:n_cry) ,  &
               r_rho_md , r_rho_g(1:n_gas) , alfa_g_in(1:n_gas) )
          
       END IF
       
       
       DO i = 1,n_gas
          
          alfa_g_in(i) = MAX( alfa_g_in(i) , 1.D-10 )
          
       END DO
       
       alfa2_in = SUM( alfa_g_in(1:n_gas) )
       alfa1_in = 1.D0 - alfa2_in
       
       x_g = alfa_g_in * r_rho_g / ( alfa1_in * r_rho_1 + alfa2_in * r_rho_2 )
       
       alfa_g_2 = alfa_g_in / alfa2_in
       
       error_iter = (MAXVAL(ABS(x_g - x_g_old)))
       
    END DO
    
    IF (iter .GE. max_iter) THEN
       
       WRITE(*,*) 'No convergence for initial gas dissolve mass fraction'
       STOP
       
    END IF
    
    ! initialize the fragmentation efficiency to zero
    
    r_frag_eff = frag_eff_in

    ! -------- define the indexes of variables and equations --------------------

    idx_p1 = 1
    idx_p2 = 2
    idx_u1 = 3
    idx_u2 = 4
    idx_T = 5
    idx_xd_first = 5+1
    idx_xd_last = 5+n_gas
    idx_alfa_first = 5+n_gas+1
    idx_alfa_last = 5+2*n_gas
    idx_beta_first = 5+2*n_gas+1
    
    idx_mix_mass_eqn = 1
    idx_vol1_eqn = 2
    idx_mix_mom_eqn = 3
    idx_rel_vel_eqn = 4
    idx_mix_engy_eqn = 5
    idx_dis_gas_eqn_first = 5+1
    idx_dis_gas_eqn_last = 5+n_gas
    idx_ex_gas_eqn_first = 5+n_gas+1
    idx_ex_gas_eqn_last = 5+2*n_gas 
    idx_cry_eqn_first = 5+2*n_gas+1

    IF ( method_of_moments_flag ) THEN

       idx_cry_eqn_last = 5 + 2*n_gas + 2*n_cry*n_mom
       idx_beta_last = 5 + 2*n_gas + 2*n_cry*n_mom

       idx_components_eqn_first = 5 + 2*n_gas + 2*n_cry*n_mom + 1
       idx_components_first = 5 + 2*n_gas + 2*n_cry*n_mom + 1

       idx_components_eqn_last = 5 + 2*n_gas + 2*n_cry*n_mom + n_components
       idx_components_last = 5 + 2*n_gas + 2*n_cry*n_mom + n_components

    ELSE
       
       idx_cry_eqn_last = 5 + 2*n_gas + n_cry
       idx_beta_last = 5 + 2*n_gas + n_cry
       
    END IF
    
    ! --------- define the vector of primitive variables ------------------------
    
    ! Frist phase pressure
    qp(idx_p1) = p1_in
    
    ! Second phase pressure
    qp(idx_p2) = p2_in

    ! Frist phase velocity
    qp(idx_u1) = r_u_1

    ! Second phase velocity
    qp(idx_u2) = r_u_2

    ! Temperature
    qp(idx_T) = T_in

    ! Dissolved gas mass fractions
    idx = idx_xd_first
    DO i = 1,n_gas

       qp(idx) = xd_md_in(i)

       idx = idx + 1

    END DO

    ! Exsolved gas volume fraction
    idx = idx_alfa_first
    DO i = 1,n_gas
              
       qp(idx) = alfa_g_in(i)

       idx = idx + 1

    END DO
    
    ! Crystal volume fractions
    idx = idx_beta_first

    IF ( method_of_moments_flag ) THEN
 
       DO i = 1,n_cry

          DO j = 0,n_mom-1

             ! k is the index for microlith and phenocryst
             DO k = 1,2
  
                IF(k .EQ. 2) THEN

                    mom_cry(i,j,k) = beta_in(i) * alfa1_in / cry_shape_factor_in(i)    &
                         * L0_cry_in(i)**j / L0_cry_in(i)**3.0 !Revisar
                
                ELSE

                    mom_cry(i,j,k) = L_nucleus_in(i) * alfa1_in / cry_shape_factor_in(i)    &
                         * L_nucleus_in(i)**j  / L_nucleus_in(i)**3.0 !Revisar

                END IF 
                
                qp(idx) = mom_cry(i,j,k)
                
                idx = idx + 1
       
             END DO

          END DO

       END DO

       idx = idx_components_first

       DO i = 1, n_components

          wt_components_init(i) = wt_tot_0  * ( 1.0 - xtot_in ) * wt_components_fit(i)  

          DO j = 1, n_cry

             wt_components_init(i) = wt_components_init(i) - cry_init_solid_solution(i,j) *   &
                ( beta_in(j) * alfa1_in * REAL(rho_c(j)) / (alfa1_in * REAL(rho_1)  +         &
                (1.D0 - alfa1_in) * REAL(rho_2) ) )

          END DO

          IF( ( wt_components_init(i) ) .LT. 0.D0 ) THEN

             WRITE(*,*) 'Initial volume of crystals is not compatible with composition'

             STOP

          ENDIF

          rhoB_components(i) =  wt_components_init(i) * REAL( (rho_1) * alfa1_in + (rho_2) * (1.D0 - alfa1_in) )

          qp(idx) =  rhoB_components(i)

          idx = idx + 1

       END DO
       
    ELSE

       DO i = 1,n_cry

          qp(idx) = beta_in(i)
          
          idx = idx + 1
          
       END DO

    END IF

  END SUBROUTINE init_steady

END MODULE init
