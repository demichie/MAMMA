!*********************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!*********************************************************************

MODULE init

  USE parameters, ONLY : n_vars , n_cry , n_gas

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

    IMPLICIT none

    REAL*8, INTENT(IN) :: u_0
    REAL*8, INTENT(OUT) :: qp(n_vars)

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

    INTEGER :: i,idx,iter,max_iter
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
       
       CALL f_beta_eq

       beta_in(1:n_cry) = REAL(beta_eq(1:n_cry))
       
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
          
          CALL f_alfa( x_ex_dis_in(1:n_gas) , xd_md_in(1:n_gas) , beta_in(1:n_cry) ,  &
               r_rho_md , r_rho_2 , alfa_g_in(1:n_gas) )
          
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
    
    ! define the vector of primitive variables
    
    idx = 0
    
    DO i = 1,n_gas
       
       idx = idx + 1
       
       qp(idx) = alfa_g_in(i)
       
    END DO
    
    ! Frist phase pressure
    idx = idx + 1
    
    qp(idx) = p1_in
    
    ! Second phase pressure
    idx = idx + 1
    
    qp(idx) = p2_in

    ! Frist phase velocity
    idx = idx + 1

    qp(idx) = r_u_1

    ! Second phase velocity
    idx = idx + 1

    qp(idx) = r_u_2

    ! Temperature
    idx = idx + 1
    qp(idx) = T_in

    ! Crystal volume fractions
    
    DO i = 1,n_cry

       idx = idx + 1

       qp(idx) = beta_in(i)

    END DO

    ! Dissolved gas mass fractions

    DO i = 1,n_gas

       idx = idx + 1

       qp(idx) = xd_md_in(i)

    END DO

  END SUBROUTINE init_steady

END MODULE init
