!********************************************************************************
!> \brief Governing equations
!> @author 
!> Mattia de' Michieli Vitturi
!> This module contains the terms defining the governing equations of the
!> system. 
!> \date 04/02/2017
!********************************************************************************
MODULE equations

  USE constitutive
  
  USE geometry, ONLY : pi
  USE parameters, ONLY : verbose_level
  USE parameters, ONLY : method_of_moments_flag
  USE parameters, ONLY : n_eqns , n_vars, n_components
  USE parameters, ONLY : n_cry , n_gas , n_mom

  USE parameters, ONLY : idx_p1 , idx_p2 , idx_u1 , idx_u2 , idx_T ,            &
       idx_xd_first , idx_xd_last , idx_alfa_first , idx_alfa_last ,            &
       idx_beta_first , idx_beta_last
  
  USE parameters, ONLY : idx_mix_mass_eqn , idx_vol1_eqn , idx_mix_mom_eqn ,    &
       idx_rel_vel_eqn , idx_mix_engy_eqn , idx_dis_gas_eqn_first ,             &
       idx_dis_gas_eqn_last , idx_ex_gas_eqn_first , idx_ex_gas_eqn_last ,      &
       idx_cry_eqn_first , idx_cry_eqn_last
  
  USE melts_fit_module, ONLY: rel_cry_components 
  
  IMPLICIT none

  !> Flag for isothermal runs:\n
  !> - isothermal = .TRUE.  => constant temperature
  !> - isothermal = .FALSE. => variable temperature (an energy equation is solved)
  !> .
  LOGICAL :: isothermal

  !> Temperature for isothermal runs
  REAL*8 :: fixed_temp

  !> Input flag for lateral degassing:\n
  !> - lateral_degassing_flag = .TRUE.  => lateral degassing is considered
  !> - lateral_degassing_flag = .FALSE. => no lateral degassing
  !> .
  LOGICAL :: lateral_degassing_flag

  !> Flag for lateral degassing:\n
  !> - lateral_degassing = .TRUE.  => lateral degassing is active
  !> - lateral_degassing = .FALSE. => no lateral degassing
  !> .
  LOGICAL :: lateral_degassing

  !> Exsolved gas volume fraction threshold for lateral degassing 
  REAL*8 :: alfa2_lat_thr

  !> Flag to activate the injection of external water:\n
  !> - ext_water_flag = .TRUE.   => interaction with external water
  !> - ext_water_flag = .FALSE.  => no interaction
  LOGICAL :: ext_water_flag
  
  !> Flag for instantaneous vaporization:\n
  !> - inst_vaporization = .TRUE.   => instantaneous vaporization
  !> - inst_vaporization = .FALSE.  => injection as liquid
  LOGICAL :: inst_vaporization

  !> Aquifer type:\n
  !> - 'confined'
  !> - 'unconfined'
  !> .
  CHARACTER*30 :: aquifer_type

  !> Total water influx
  REAL*8 :: total_water_influx

  !> Minimum depth of influx
  REAL*8 :: min_z_influx

  REAL*8 :: delta_z_influx

  !> Maximum depth of influx
  REAL*8 :: max_z_influx

  !> Water temperature in the aquifer
  REAL*8 :: T_w

  !> Boiling temperature
  COMPLEX*16 :: T_boiling

  !> Latent heat of vaporization
  REAL*8 :: lambda_w

  !> Heat capacity of water
  REAL*8 :: cv_w

  !> Water viscosity
  REAL*8 :: visc_w

CONTAINS


  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Physical variables
  !
  !> This subroutine evaluates, from the array qp,
  !> the local physical variables  (\f$\alpha_i, u_i, \rho_i, x_i, p_i \f$)
  !> for the two phases and the mixture denisity \f$ \rho_{mix} \f$.
  !> Also the temperature, the crystals volume fractions and the dissolved gas
  !> mass fractions are computed.
  !> In addition, the volume fractions of the gas components in the mixture and 
  !> in the gas phase are evaluated (\f$\alpha_g(i), \alpha_g2(i), \f$).
  !> \param[in]    r_qp     real array variables 
  !> \param[in]    c_qp     complex array variables 
  !> \date 02/10/2012
  !******************************************************************************

  SUBROUTINE phys_var_qp(r_qp,c_qp)

    USE complexify 
    IMPLICIT none

    REAL*8, INTENT(IN), OPTIONAL :: r_qp(n_eqns)
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qp(n_eqns)

    COMPLEX*16 :: qp(n_eqns)

    INTEGER :: i , j , k

    IF ( present(c_qp) ) THEN

       qp = c_qp

    ELSE

       qp = DCMPLX( r_qp , 0.D0 )

    END IF

    p_1 = qp(idx_p1)
    p_2 = qp(idx_p2)
    u_1 = qp(idx_u1)
    u_2 = qp(idx_u2)
    T = qp(idx_T)

    ! dissolved and axsolved gas variables
    x_d_md(1:n_gas) = qp(idx_xd_first:idx_xd_last)

    alfa_g(1:n_gas) = qp(idx_alfa_first:idx_alfa_last)

    alfa_2 = SUM( alfa_g(1:n_gas) )

    alfa_g_2(1:n_gas) = alfa_g(1:n_gas) / alfa_2

    alfa_1 = 1.D0 - alfa_2

    ! crystal variables
    IF ( method_of_moments_flag ) THEN
       
       DO i = 1,n_cry

          DO k = 1,2
          
             DO j = 0,n_mom-1
             
                mom_cry(i,j,k) = qp(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)+j) 
                
             END DO

          END DO
          
          beta(i) = cry_shape_factor(i) * SUM( mom_cry( i , 3 , 1:2 ) )  / alfa_1          
          
       END DO

       DO i = 1,n_components

          rhoB_components(i) = qp(idx_cry_eqn_first+2*n_mom*(n_cry)-1 + i) 

       ENDDO

       DO i = 1,n_cry 

          sum_rhoB_components(i) = 0.0

          DO j = 1,n_components

             sum_rhoB_components(i) = sum_rhoB_components(i) + rhoB_components(j) * rel_cry_components(j,i) 

          END DO

       END DO !CAMBIO

    ELSE

       beta(1:n_cry) = qp(idx_beta_first:idx_beta_last)

    END IF
       
    ! eval_densities requires: beta, x_d_md , p_1 , p_2 , T
    CALL eval_densities

    alfarho_2 = alfa_2 * rho_2

    x_c_1(1:n_cry) = rho_c(1:n_cry) * beta(1:n_cry) / rho_1
    x_md_1 = DCMPLX(1.D0,0.D0) - SUM(x_c_1)

    x_d_1(1:n_gas) = x_d_md(1:n_gas) * x_md_1
    x_m_1 = DCMPLX(1.D0,0.D0) - SUM(x_d_1) - SUM(x_c_1)

    alfa_d_1(1:n_gas) = rho_1 * x_d_1(1:n_gas) / rho_d(1:n_gas) 
    alfa_m_1 = rho_1 * x_m_1 / rho_m

    rho_mix = alfa_1 * rho_1 + alfa_2 * rho_2

    x_1 = alfa_1 * rho_1 / rho_mix
    x_2 = alfa_2 * rho_2 / rho_mix

    x_g(1:n_gas) = alfa_g(1:n_gas) * rho_g(1:n_gas) / rho_mix

    x_g_2(1:n_gas) = x_g(1:n_gas) / x_2

    x_d(1:n_gas) = x_d_1(1:n_gas) * x_1
    x_m = x_m_1 * x_1
    x_c(1:n_cry) = x_c_1(1:n_cry) * x_1

    u_mix = x_1 * u_1 + x_2 * u_2

    rhoB_m = rho_mix * x_m
    rhoB_c(1:n_cry) = rho_mix * x_c(1:n_cry)

    cv_mix = x_m * cv_m + SUM( x_c(1:n_cry) * cv_c(1:n_cry) )                   &
         + SUM( x_d(1:n_gas) * cv_d(1:n_gas) )                                  &
         + SUM( x_g(1:n_gas) * cv_g(1:n_gas) )

  END SUBROUTINE phys_var_qp

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the physical variables qp the 
  !> array of additional physical variables qp2:\n
  !> - qp2(1) = \f$ \rho_1 \f$
  !> - qp2(2:1+ngas) = \f$ \rho_g(1:ngas) \f$
  !> - qp2(2+ngas:1+ngas+ncry) = \f$ \beta^{eq}(1:ncry) \f$
  !> - qp2(2+ngas+ncry) = \f$ x_{d,md}^eq(1:ngas) \f$
  !> .
  !> \param[in]     r_qp   real physical variables 
  !> \param[out]    qp2    real additional physical variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE eval_qp2( r_qp , qp2 )

    IMPLICIT none

    REAL*8, INTENT(IN) :: r_qp(n_vars)
    REAL*8, INTENT(OUT) :: qp2(1+n_cry+n_gas+n_gas+4)

    COMPLEX*16 :: qp(n_vars)

    INTEGER :: i

    DO i = 1,n_vars

       qp(i) = DCMPLX( r_qp(i), 0.D0 )

    END DO

    CALL phys_var_qp( c_qp = qp )
    CALL eos

    CALL f_beta_eq

    CALL f_xdis_eq

    qp2(1) = REAL(rho_1)
    qp2(1+1:1+n_gas) = REAL(rho_g(1:n_gas))
    qp2(1+n_gas+1:1+n_gas+n_cry) = REAL(beta_eq(1:n_cry))
    qp2(1+n_gas+n_cry+1:1+n_gas+n_cry+n_gas) = REAL( x_d_md_eq(1:n_gas) )

    CALL mixture_viscosity

    qp2(1+n_gas+n_cry+n_gas+1) = REAL( visc_mix )

    CALL f_viscmelt

    qp2(1+n_gas+n_cry+n_gas+2) = REAL( visc_melt )

    CALL f_theta

    qp2(1+n_gas+n_cry+n_gas+3) = REAL( theta )

    CALL f_bubbles

    qp2(1+n_gas+n_cry+n_gas+4) = REAL( visc_rel_bubbles )

  END SUBROUTINE eval_qp2

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the physical
  !> variables qp.
  !> \date 15/08/2011
  !> \param[in]     c_qp     complex pysical variables 
  !> \param[in]     r_qp     real physical variables 
  !> \param[out]    c_flux   complex analytical fluxes    
  !> \param[out]    r_flux   real analytical fluxes    
  !******************************************************************************

  SUBROUTINE eval_fluxes_qp(c_qp,r_qp,c_flux,r_flux)

    USE complexify
    USE geometry, ONLY : radius

    IMPLICIT none

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qp(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_flux(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qp(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_flux(n_eqns)

    COMPLEX*16 :: qp(n_eqns)
    COMPLEX*16 :: flux(n_eqns)

    INTEGER :: i , j , k

    IF ( present(c_qp) .AND. present(c_flux) ) THEN

       qp = c_qp

    ELSEIF ( present(r_qp) .AND. present(r_flux) ) THEN

       qp = DCMPLX( r_qp , 0.D0 )

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    CALL phys_var_qp( c_qp = qp )
    CALL eos

    flux(1:n_eqns) = DCMPLX(0.D0,0.D0)

    !---- Mixture Density -------------------------------------------------------
    flux(idx_mix_mass_eqn) = rho_mix * u_mix * radius**2

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Mixture mass flux'
       WRITE(*,*) flux(idx_mix_mass_eqn)

    END IF
       
    !---- Volumetric Fraction First Phase ---------------------------------------
    IF ( p_relax_model .EQ. 'single' ) THEN

       flux(idx_vol1_eqn) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       flux(idx_vol1_eqn) = rho_mix * u_mix * alfa_1 * radius**2

    END IF

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Volume Fraction1 flux'
       WRITE(*,*) flux(idx_vol1_eqn)

    END IF
    
    !---- Mixture Momentum ------------------------------------------------------
    flux(idx_mix_mom_eqn) = ( alfa_1 * rho_1 * u_1**2 + alfa_2 * rho_2 * u_2**2 &
         + alfa_1 * p_1 + alfa_2 * p_2 ) * radius**2

    IF ( verbose_level .GE. 3 ) THEN
       
       WRITE(*,*) 'Mixture momentum flux'
       WRITE(*,*) flux(idx_mix_mom_eqn)
       
    END IF
           
    !---- Relative Velocity -----------------------------------------------------
    IF ( drag_funct_model .EQ. 'single_velocity' ) THEN

       flux(idx_rel_vel_eqn) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       flux(idx_rel_vel_eqn) = ( 0.5D0 * ( u_1 * u_1 - u_2 * u_2 )              &
            + ( mu_1 - mu_2 ) ) * radius**2

    END IF

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Relative velocity flux'
       WRITE(*,*) flux(idx_rel_vel_eqn)

    END IF
       
    !---- Total Energy ----------------------------------------------------------
    IF ( isothermal ) THEN

       flux(idx_mix_engy_eqn) = DCMPLX( 0.D0 , 0.D0 )

    ELSE

       flux(idx_mix_engy_eqn) = ( alfa_1 * rho_1 * u_1 * ( e_1 + p_1/rho_1      &
            + 0.5D0 * u_1 * u_1 ) + alfa_2 * rho_2 * u_2 * ( e_2 + p_2/rho_2    &
            + 0.5D0 * u_2* u_2 ) - rho_mix * x_1 * x_2 * ( u_1 - u_2 )          &
            * ( s_1 - s_2 ) * T ) * radius**2
    END IF

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Mixture energy flux'
       WRITE(*,*) flux(idx_mix_engy_eqn)

    END IF

    !----- Dissolved Gas Phases -------------------------------------------------
    flux(idx_dis_gas_eqn_first:idx_dis_gas_eqn_last) =                          &
         ( x_d_md(1:n_gas)*alfa_1 *                                             &
         ( rho_1 - SUM( beta(1:n_cry) * rho_c(1:n_cry) ) ) * u_1 ) * radius**2

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Dissolved gas mass fraction fluxes'
       WRITE(*,*) flux(idx_ex_gas_eqn_first:idx_ex_gas_eqn_last)

    END IF

    !---- Mass Fraction Exsolved Gas Phases -------------------------------------
    flux(idx_ex_gas_eqn_first:idx_ex_gas_eqn_last) = alfa_g(1:n_gas) *          &
         rho_g(1:n_gas) * u_2 * radius**2

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'Mass fraction exsolved gas'
       WRITE(*,*) flux(idx_ex_gas_eqn_first:idx_ex_gas_eqn_last)

    END IF

    !----- Crystal Phases -------------------------------------------------------

    IF ( method_of_moments_flag ) THEN
       
       DO i=1,n_cry
          
          DO k=1,2
          
             DO j=0,n_mom-1
                
                flux(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)+j) = mom_cry(i,j,k) * u_1    &
                     * radius**2

             END DO

          END DO

       END DO

       DO i = 1,n_components

	        flux(idx_cry_eqn_first + 2*n_mom*(n_cry) - 1 + i) = rhoB_components(i) * u_1 * radius ** 2

       ENDDO

    ELSE

       flux(idx_cry_eqn_first:idx_cry_eqn_last) = alfa_1 * rho_c(1:n_cry) *     &
            beta(1:n_cry) * u_1 * radius**2

       IF ( verbose_level .GE. 3 ) THEN
       
          WRITE(*,*) 'Crystal volume fraction'
          WRITE(*,*) flux(idx_cry_eqn_first:idx_cry_eqn_last)
          
       END IF
       
    END IF
       
    IF ( present(c_qp) .AND. present(c_flux) ) THEN

       c_flux = flux

    ELSEIF ( present(r_qp) .AND. present(r_flux) ) THEN

       r_flux = REAL( flux )

    END IF

  END SUBROUTINE eval_fluxes_qp

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qp            complex physical variables 
  !> \param[in]     r_qp            real physical variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms_qp( c_qp , c_nh_term_impl , r_qp ,        &
       r_nh_term_impl )

    USE COMPLEXIFY 
    USE geometry, ONLY : radius
    IMPLICIT NONE

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qp(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qp(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX*16 :: qp(n_eqns)

    COMPLEX*16 :: nh_term_impl(n_eqns)

    COMPLEX*16 :: relaxation_term(n_eqns)
    COMPLEX*16 :: forces_term(n_eqns)
    COMPLEX*16 :: source_term(n_eqns)

    INTEGER :: i

    IF ( present(c_qp) .AND. present(c_nh_term_impl) ) THEN

       qp = c_qp

    ELSEIF ( present(r_qp) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_eqns

          qp(i) = DCMPLX( r_qp(i) , 0.D0 )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    CALL phys_var_qp( c_qp = qp )
    CALL eos

    !--------------- Evaluate the relaxation terms -----------------------------
    CALL eval_relaxation_terms( relaxation_term )

    !--------------- Evaluate the forces terms ---------------------------------
    CALL eval_forces_terms( forces_term )

    !--------------- Evaluate the forces terms ---------------------------------
    CALL eval_source_terms( source_term )

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'relaxation_term'
       WRITE(*,*) relaxation_term / ( radius**2 )
    
       WRITE(*,*) 'forces_term'
       WRITE(*,*) forces_term / ( radius**2 )
    
       WRITE(*,*) 'source_term'
       WRITE(*,*) source_term / ( radius**2 )
    
    END IF

    nh_term_impl = relaxation_term + forces_term + source_term

    IF ( present(c_qp) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term_impl

    ELSEIF ( present(r_qp) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term_impl )

    END IF

  END SUBROUTINE eval_nonhyperbolic_terms_qp

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Relaxation terms
  !
  !> This subroutine evaluate the relaxation terms for pressures, velocities, 
  !> exsolution, crystallization and fragmentation.
  !> \date 20/06/2013
  !> \param[out]    relaxation_term        complex relaxation terms
  !******************************************************************************

  SUBROUTINE eval_relaxation_terms( relaxation_term )

    USE complexify 
    USE geometry, ONLY : radius
    IMPLICIT none

    COMPLEX*16, INTENT(OUT) :: relaxation_term(n_eqns)

    COMPLEX*16 :: pressure_relaxation 
    COMPLEX*16 :: velocity_relaxation
    !COMPLEX*16 :: frag_relaxation

    INTEGER :: i , j , k
    INTEGER :: idx

    idx = 0

    relaxation_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! relaxation term for the mixture density -----------------------------------
    relaxation_term(idx_mix_mass_eqn) = DCMPLX(0.D0,0.D0)

    ! relaxation term for the volume fraction equation --------------------------
    CALL press_relax_term( pressure_relaxation )
    relaxation_term(idx_vol1_eqn) = pressure_relaxation * radius**2

    ! relaxation term for the mixture momentum ----------------------------------
    relaxation_term(idx_mix_mom_eqn) = DCMPLX(0.D0,0.D0)

    ! relaxation term for relative velocity -------------------------------------
    CALL vel_relax_term( velocity_relaxation )
    relaxation_term(idx_rel_vel_eqn) = velocity_relaxation * radius**2

    ! relaxation term for the mixture energy ------------------------------------
    relaxation_term(idx_mix_engy_eqn) = DCMPLX(0.D0,0.D0)

    ! relaxation term for dissolved gas -----------------------------------------
    CALL f_xdis_eq
    relaxation_term(idx_dis_gas_eqn_first:idx_dis_gas_eqn_last) =               &
         - ( x_d_md(1:n_gas) - x_d_md_eq(1:n_gas) ) * ( 1.D0 - alfa_2 )         &
         * ( rho_1 - SUM( beta(1:n_cry) * rho_c(1:n_cry) ) ) / tau_d(1:n_gas)   &
         * radius ** 2
    
    DO i=1,n_gas

       IF ( REAL( x_d_md(i) ) .LT. REAL( x_d_md_eq(i) ) ) THEN

          relaxation_term(idx_dis_gas_eqn_first+i-1) = DCMPLX(0.D0,0.D0)

       END IF

    END DO

    ! relaxation term for exsolved gas ------------------------------------------
    CALL f_xdis_eq

    relaxation_term(idx_ex_gas_eqn_first:idx_ex_gas_eqn_last) =                 &
         ( x_d_md(1:n_gas) - x_d_md_eq(1:n_gas) ) * ( 1.D0 - alfa_2)            &
         * ( rho_1 - SUM( beta(1:n_cry) * rho_c(1:n_cry) ) ) / tau_d(1:n_gas)   &
         * radius **2
    
    DO i=1,n_gas

       IF ( REAL( x_d_md(i) ) .LT. REAL( x_d_md_eq(i) ) ) THEN

          relaxation_term(idx_ex_gas_eqn_first+i-1) = DCMPLX(0.D0,0.D0)

       END IF

    END DO

    ! relaxation term for crystallization ---------------------------------------

    IF ( method_of_moments_flag ) THEN

       CALL update_kinetics
      
       DO i = 1,n_cry

          DO k = 1,2
          
             IF(k == 1) THEN 

	        relaxation_term(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)) =                 &
		   nucleation_rate(i) * L_nucleus(i)**j * radius**2 * sum_rhoB_components(i) 

             ELSE

	        relaxation_term(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)+j) =                 &
		   0.0

             END IF
          
             DO j = 1,n_mom-1

		IF(k == 1) THEN 

	           relaxation_term(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)+j) =                 &
		      nucleation_rate(i) * L_nucleus(i)**j * radius**2 * sum_rhoB_components(i) +   &
	              radius**2 * sum_rhoB_components(i) * j * growth_rate(i) * mom_cry(i,j-1,k)

		ELSE

	           relaxation_term(idx_cry_eqn_first+2*n_mom*(i-1)+n_mom*(k-1)+j) =                 &
		      radius**2 * sum_rhoB_components(i) * j * growth_rate(i) * mom_cry(i,j-1,k)

		END IF

             END DO
                
          END DO

       END DO
       
       DO i = 1,n_components

          relaxation_term(idx_cry_eqn_first + 2*n_mom*(n_cry) - 1 + i) = 0.0

          DO j = 1,n_cry

             DO k = 1,2

                relaxation_term(idx_cry_eqn_first + 2*n_mom*(n_cry) - 1 + i) =                         &
                   relaxation_term(idx_cry_eqn_first + 2*n_mom*(n_cry) - 1 + i) - cry_shape_factor(j)* &
                   relaxation_term(idx_cry_eqn_first + 2*n_mom*(j-1) + n_mom*(k-1) + 3) *              &
		   cry_current_solid_solution(j,i) * rho_c(j)

             ENDDO   

          ENDDO

       ENDDO

    ELSE

       CALL f_beta_eq

       relaxation_term(idx_cry_eqn_first:idx_cry_eqn_last) =                    &
            - ( 1.d0 - alfa_2 ) * rho_c(1:n_cry)                                &
            * ( beta(1:n_cry) - beta_eq(1:n_cry) ) / tau_c(1:n_cry) * radius ** 2

    END IF
       
  END SUBROUTINE eval_relaxation_terms

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Force terms
  !
  !> This subroutine evaluates the forces terms due to gravity and viscosity.
  !> \date 20/06/2013
  !> \param[out]    force_term        complex force terms
  !******************************************************************************

  SUBROUTINE eval_forces_terms( force_term )

    USE geometry, ONLY : radius, f_eccen_b
    USE complexify 
    IMPLICIT none

    COMPLEX*16, INTENT(OUT) :: force_term(n_eqns)

    COMPLEX*16 :: visc_force_1 , visc_force_2
    COMPLEX*16 :: visc_force_1_rel , visc_force_2_rel

    REAL*8 :: gas_wall_drag

    INTEGER :: i

    INTEGER :: idx

    idx = 0

    force_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! Term for the mixture density ----------------------------------------------
    force_term(idx_mix_mass_eqn) = DCMPLX(0.D0,0.D0)

    ! Term for the first phase volume fraction ----------------------------------
    force_term(idx_vol1_eqn) = DCMPLX(0.D0,0.D0)

    ! Term for the mixture momentum equation ------------------------------------
    force_term(idx_mix_mom_eqn) = - rho_mix * grav * radius**2.D0

    CALL mixture_viscosity

    visc_force_1 = - 8.D0 * visc_mix * u_1 / f_eccen_b / f_eccen_b

    ! Turbulent gas-wall friction (Degruyter et al. 2012)

    gas_wall_drag = 0.03D0

    visc_force_2 = - gas_wall_drag / 4.D0 * radius * rho_2 * CDABS(u_2) * u_2   &
         / f_eccen_b / f_eccen_b

    visc_force_1 = visc_force_1 * ( 1.D0 - frag_eff )
    visc_force_2 = visc_force_2 * frag_eff

    !visc_force_1 =  visc_force_1 * t
    !visc_force_2 =  visc_force_2 * ( 1.0 - t )

    force_term(idx_mix_mom_eqn) = force_term(idx_mix_mom_eqn) + visc_force_2    &
         + visc_force_1 

    ! Term for the relative velocity equation -----------------------------------
    force_term(idx_rel_vel_eqn) = DCMPLX(0.0,0.0)

    IF ( drag_funct_model .EQ. 'single_velocity' ) THEN

       force_term(idx_rel_vel_eqn) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       visc_force_1_rel = visc_force_1 / ( ( 1.D0 - alfa_2 ) * rho_1 ) 

       visc_force_2_rel = - visc_force_2 / ( alfa_2 * rho_2 )

       force_term(idx_rel_vel_eqn) = force_term(idx_rel_vel_eqn)                &
            + visc_force_2_rel + visc_force_1_rel

    END IF

    ! Term for the mixture energy -----------------------------------------------
    IF ( isothermal ) THEN

       force_term(idx_mix_engy_eqn) = DCMPLX(0.d0,0.D0)

    ELSE

       force_term(idx_mix_engy_eqn) = - rho_mix * u_mix * grav * radius**2

       force_term(idx_mix_engy_eqn) = force_term(idx_mix_engy_eqn)              &
            + visc_force_2 * u_2 + visc_force_1 * u_1

    END IF

    ! Terms for the exsolved gas phases -----------------------------------------
    DO i=idx_dis_gas_eqn_first,idx_dis_gas_eqn_last

       force_term(i) = DCMPLX(0.D0,0.D0)

    END DO
    
    ! Terms for the dissolved gas phases ----------------------------------------
    DO i = idx_ex_gas_eqn_first,idx_ex_gas_eqn_last

       force_term(i) = DCMPLX(0.d0,0.D0)

    END DO

    ! Terms for the crystal phases ----------------------------------------------
    DO i= idx_cry_eqn_first,idx_cry_eqn_last

       force_term(i) = DCMPLX(0.d0,0.D0)

    END DO

  END SUBROUTINE eval_forces_terms

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Source terms
  !
  !> This subroutine evaluates the loss of mass due to lateral degassing or the
  !> inlet of mass due to waterfrom a surrounding aquifer and the corresponding
  !> terms to add in the system of equations. The terms are added alson in the
  !> energy equation.
  !> \date 20/06/2013
  !> \param[out]    source_term        complex source terms
  !******************************************************************************

  SUBROUTINE eval_source_terms( source_term )

    USE geometry, ONLY : radius, f_eccen_a, f_eccen_b
    USE complexify 
    IMPLICIT none

    COMPLEX*16, INTENT(OUT) :: source_term(n_eqns)

    COMPLEX*16 :: q_lat

    COMPLEX*16 :: water_mass_flux

    COMPLEX*16 :: heat_flux

    INTEGER :: i

    INTEGER :: idx

    idx = 0

    q_lat = DCMPLX(0.D0,0.D0)

    IF ( lateral_degassing ) THEN
       
       CALL lithostatic_pressure
       
       IF ( p_2 .GE. p_lith ) THEN
          
          q_lat = ( rho_2 * alfa_2 * k_cr * ( p_2 - p_lith ) ) /                &
               ( visc_2 * radius * f_eccen_b )
          
       ELSE
          
          q_lat = DCMPLX(0.D0,0.D0)
          
       END IF
       
    END IF

    water_mass_flux = DCMPLX(0.D0,0.D0)

    source_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! --- TOTAL MASS SOURCE TERM ------------------------------------------------
    IF ( lateral_degassing ) THEN

       source_term(idx_mix_mass_eqn) = - 2.D0 * q_lat * radius * f_eccen_a

    ELSE

       source_term(idx_mix_mass_eqn) = DCMPLX(0.D0,0.D0)

    END IF

    IF ( ext_water_flag ) THEN

       rho_w = 971.80  ! 80C. For shallow aquifers, the effect of P is expected
                       ! to be less important
       
       IF ( ( zeta_lith .GT. min_z_influx ) .AND.                               &
            ( zeta_lith .LT. min_z_influx + delta_z_influx ) ) THEN
          
          IF ( total_water_influx .GT. 0.D0 ) THEN
             
             water_mass_flux = DCMPLX( total_water_influx/delta_z_influx , 0.D0)
             
          ELSE
             
             IF ( aquifer_type .EQ. 'unconfined' ) THEN
                
                CALL hydrostatic_pressure
  
                IF ( p_hydro .GE. REAL(p_1) ) THEN
                   
                   visc_w = 2.414D-5 * 10.D0 ** ( 247.8D0 / ( T_w - 140.D0 ) )
                   
                   water_mass_flux = rho_w * k_cr / visc_w * ( p_hydro - p_1 )  &
                        / radius / f_eccen_b
                   
                ELSE
                   
                   water_mass_flux = DCMPLX( 0.D0 , 0.D0 )
                   
                END IF

             ELSEIF ( aquifer_type .EQ. 'confined') THEN

                CALL lithostatic_pressure
         
                IF ( p_lith .GE. REAL(p_1) ) THEN
                   
                   visc_w = 2.414D-5 * 10.D0 ** ( 247.8D0 / ( T_w - 140.D0 ) )
                   
                   water_mass_flux = rho_w * k_cr / visc_w * ( p_lith - p_1 )   &
                        / radius / f_eccen_b
                   
                ELSE
                   
                   water_mass_flux = DCMPLX( 0.D0 , 0.D0 )
                   
                END IF
                
             END IF
             
          END IF
          
       ELSE
          
          water_mass_flux = DCMPLX(0.D0,0.D0)
          
       END IF
       
       source_term(idx_mix_mass_eqn) = source_term(idx_mix_mass_eqn)            &
            + 2.D0 * radius * water_mass_flux * f_eccen_a
       
    END IF
    
    ! --- FIRST PHASE VOLUME FRACTION -------------------------------------------
    source_term(idx_vol1_eqn) = DCMPLX(0.D0,0.D0)

    ! --- Mixture Momentum ------------------------------------------------------
    IF ( lateral_degassing ) THEN

       source_term(idx_mix_mom_eqn) = - 2.D0 * q_lat * radius * u_2 * f_eccen_a

    ELSE

       source_term(idx_mix_mom_eqn) = DCMPLX(0.D0,0.D0)

    END IF

    ! --- Relative Velocity -----------------------------------------------------
    source_term(idx_rel_vel_eqn) = DCMPLX(0.D0,0.D0)

    ! --- MIXTURE ENERGY --------------------------------------------------------
    IF ( isothermal ) THEN

       source_term(idx_mix_engy_eqn) = DCMPLX(0.D0,0.D0)

    ELSE

       IF ( lateral_degassing ) THEN

          source_term(idx_mix_engy_eqn) = - 2.D0 * q_lat * radius * f_eccen_a * &
               ( T * cv_2 + 0.5D0 * u_2*u_2 )

       ELSE

          source_term(idx_mix_engy_eqn) = DCMPLX(0.D0,0.D0)

       END IF

       IF ( ext_water_flag ) THEN

          IF ( ( zeta_lith .GT. min_z_influx ) .AND.                            &
               ( zeta_lith .LT. min_z_influx + delta_z_influx ) ) THEN

             IF ( inst_vaporization ) THEN

                heat_flux = - 2.D0 * radius * f_eccen_a * water_mass_flux * &
                     ( cv_d(1) * ( T_boiling-T_w ) + cv_2*( -T_boiling ) &
                     + lambda_w )

             ELSE

                heat_flux = 2.D0 * radius * f_eccen_a * water_mass_flux *       &
                     cv_d(1) * T_w

             END IF

             source_term(idx_mix_engy_eqn) = source_term(idx_mix_engy_eqn)      &
                  + heat_flux

          END IF

       END IF

    END IF

    ! --- DISSOLVED GAS BULK DENSITY --------------------------------------------

    ! H2O source due to inlet
    IF ( ext_water_flag .AND. .NOT.inst_vaporization ) THEN

       source_term(idx_dis_gas_eqn_first) = 2.D0 * radius * f_eccen_a           &
            * water_mass_flux

    ELSE

       source_term(idx_dis_gas_eqn_first) = DCMPLX(0.D0,0.D0)

    END IF

    ! Other gas phases (from 2 to n_gas)
    DO i = idx_dis_gas_eqn_first+1,idx_dis_gas_eqn_last

       source_term(i) = DCMPLX(0.D0,0.D0)

    END DO

    ! --- EXSOLVED GAS MASS FRACTIONS -------------------------------------------

    DO i = idx_ex_gas_eqn_first,idx_ex_gas_eqn_last

       IF ( lateral_degassing ) THEN

          source_term(i) = - 2.D0 * q_lat * radius * f_eccen_a

       ELSE

          source_term(i) = DCMPLX(0.D0,0.D0)

       END IF

    END DO

    ! --- CRYSTALS BULK DENSITY -------------------------------------------------

    DO i = idx_cry_eqn_first,idx_cry_eqn_last

       source_term(i) = DCMPLX(0.D0,0.D0)

    END DO

  END SUBROUTINE eval_source_terms

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Explicit Forces terms
  !
  !> This subroutine evaluates the forces to be treated explicitely.
  !> \date 01/06/2012
  !> \param[out]    expl_forces_term   forces term
  !******************************************************************************

  SUBROUTINE eval_explicit_forces( expl_forces_term )

    USE geometry, ONLY : radius

    IMPLICIT NONE

    REAL*8, INTENT(OUT) :: expl_forces_term(n_eqns)  !< explicit forces 

    INTEGER :: i 

    INTEGER :: idx

    idx = 0

    expl_forces_term(1:n_eqns) = 0.D0

    ! --- TOTAL MASS TERM -------------------------------------------------------
    expl_forces_term(idx_mix_mass_eqn) = 0.D0

    ! --- FIRST PHASE VOLUME FRACTION -------------------------------------------
    expl_forces_term(idx_vol1_eqn) = 0.D0

    ! --- Mixture Momentum ------------------------------------------------------
    expl_forces_term(idx_mix_mom_eqn) = DREAL( rho_mix * grav * radius ** 2 )

    ! --- Relative Velocity -----------------------------------------------------
    expl_forces_term(idx_rel_vel_eqn) = 0.D0

    ! --- MIXTURE ENERGY --------------------------------------------------------
    IF ( isothermal ) THEN

       expl_forces_term(idx_mix_engy_eqn) = 0.D0

    ELSE

       expl_forces_term(idx_mix_engy_eqn) = DREAL( rho_mix * u_mix * grav *     &
            radius **2 )

    END IF

    ! --- DISSOLVED GAS BULK DENSITY --------------------------------------------
    DO i = idx_dis_gas_eqn_first,idx_dis_gas_eqn_last

       expl_forces_term(i) = 0.D0

    END DO

    ! --- EXSOLVED GAS MASS FRACTIONS -------------------------------------------
    DO i = idx_ex_gas_eqn_first,idx_ex_gas_eqn_last

       expl_forces_term(i) = 0.D0

    END DO

    ! --- CRYSTALS BULK DENSITY ------------------------------------------------
    DO i = idx_cry_eqn_first,idx_cry_eqn_last

       expl_forces_term(i) = 0.D0

    END DO

  END SUBROUTINE eval_explicit_forces

END MODULE equations
