
MODULE steady_solver
  ! external variables
  USE geometry, ONLY : z0 , zN, zeta_exit , radius

  USE parameters, ONLY : n_eqns , n_vars
  USE parameters, ONLY : max_nl_iter
  USE parameters, ONLY : verbose_level

  USE parameters, ONLY : idx_p1 , idx_p2 , idx_u1 , idx_u2 , idx_T ,            &
       idx_xd_first , idx_xd_last , idx_alfa_first , idx_alfa_last ,            &
       idx_beta_first , idx_beta_last
  
  USE parameters, ONLY : idx_mix_mass_eqn , idx_vol1_eqn , idx_mix_mom_eqn ,    &
       idx_rel_vel_eqn , idx_mix_engy_eqn , idx_dis_gas_eqn_first ,             &
       idx_dis_gas_eqn_last , idx_ex_gas_eqn_first , idx_ex_gas_eqn_last ,      &
       idx_cry_eqn_first , idx_cry_eqn_last
  
  USE equations, ONLY : lateral_degassing_flag
  USE equations, ONLY : lateral_degassing
  USE equations, ONLY : alfa2_lat_thr

  USE constitutive, ONLY : explosive

  USE init, ONLY : p_out

  USE equations, ONLY : isothermal

  IMPLICIT NONE

  REAL*8 :: zeta
  REAL*8, ALLOCATABLE :: fluxes_old(:)
  REAL*8, ALLOCATABLE :: nh_terms_old(:)

  REAL*8 :: zeta_old , dz_max

  LOGICAL :: increase_flow_rate

  INTEGER :: nl_iter

CONTAINS

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> @brief Shooting Method
  !
  !> This subroutine search for the steady solution using a Shooting Method for 
  !> Two-Point Boundary Value Problems. A bisection method is used to find the
  !> value of inlet velocity, in order to reach the desired prescribed boundary
  !> condition at the exit.
  !> \date 13/03/2013
  !******************************************************************************

  SUBROUTINE steady_shooting

    ! external procedures
    USE constitutive, ONLY : eos
    USE constitutive, ONLY : eval_densities
    USE equations, ONLY : phys_var_qp
    USE init, ONLY : init_steady

    ! external variables
    USE constitutive, ONLY : rho_m, rho_c, rho_mix, u_mix, u_1, u_2
    USE constitutive, ONLY : p_1 , p_2 , T , s_1, s_2, e_1, e_2
    USE constitutive, ONLY : rho_2 , rho_g, rho_1, rho_c, rho_m
    USE constitutive, ONLY : mu_1, mu_2, s_c, mu_c, s_g, mu_g
    USE constitutive, ONLY : e_c, e_g, e_m, mu_m, s_m
    USE init, ONLY : u1_in
    USE parameters, ONLY : shooting
    USE parameters, ONLY : eps_conv
    USE geometry, ONLY : update_radius

    IMPLICIT NONE

    REAL*8, ALLOCATABLE :: qp(:)

    REAL*8 :: V_0 , V_1, V_2
    REAL*8 :: u_inlet

    LOGICAL :: check_flow_rate

    INTEGER :: iter

    REAL*8 :: extrap_z , extrap_z_p , extrap_z_mach
    INTEGER :: extrap_flag
    REAL*8 :: r_p_1 , r_p_2 , r_p
    REAL*8 :: mach

    REAL*8 :: extrap_z_0 , extrap_z_2

    REAL*8 :: V_temp , V_coeff
    LOGICAL :: check_interval

    LOGICAL :: flag_output , increase_temp

    REAL*8 :: initial_mass_flow_rate, final_mass_flow_rate


    ALLOCATE( qp(n_vars) ) 
    ALLOCATE( fluxes_old(n_eqns) )
    ALLOCATE( nh_terms_old(n_eqns) )


    IF ( shooting ) THEN

       IF ( explosive ) THEN

          V_temp = u1_in

       ELSE

          V_temp = u1_in

       END IF

       flag_output = .FALSE.

    ELSE

       V_temp = u1_in

       flag_output = .FALSE.

    END IF

    ! ------------------- SOLVE with u_inlet = V_0 -----------------------------

    u_inlet = V_temp

    CALL init_steady( u_inlet , qp )

    CALL phys_var_qp(qp)
    CALL eos
    CALL eval_densities

    IF ( verbose_level .GE. 1 ) THEN
       p_1 = DCMPLX( qp(idx_p1) , 0.D0 ) 
       p_2 = DCMPLX( qp(idx_p2) , 0.D0 )
       u_1 = DCMPLX( qp(idx_u1) , 0.D0 ) 
       u_2 = DCMPLX( qp(idx_u2) , 0.D0 ) 
       T = DCMPLX( qp(idx_T) , 0.D0 )

       WRITE(*,*) ''
       WRITE(*,*) '<<<<<<<<<<<< ThermoPhysical quantities >>>>>>>>>>>'
       WRITE(*,*) 'P_1 =', REAL(p_1)
       WRITE(*,*) 'P_2 =', REAL(p_2)
       WRITE(*,*) 'T =', REAL(T)
       WRITE(*,*) ''

       WRITE(*,*) 'rho_c =', REAL(rho_c)
       WRITE(*,*) 's_c =', REAL(s_c)
       WRITE(*,*) 'e_c =', REAL(e_c)
       WRITE(*,*) 'enth_c =', REAL(e_c + p_1 / rho_c)
       WRITE(*,*) 'gibbs_c =', REAL(mu_c)
       WRITE(*,*) ''

       WRITE(*,*) 'rho_g =', REAL(rho_g)
       WRITE(*,*) 's_g =', REAL(s_g)
       WRITE(*,*) 'e_g =', REAL(e_g)
       WRITE(*,*) 'enth_g =', REAL(e_g + p_2 / rho_g)
       WRITE(*,*) 'gibbs_g =', REAL(mu_g)
       WRITE(*,*) ''

       WRITE(*,*) 'rho_m =', REAL(rho_m)
       WRITE(*,*) 's_m =', REAL(s_m)
       WRITE(*,*) 'e_m =', REAL(e_m)
       WRITE(*,*) 'enth_m =', REAL(e_m + p_1 / rho_m)
       WRITE(*,*) 'gibbs_m =', REAL(mu_m)
       WRITE(*,*) ''


       WRITE(*,*) 'rho_1 =', REAL(rho_1)
       WRITE(*,*) 's_1 =', REAL(s_1)
       WRITE(*,*) 'e_1 =', REAL(e_1)
       WRITE(*,*) 'enth_1 =', REAL(e_1 + p_1 / rho_1)
       WRITE(*,*) 'gibbs_1 =', REAL(mu_1)
       WRITE(*,*) ''

       WRITE(*,*) 'rho_2 =', REAL(rho_2)
       WRITE(*,*) 's_2 =', REAL(s_2)
       WRITE(*,*) 'e_2 =', REAL(e_2)
       WRITE(*,*) 'enth_2 =', REAL(e_2 + p_2 / rho_2)
       WRITE(*,*) 'gibbs_2 =', REAL(mu_2)
       WRITE(*,*) ''
       WRITE(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>'

    END IF

    WRITE(*,*) ''
    WRITE(*,*) '********* V_inlet = ',V_temp,'**********'

    CALL update_radius(zeta)

    initial_mass_flow_rate = REAL(rho_mix * u_mix)*3.1415*radius*radius

    WRITE(*,*) 'Initial mass flow rate = ',initial_mass_flow_rate,' kg/s'

    CALL integrate_equations(qp , flag_output , extrap_z , extrap_z_p ,         &
         extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

    IF ( .NOT. shooting ) THEN

       flag_output = .TRUE.
       zeta_exit = zeta   
       u_inlet = V_temp

       CALL init_steady( u_inlet , qp )

       CALL phys_var_qp(qp)

       WRITE(*,*) ''
       WRITE(*,*) 'Saving solution for V at the inlet',V_temp

       CALL integrate_equations(qp , flag_output , extrap_z , extrap_z_p ,      &
            extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

    END IF

    CALL update_radius(zeta)

    final_mass_flow_rate = REAL(rho_mix * u_mix)*3.1415*radius*radius

    WRITE(*,*) 'Initial mass flow rate = ',initial_mass_flow_rate,' kg/s'
    WRITE(*,*) 'Final mass flow rate   = ',final_mass_flow_rate,' kg/s'

    extrap_z_0 = extrap_z

    IF ( increase_flow_rate ) THEN

       WRITE(*,*) 'Pressure 1 at exit:', r_p_1
       WRITE(*,*) 'Pressure 2 at exit:', r_p_2
       WRITE(*,*) 'Mach at exit:', mach

       IF ( extrap_flag .EQ. -3 ) THEN

          WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

       ELSE

          WRITE(*,*) 'Exit pressure is reached at:',extrap_z

       END IF

       V_0 = V_temp
       V_coeff = 1.D1
       V_2 = V_temp * V_coeff

    ELSE

       V_2 = V_temp
       V_coeff = 5.D-1
       V_0 = V_temp * V_coeff

       WRITE(*,*) 'At zeta ',zeta
       WRITE(*,*) 'Pressure 1 is: ',r_p_1
       WRITE(*,*) 'Pressure 2 is: ',r_p_2
       WRITE(*,*) 'Mach is:', mach

       IF ( extrap_flag .EQ. 3 ) THEN

          WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

       ELSE

          WRITE(*,*) 'Extrap_flag:',extrap_flag
          WRITE(*,*) 'Exit pressure is reached at:',extrap_z

       END IF

    END IF

    IF ( .NOT.shooting ) THEN

       RETURN

    END IF

    IF ( verbose_level .GE. 1 ) READ(*,*)


    ! -------------------- SOLVE with u_inlet = V_2 ----------------------------

    check_interval = .FALSE.
    increase_temp = increase_flow_rate

    DO WHILE ( .NOT. check_interval )

       IF (V_coeff .LT. 1.0) THEN

          V_2 = V_temp

       ELSE

          V_0 = V_temp

       END IF

       V_temp = V_temp * V_coeff

       WRITE(*,*) '********* V_inlet = ',V_temp,'**********'

       u_inlet = V_temp

       CALL init_steady( u_inlet , qp )

       CALL phys_var_qp(qp)

       flag_output = .FALSE.

       zeta = z0
       CALL update_radius(zeta)
       initial_mass_flow_rate = REAL(rho_mix * u_mix)*3.1415*radius*radius

       CALL integrate_equations(qp , flag_output , extrap_z , extrap_z_p ,      &
            extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

       CALL update_radius(zeta)
       final_mass_flow_rate = REAL(rho_mix * u_mix)*3.1415*radius*radius

       WRITE(*,*) 'Initial mass flow rate = ',initial_mass_flow_rate,' kg/s'
       WRITE(*,*) 'Final mass flow rate   = ',final_mass_flow_rate,' kg/s'

       extrap_z_2 = extrap_z


       IF ( increase_flow_rate ) THEN

          WRITE(*,*) 'Pressure 1 at exit:', r_p_1
          WRITE(*,*) 'Pressure 2 at exit:', r_p_2
          WRITE(*,*) 'Mach at exit:', mach

          IF ( extrap_flag .EQ. -3 ) THEN

             WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

          ELSE

             WRITE(*,*) 'Exit pressure is reached at:',extrap_z

          END IF


          IF ( .NOT. increase_temp ) THEN

             check_interval = .TRUE.
             V_0 = V_temp

          END IF

       ELSE

          WRITE(*,*) 'At zeta ',zeta
          WRITE(*,*) 'Pressure 1 is: ',r_p_1
          WRITE(*,*) 'Pressure 2 is: ',r_p_2
          WRITE(*,*) 'Mach is:', mach

          IF ( extrap_flag .EQ. 3 ) THEN

             WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

          ELSE

             WRITE(*,*) 'Extrap_flag:',extrap_flag
             WRITE(*,*) 'Exit pressure is reached at:',extrap_z

          END IF

          IF ( increase_temp ) THEN

             check_interval = .TRUE.
             V_2 = V_temp

          END IF

       END IF

    END DO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    ! ----------------------START LOOP FOR SOLUTION -----------------------------

    check_flow_rate = .FALSE.

    iter = 0

    DO WHILE ( .NOT. check_flow_rate )

       iter = iter + 1

       V_1 = 0.5D0 * ( V_0 + V_2 )

       WRITE(*,*) '********* V_inlet = ',V_1,'********** iter = ', iter
       WRITE(*,*) 'V_0 = ',V_0,' V_1 = ',V_1,' V_2 = ',V_2
       WRITE(*,*) 'Delta V = ', 0.5D0 * ( V_2 - V_0 )

       u_inlet = V_1

       CALL init_steady( u_inlet , qp )

       CALL phys_var_qp(qp)

       flag_output = .FALSE.

       zeta = z0
       CALL update_radius(zeta)
       initial_mass_flow_rate = REAL(rho_mix * u_mix)*3.14*radius*radius
       
       CALL integrate_equations( qp , flag_output , extrap_z , extrap_z_p ,     &
            extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

       CALL update_radius(zeta)
       final_mass_flow_rate = REAL(rho_mix * u_mix)*3.14*radius*radius

       WRITE(*,*) 'Initial mass flow rate = ',initial_mass_flow_rate,' kg/s'
       WRITE(*,*) 'Final mass flow rate   = ',final_mass_flow_rate,' kg/s'

       WRITE(*,*) 'Extrap_flag:',extrap_flag
       WRITE(*,*) 'Extrap_z_p',extrap_z_p
       WRITE(*,*) 'Extrap_z_mach',extrap_z_mach

       r_p = MIN( r_p_1 , r_p_2 )

       IF ( increase_flow_rate ) THEN

          WRITE(*,*) 'Pressure 1 at exit:', r_p_1
          WRITE(*,*) 'Pressure 2 at exit:', r_p_2
          WRITE(*,*) 'Mach at exit:', mach


          IF ( extrap_z_mach .LE. extrap_z_p ) THEN

             WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

          ELSE

             WRITE(*,*) 'Extrap_flag:',extrap_flag
             WRITE(*,*) 'Exit pressure is reached at:',extrap_z

          END IF

          IF ( ABS(r_p - p_out) .LT. eps_conv * p_out ) THEN

             check_flow_rate = .TRUE.

          END IF

       ELSE 

          IF ( ( ABS( mach - 1.D0 ) .LT. eps_conv ) .AND. ( zeta .EQ. zN ) ) THEN

             WRITE(*,*) 'Mach at exit:', mach
             WRITE(*,*) 'Pressure 1 at exit:', r_p_1
             WRITE(*,*) 'Pressure 2 at exit:', r_p_2

             check_flow_rate = .TRUE.

          ELSE

             IF ( zeta .GE. zN ) THEN

                WRITE(*,*) 'Mach at exit:', mach
                WRITE(*,*) 'Pressure 1 at exit:', r_p_1
                WRITE(*,*) 'Pressure 2 at exit:', r_p_2

             ELSE

                WRITE(*,*) 'At zeta ',zeta
                WRITE(*,*) 'Pressure 1 is: ',r_p_1
                WRITE(*,*) 'Pressure 2 is: ',r_p_2
                WRITE(*,*) 'Mach is:', mach

                IF ( extrap_flag .EQ. 3 ) THEN

                   WRITE(*,*) 'Mach = 1 is reached at:',extrap_z

                ELSE

                   WRITE(*,*) 'Extrap_flag:',extrap_flag
                   WRITE(*,*) 'Exit pressure is reached at:',extrap_z

                END IF

             END IF

          END IF

       END IF

       IF ( increase_flow_rate ) THEN

          V_0 = V_1


          extrap_z_0 = extrap_z

       ELSE

          V_2 = V_1

          extrap_z_2 = extrap_z

       END IF


       IF ( ABS(r_p - p_out) .GT. p_out ) THEN

          check_flow_rate = .FALSE.

       END IF


       IF ( iter .GT. 20 ) THEN

          IF ( ( ( V_2 - V_0 ) / V_0 .LT. eps_conv ) .AND.                      &
               ( zeta .GE. zN - (iter-20)/100.D0 ) .AND. ( extrap_flag .NE. 4 ) &
               .AND. ( extrap_flag .NE. -3 ) ) THEN

             WRITE(*,*) 'Relative change in flow rate', ( V_2 - V_0 ) / V_0

             check_flow_rate = .TRUE.

          END IF

          IF ( ( ( V_2 - V_0 ) / V_0 .LT. eps_conv ) .AND.                      &
               ( zeta .GE. zN - (iter-20)/100.D0 ) .AND.                        &
               ( ABS(mach - 1.0) .LT. 0.01) ) THEN

             WRITE(*,*) 'Relative change in flow rate', ( V_2 - V_0 ) / V_0

             check_flow_rate = .TRUE.

          END IF

       END IF

       IF ( iter .GT. 100 ) THEN

          IF ( ( zeta .GE. zN - (iter-20)/100.D0 ) .AND. ( extrap_flag .NE. 4 ) &
               .AND. ( extrap_flag .NE. -3 ) ) THEN

             WRITE(*,*) 'Relative change in flow rate', ( V_2 - V_0 ) / V_0

             check_flow_rate = .TRUE.

          END IF

          IF ( ( zeta .GE. zN - (iter-20)/100.D0 ) .AND. ( ABS(r_p - p_out)     &
               .LT. (1.0D0 + (iter-50)/10.D0) * p_out) ) THEN

             WRITE(*,*) 'Relative change in flow rate', ( V_2 - V_0 ) / V_0

             check_flow_rate = .TRUE.

          END IF

       END IF

       IF ( ( ( V_2 - V_0 ) / V_0 .LT. eps_conv ) .AND.                         &
            ( zeta .GE. zN ) .AND. ( extrap_flag .NE. 4 )                       &
            ! .AND. ( extrap_flag .NE. -3 )                                     & 
          ) THEN

          WRITE(*,*) 'Relative change in flow rate', ( V_2 - V_0 ) / V_0

          check_flow_rate = .TRUE.

       END IF

       IF ( verbose_level .GE. 1 ) READ(*,*)

    END DO

    ! --- Integrate again with the found value of u_inlet saving the output -----

    flag_output = .TRUE.
    zeta_exit = zeta   
    u_inlet = V_1

    CALL init_steady( u_inlet , qp )

    CALL phys_var_qp(qp)

    WRITE(*,*) 'Number of iterations',iter
    WRITE(*,*) 'Saving solution for V at the inlet',V_1

    CALL integrate_equations(qp , flag_output , extrap_z , extrap_z_p ,         &
         extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

  END SUBROUTINE steady_shooting

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> @brief Steady state integration
  !
  !> This subroutine integrate the steady equations for magma ascent along the 
  !> conduit, givent the state qp at the base of the conduit. The integration
  !> ends if one of the top boundary condition is reached befor the top or if the
  !> top of the conduit is reached.
  !> \param[in,out] qp           physical variables at the inlet
  !> \param[in]     flag_output  flag to determine if the output has to be saved
  !> \param[out]    extrap_z     depht at which the boundary condition is reached
  !> \param[out]    extrap_z_p   depht at which the boundary condition is reached
  !> \param[out]   extrap_z_mach depht at which the boundary condition is reached
  !> \param[out]    extrap_flag  flag for the boundary condition reached
  !> \param[out]    r_p_1        exit phase 1 pressure
  !> \param[out]    r_p_2        exit phase 2 (exsolved gas) pressure
  !> \param[out]    mach         Mach at the exit 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE integrate_equations( qp , flag_output , extrap_z , extrap_z_p ,    &
       extrap_z_mach , extrap_flag , r_p_1 , r_p_2 , mach )

    ! external procedures
    USE constitutive, ONLY : eos
    USE constitutive, ONLY : sound_speeds

    USE constitutive, ONLY : eval_densities
    USE constitutive, ONLY : f_alfa3, f_alfa

    USE equations, ONLY : eval_nonhyperbolic_terms_qp
    USE equations, ONLY : eval_fluxes_qp
    USE equations, ONLY : phys_var_qp

    
    USE geometry, ONLY : update_radius

    ! external variables
    USE constitutive, ONLY : frag_thr, frag_eff
    USE constitutive, ONLY : rho_mix, u_mix

    USE inpout, ONLY : output_steady

    ! external variables
    USE constitutive, ONLY : zeta_lith

    USE geometry, ONLY : z_comp , comp_cells, radius

    USE parameters, ONLY : shooting

    USE parameters, ONLY : tol_abs , tol_rel, n_gas

    IMPLICIT NONE

    REAL*8, INTENT(INOUT) :: qp(n_eqns)
    LOGICAL, INTENT(IN) :: flag_output
    REAL*8,INTENT(OUT) :: extrap_z , extrap_z_p , extrap_z_mach
    INTEGER,INTENT(OUT) :: extrap_flag
    REAL*8, INTENT(OUT) :: r_p_1 , r_p_2
    REAL*8, INTENT(OUT) :: mach

    !> Integration step
    REAL*8 :: dz

    !> Sound speed of the mixture
    REAL*8 :: C_mix

    !> Variables used to define the new value of the integration step
    REAL*8 :: check_error , check_error_old, coeff_z

    REAL*8 :: r_alfa_2

    LOGICAL :: fragmentation

    INTEGER :: idx_zeta

    REAL*8 :: coeff_zeta

    !> Solution at the computational grid point
    REAL*8 :: qp_comp(n_eqns)

    !> Solution at the previous integration step
    REAL*8 :: qp_old(n_eqns)

    !> Solution obtained integrating with full step 
    REAL*8 :: qp_full(n_eqns)

    !> Solution obtained integrating with one half step
    REAL*8 :: qp_half(n_eqns)

    !> Solution obtained integrating with two half steps
    REAL*8 :: qp_half2(n_eqns)

    REAL*8 :: fluxes_temp(n_vars) , nh_terms_temp(n_vars)

    REAL*8 :: dqp_dz(n_vars)
    !    REAL*8 :: dqp_dz_half(n_vars)

    LOGICAL :: check_convergence

    REAL*8 :: qp_guess(n_vars)

    REAL*8 :: delta_full , delta_half2

    LOGICAL :: fragmentation_half2, fragmentation_full

    REAL*8 :: alfa_2_half2, alfa_2_full, alfa_2_qp

    REAL*8 :: strain_rate_qp

    REAL*8 :: u_1_old

    INTEGER :: counter

    delta_full = 0.D0
    delta_half2 = 0.D0


    zeta = z0

    idx_zeta = 1

    IF ( flag_output) shooting = .FALSE.

    CALL update_radius(zeta)

    IF ( ( .NOT. shooting ) .OR. ( flag_output ) ) CALL output_steady( zeta , &
         qp , radius )

    dz_max = ( zN - z0 ) / comp_cells

    dz = dz_max

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'qp = '
       WRITE(*,*) qp

       READ(*,*)

    ENDIF

    lateral_degassing = .FALSE.
    fragmentation = .FALSE.

    check_error = 1.D0
    check_error_old = 1.D0

    dqp_dz(1:n_vars) = 0.D0 

    counter = 1

    frag_eff = 0.D0

    u_1_old = qp(idx_u1) 
    strain_rate_qp = 0.0

    zeta_integration:DO

       counter = counter + 1

       IF ( counter .GT. 10 * comp_cells ) THEN

          WRITE(*,*)'Convergence error: Too many iterations'
          WRITE(*,*)'counter = ', counter
          STOP

       END IF

       zeta_old = zeta

       qp_old = qp

       u_1_old = qp_old(idx_u1) 

       zeta_lith = zeta

       CALL update_radius(zeta)
       
       IF ( verbose_level .GE. 1 ) THEN 
       
          WRITE(*,*) 'zeta, radius, counter',zeta, radius, counter

       END IF

       CALL eval_fluxes_qp( r_qp = qp_old , r_flux = fluxes_old )

       CALL eval_nonhyperbolic_terms_qp( r_qp = qp_old , r_nh_term_impl =       &
            nh_terms_old )

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'Non hyperbolic term = '
          WRITE(*,*) nh_terms_old/(radius**2)
          READ(*,*)

       END IF

       IF ( isothermal ) THEN

          fluxes_old(idx_mix_engy_eqn) = qp_old(idx_T)
          nh_terms_old(idx_mix_engy_eqn) = 0.D0

       END IF


       dz = MIN( dz , zN - zeta_old )

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'z = ',zeta,'dz = ',dz

       find_step_loop:DO 

          zeta = zeta_old + dz
          zeta_lith = zeta 
          CALL update_radius(zeta)

          ! ------------integrate the whole step --------------------------------

          IF ( verbose_level .GE. 1 ) THEN
             WRITE(*,*)''
             WRITE(*,*) '/--------full step---------/'
          END IF

          qp_full = qp_old

          qp_guess = qp_full

          CALL perturbe_qp(qp_full)

          CALL advance_dz( qp_full , dz , check_convergence )

          IF ( .NOT. check_convergence ) THEN

             qp_full = qp_old + dqp_dz * dz

             CALL perturbe_qp(qp_full)

             CALL advance_dz( qp_full , dz , check_convergence )

          END IF

          fluxes_temp = fluxes_old
          nh_terms_temp = nh_terms_old

          IF ( check_convergence ) THEN

             ! ------------integrate the first half step ------------------------

             zeta = zeta_old + 0.5D0 * dz
             zeta_lith = zeta 
             CALL update_radius(zeta)


             IF ( verbose_level .GE. 1 ) THEN
                WRITE(*,*)''
                WRITE(*,*) '/--------fist half step---------/'
             END IF

             qp_half(1:n_vars) = qp_old(1:n_vars)

             CALL perturbe_qp(qp_half)

             CALL advance_dz( qp_half , 0.5 * dz , check_convergence )

             ! if the integration failed try with a different initial guess:
             ! use an average between qp_old and qp_full
             IF ( .NOT. check_convergence ) THEN

                qp_half(1:n_vars) = 0.5D0 * (qp_old(1:n_vars) +                 &
                     qp_full(1:n_vars) )

                CALL perturbe_qp(qp_half)

                CALL advance_dz( qp_half , 0.5 * dz , check_convergence )

             END IF

             ! if the integration failed try with a different initial guess:
             ! use qp_full
             IF ( .NOT. check_convergence ) THEN

                qp_half = qp_full

                CALL perturbe_qp(qp_half)

                CALL advance_dz( qp_half , 0.5 * dz , check_convergence )

             END IF


          END IF

          IF ( check_convergence ) THEN

             ! ------------integrate the second half step -----------------------
             
             zeta = zeta_old + dz

             IF ( verbose_level .GE. 1 ) THEN
                WRITE(*,*)''
                WRITE(*,*) '/--------second half step---------/'
             END IF

             ! use the values from the first step as old values
             CALL update_radius(zeta_old + 0.5D0 * dz)
             CALL eval_fluxes_qp( r_qp = qp_half , r_flux = fluxes_old )

             CALL eval_nonhyperbolic_terms_qp( r_qp = qp_half , r_nh_term_impl =&
                  nh_terms_old )

             IF ( isothermal ) THEN

                fluxes_old(idx_mix_engy_eqn) = qp_old(idx_T)
                nh_terms_old(idx_mix_engy_eqn) = 0.D0

             END IF

             zeta_lith = zeta 
             CALL update_radius(zeta)

             ! use the solution of the whole step as initial guess
             qp_half2 = qp_half

             CALL perturbe_qp(qp_half2)

             CALL advance_dz( qp_half2 , 0.5 * dz , check_convergence )

             ! if the integration failed try with a different initial guess:
             ! use the solution from the whole step
             IF ( .NOT. check_convergence ) THEN

                qp_half2 = qp_full

                CALL perturbe_qp(qp_half2)

                CALL advance_dz( qp_half2 , 0.5 * dz , check_convergence )

             END IF

             ! if the integration failed again try with another initial guess:
             ! use a solution extrapolated from the first half step 
             IF ( .NOT. check_convergence ) THEN

                qp_half2 = qp_old + 0.50 * dz * ( qp_half - qp_old )

                CALL perturbe_qp(qp_half2)

                CALL advance_dz( qp_half2 , 0.5 * dz , check_convergence )

             END IF

          END IF

          IF ( verbose_level .GE. 1 ) THEN
             WRITE(*,*)''
             WRITE(*,*) '<<<<<<<<< check_convergence >>>>>>>> ',check_convergence
             READ(*,*)
          END IF
          
          IF ( check_convergence ) THEN
             
             ! -------- compare the solutions obtained with dz and twice dz/2 ---
             
             check_error = MAXVAL(ABS( qp_full(1:n_vars) -                      &
                  qp_half2(1:n_vars) ) / ( tol_abs + tol_rel *                  &
                  MAX( ABS( qp_full(1:n_vars) ) ,                               &
                  ABS( qp_half2(1:n_vars) ) ) ) ) 

             delta_full = MAXVAL(ABS( qp_full(1:n_vars) -                       &
                  qp_guess(1:n_vars) ) / ( tol_abs + tol_rel *                  &
                  MAX( ABS( qp_full(1:n_vars) ) ,                               &
                  ABS( qp_guess(1:n_vars) ) ) ) ) 

             delta_half2 = MAXVAL(ABS( qp_half2(1:n_vars) -                     &
                  qp_guess(1:n_vars) ) / ( tol_abs + tol_rel *                  &
                  MAX( ABS( qp_half2(1:n_vars) ) ,                              &
                  ABS( qp_guess(1:n_vars) ) ) ) ) 

             IF ( verbose_level .GE. 1 ) THEN

                WRITE(*,*) 'check dz and dz/2 error',check_error

             END IF

             IF ( MAX( qp_half2(idx_p1) , qp_half2(idx_p2) )  &
                  .LT. p_out) THEN

                check_convergence = .FALSE.
                IF ( verbose_level .GE. 1 ) WRITE(*,*) 'pressure',              &
                     MAX( qp_half2(idx_p1) , qp_half2(idx_p2) )

             ELSE

                IF ( verbose_level .GE. 2 ) THEN
                   WRITE(*,*)''
                   WRITE(*,*) 'qp_half2'
                   WRITE(*,*) qp_half2
                   WRITE(*,*)''
                   WRITE(*,*) 'qp_full'
                   WRITE(*,*) qp_full
                   WRITE(*,*)''
                   WRITE(*,*) 'alfa_2',SUM(qp_half2(1:n_gas))
                   READ(*,*)

                END IF

             END IF

          END IF


          ! ---- Check if the fragmentation threshold is reached ---------------

          IF ( EXPLOSIVE ) THEN

             alfa_2_half2 = SUM(qp_half2(idx_alfa_first:idx_alfa_last)) 
             alfa_2_full = SUM(qp_full(idx_alfa_first:idx_alfa_last)) 
             alfa_2_qp = SUM(qp(idx_alfa_first:idx_alfa_last))

             IF ( alfa_2_half2 .GT. frag_thr ) THEN

                fragmentation_half2 = .TRUE.

             ELSE

                fragmentation_half2 = .FALSE.

             END IF

             IF ( alfa_2_full .GT. frag_thr ) THEN

                fragmentation_full = .TRUE.

             ELSE

                fragmentation_full = .FALSE.

             END IF


             IF ( ( .NOT.fragmentation )                                        &
                  .AND. ( fragmentation_half2 .OR. fragmentation_full )         &
                  .AND. ( frag_thr - alfa_2_qp .GT. 1.D-4 )                     &
!                 .AND. ( MAX(alfa_2_half2,alfa_2_full) - frag_thr .GT. 1.D-4 ) &
                ) THEN

                check_convergence = .FALSE.

             END IF

          END IF

          IF ( ( check_error .LT. 1.D0 ) .AND. ( check_convergence ) ) THEN

             ! --- if the error is small accept the solution obtained with dz/2
             ! --- and exit from the search loop

             coeff_z = MAX( 1.05D0 , MIN( 1.50D0 , check_error ** ( -0.35D0 ) * &
                  check_error_old ** ( 0.2D0 ) ) )

             check_error_old = check_error

             ! evaluate the slope of the solution for the guess at the next step

             IF ( delta_half2 .LT. delta_full ) THEN

                qp = qp_half2
                dqp_dz = ( qp_half2 - qp_old ) / dz

             ELSE

                qp = qp_full
                dqp_dz = ( qp_full - qp_old ) / dz

             END IF

             CALL phys_var_qp(qp)

             IF ( verbose_level .GE. 1 ) THEN

                WRITE(*,*) '::::::::: checks ok :::::::::'
                WRITE(*,*) 'zeta = ',zeta,' check_error = ',check_error,        &
                     ' coeff_z = ',coeff_z
                WRITE(*,*) ''
                WRITE(*,*) 'qp ='
                WRITE(*,*) qp
                WRITE(*,*) ''
                WRITE(*,*) 'alfa_2 =',SUM(qp(idx_alfa_first:idx_alfa_last)),    &
                     ' fragmentation = ', frag_eff, 'Mass flow rate = ',        &
                     REAL(rho_mix * u_mix)*3.14*radius*radius
                WRITE(*,*)''
                READ(*,*)
             END IF

             EXIT find_step_loop

          ELSE


             ! --- if the error is big repeat the step with smaller dz

             dz = 0.950D0 * dz

             fluxes_old = fluxes_temp
             nh_terms_old = nh_terms_temp

             IF ( verbose_level .GE. 1 ) WRITE(*,*) 'zeta_old',zeta_old,        &
                  'dz = ',dz

             IF ( dz .LT. 1E-20 ) THEN

                WRITE(*,*)'Convergence Error: dz too small'
                WRITE(*,*)'dz =', dz
                STOP

             END IF

          END IF

       END DO find_step_loop

       ! ---------- Write the output on file when the solution is found ---------

       IF ( flag_output ) THEN

          IF ( shooting ) THEN

             DO WHILE ( ( z_comp(idx_zeta) .LT. zeta ) .AND.                    &
                  ( idx_zeta .LE. comp_cells ) )

                coeff_zeta = (  z_comp(idx_zeta) - ( zeta - dz ) ) / dz 

                qp_comp = coeff_zeta * qp + ( 1.D0 - coeff_zeta ) * qp_old

                CALL output_steady(z_comp(idx_zeta),qp_comp,radius)

                idx_zeta = idx_zeta + 1

             END DO

             IF ( zeta .EQ. zN ) CALL output_steady(zeta,qp,radius)

          ELSE

             CALL output_steady(zeta,qp,radius)

          END IF

       END IF

       ! ----- Update the integration step with the coefficient coeff_z ---------

       dz = MIN( coeff_z * dz , dz_max )

       ! ----- Evaluate some phys. variables to pass out of the subroutine ------

       r_alfa_2 = SUM(qp_half2(idx_alfa_first:idx_alfa_last))

       r_p_1 = qp(idx_p1)
       r_p_2 = qp(idx_p2)

       CALL phys_var_qp( r_qp = qp )
       CALL eos
       CALL sound_speeds( C_mix , mach ) 

       ! ------ Check if the extrapolated solution or the solution at zeta ------
       ! ------ reach the boundary conditions -----------------------------------

       CALL linear_extrapolation(zeta_old,qp_old,zeta,qp,extrap_z,extrap_z_p,   &
            extrap_z_mach,extrap_flag)

       IF ( extrap_flag .GT. 0 ) THEN

          increase_flow_rate = .FALSE.

          RETURN

       END IF

       ! ----- Check if the top is reached --------------------------------------

       IF ( zeta .GE. zN ) THEN

          increase_flow_rate = .TRUE.

          RETURN

       END IF

       ! ---- Check if the fragmentation threshold is reached ---------------

       IF ( EXPLOSIVE ) THEN

          alfa_2_qp = SUM(qp(idx_alfa_first:idx_alfa_last))

          IF ( ( alfa_2_qp .GT. frag_thr) .AND.                                 &	 
               ( .NOT. fragmentation ) )THEN

             frag_eff = 1.0D0
             fragmentation = .TRUE.

             WRITE(*,*) 'Fragmentation at z = ',zeta
             ! verbose_level = 3

          END IF

       END IF

       ! ---- Check if the lateral degassing threshold is reached ---------------

       IF ( ( r_alfa_2 .GE. alfa2_lat_thr ) .AND.                               &
            ( lateral_degassing_flag ) .AND. ( .NOT.lateral_degassing ) ) THEN

          lateral_degassing = .TRUE.

          WRITE(*,*) 'Lateral degassing from z = ',zeta

       END IF

    END DO zeta_integration

    !WRITE(*,*) counter

  END SUBROUTINE integrate_equations

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> @brief Solution preturbation
  !
  !> This subroutine perturbes the initial guess for the Newthon-Raphson method
  !> in order to improve the convergence.
  !> \param[in,out] qp           physical variables 
  !> \date 08/07/2013
  !******************************************************************************

  SUBROUTINE perturbe_qp(qp)

    IMPLICIT NONE

    REAL*8, INTENT(INOUT) :: qp(n_vars)

    !Pressure phase 2 
    qp(idx_p2) = qp(idx_p2) + 1.D-2

    !Velocity phase 1 
    IF ( qp(idx_u1) - qp(idx_u2) .GT. -1D-7 )  THEN

       qp(idx_u2) = qp(idx_u2) * ( 1.0001D0 )

    END IF
    
  END SUBROUTINE perturbe_qp

  !******************************************************************************
  !> \brief Solution advance in space
  !
  !> This subroutine integrates the conservative variables in space from zeta to
  !> zeta+dz.
  !> \param[in,out]  qp                  actual physical variables
  !> \param[in]      dz                  integration step
  !> \param[out]     check_convergence   logical for convergence check
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE advance_dz( qp , dz , check_convergence )

    ! external variables
    USE parameters, ONLY : alfa_impl

    ! external subroutines
    USE constitutive, ONLY : eos
    USE equations, ONLY : phys_var_qp

    IMPLICIT none

    REAL*8, INTENT(INOUT) :: qp(n_eqns)
    REAL*8, INTENT(IN) :: dz
    LOGICAL, INTENT(OUT) :: check_convergence

    REAL*8 :: qp_rel(n_vars)
    REAL*8 :: qp_org(n_vars)

    REAL*8 :: left_matrix(n_eqns,n_eqns)
    REAL*8 :: right_term(n_eqns)

    REAL*8 :: scal_f , scal_f_old

    INTEGER :: pivot(n_eqns)
    INTEGER :: ok

    INTEGER :: idx_max(1)

    REAL*8 :: delta_qp_rel(n_eqns)
    REAL*8 :: grad_f(n_eqns)

    REAL*8, DIMENSION(size(qp)) :: qp_NR_old , desc_dir

    REAL*8 :: qp_rel_NR_old(n_vars)

    REAL*8 :: check_NR_error

    REAL*8 :: scal_f_init

    LOGICAL :: opt_search_NL
    REAL*8, PARAMETER :: STPMX=100.D0
    REAL*8 :: stpmax
    LOGICAL :: check 
    !    LOGICAL :: check_phys_var

    REAL*8, PARAMETER :: tol_rel_NR = 1.D-20 , tol_abs_NR = 1.D-20


    LOGICAL :: normalize_qp
    LOGICAL :: normalize_f

    REAL*8 :: coeff_f(n_eqns)

    REAL*8 :: arg_check(n_eqns)

    INTEGER :: i,j

    check_convergence = .FALSE.

    normalize_qp = .TRUE.
    normalize_f = .FALSE.

    IF ( normalize_qp ) THEN

       DO i = 1,n_vars

          IF ( qp(i) .EQ. 0.D0 ) THEN

             qp_org(i) = 1.D0

          ELSE

             qp_org(i) = qp(i)

          END IF

       END DO

    ELSE 

       qp_org(1:n_vars) = 1.0D0

    END IF

    qp_rel = qp / qp_org

    IF ( normalize_qp ) THEN

       qp_rel = qp_rel * 1.00001D0

    END IF


    coeff_f(1:n_eqns) = 1.D0

    IF ( normalize_f ) THEN

       CALL eval_f( qp_rel , qp_org , dz , coeff_f , right_term , scal_f )

       DO i = 1,n_vars

          IF ( ABS(right_term(i)) .GE. 1.D-5 ) THEN

             coeff_f(i) = 1.D0 / right_term(i)

          END IF

       END DO

    END IF

    IF ( verbose_level .GE. 3 ) THEN

       CALL eval_f( qp_rel , qp_org , dz , coeff_f , right_term , scal_f )
       WRITE(*,*) 'before iter'
       WRITE(*,*) 'right_term',right_term
       WRITE(*,*) 'qp'
       WRITE(*,*) qp

    END IF

    DO nl_iter = 1,max_nl_iter

       CALL eval_jacobian( qp_rel , qp_org , dz , coeff_f , left_matrix )

       CALL eval_f( qp_rel , qp_org , dz , coeff_f , right_term , scal_f )

       IF ( nl_iter .EQ. 1 ) scal_f_init = scal_f

       delta_qp_rel = right_term

       call DGESV(n_eqns, 1, left_matrix , n_eqns, pivot, delta_qp_rel , n_eqns,&
            ok)

       IF ( ok .EQ. 0 ) THEN

          qp_rel_NR_old = qp_rel
          qp_NR_old = qp_rel * qp_org
          scal_f_old = scal_f
          desc_dir = - delta_qp_rel

          opt_search_NL = .TRUE.

          IF ( ( opt_search_NL ) .AND. ( nl_iter .GT. 0 ) ) THEN

             stpmax = STPMX * MAX( DSQRT( DOT_PRODUCT(qp_rel,qp_rel) ) ,        &
                  DBLE(SIZE(qp_rel)) )

             grad_f = MATMUL( right_term , left_matrix )

             CALL steady_lnsrch( qp_rel_NR_old , qp_org , scal_f_old , grad_f , &
                  desc_dir , dz , coeff_f , qp_rel , scal_f , right_term ,      &
                  stpmax , check , eval_f )

             IF ( check ) THEN

                desc_dir = - 0.5D0 * delta_qp_rel

                qp_rel = qp_rel_NR_old + desc_dir

                CALL eval_f( qp_rel , qp_org , dz , coeff_f , right_term ,      &
                     scal_f )

             END IF

             DO i=idx_xd_first,idx_xd_last

                qp_rel(i) = MAX(qp_rel(i),0.D0)

             END DO

             DO i=idx_beta_first,idx_beta_last

                qp_rel(i) = MAX(qp_rel(i),0.D0)

             END DO

          ELSE

             qp_rel = qp_rel_NR_old + desc_dir

             DO i=idx_xd_first,idx_xd_last

                qp_rel(i) = MAX(qp_rel(i),0.D0)

             END DO

             DO i=idx_beta_first,idx_beta_last

                qp_rel(i) = MAX(qp_rel(i),0.D0)

             END DO

             CALL eval_f( qp_rel , qp_org , dz , coeff_f , right_term , scal_f )

          END IF

          ! Sometimes it happens that the crystal content becomes negative 
          ! (even if very small, for example -1.e-140). Thus I force the
          ! crystal content and the dissolved gas to be non-negative.


          DO i=idx_xd_first,idx_xd_last
             
             qp_rel(i) = MAX(qp_rel(i),0.D0)
             
          END DO
          
          DO i=idx_beta_first,idx_beta_last
             
             qp_rel(i) = MAX(qp_rel(i),0.D0)
             
          END DO
          
          qp = qp_rel * qp_org

          arg_check = ABS( qp_rel(1:n_vars) - qp_rel_NR_old(1:n_vars) ) /       &
               ( tol_abs_NR + tol_rel_NR * MAX( ABS(qp_rel(1:n_vars)) ,         &
               ABS(qp_rel_NR_old(1:n_vars)) ) )

          check_NR_error = MAXVAL( arg_check )

          idx_max = MAXLOC( arg_check )

          IF ( verbose_level .GE. 3 ) THEN

             WRITE(*,*) 'iter = ',nl_iter                    
             WRITE(*,*) 'right_term'
             WRITE(*,*) right_term
             WRITE(*,*) ''
             WRITE(*,*) 'desc_dir'
             WRITE(*,*) desc_dir
             WRITE(*,*) ''
             WRITE(*,*) 'qp_NR_old'
             WRITE(*,*) qp_NR_old
             WRITE(*,*) ''
             WRITE(*,*) 'qp'
             WRITE(*,*) qp
             WRITE(*,*) ''

             WRITE(*,*) 'scal_f = '
             WRITE(*,*) scal_f_init , scal_f , scal_f/scal_f_old ,              &
                  scal_f/scal_f_init
             WRITE(*,*) ''

             WRITE(*,*) 'check_NR_error = '
             WRITE(*,*) check_NR_error, idx_max,                                &
                  qp_rel(idx_max) * qp_org(idx_max),                            &
                  qp_rel_NR_old(idx_max) * qp_org(idx_max)
             WRITE(*,*) ''
             WRITE(*,*) 'zeta = '
             WRITE(*,*) zeta


             READ(*,*)

          END IF

       ELSE 

          IF ( verbose_level .GE. 3 ) THEN

             WRITE(*,*) 'advance_dz : num_cond too small'
             WRITE(*,*) 'left_matrix'

             DO j=1,n_eqns
                WRITE(*,*) left_matrix(j,:)
             END DO

             READ(*,*)

          END IF

          EXIT

       END IF


       IF ( qp_rel(1) .LT. 0.0D0 ) THEN

          check_convergence = .FALSE.

          EXIT

       END IF


       IF ( check_NR_error .LT. 0.10D0 ) THEN

          IF ( scal_f / scal_f_init .LT. 1.d-10 ) check_convergence = .TRUE. 

          EXIT

       END IF

       IF ( scal_f / scal_f_init .LT. 1.d-14 ) THEN  

          check_convergence = .TRUE.

          EXIT

       END IF

    END DO

    IF ( scal_f / scal_f_init .LT. 1.d-11 ) check_convergence = .TRUE. 

    IF ( SUM(qp_rel(idx_alfa_first:idx_alfa_last)) .LT. 0.0D0 ) THEN

       check_convergence = .FALSE.

    END IF

    IF ( SUM(qp_rel(idx_alfa_first:idx_alfa_last)) .GE. 1.0D0 ) THEN

       check_convergence = .FALSE.

    END IF

    
    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*)''
       WRITE(*,*) 'scal_f = ' , nl_iter , idx_max ,    &
            scal_f/scal_f_init,check_NR_error, idx_max, ok
       WRITE(*,*)''
       WRITE(*,*)'qp_new = '
       WRITE(*,*) qp_rel * qp_org 

    END IF

  END SUBROUTINE advance_dz

  !******************************************************************************
  !> \brief Nonlinear function evaluation
  !
  !> This subroutine evaluates the residual of the nonlinear functions obtained 
  !> from the integration in space of the steady equations.
  !> \param[in]  qp_rel       actual physical variables
  !> \param[in]  qp_org       actual physical variables
  !> \param[in]  dz           integration step
  !> \param[in]  coeff_f      weighting coefficients
  !> \param[out] right_term   array of residuals
  !> \param[out] scal_f       magnitude of the array of residuals
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE eval_f( qp_rel , qp_org , dz , coeff_f , right_term  , scal_f )

    USE equations, ONLY : eval_fluxes_qp
    USE equations, ONLY : eval_nonhyperbolic_terms_qp

    USE parameters, ONLY : alfa_impl

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: qp_rel(:)
    REAL*8, INTENT(IN) :: qp_org(:)
    REAL*8, INTENT(IN) :: dz
    REAL*8, INTENT(IN) :: coeff_f(:)
    REAL*8, INTENT(OUT) :: right_term(:)
    REAL*8, INTENT(OUT) :: scal_f

    REAL*8 :: qp(n_vars)
    REAL*8 :: fluxes(n_eqns)
    REAL*8 :: nh_terms(n_eqns)


    qp = qp_rel * qp_org

    CALL eval_fluxes_qp( r_qp=qp , r_flux=fluxes)

    CALL eval_nonhyperbolic_terms_qp( r_qp = qp , r_nh_term_impl = nh_terms )
    ! correction for isothermal model


    !WRITE(*,*) 'fluxes'
    !WRITE(*,*) fluxes

    !WRITE(*,*) 'nh_terms'
    !WRITE(*,*) nh_terms
    !READ(*,*)

    IF ( isothermal ) THEN

       fluxes(idx_mix_engy_eqn) = qp(idx_T)
       nh_terms(idx_mix_engy_eqn) = 0.D0
       right_term(idx_mix_engy_eqn) = fluxes(idx_mix_engy_eqn)                  &
            - fluxes_old(idx_mix_engy_eqn)

    END IF

    right_term = ( fluxes - fluxes_old ) - dz * ( alfa_impl * nh_terms          &
         + ( 1.d0 - alfa_impl ) * nh_terms_old )     
    
    right_term = right_term * coeff_f

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*) 'eval_f'
       WRITE(*,*) 'fluxes',REAL(fluxes)
       WRITE(*,*) 'fluxes_old',REAL(fluxes_old)
       WRITE(*,*) 'nh_terms',REAL(nh_terms)
       WRITE(*,*) 'nh_terms_old',REAL(nh_terms_old)
       READ(*,*)

    END IF
    
    scal_f = 0.5D0 * DOT_PRODUCT( right_term , right_term )

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Jacobian evaluation
  !
  !> This subroutine evaluates with a complex step derivative procedure the 
  !> Jacobian of the nonlinear system obtained from the integration in space of 
  !> the steady equations.
  !> \param[in]  qp_rel       actual physical variables
  !> \param[in]  qp_org       actual physical variables
  !> \param[in]  dz           integration step
  !> \param[in]  coeff_f      weighting coefficients
  !> \param[out] left_matrix  Jacobian matrix
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE eval_jacobian( qp_rel , qp_org , dz , coeff_f , left_matrix)

    ! external procedures
    USE equations, ONLY : eval_fluxes_qp
    USE equations, ONLY : eval_nonhyperbolic_terms_qp

    USE parameters, ONLY : alfa_impl

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: qp_rel(n_vars)
    REAL*8, INTENT(IN) :: qp_org(n_vars)
    REAL*8, INTENT(IN) :: dz
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: left_matrix(n_eqns,n_eqns)

    REAL*8 :: h

    REAL*8 :: qp(n_vars)
    COMPLEX*16 :: qp_cmplx(n_eqns)
    COMPLEX*16 :: fluxes_cmplx(n_eqns)
    COMPLEX*16 :: nh_terms_cmplx(n_eqns)

    REAL*8 :: Jacob_fluxes(n_eqns,n_eqns)
    REAL*8 :: Jacob_nh(n_eqns,n_eqns)

    INTEGER :: i

    h = n_eqns * epsilon(1.d0)

    qp = qp_rel * qp_org

    qp_cmplx(1:n_eqns) = DCMPLX( qp(1:n_eqns) , 0.D0 )

    DO i=1,n_vars

       qp_cmplx(i) = DCMPLX( qp_rel(i) , h ) * qp_org(i)

       CALL eval_fluxes_qp( c_qp=qp_cmplx , c_flux=fluxes_cmplx )
       CALL eval_nonhyperbolic_terms_qp( c_qp = qp_cmplx , c_nh_term_impl =     &
            nh_terms_cmplx )

       IF ( isothermal ) THEN

          fluxes_cmplx(idx_mix_engy_eqn) = qp_cmplx(idx_T)
          nh_terms_cmplx(idx_mix_engy_eqn) = DCMPLX( 0.D0 , 0.D0)

       END IF

       Jacob_fluxes(1:n_eqns,i) = DIMAG(fluxes_cmplx) / h
       Jacob_nh(1:n_eqns,i) = DIMAG(nh_terms_cmplx) / h

       Jacob_fluxes(1:n_eqns,i) = Jacob_fluxes(1:n_eqns,i) * coeff_f
       Jacob_nh(1:n_eqns,i) = Jacob_nh(1:n_eqns,i) * coeff_f

       qp_cmplx(i) = DCMPLX( qp(i) , 0.D0 )

    END DO

    left_matrix = Jacob_fluxes - dz * alfa_impl * Jacob_nh


  END SUBROUTINE eval_jacobian

  !******************************************************************************
  !> \brief Search the descent stepsize
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> \param[in]     x_rel_init      coefficients for initial guess
  !> \param[in]     x_org           initial guess
  !> \param[in]     scal_f_old      old value of the function f
  !> \param[in]     grad_f          gradient of the f
  !> \param[in,out] desc_dir        descent direction
  !> \param[in]     dz              integration step
  !> \param[in]     coeff_f      weighting coefficients
  !> \param[out]    x_rel_new       optimal solution
  !> \param[out]    scal_f          new value of the function f (at x_new)
  !> \param[out]    f_nl            values of the vector function f_nl
  !> \param[in]     stpmax          max number of steps
  !> \param[out]    check           logical for convergence
  !> \param[in]     callf           name of the subroutine for the function f
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE steady_lnsrch( x_rel_init , x_org , scal_f_old , grad_f , desc_dir &
       , dz , coeff_f , x_rel_new , scal_f , f_nl , stpmax , check , callf )

    IMPLICIT NONE

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: x_rel_init

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: x_org

    !> Gradient at xold
    REAL*8, DIMENSION(:), INTENT(IN) :: grad_f

    !> Value of the function at xold
    REAL*8, INTENT(IN) :: scal_f_old

    !> Descent direction (usually Newton direction)
    REAL*8, DIMENSION(:), INTENT(INOUT) :: desc_dir

    REAL*8, INTENT(IN) :: dz

    REAL*8, INTENT(IN) :: stpmax

    REAL*8, DIMENSION(:), INTENT(IN) :: coeff_f

    !> Updated solution
    REAL*8, DIMENSION(:), INTENT(OUT) :: x_rel_new

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: scal_f

    !> Values of the nonlinear functions at x
    REAL*8, DIMENSION(:), INTENT(OUT) :: f_nl

    !> Output quantity check is false on a normal exit 
    LOGICAL, INTENT(OUT) :: check

    INTERFACE
       SUBROUTINE callf( x_rel , x_org , dz , coeff_f , f_nl , scal_f )

         IMPLICIT NONE

         REAL*8, INTENT(IN) :: x_rel(:)
         REAL*8, INTENT(IN) :: x_org(:)
         REAL*8, INTENT(IN) :: dz
         REAL*8, INTENT(IN) :: coeff_f(:)
         REAL*8, INTENT(OUT) :: f_nl(:)
         REAL*8, INTENT(OUT) :: scal_f

       END SUBROUTINE callf
    END INTERFACE

    REAL*8, PARAMETER :: TOLX=epsilon(x_rel_init)

    INTEGER, DIMENSION(1) :: ndum
    REAL*8 :: ALF , a,alam,alam2,alamin,b,disc
    REAL*8 :: scal_f2
    REAL*8 :: desc_dir_abs
    REAL*8 :: rhs1 , rhs2 , slope, tmplam

    ALF = 1.0d-4

    IF ( size(grad_f) == size(desc_dir) .AND. size(grad_f) == size(x_rel_new)   &
         .AND. size(x_rel_new) == size(x_rel_init) ) THEN

       ndum = size(grad_f)

    ELSE

       WRITE(*,*) 'nrerror: an assert_eq failed with this tag:', 'lnsrch'
       STOP 'program terminated by assert_eq4'

    END IF

    check = .FALSE.

    desc_dir_abs = DSQRT( DOT_PRODUCT(desc_dir,desc_dir) )

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax /           &
         desc_dir_abs  

    slope = DOT_PRODUCT(grad_f,desc_dir)

    alamin = TOLX / MAXVAL( ABS( desc_dir(:)) / MAX( ABS( x_rel_init(:)) ,      &
         1.D0 ) )

    alamin = 1.D-20

    IF ( alamin .EQ. 0.d0) THEN

       x_rel_new(:) = x_rel_init(:)

       WRITE(*,*) 'alam 0'

       RETURN

    END IF

    alam = 1.0D0
    alam2 = alam

    optimal_step_search: DO

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'alam',alam

       END IF

       x_rel_new = x_rel_init + alam * desc_dir

       CALL callf( x_rel_new , x_org , dz , coeff_f , f_nl , scal_f )

       scal_f2 = scal_f

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'x',x_rel_new,x_rel_init
          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          READ(*,*)

       END IF

       IF ( alam < alamin ) THEN   
          ! convergence on Delta_x

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) ' convergence on Delta_x',alam,alamin

          END IF

          x_rel_new(:) = x_rel_init(:)
          scal_f = scal_f_old
          check = .TRUE.

          EXIT optimal_step_search

       ELSE IF ( scal_f .LE. scal_f_old + ALF * alam * slope ) THEN   

          EXIT optimal_step_search

       ELSE  

          IF ( alam .EQ. 1.D0 ) THEN

             tmplam = - slope / ( 2.0D0 * ( scal_f - scal_f_old - slope ) )

          ELSE

             rhs1 = scal_f - scal_f_old - alam*slope
             rhs2 = scal_f2 - scal_f_old - alam2*slope

             a = ( rhs1/alam**2.D0 - rhs2/alam2**2.D0 ) / ( alam - alam2 )
             b = ( -alam2*rhs1/alam**2 + alam*rhs2/alam2**2 ) / ( alam - alam2 )

             IF ( a .EQ. 0.D0 ) THEN

                tmplam = - slope / ( 2.0D0 * b )

             ELSE

                disc = b*b - 3.0D0*a*slope

                IF ( disc .LT. 0.D0 ) THEN

                   tmplam = 0.5D0 * alam

                ELSE IF ( b .LE. 0.D0 ) THEN

                   tmplam = ( - b + DSQRT(disc) ) / ( 3.D0 * a )

                ELSE

                   tmplam = - slope / ( b + DSQRT(disc) )

                ENDIF

             END IF

             IF ( tmplam .GT. 0.5D0*alam ) tmplam = 0.5D0 * alam

          END IF

       END IF

       alam2 = alam
       scal_f2 = scal_f
       alam = MAX( tmplam , 0.1D0*alam )

    END DO optimal_step_search

  END SUBROUTINE steady_lnsrch

  !******************************************************************************
  !> \brief Linear extrapolation of the solution
  !
  !> This subroutine extrapolate linearly some primitive variables of the 
  !> solution to find where the boundary conditions are reached (below or above
  !> the vent).
  !> \param[in]  zeta_old     zeta at the previous step
  !> \param[in]  qp_old       physical variables at the previous step
  !> \param[in]  zeta         actual zeta
  !> \param[in]  qp           physical conservative variables
  !> \param[out] extrap_z     zeta at which the boundary condition is reached
  !> \param[out] extrap_z_p     zeta at which the boundary condition is reached
  !> \param[out] extrap_z_mach     zeta at which the boundary condition is reached
  !> \param[out] extrap_flag  flag for the boundary condition reached first 
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE linear_extrapolation(zeta_old,qp_old,zeta,qp,extrap_z,extrap_z_p,  &
       extrap_z_mach,extrap_flag)

    ! external procedures
    USE constitutive, ONLY : eos , alfa_2
    USE constitutive, ONLY : sound_speeds
    USE equations, ONLY : phys_var_qp
    
    ! external variables
    USE constitutive, ONLY : x_d_md, x_ex_dis_in
    USE init, ONLY : p1_in , p2_in

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: zeta_old
    REAL*8, INTENT(IN) :: qp_old(n_eqns)
    REAL*8, INTENT(IN) :: zeta
    REAL*8, INTENT(IN) :: qp(n_eqns)

    REAL*8, INTENT(OUT) :: extrap_z
    REAL*8, INTENT(OUT) :: extrap_z_p
    REAL*8, INTENT(OUT) :: extrap_z_mach
    INTEGER, INTENT(OUT) :: extrap_flag

    REAL*8 :: r_p_2 , r_p_2_old
    REAL*8 :: r_p_1 , r_p_1_old
    REAL*8 :: mach , mach_old

    REAL*8 :: C_mix

    REAL*8 :: p_check

    REAL*8 :: extrap_z_2 
    REAL*8 :: extrap_z_1

    REAL*8 :: coeff_p_1 , coeff_p_2

    REAL*8 :: p_out_1 , p_out_2

    REAL*8 :: eps_extrap

    REAL*8 :: r_alfa_2

    eps_extrap = 1.D-10

    extrap_flag = 0

    !--- check on small values of the phisical variables ------------------------


    r_p_1 = qp(idx_p1)
    r_p_2 = qp(idx_p2)

    IF ( qp(n_vars) .GT. 0.D0 ) THEN

       p_out_1 = p_out
       p_out_2 = p_out

    ELSE

       p_out_1 = p_out
       p_out_2 = p_out

    END IF

    CALL phys_var_qp( r_qp = qp )
    CALL eos
    CALL sound_speeds(C_mix,mach) 



    IF ( verbose_level .GT. 1 ) THEN

       WRITE(*,*) 'p1,p2,mach',r_p_1,r_p_2,mach

    END IF

    ! ---- Check if the solution at zeta has already reached a boundary condition

    p_check = MIN(r_p_2,r_p_1)

    r_alfa_2 = REAL(alfa_2)

    IF ( ( r_p_1 .LT. p_out_1 ) .OR. ( r_p_2 .LT. p_out_2 ) .OR.                &
         ( mach .GT. 1.d0 ) ) THEN

       IF ( verbose_level .GE. -1 ) THEN

          WRITE(*,*) 'Boundary conditions before the top' 
          WRITE(*,*) 'z = ',zeta,'Pressures = ',r_p_1,r_p_2
          WRITE(*,*) 'Mach =',mach

       END IF

       extrap_z = zeta
       extrap_flag = 4

       RETURN

    END IF

    IF ( r_alfa_2 .GE. 1.D0 ) THEN
       
       IF ( verbose_level .GE. -1 ) THEN
          
          WRITE(*,*) 'Boundary conditions before the top' 
          WRITE(*,*) 'z = ',zeta,'alfa gas = ',r_alfa_2
                    
       END IF
       
       extrap_z = zeta
       extrap_flag = 4
       
       RETURN
       
    END IF

    IF ( (ABS(zeta - zeta_old) .LT. 1E-11) .AND.                                &
         (r_p_1 .LT. 1.0E+6) ) THEN

       IF ( verbose_level .GE. -1 ) THEN

          WRITE(*,*) 'Pressure Conditions reached before the exit' 
          !READ(*,*)

       END IF

       extrap_z = zeta
       extrap_flag = 5

       RETURN

    END IF

    IF ( (ABS(zeta - zeta_old) .LT. 1E-6) .AND.        &
         ( SUM(REAL(x_d_md)) .LT. SUM(x_ex_dis_in) ) ) THEN

       IF ( verbose_level .GE. -1 ) THEN

          WRITE(*,*) 'Dissolved gas Conditions reached before the exit' 
          !READ(*,*)

       END IF

       extrap_z = zeta
       extrap_flag = 6

       RETURN

    END IF



    ! --------------------------------------------------------------------------

    r_p_1_old = qp_old(idx_p1)
    r_p_2_old = qp_old(idx_p2)

    CALL phys_var_qp( r_qp = qp_old )
    CALL eos
    CALL sound_speeds(C_mix,mach_old) 

    ! ---- gas pressure extrapolation -------------------------------------------
    IF ( r_p_2 - r_p_2_old .LT. 0.d0 ) THEN

       extrap_z_2 = zeta_old - ( zeta - zeta_old ) / ( r_p_2 - r_p_2_old) *     &
            ( r_p_2_old )

       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,*) 'gas pressure extrapolation'
          WRITE(*,*) 'extrap_z_2 =',extrap_z_2,r_p_2

       END IF

    ELSE

       extrap_z_2 = zN

    END IF

    ! --- liquid pressure extrapolation -----------------------------------------
    IF ( r_p_1 - r_p_1_old .LT. 0.d0 ) THEN

       extrap_z_1 = zeta_old + ( zeta - zeta_old ) / ( r_p_1 - r_p_1_old) *     &
            ( p_out_1 - r_p_1_old )

       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,*) 'liquid pressure extrapolation'
          WRITE(*,*) 'extrap_z_1 =',extrap_z_1,r_p_1

       END IF

    ELSE

       extrap_z_1 = zN

    END IF

    extrap_z_p = MIN( extrap_z_1 , extrap_z_2 )

    ! --- Mach number extrapolation ---------------------------------------------
    IF ( mach - mach_old .GT. 0.D0 ) THEN

       extrap_z_mach = zeta_old + ( zeta - zeta_old ) / ( mach - mach_old ) *   &
            ( 1.d0 - mach_old )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'Mach number extrapolation'
          WRITE(*,*) 'extrap_z_mach',extrap_z_mach,mach

       END IF

    ELSE

       extrap_z_mach = zN

    END IF

    ! ---------------------------------------------------------------------------
    ! If the top is reached without fulfilling the boundary conditions then save
    ! the extrapolated point at which the condition is satisfied

    IF ( zeta .EQ. zN ) THEN

       IF ( extrap_z_1 .GT. zN ) extrap_z_p = extrap_z_1

       IF ( extrap_z_2 .GT. zN ) extrap_z_p = MIN( extrap_z_p , extrap_z_2 )

       extrap_z = extrap_z_p

       IF ( extrap_z_mach .GT. zN ) THEN

          extrap_z = MIN( extrap_z_p , extrap_z_mach )
          extrap_flag = -3

       END IF

       RETURN

    END IF

    ! --- Check on the extrapolated values --------------------------------------

    IF ( ( extrap_z_mach - zeta  ) * ( 1.D0 - mach ) .LT. eps_extrap ) THEN

       IF ( extrap_z_mach .LT. zN ) THEN 

          extrap_z = extrap_z_mach
          extrap_flag = 3

          RETURN

       END IF

    END IF

    coeff_p_1 = ( p1_in - r_p_1 ) / ( p1_in - p_out_1 )

    IF ( ( extrap_z_1 - zeta  ) * ( 1.D0 - coeff_p_1 ) .LT. eps_extrap ) THEN

       IF ( extrap_z_1 .LT. zN ) THEN

          extrap_z = extrap_z_1
          extrap_flag = 1

          RETURN

       END IF

    END IF

    coeff_p_2 = ( r_p_2 - p2_in ) / ( p_out_2 - p2_in )

    IF ( ( extrap_z_2 - zeta  ) * ( 1.D0 - coeff_p_2 ) .LT. eps_extrap ) THEN

       IF ( extrap_z_2 .LT. zN ) THEN

          extrap_z = extrap_z_2
          extrap_flag = 2

          RETURN

       END IF

    END IF

  END SUBROUTINE linear_extrapolation

END MODULE steady_solver
