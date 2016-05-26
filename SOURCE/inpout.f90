!********************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!> \date 15/08/2011
!********************************************************************

MODULE inpout

  ! -- Variables for the namelist TRANSIENT_PARAMETERS
  USE parameters, ONLY : verbose_level
  USE parameters, ONLY : n_cry , n_gas , n_eqns , n_vars

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry, ONLY : z0 , zN , radius_fixed, radius_min, radius_max,          &
                       radius_z, radius_z_sig, radius_model , comp_cells ,      &
                       update_radius

  ! -- Variables for the namelist STEADY_BOUNDARY_CONDITIONS
  USE init, ONLY : T_in , p1_in , delta_p_in , p_out , u1_in
  USE constitutive, ONLY : x_ex_dis_in
  USE parameters, ONLY : shooting
  USE parameters, ONLY : eps_conv

  ! -- Variables for the namelist EXSOLVED_GAS_PARAMETERS
  USE constitutive, ONLY : rho0_g , cv_g , gamma_g , T0_g ,                     &
       bar_e_g , visc_2 , lateral_degassing_flag ,                              &
       alfa2_lat_thr , perm0, Pc_g, Tc_g, a_g, b_g, s0_g, gas_law

  ! -- Variables for the namelist DISSOLVED_GAS_PARAMETERS
  USE constitutive, ONLY : rho0_d , C0_d , cv_d , gamma_d , p0_d , T0_d ,       &
       bar_e_d , bar_p_d , s0_d , exsol_model , solub , solub_exp

  ! -- Variables for the namelist CRYSTALS_PARAMETERS
  USE constitutive, ONLY : rho0_c , C0_c , cv_c , gamma_c , p0_c , T0_c ,       &
       bar_e_c , bar_p_c , s0_c , beta0 , beta_max, crystallization_model

  ! -- Variables for the namelist MELT_PARAMETERS
  USE constitutive, ONLY : rho0_m , C0_m , cv_m , gamma_m , p0_m , T0_m ,       &
       bar_e_m , bar_p_m , s0_m

  ! -- Variables for the namelist VISCOSITY_PARAMETERS
  USE constitutive, ONLY : visc_melt_model , bubbles_model , theta_model ,      &
       theta_fixed

  ! -- Variables for the namelist TEMPERATURE_PARAMETERS
  USE constitutive, ONLY : isothermal , fixed_temp

  ! -- Variables for the namelist FRAGMENTATION_PARAMETERS
  USE constitutive, ONLY : explosive , fragmentation_model , frag_thr

  ! -- Variables for the namelist EXTERNAL_WATER_PARAMETERS
  USE constitutive, ONLY : ext_water , total_water_influx , min_z_influx ,      &
       delta_z_influx , T_w , inst_vaporization

  ! -- Variables for the namelist SOURCE_PARAMETERS
  USE constitutive, ONLY : grav

  ! -- Variables for the namelist RELAXATION_PARAMETERS
  USE constitutive, ONLY : drag_funct_model , drag_funct_coeff , p_relax_model ,&
       tau_p_coeff , tau_c , tau_d

  ! -- Variables for the namelist FORCHHEIMER_PARAMETERS
  USE constitutive, ONLY : bubble_number_density , tortuosity_factor ,          &
       throat_bubble_ratio , friction_coefficient , C_D , r_a, xa, xb, xc,       &
       log10_bubble_number_density 

  !  USE parameters, ONLY : atmospheric_pressure, chocked_flow

  IMPLICIT NONE

  CHARACTER(LEN=50) :: run_name           !< Name of the run
  CHARACTER(LEN=50) :: bak_name           !< Name of the backup file for the parameters
  CHARACTER(LEN=50) :: input_file         !< Name of the file with the run parameters
  CHARACTER(LEN=50) :: output_q_file      !< Name of the output files
  CHARACTER(LEN=50) :: output_p_file      !< Name of the output files
  CHARACTER(LEN=50) :: steady_p_file      !< Name of the steady output file 
  CHARACTER(LEN=50) :: steady_q_file      !< Name of the steady output file 
  CHARACTER(LEN=50) :: exit_file          !< Name of the steady output file 

  INTEGER, PARAMETER :: input_unit = 7            !< Input data unit
  INTEGER, PARAMETER :: backup_unit = 8           !< Backup input data unit
  INTEGER, PARAMETER :: output_q_unit = 10        !< Output data unit
  INTEGER, PARAMETER :: output_p_unit = 11        !< Output data unit
  INTEGER, PARAMETER :: steady_p_output_unit = 13 !< Output unit
  INTEGER, PARAMETER :: steady_q_output_unit = 14 !< Output unit
  INTEGER, PARAMETER :: dakota_unit = 15          !< Dakota Output unit
  INTEGER, PARAMETER :: dakota_unit2 = 17          !< Dakota Output unit
  INTEGER, PARAMETER :: exit_unit = 16            !< Exit Output unit


  LOGICAL :: close_units

  ! -- Variables for the namelist RELAXATION_PARAMETERS
  REAL*8 :: log10_drag_funct_coeff , log10_tau_p_coeff

  REAL*8, ALLOCATABLE :: log10_tau_c(:)
  REAL*8, ALLOCATABLE :: log10_tau_d(:)

  NAMELIST / run_parameters / run_name , verbose_level

  NAMELIST / geometry_parameters / z0 , zN , radius_model , radius_fixed ,      &
       radius_min, radius_max, radius_z, radius_z_sig,  comp_cells

  NAMELIST / phases_parameters / n_gas , n_cry

  NAMELIST / steady_boundary_conditions / T_in , p1_in , delta_p_in ,           &
       x_ex_dis_in , p_out , u1_in , eps_conv, shooting

  NAMELIST / exsolved_gas_parameters / gas_law, Pc_g , Tc_g , cv_g , gamma_g ,  &
       rho0_g , T0_g , bar_e_g , s0_g, visc_2 , lateral_degassing_flag ,        &
       alfa2_lat_thr , perm0 

  NAMELIST / dissolved_gas_parameters / rho0_d , C0_d , cv_d , gamma_d , p0_d , &
       T0_d , bar_e_d , bar_p_d , s0_d , exsol_model , solub , solub_exp

  NAMELIST / crystals_parameters / rho0_c , C0_c , cv_c , gamma_c , p0_c ,      &
       T0_c , bar_e_c , bar_p_c , s0_c , beta0 , beta_max, crystallization_model

  NAMELIST / melt_parameters / rho0_m , C0_m , cv_m , gamma_m , p0_m , T0_m ,   &
       bar_e_m , bar_p_m , s0_m

  NAMELIST / viscosity_parameters / visc_melt_model , bubbles_model ,           &
       theta_model , theta_fixed

  NAMELIST / temperature_parameters / isothermal , fixed_temp

  NAMELIST / fragmentation_parameters / explosive , fragmentation_model ,       &
       frag_thr

  NAMELIST / external_water_parameters / ext_water , total_water_influx ,       &
       min_z_influx , delta_z_influx , T_w , inst_vaporization

  NAMELIST / source_parameters /  grav

  NAMELIST / relaxation_parameters / drag_funct_model , log10_drag_funct_coeff ,&
       p_relax_model , log10_tau_p_coeff, log10_tau_c , log10_tau_d

  NAMELIST / bubbles_parameters / bubble_number_density

  NAMELIST / forchheimer_parameters / log10_bubble_number_density ,             &
       tortuosity_factor , throat_bubble_ratio , friction_coefficient , C_D ,   &
       r_a

  NAMELIST / permeability_parameters / xa, xb, xc


CONTAINS

  !*********************************************************************
  !> \brief Initialization of the variables read from the input file
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !> \date 26/08/2011
  !*********************************************************************

  SUBROUTINE init_param

    IMPLICIT none

    LOGICAL :: lexist

    INTEGER :: n_gas_init

    INTEGER :: n_cry_init

    ! -- Variables for the namelist DISSOLVED_GAS_PARAMETERS
    REAL*8 :: rho0_d_init, C0_d_init , cv_d_init , gamma_d_init , p0_d_init , &
         T0_d_init , bar_e_d_init , bar_p_d_init , s0_d_init , solub_init ,   &
         solub_exp_init

    ! -- Variables for the namelist CRYSTALS_PARAMETERS
    REAL*8 :: rho0_c_init , C0_c_init , cv_c_init , gamma_c_init , p0_c_init ,&
         T0_c_init , bar_e_c_init , bar_p_c_init , s0_c_init , beta0_init ,   &
         beta_max_init

    REAL*8 :: log10_tau_c_init , log10_tau_d_init

    REAL*8 :: x_ex_dis_in_init, xa, xb, xc

    REAL*8 :: Pc_g_init , Tc_g_init , cv_g_init , gamma_g_init , rho0_g_init ,&
         T0_g_init , bar_e_g_init, s0_g_init

    CHARACTER*20 :: gas_law_init 

    NAMELIST / phases_parameters_init / n_gas_init , n_cry_init

    NAMELIST / steady_boundary_conditions_init / T_in , p1_in , delta_p_in ,  &
         x_ex_dis_in_init , p_out , u1_in , eps_conv, shooting

    NAMELIST / exsolved_gas_parameters_init / gas_law_init, Pc_g_init , Tc_g_init ,         &
         cv_g_init , gamma_g_init , rho0_g_init , T0_g_init , bar_e_g_init ,  &
         s0_g_init, visc_2 , lateral_degassing_flag , alfa2_lat_thr , perm0 

    NAMELIST / dissolved_gas_parameters_init / rho0_d_init ,                  &
         C0_d_init , cv_d_init , gamma_d_init , p0_d_init , T0_d_init ,       &
         bar_e_d_init , bar_p_d_init , s0_d_init , exsol_model ,              &
         solub_init , solub_exp_init

    NAMELIST / crystals_parameters_init / rho0_c_init ,                       &
         C0_c_init , cv_c_init , gamma_c_init , p0_c_init , T0_c_init ,       &
         bar_e_c_init , bar_p_c_init , s0_c_init , beta0_init , beta_max_init


    NAMELIST / relaxation_parameters_init / drag_funct_model ,                &
         log10_drag_funct_coeff , p_relax_model , log10_tau_p_coeff ,         &
         log10_tau_c_init , log10_tau_d_init


    !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
    run_name = 'test'
    verbose_level = 4 

    !-- Inizialization of the Variables for the namelist geometry_parameters
    z0 = 0.D0
    zN = 5200.D0
    radius_fixed = 0.D0
    radius_min = 0.D0
    radius_max = 0.D0
    radius_z = 0.D0
    radius_z_sig = 0.D0
    radius_model = 'fixed'
    comp_cells = 500

    !-- Inizialization of the Variables for the namelist 
    !-- steady_boundary_condition_parameters
    T_in = 1323.D0
    p1_in = 135000000.D0
    delta_p_in = 0.D0
    x_ex_dis_in_init = 5.0D-2
    p_out = 101300
    u1_in = 0.D0
    shooting = .TRUE.
    eps_conv = 1.D-5

    ! exsolved gas parameters
    gas_law_init = 'VDW'
    Pc_g_init = 22064000.D0
    Tc_g_init = 647.D0
    cv_g_init = 1571.D0
    gamma_g_init = 1.324D0
    rho0_g_init = 0.52924425D0
    T0_g_init = 373.0
    s0_g_init = 2373.0
    visc_2 = 1.5D-2
    lateral_degassing_flag = .FALSE.
    alfa2_lat_thr = 1.1D0
    perm0 = 5.0D-3

    n_gas_init = 1
    n_cry_init = 1

    n_vars = 5 + 2 * n_gas_init + n_cry_init
    n_eqns = n_vars

    ! dissolved gas parameters

    rho0_d_init = 1000.D0
    C0_d_init = 2000.D0
    cv_d_init = 1200.D0
    gamma_d_init = 2.3D0
    p0_d_init = 1.0d8
    s0_d_init = 0.d0
    exsol_model = 'Henry'
    solub_init = 4.11D-06 
    solub_exp_init = 0.5D0

    ! crystals parameters


    rho0_c_init = 2600.D0
    C0_c_init = 2000.D0
    cv_c_init = 1200.D0
    gamma_c_init = 2.3D0
    p0_c_init = 1.0d8
    s0_c_init = 0.d0
    beta0_init = 0.0D0
    beta_max_init = 0.6D0
    crystallization_model = 'Vitturi2010'

    ! MELT_PARAMETERS
    rho0_M = 2300.0000000000000  
    c0_m = 2000.0000000000000    
    cv_m = 1200.0000000000000    
    gamma_m = 2.3    
    p0_m = 100000000.00000000  
    s0_m = 0.0000000000000000  

    !-- Inizialization of the Variables for the namelist viscosity_parameters
    visc_melt_model = 'Hess_and_Dingwell1996'
    bubbles_model = 'none'
    theta_model = 'Fixed_value'
    theta_fixed = 50.d0

    !-- Inizialization of the Variables for the namelist temperature_parameters
    isothermal = .FALSE.
    fixed_temp = 1323.D0

    !-- Inizialization of the Variables for the namelist 
    !-- fragmentation_parameters
    explosive = .TRUE.
    fragmentation_model = 1
    frag_thr = 0.60D+0
    !tau_frag_coeff = 1.D0
    !tau_frag_exp = 5.0D0

    !-- Inizialization of the Variables for the namelist 
    !-- external_water_parameters 
    ext_water = .FALSE. 
    total_water_influx = 0.D0
    min_z_influx = 0.D0
    delta_z_influx = 0.D0
    T_w = 0.D0
    inst_vaporization = .FALSE.

    !-- Inizialization of the Variables for the namelist source_parameters
    grav = 9.81D0

    !-- Initialization of the variables for the namelist relaxation_parameters
    drag_funct_model = 'forchheimer'
    log10_drag_funct_coeff = 1.0D0
    p_relax_model = 'constant'
    log10_tau_p_coeff = -8.D0
    log10_tau_d_init = -4.0D0
    log10_tau_c_init = -4.0D0

    drag_funct_coeff = 10.D0 ** log10_drag_funct_coeff
    tau_p_coeff = 10.D0 ** log10_tau_p_coeff
    tau_d(1:n_gas) = 10.D0 ** log10_tau_d
    tau_c(1:n_cry) = 10.D0 ** log10_tau_c(1:n_cry)

    ! Forchheimer (Eq. 16 Degruyter et al. 2012)

    log10_bubble_number_density = 15.0D0   
    bubble_number_density = 10.0D0 ** log10_bubble_number_density
    tortuosity_factor = 3.5 
    throat_bubble_ratio = 0.1 
    friction_coefficient = 10.D0

    C_D = 0.8D0
    r_a = 1.D-3

    xa = 1.0
    xb = 1.0
    xc = 1.0

    input_file = 'conduit_solver.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF ( .NOT. lexist ) THEN

       OPEN(input_unit,FILE=input_file,STATUS='NEW')

       WRITE(input_unit, run_parameters )

       WRITE(input_unit, geometry_parameters )

       WRITE(input_unit, phases_parameters_init )

       WRITE(input_unit, steady_boundary_conditions_init )

       WRITE(input_unit, exsolved_gas_parameters_init )

       WRITE(input_unit, dissolved_gas_parameters_init )

       WRITE(input_unit, crystals_parameters_init )

       WRITE(input_unit, melt_parameters )

       WRITE(input_unit, viscosity_parameters )

       WRITE(input_unit, temperature_parameters )

       WRITE(input_unit, fragmentation_parameters )

       WRITE(input_unit, external_water_parameters )

       WRITE(input_unit, source_parameters )

       WRITE(input_unit, relaxation_parameters_init )

       IF ( drag_funct_model .EQ. 'eval' ) THEN

          WRITE(input_unit, bubbles_parameters )

       ELSEIF ( ( drag_funct_model .EQ. 'darcy' ) .OR.                          &
            ( drag_funct_model .EQ. 'darcy_Bai2011' ) .OR.                      &
            ( drag_funct_model .EQ. 'darcy_Bai2010_LB' ) .OR.                   &
            ( drag_funct_model .EQ. 'darcy_Bai2010_meas' ) .OR.                 &
            ( drag_funct_model .EQ. 'darcy_Polacci2009' ) .OR.                  &
            ( drag_funct_model .EQ. 'forchheimer_Bai2011' ) .OR.                &
            ( drag_funct_model .EQ. 'forchheimer_Bai2010_LB' ) .OR.             &
            ( drag_funct_model .EQ. 'forchheimer_Bai2010_meas' ) .OR.           &
            ( drag_funct_model .EQ. 'forchheimer_Polacci2009' ) .OR.            &
            drag_funct_model .EQ. 'forchheimer' ) THEN

          WRITE(input_unit, forchheimer_parameters )

       ELSEIF ( (drag_funct_model .EQ. 'forchheimer_mod') .OR.                  & 
            (drag_funct_model .EQ. 'forchheimer_mod2') .OR.                     &  
            (drag_funct_model .EQ. 'forchheimer_mod3') ) THEN

          WRITE(input_unit, permeability_parameters )

       END IF

       CLOSE(input_unit)

       WRITE(*,*) 'Input file not found'
       WRITE(*,*) 'A new one with default values has been created'
       STOP

    ELSE


    END IF

  END SUBROUTINE init_param

  !*********************************************************************
  !> \brief Read the input file
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !> \date 26/08/2011
  !*********************************************************************

  SUBROUTINE read_param

    USE constitutive, ONLY : T0_c , bar_p_c !, bar_e_c
    USE constitutive, ONLY : T0_m , bar_p_m !, bar_e_m

    USE init, ONLY : beta_in , xd_md_in

    USE constitutive, ONLY : allocate_phases_parameters

    IMPLICIT none

    OPEN(input_unit,FILE=input_file,STATUS='old')


    ! ------- READ run_parameters NAMELIST -----------------------------------

    READ(input_unit, run_parameters )

    READ(input_unit,geometry_parameters)

    IF ( radius_model == 'fixed' ) THEN

       IF ( radius_min .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_min not used'
       IF ( radius_max .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_max not used'
       IF ( radius_z .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_z not used'
       IF ( radius_z_sig .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_z_sig not used'

    ELSEIF ( radius_model == 'linear' ) THEN

       IF ( radius_fixed .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_fixed not used'
       IF ( radius_z .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_z not used'
       IF ( radius_z_sig .NE. 0.D0 ) WRITE(*,*) 'WARNING: radius_z_sig not used'

    END IF
              

    ! ------- READ phases_parameters NAMELIST --------------------------
    READ(input_unit, phases_parameters )

    n_vars = 5 + 2 * n_gas + n_cry
    n_eqns = n_vars

    ALLOCATE( beta_in(n_cry) )
    ALLOCATE( xd_md_in(n_gas) )

    CALL allocate_phases_parameters

    READ(input_unit,steady_boundary_conditions)

    IF ( u1_in .GT. 0.D0 ) THEN

       IF ( shooting .EQV. .TRUE.) THEN
          WRITE(*,*) 'Shooting technique to search for inlet velocity'
       ELSE	
          WRITE(*,*) 'Single run: no shooting'
       END IF

    ELSE

       shooting = .TRUE.
       WRITE(*,*) 'Shooting technique to search for inlet velocity'

    END IF


    ! ------- READ exsolved_gas_parameters NAMELIST --------------------------
    READ(input_unit, exsolved_gas_parameters )

    IF ( gas_law .EQ. 'IDEAL' )  THEN

       a_g(1:n_gas) = 0.D0
       b_g(1:n_gas) = 0.D0

    ELSEIF (  gas_law .EQ. 'VDW' ) THEN

       a_g(1:n_gas) = 27.0 / 64.0 * ( cv_g(1:n_gas)*(gamma_g(1:n_gas) - 1)      &
            * Tc_g(1:n_gas))**2.0 / Pc_g(1:n_gas)

       b_g(1:n_gas) = 1.0 / 8.0 * ( cv_g(1:n_gas)*(gamma_g(1:n_gas) - 1)        &
            * Tc_g(1:n_gas)) / Pc_g(1:n_gas)

    ELSE 

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong gas law chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'IDEAL'
       WRITE(*,*) 'VDW'
       WRITE(*,*) ''
       CALL abort

    END IF


    ! ------- READ dissolved_gas_parameters NAMELIST --------------------------
    READ(input_unit, dissolved_gas_parameters )

    T0_d(1:n_gas) = C0_d(1:n_gas) **2.D0 / ( cv_d(1:n_gas) * gamma_d(1:n_gas)   &
         * ( gamma_d(1:n_gas) - 1.D0 ) )

    bar_p_d(1:n_gas) = ( rho0_d(1:n_gas) * C0_d(1:n_gas)**2.d0 -                &
         gamma_d(1:n_gas) * p0_d(1:n_gas) ) / gamma_d(1:n_gas)


    ! ------- READ crystals_parameters NAMELIST ------------------------------
    READ(input_unit, crystals_parameters )

    T0_c(1:n_cry) = C0_c(1:n_cry) **2.D0 / ( cv_c(1:n_cry) * gamma_c(1:n_cry)   &
         * ( gamma_c(1:n_cry) - 1.D0 ) )

    bar_p_c(1:n_cry) = ( rho0_c(1:n_cry) * C0_c(1:n_cry) ** 2.d0 -              &
         gamma_c(1:n_cry) * p0_c(1:n_cry) ) / gamma_c(1:n_cry)


    IF ( (.NOT. (crystallization_model .EQ. 'Vitturi2010' ) ) .AND.             & 
         (.NOT. (crystallization_model .EQ. 'Stromboli_ST133s' ) ) .AND.        &
         (.NOT. (crystallization_model .EQ. 'Stromboli_ST130p' ) ) .AND.        &
         (.NOT. (crystallization_model .EQ. 'Stromboli_STR180307' ) ) .AND.     &
         (.NOT. (crystallization_model .EQ. 'Etna_240701D' ) ) .AND.            &
         (.NOT. (crystallization_model .EQ. 'Etna_260701C' ) ) .AND.            &
         (.NOT. (crystallization_model .EQ. 'Kilauea_07-Jul-92' ) ) ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong crystallization model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'Vitturi2010'
       WRITE(*,*) 'Stromboli_ST133s'
       WRITE(*,*) 'Stromboli_ST130p'
       WRITE(*,*) 'Stromboli_STR180307'
       WRITE(*,*) 'Etna_240701D'
       WRITE(*,*) 'Etna_260701C'
       WRITE(*,*) 'Kilauea_07-Jul-92'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Vitturi2010' ) .AND.                      & 
         (.NOT.(n_cry .EQ. 1.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Vitturi2010, n_cry has to be 1'
       WRITE(*,*) ''

       CALL ABORT

    END IF


    IF ( (crystallization_model .EQ. 'Stromboli_ST133s' ) .AND.                 & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Stromboli_ST133s, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Stromboli_ST130p' ) .AND.                 & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Stromboli_ST130p, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Stromboli_STR180307' ) .AND.              & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Stromboli_STR180307, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Etna_240701D' ) .AND.                     & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Etna_240701D, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Etna_260701C' ) .AND.                     & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Etna_260701C, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Kilauea_07-Jul-92' ) .AND.                & 
         (.NOT.(n_cry .EQ. 3.0))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Kilauea_07-Jul-92, n_cry has to be 3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    ! ------- READ melt_parameters NAMELIST -----------------------------------
    READ(input_unit, melt_parameters )

    T0_m = C0_m **2.D0 / ( cv_m * gamma_m * ( gamma_m - 1.D0 ) )
    bar_p_m = ( rho0_m * C0_m ** 2.d0 - gamma_m * p0_m ) / gamma_m

    ! ------- READ viscosity_parameters NAMELIST ------------------------------
    READ(input_unit, viscosity_parameters)


    IF ( (.NOT. (visc_melt_model .EQ. 'Hess_and_Dingwell1996' ) ) .AND.         & 
         (.NOT. (visc_melt_model .EQ. 'Giordano_et_al2008' ) ) .AND.            & 
         (.NOT. (visc_melt_model .EQ. 'Di_Genova_et_al2013_eqn_3,5' ) ) .AND.   & 
         (.NOT. (visc_melt_model .EQ. 'Di_Genova_et_al2013_eqn_4,5' ) ) .AND.   & 
         (.NOT. (visc_melt_model .EQ. 'Giordano_et_al2009' ) ) .AND.            & 
         (.NOT. (visc_melt_model .EQ. 'Romano_et_al2003' ) ) ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong melt viscosity model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'Hess_and_Dingwell1996'
       WRITE(*,*) 'Romano_et_al2003'
       WRITE(*,*) 'Giordano_et_al2008'
       WRITE(*,*) 'Giordano_et_al2009'
       WRITE(*,*) 'Di_Genova_et_al2013_eqn_3,5'
       WRITE(*,*) 'Di_Genova_et_al2013_eqn_4,5'
       WRITE(*,*) ''

       CALL ABORT

    END IF


    IF ( (.NOT. (theta_model .EQ. 'Lejeune_and_Richet1995' ) ) .AND.            & 
         (.NOT. (theta_model .EQ. 'Dingwell1993' ) ) .AND.                      & 
         (.NOT. (theta_model .EQ. 'Melnik_and_Sparks1999' ) ) .AND.             & 
         (.NOT. (theta_model .EQ. 'Costa2005' ) ) .AND.                         & 
         (.NOT. (theta_model .EQ. 'Melnik_and_Sparks2005' ) ) .AND.             & 
         (.NOT. (theta_model .EQ. 'Vona_et_al2011' ) ) .AND.                    & 
         (.NOT. (theta_model .EQ. 'Vona_et_al2011_mod' ) ) .AND.                & 
         (.NOT. (theta_model .EQ. 'Fixed_value' ) ) ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong theta model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'Lejeune_and_Richet1995'
       WRITE(*,*) 'Dingwell1993'
       WRITE(*,*) 'Melnik_and_Sparks1999'
       WRITE(*,*) 'Costa2005'
       WRITE(*,*) 'Melnik_and_Sparks2005'
       WRITE(*,*) 'Vona_et_al2011'
       WRITE(*,*) 'Vona_et_al2011_mod'
       WRITE(*,*) 'Fixed_value'
       WRITE(*,*) ''

       CALL ABORT

    END IF


    ! ------- READ temperature_parameters NAMELIST ----------------------------
    READ(input_unit, temperature_parameters)

    IF ( isothermal ) THEN

       T_in = fixed_temp

    END IF

    ! ------- READ fragmentation_parameters NAMELIST --------------------------
    READ(input_unit, fragmentation_parameters)

    IF (fragmentation_model .NE. 1 ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong fragmentation model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) '1 (Volume_Fraction)'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    ! ------- READ external_water_parameters NAMELIST --------------------------
    READ(input_unit, external_water_parameters)

    ! ------- READ source_parameters NAMELIST ---------------------------------
    READ(input_unit, source_parameters)

    ! ------- READ relaxation_parameters NAMELIST -----------------------------
    ALLOCATE( log10_tau_d(n_gas) )
    ALLOCATE( log10_tau_c(n_cry) )

    READ(input_unit, relaxation_parameters )

    drag_funct_coeff = 10.D0 ** log10_drag_funct_coeff
    tau_p_coeff = 10.D0 ** log10_tau_p_coeff
    tau_d(1:n_gas) = 10.D0 ** log10_tau_d(1:n_gas)
    tau_c(1:n_cry) = 10.D0 ** log10_tau_c(1:n_cry)

    IF ( (drag_funct_model .EQ. 'eval' ) .OR. (drag_funct_model .EQ. 'drag' ) ) &
         THEN

       READ(input_unit, bubbles_parameters )
       
    ELSEIF ( ( drag_funct_model .EQ. 'darcy' ) .OR.                             & 
         ( drag_funct_model .EQ. 'darcy_Bai2011' ) .OR.                         &
         ( drag_funct_model .EQ. 'darcy_Bai2010_LB' ) .OR.                      &
         ( drag_funct_model .EQ. 'darcy_Bai2010_meas' ) .OR.                    &
         ( drag_funct_model .EQ. 'darcy_Polacci2009' ) .OR.                     &
         ( drag_funct_model .EQ. 'forchheimer_Bai2011' ) .OR.                   &
         ( drag_funct_model .EQ. 'forchheimer_Bai2010_LB' ) .OR.                &
         ( drag_funct_model .EQ. 'forchheimer_Bai2010_meas' ) .OR.              &
         ( drag_funct_model .EQ. 'forchheimer_Polacci2009' ) .OR.               &
         ( drag_funct_model .EQ. 'forchheimer') ) THEN

       READ(input_unit, forchheimer_parameters )

       bubble_number_density = 10.0D0 ** log10_bubble_number_density

    ELSEIF ( (drag_funct_model .EQ. 'forchheimer_mod') .OR.                     &
         (drag_funct_model .EQ. 'forchheimer_mod2') .OR.                        &
         (drag_funct_model .EQ. 'forchheimer_mod3') ) THEN

       READ(input_unit, permeability_parameters )


    ELSE

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong drag_funct_model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'eval'
       WRITE(*,*) 'drag'
       WRITE(*,*) 'darcy'
       WRITE(*,*) 'darcy_Bai2011'
       WRITE(*,*) 'darcy_Bai2010_LB'
       WRITE(*,*) 'darcy_Bai2010_meas'
       WRITE(*,*) 'darcy_Polacci2009'
       WRITE(*,*) 'forchheimer'
       WRITE(*,*) 'forchheimer_Bai2011'
       WRITE(*,*) 'forchheimer_Bai2010_LB'
       WRITE(*,*) 'forchheimer_Bai2010_meas'
       WRITE(*,*) 'forchheimer_Polacci2009'
       WRITE(*,*) 'forchheimer_mod'
       WRITE(*,*) 'forchheimer_mod2'
       WRITE(*,*) 'forchheimer_mod3'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    CLOSE( input_unit )

    bak_name = TRIM(run_name)//'.bak'

    OPEN(backup_unit,file=bak_name,status='unknown')

    WRITE(backup_unit, run_parameters )

    WRITE(backup_unit, geometry_parameters )

    WRITE(backup_unit, phases_parameters )

    WRITE(backup_unit, steady_boundary_conditions )

    WRITE(backup_unit, exsolved_gas_parameters )

    WRITE(backup_unit, dissolved_gas_parameters )

    WRITE(backup_unit, crystals_parameters )

    WRITE(backup_unit, melt_parameters )

    WRITE(backup_unit, viscosity_parameters )

    WRITE(backup_unit, temperature_parameters )

    WRITE(backup_unit, fragmentation_parameters )

    WRITE(backup_unit, external_water_parameters )

    WRITE(backup_unit, source_parameters )

    WRITE(backup_unit, relaxation_parameters )

    IF ( drag_funct_model .EQ. 'eval' ) THEN

       WRITE(backup_unit, bubbles_parameters )

    ELSEIF ( (drag_funct_model .EQ. 'darcy') .OR. 				 &
         ( drag_funct_model .EQ. 'darcy_Bai2011' ) .OR.                          &
         ( drag_funct_model .EQ. 'darcy_Bai2010_LB' ) .OR.                       &
         ( drag_funct_model .EQ. 'darcy_Bai2010_meas' ) .OR.                     &
         ( drag_funct_model .EQ. 'darcy_Polacci2009' ) .OR.                      &
         ( drag_funct_model .EQ. 'forchheimer_Bai2011' ) .OR.                    &
         ( drag_funct_model .EQ. 'forchheimer_Bai2010_LB' ) .OR.                 &
         ( drag_funct_model .EQ. 'forchheimer_Bai2010_meas' ) .OR.               &
         ( drag_funct_model .EQ. 'forchheimer_Polacci2009' ) .OR.                &
         drag_funct_model .EQ. 'forchheimer' ) THEN

       WRITE(backup_unit, forchheimer_parameters )

    ELSEIF ( (drag_funct_model .EQ. 'forchheimer_mod') .OR.                      &
         (drag_funct_model .EQ. 'forchheimer_mod2') .OR.                         &
         (drag_funct_model .EQ. 'forchheimer_mod3') ) THEN

       WRITE(backup_unit, permeability_parameters )

    END IF

    CLOSE(backup_unit)

  END SUBROUTINE read_param

  !*********************************************************************
  !> \brief Write the steady solution on the output unit
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !> \param[in]   zeta    vertical coordinate
  !> \param[in]   qp      array of physical variables
  !> \date 25/09/2012
  !*********************************************************************

  SUBROUTINE output_steady(zeta,qp,radius)

    USE parameters, ONLY : n_vars , n_cry , n_gas
    USE geometry, ONLY : zeta_exit

    USE constitutive, ONLY : eval_qp2
    USE constitutive, ONLY : sound_speeds
    USE constitutive, ONLY : frag_eff

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: zeta
    REAL*8 :: qp(n_vars)
    REAL*8, INTENT(IN) :: radius

    REAL*8 :: qp2(1+n_cry+n_gas+n_gas+4)
    REAL*8 :: mass_flow_rate
    REAL*8 :: pi

    REAL*8 :: C_mix, mach, mix_velocity

    REAL*8 :: rho1, rho2(n_gas), rho_mix

    REAL*8 :: gasTotVolFrac

    INTEGER :: i,j

    CALL eval_qp2( qp, qp2 )

    IF ( zeta .EQ. z0 ) THEN

       steady_p_file = TRIM(run_name)//'_p.std'

       OPEN(steady_p_output_unit,FILE=steady_p_file,STATUS='UNKNOWN')
       !WRITE(steady_p_output_unit,*) radius_fixed
       !WRITE(steady_p_output_unit,*) radius_min
       !WRITE(steady_p_output_unit,*) radius_max
       !WRITE(steady_p_output_unit,*) radius_z
       !WRITE(steady_p_output_unit,*) radius_z_sig
       !WRITE(steady_p_output_unit,*) radius_model
       !WRITE(steady_p_output_unit,*) comp_cells

    END IF



    DO i = 1,n_vars

       IF ( ABS(qp(i)) .LT. 1D-20 ) THEN
          qp(i) = 0.D0
       END IF

    END DO

    DO i = 1,1+n_cry+n_gas+n_gas+4

       IF ( ABS(qp2(i)) .LT. 1D-20 ) THEN
          qp2(i) = 0.D0
       END IF

    END DO



    WRITE(steady_p_output_unit, *) (zeta)
    WRITE(steady_p_output_unit, 1006) (qp(i), i=1,n_vars)
    WRITE(steady_p_output_unit, 1006) (qp2(i), i=1,1+n_cry+n_gas+n_gas+4)
    WRITE(steady_p_output_unit, *) radius

1006 FORMAT(4e20.12) 

    IF ( zeta .EQ. zeta_exit ) THEN

       rho1 = qp2(1)

       gasTotVolFrac = 0.0

       DO j=1,n_gas

          rho2(j) = qp2(1+j)
          gasTotVolFrac = gasTotVolFrac + qp(j)

       END DO

       rho_mix = (1.0D0 - gasTotVolFrac) * rho1

       DO j=1,n_gas

          rho_mix = rho_mix + qp(j) * rho2(j)

       END DO

       mix_velocity = (1.0D0 - gasTotVolFrac) * rho1 / rho_mix * qp(n_gas + 3)

       DO j=1,n_gas

          mix_velocity = mix_velocity + qp(j) * rho2(j) / rho_mix * qp(n_gas + 4)

       END DO


       CLOSE(steady_p_output_unit)

       OPEN( dakota_unit,FILE='MAMMA.out',STATUS='UNKNOWN')

       pi = 4.D0 * ATAN(1.D0)

       WRITE(dakota_unit,*) 'Total Gas volume fraction', gasTotVolFrac
       WRITE(dakota_unit,*) 'Pressure 1',qp(n_gas + 1)
       WRITE(dakota_unit,*) 'Pressure 2',qp(n_gas + 2)
       WRITE(dakota_unit,*) 'Liquid/particles velocity',qp(n_gas + 3 )
       WRITE(dakota_unit,*) 'Gas velocity',qp(n_gas + 4)
       WRITE(dakota_unit,*) 'Exit Temperature',qp(n_gas + 5)
       WRITE(dakota_unit,*) 'Plagioclase volume fraction',qp(n_gas + 6)
       WRITE(dakota_unit,*) 'Pyroxene volume fraction',qp(n_gas + 7)
       WRITE(dakota_unit,*) 'Olivine volume fraction',qp(n_gas + 8)

       CALL sound_speeds(C_mix,mach) 
       CALL update_radius(zeta)

       WRITE(dakota_unit,*) 'Mach number',mach
       WRITE(dakota_unit,*) 'Exit mixture velocity', mix_velocity
       WRITE(dakota_unit,*) 'Mass flow rate', pi * radius * radius &
            * mix_velocity * rho_mix
       WRITE(dakota_unit,*) 'Fragmentation', REAL(frag_eff)

       CLOSE( dakota_unit )

       OPEN( dakota_unit2,FILE='dakota.out',STATUS='UNKNOWN')

       pi = 4.D0 * ATAN(1.D0)

       WRITE(dakota_unit2,*) gasTotVolFrac
       WRITE(dakota_unit2,*) qp(n_gas + 1)
       WRITE(dakota_unit2,*) qp(n_gas + 2)
       WRITE(dakota_unit2,*) qp(n_gas + 3)
       WRITE(dakota_unit2,*) qp(n_gas + 4)
       WRITE(dakota_unit2,*) qp(n_gas + 5)
       WRITE(dakota_unit2,*) qp(n_gas + 6)
       WRITE(dakota_unit2,*) qp(n_gas + 7)
       WRITE(dakota_unit2,*) mach
       WRITE(dakota_unit2,*) mix_velocity
       WRITE(dakota_unit2,*) pi * radius * radius * mix_velocity * rho_mix
       WRITE(dakota_unit2,*) REAL(frag_eff)

       CLOSE( dakota_unit2 )

       exit_file = TRIM(run_name)//'_exit.std'

       OPEN( exit_unit,FILE=exit_file,STATUS='UNKNOWN')

       WRITE(exit_unit,*) 'Mass_flow_rate',mass_flow_rate
       WRITE(exit_unit,*) 'Gas volume fraction',1.D0-qp(1)
       WRITE(exit_unit,*) 'Pressure 1',qp(2)
       WRITE(exit_unit,*) 'Pressure 2',qp(3)
       WRITE(exit_unit,*) 'Crystals volume fraction',qp(7)
       WRITE(exit_unit,*) 'Liquid/particles velocity',qp(4)
       WRITE(exit_unit,*) 'Gas velocity',qp(5)
       WRITE(exit_unit,*) 'Mach number',mach

       CLOSE( exit_unit )


    END IF


  END SUBROUTINE output_steady



  !----------------------------------------------------------------------
  !> \brief Numeric to String conversion
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !> \date 27/20/2009
  !> \param   k      integer to convert             (\b input)           
  !----------------------------------------------------------------------

  CHARACTER*4 FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera


END MODULE inpout

