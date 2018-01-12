!********************************************************************
!> \brief Input/Output module
!> @author 
!> Mattia de' Michieli Vitturi
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!> \date 15/08/2011
!********************************************************************

MODULE inpout

  ! -- Variables for the namelist RUN_PARAMETERS
  USE constitutive, ONLY : explosive_flag
  USE equations, ONLY : ext_water_flag , method_of_moments_flag
  USE equations, ONLY : lateral_degassing_flag

  ! -- Variables for the namelist METHOD_OF_MOMENTS_PARAMETERS
  USE parameters, ONLY : n_mom
  USE constitutive, ONLY : T_u, U_m 
  USE constitutive, ONLY : cry_shape_factor , L0_cry
  
  ! -- Variables for the namelist TRANSIENT_PARAMETERS
  USE parameters, ONLY : verbose_level
  USE parameters, ONLY : n_cry , n_gas , n_eqns , n_vars 

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry, ONLY : z0 , zN , radius_fixed, radius_min, radius_max,          &
                       radius_z, radius_z_sig, radius_model , eccen_fixed,	&
		       eccen_base, eccen_top, eccen_z_base, eccen_z_top,	&
		       eccen_axis_b, comp_cells, update_radius

  ! -- Variables for the namelist STEADY_BOUNDARY_CONDITIONS
  USE init, ONLY : T_in , p1_in , delta_p_in , p_out , u1_in
  USE constitutive, ONLY : x_ex_dis_in
  USE parameters, ONLY : shooting
  USE parameters, ONLY : eps_conv

  ! -- Variables for the namelist EXSOLVED_GAS_PARAMETERS
  USE constitutive, ONLY : rho0_g , cv_g , gamma_g , T0_g ,                     &
       bar_e_g , visc_2 , perm0, Pc_g, &
       Tc_g, a_g, b_g, s0_g, gas_law

  USE equations, ONLY : alfa2_lat_thr

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
  USE equations, ONLY : isothermal , fixed_temp

  ! -- Variables for the namelist FRAGMENTATION_PARAMETERS
  USE constitutive, ONLY : fragmentation_model , frag_thr

  ! -- Variables for the namelist EXTERNAL_WATER_PARAMETERS
  USE equations, ONLY : total_water_influx , min_z_influx , delta_z_influx ,    &
       T_w , inst_vaporization , aquifer_type

  ! -- Variables for the namelist SOURCE_PARAMETERS
  USE constitutive, ONLY : grav

  ! -- Variables for the namelist COUNTRY_ROCK_PARAMETERS
  USE constitutive, ONLY : rho_cr, k_cr, log10_k_cr

  ! -- Variables for the namelist RELAXATION_PARAMETERS
  USE constitutive, ONLY : drag_funct_model , drag_funct_coeff , p_relax_model ,&
       tau_p_coeff , tau_c , tau_d

  ! -- Variables for the namelist FORCHHEIMER_PARAMETERS
  USE constitutive, ONLY : bubble_number_density , tortuosity_factor ,          &
       throat_bubble_ratio , friction_coefficient , C_D , r_a, xa, xb, xc,      &
       log10_bubble_number_density 

  !  USE parameters, ONLY : atmospheric_pressure, chocked_flow

  ! -- Variables for the card meltcomposition
  USE constitutive, ONLY : wt_init

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

  NAMELIST / run_parameters / run_name , verbose_level , ext_water_flag ,       &
       lateral_degassing_flag , explosive_flag , method_of_moments_flag

  NAMELIST / geometry_parameters / z0 , zN , radius_model , radius_fixed ,      &
       radius_min, radius_max, radius_z, radius_z_sig, eccen_fixed,		&
       eccen_base, eccen_top, eccen_z_base, eccen_z_top, eccen_axis_b, comp_cells		

  NAMELIST / phases_parameters / n_gas , n_cry

  NAMELIST / steady_boundary_conditions / T_in , p1_in , delta_p_in ,           &
       x_ex_dis_in , p_out , u1_in , eps_conv, shooting

  NAMELIST / exsolved_gas_parameters / gas_law, Pc_g , Tc_g , cv_g , gamma_g ,  &
       rho0_g , T0_g , bar_e_g , s0_g, visc_2 , alfa2_lat_thr , perm0
       
  NAMELIST / dissolved_gas_parameters / rho0_d , C0_d , cv_d , gamma_d , p0_d , &
       T0_d , bar_e_d , bar_p_d , s0_d , exsol_model , solub , solub_exp

  NAMELIST / crystals_parameters / rho0_c , C0_c , cv_c , gamma_c , p0_c ,      &
       T0_c , bar_e_c , bar_p_c , s0_c , beta0 , beta_max, crystallization_model

  NAMELIST / melt_parameters / rho0_m , C0_m , cv_m , gamma_m , p0_m , T0_m ,   &
       bar_e_m , bar_p_m , s0_m

  NAMELIST / viscosity_parameters / visc_melt_model , bubbles_model ,           &
       theta_model , theta_fixed

  NAMELIST / temperature_parameters / isothermal , fixed_temp

  NAMELIST / fragmentation_parameters / fragmentation_model , frag_thr

  NAMELIST / external_water_parameters / total_water_influx ,                   &
       min_z_influx , delta_z_influx , T_w , inst_vaporization , aquifer_type

  NAMELIST / source_parameters /  grav

  NAMELIST / country_rock_parameters / rho_cr , log10_k_cr
 
  NAMELIST / relaxation_parameters / drag_funct_model , log10_drag_funct_coeff ,&
       p_relax_model , log10_tau_p_coeff, log10_tau_c , log10_tau_d

  NAMELIST / bubbles_parameters / bubble_number_density

  NAMELIST / forchheimer_parameters / log10_bubble_number_density ,             &
       tortuosity_factor , throat_bubble_ratio , friction_coefficient , C_D ,   &
       r_a

  NAMELIST / permeability_parameters / xa, xb, xc

  NAMELIST / method_of_moments_parameters / n_mom , T_u , U_m , L0_cry ,        &
       cry_shape_factor

CONTAINS

  !******************************************************************************
  !> \brief Initialization of the variables read from the input file
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !> \date 26/08/2011
  !******************************************************************************

  SUBROUTINE init_param

    ! external subroutines
    USE constitutive, ONLY : allocate_phases_parameters
    
    IMPLICIT none

    LOGICAL :: lexist

    INTEGER :: n_gas_init

    INTEGER :: n_cry_init

    !  Variables for the namelist DISSOLVED_GAS_PARAMETERS
    REAL*8 :: rho0_d_init, C0_d_init , cv_d_init , gamma_d_init , p0_d_init , &
         T0_d_init , bar_e_d_init , bar_p_d_init , s0_d_init , solub_init ,   &
         solub_exp_init

    !  Variables for the namelist CRYSTALS_PARAMETERS
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

    NAMELIST / exsolved_gas_parameters_init / gas_law_init, Pc_g_init ,       &
         Tc_g_init , cv_g_init , gamma_g_init , rho0_g_init , T0_g_init ,     &
         bar_e_g_init , s0_g_init, visc_2 , lateral_degassing_flag ,          &
         alfa2_lat_thr , perm0

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

    ! Inizialization of the Variables for the namelist RUN_PARAMETERS
    run_name = 'test'
    ext_water_flag = .FALSE.
    lateral_degassing_flag = .FALSE.
    explosive_flag = .TRUE.
    method_of_moments_flag = .FALSE.
    verbose_level = 4 

    ! Inizialization of the Variables for the namelist geometry_parameters
    z0 = 0.D0
    zN = 5200.D0
    radius_fixed = 0.D0
    radius_min = 0.D0
    radius_max = 0.D0
    radius_z = 0.D0
    radius_z_sig = 0.D0
    radius_model = 'fixed'
    eccen_fixed = 0.D0
    eccen_base = 0.D0
    eccen_top = 0.D0
    eccen_z_base = 0.D0
    eccen_z_top = 0.D0
    eccen_axis_b = 0.D0
    comp_cells = 500

    ! Inizialization of the Variables for the namelist 
    ! steady_boundary_condition_parameters
    T_in = 1323.D0
    p1_in = 135000000.D0
    delta_p_in = 0.D0
    x_ex_dis_in_init = 5.0D2
    p_out = 101300
    u1_in = 0.D0
    shooting = .TRUE.
    eps_conv = 1.D5

    
    ! exsolved gas parameters
    gas_law_init = 'VDW'
    Pc_g_init = 22064000.D0
    Tc_g_init = 647.D0
    cv_g_init = 1571.D0
    gamma_g_init = 1.324D0
    rho0_g_init = 0.52924425D0
    T0_g_init = 373.0
    s0_g_init = 2373.0
    visc_2 = 1.5D2
    alfa2_lat_thr = 1.1D0
    perm0 = 5.0D3

    n_gas_init = 1
    n_cry_init = 1

    ! dissolved gas parameters
    rho0_d_init = 1000.D0
    C0_d_init = 2000.D0
    cv_d_init = 1200.D0
    gamma_d_init = 2.3D0
    p0_d_init = 1.0d8
    s0_d_init = 0.d0
    exsol_model = 'Henry'
    solub_init = 4.11D06 
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

    ! Inizialization of the Variables for the namelist viscosity_parameters
    visc_melt_model = 'Hess_and_Dingwell1996'
    bubbles_model = 'none'
    theta_model = 'Fixed_value'
    theta_fixed = 50.d0

    ! Inizialization of the Variables for the namelist temperature_parameters
    isothermal = .FALSE.
    fixed_temp = 1323.D0

    ! Inizialization of the Variables for the namelist 
    ! fragmentation_parameters
    fragmentation_model = 1
    frag_thr = 0.60D+0
    !tau_frag_coeff = 1.D0
    !tau_frag_exp = 5.0D0

    ! Inizialization of the Variables for the namelist 
    ! external_water_parameters 
    total_water_influx = 0.D0
    min_z_influx = 0.D0
    delta_z_influx = 0.D0
    T_w = 0.D0
    inst_vaporization = .FALSE.
    aquifer_type = "unconfined"
    
    ! Inizialization of the Variables for the namelist source_parameters
    grav = 9.81D0

    ! Inizialization of the Variables for the namelist country_rock_parameters
    rho_cr = 2600.D0
    log10_k_cr = -12.D0
    k_cr = 10.D0 ** log10_k_cr

    ! Initialization of the variables for the namelist relaxation_parameters
    drag_funct_model = 'forchheimer'
    log10_drag_funct_coeff = 1.0D0
    p_relax_model = 'constant'
    log10_tau_p_coeff = 8.D0
    log10_tau_d_init = 4.0D0
    log10_tau_c_init = 4.0D0

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
    r_a = 1.D3

    xa = 1.0
    xb = 1.0
    xc = 1.0

    input_file = 'conduit_solver.inp'

    INQUIRE (FILE=input_file,exist=lexist)
    
    IF ( .NOT. lexist ) THEN
       
       !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
       RUN_NAME = "MSH_Degruyter_2012"
       ext_water_flag = .FALSE.
       lateral_degassing_flag = .FALSE.
       explosive_flag = .TRUE. 
       method_of_moments_flag = .FALSE.
       VERBOSE_LEVEL = 0

       !-- Inizialization of the Variables for the namelist GEOMETRY_PARAMETERS
       Z0 = 0.D0
       ZN = 5291.D0
       RADIUS_MODEL = "fixed"
       RADIUS_FIXED = 30.D0
       COMP_CELLS = 1000

       !-- Initialization of the variables for the namelist PHASES_PARAMETERS
       N_GAS = 1
       N_CRY = 1

       CALL allocate_phases_parameters

       !-- Inizialization of the Variables for the namelist 
       !-- STEADY_BOUNDARY_CONDITIONS
       T_in = 1159.0D0
       p1_in = 140000000.0D0
       delta_p_in = 0.D0
       x_ex_dis_in = 4.6D-2
       p_out = 101300
       u1_in = 1.9D0
       shooting = .TRUE.
       eps_conv = 1.D-5

       !-- Inizialization of the Variables for the namelist 
       !-- EXSOLVED_GAS_CONDITIONS
       gas_law = 'IDEAL'
       Pc_g = 22064000.D0
       Tc_g = 647.D0
       cv_g = 1571.D0
       gamma_g = 1.29D0
       rho0_g = 0.588460051D0
       T0_g = 373.D0
       s0_g = 0.D0
       visc_2 = 1.5D-2
       alfa2_lat_thr = 1.1D0
       perm0 = 5.0D-3

       !-- Inizialization of the Variables for the namelist 
       !-- DISSOLVED_GAS_PARAMETERS
       rho0_d = 1000.D0
       C0_d = 407.02249D0
       cv_d = 3637.5787820D0
       gamma_d = 1.11D0
       p0_d = 1.0d8
       T0_D = 372.99994092673325
       BAR_E_D = 0.0000000000000000
       BAR_P_D = 49249833.789414391 
       s0_d = 0.d0
       exsol_model = 'Henry'
       solub = 4.11D-06 
       solub_exp = 0.5D0

       ! crystals parameters
       RHO0_C = 2800.0000000000000     
       C0_C = 2000.0000000000000     
       CV_C = 360.00000000000000     
       GAMMA_C = 3.3999999999999999     
       P0_C = 250000000.00000000     
       T0_C = 1361.6557734204794     
       BAR_E_C = 0.0000000000000000     
       BAR_P_C = 3044117647.0588236     
       S0_C = 0.0000000000000000     
       BETA0 = 0.40000000000000002     
       BETA_MAX = 0.40000000000000002     
       CRYSTALLIZATION_MODEL="Vitturi2010"

       ! MELT_PARAMETERS
       RHO0_M = 2500.0000000000000     
       C0_M = 1366.2740410693600     
       CV_M = 707.00000000000000     
       GAMMA_M = 2.0899999999999999     
       P0_M = 140000000.00000000     
       T0_M = 1158.9999999999998     
       BAR_E_M = 0.0000000000000000     
       BAR_P_M = 2092900424.9999993     
       S0_M = 0.0000000000000000

       !-- Inizialization of the Variables for the namelist viscosity_parameters
       visc_melt_model = 'Hess_and_Dingwell1996'
       bubbles_model = 'none'
       theta_model = 'Costa2005'
       theta_fixed = 1.d0

       !-- Inizialization of the Variables for the namelist temperature_parameters
       isothermal = .TRUE.
       fixed_temp = 1159.D0

       !-- Inizialization of the Variables for the namelist 
       !-- fragmentation_parameters
       fragmentation_model = 1
       frag_thr = 0.80D+0

       !-- Inizialization of the Variables for the namelist source_parameters
       grav = 9.81D0

       !-- Inizialization of the Variables for the namelist country_rock_parameters
       rho_cr =  2600.000000000000                                                                  
       log10_k_cr = -12.000000000000  
       
       !-- Initialization of the variables for the namelist relaxation_parameters

       ! ------- READ relaxation_parameters NAMELIST -------------------------------
       ALLOCATE( log10_tau_d(n_gas) )
       ALLOCATE( log10_tau_c(n_cry) )


       drag_funct_model = 'forchheimer'
       log10_drag_funct_coeff = 0.0D0
       p_relax_model = 'single'
       log10_tau_p_coeff = 1.D0
       log10_tau_d = -8.0D0
       log10_tau_c = -8.0D0


       !-- Forchheimer (Eq. 16 Degruyter et al. 2012)
       log10_bubble_number_density = 15.0D0   
       tortuosity_factor = 3.5 
       throat_bubble_ratio = 0.1 
       friction_coefficient = 10.D0
       C_D = 0.8D0
       r_a = 1.D-3


       OPEN(input_unit,FILE=input_file,STATUS='NEW')

       WRITE(input_unit, run_parameters )

       WRITE(input_unit, geometry_parameters )

       WRITE(input_unit, phases_parameters )

       WRITE(input_unit, steady_boundary_conditions )

       WRITE(input_unit, exsolved_gas_parameters )

       WRITE(input_unit, dissolved_gas_parameters )

       WRITE(input_unit, crystals_parameters )

       WRITE(input_unit, melt_parameters )

       WRITE(input_unit, viscosity_parameters )

       WRITE(input_unit, temperature_parameters )

       IF ( explosive_flag ) THEN
       
          WRITE(input_unit, fragmentation_parameters )

       END IF
          
       WRITE(input_unit, source_parameters )

       IF ( ext_water_flag .OR. lateral_degassing_flag ) THEN
       
          WRITE(input_unit, country_rock_parameters )

       END IF
          
       WRITE(input_unit, relaxation_parameters )

       WRITE(input_unit, forchheimer_parameters )

       CLOSE(input_unit)

       WRITE(*,*) 'Input file not found'
       WRITE(*,*) 'A new one with default values has been created'
       WRITE(*,*) 'with condition to reproduce the solution for'
       WRITE(*,*) 'MSH1980 (see Degruyter et al. 2012)'
       STOP
       
    END IF

  END SUBROUTINE init_param

  !******************************************************************************
  !> \brief Read the input file
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !> \date 26/08/2011
  !******************************************************************************

  SUBROUTINE read_param

    ! external subroutines
    USE constitutive, ONLY : allocate_phases_parameters

    ! external variables
    USE constitutive, ONLY : n_drag_models , available_drag_models
    USE constitutive, ONLY : n_theta_models , available_theta_models
    USE constitutive, ONLY : n_bubble_models, available_bubble_models
    USE constitutive, ONLY : n_visc_melt_models , available_visc_melt_models
    USE constitutive, ONLY : T0_c , bar_p_c !, bar_e_c
    USE constitutive, ONLY : T0_m , bar_p_m !, bar_e_m

    USE constitutive, ONLY : T_m , mom_cry, growth_mom 
    
    USE init, ONLY : beta_in , xd_md_in
    
    IMPLICIT none

    INTEGER :: i
    LOGICAL :: tend1
    CHARACTER(LEN=80) :: card

    LOGICAL :: check_model
    
    INTEGER :: ios

    OPEN(input_unit,FILE=input_file,STATUS='old')

    ! ------- READ run_parameters NAMELIST --------------------------------------

    READ(input_unit, run_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF
    
    READ(input_unit,geometry_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist GEOMETRY_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF
    
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

    ! ------- READ phases_parameters NAMELIST -----------------------------------
    READ(input_unit, phases_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist PHASES_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    ALLOCATE( beta_in(n_cry) )
    ALLOCATE( xd_md_in(n_gas) )

    CALL allocate_phases_parameters
   
    ! ------- READ method_of_moments_parameters NAMELIST ------------------------
    IF ( method_of_moments_flag ) THEN

       ALLOCATE( T_m(n_cry) , T_u(n_cry) , U_m(n_cry) )
       
       READ(input_unit, method_of_moments_parameters , IOSTAT = ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist METHOD_OF_MOMENTS_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE

          ALLOCATE( mom_cry(1:n_cry,0:n_mom-1) )
          ALLOCATE( growth_mom(1:n_cry,0:n_mom-1) )
          WRITE(*,*) 'Solving for ',n_mom,' moments for each crystal phase'          
          REWIND(input_unit)
          
       END IF


    ELSE

       WRITE(*,*) 'Solving for crystal volume fraction only'
       
    END IF
    

    IF ( method_of_moments_flag ) THEN
    
       n_vars = 5 + 2 * n_gas + n_cry * n_mom

    ELSE

       n_vars = 5 + 2 * n_gas + n_cry

    END IF
       
    n_eqns = n_vars


    READ(input_unit,steady_boundary_conditions , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist STEADY_BOUNDARY_CONDITIONS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
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


    ! ------- READ exsolved_gas_parameters NAMELIST ----------------------------
    READ(input_unit, exsolved_gas_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist EXSOLVED_GAS_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
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

    ! ------- READ dissolved_gas_parameters NAMELIST ---------------------------
    READ(input_unit, dissolved_gas_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist DISSOLVED_GAS_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    T0_d(1:n_gas) = C0_d(1:n_gas) **2.D0 / ( cv_d(1:n_gas) * gamma_d(1:n_gas)   &
         * ( gamma_d(1:n_gas) - 1.D0 ) )

    bar_p_d(1:n_gas) = ( rho0_d(1:n_gas) * C0_d(1:n_gas)**2.d0 -                &
         gamma_d(1:n_gas) * p0_d(1:n_gas) ) / gamma_d(1:n_gas)


    ! ------- READ crystals_parameters NAMELIST ---------------------------------
    READ(input_unit, crystals_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist CRYSTALS_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    T0_c(1:n_cry) = C0_c(1:n_cry) **2.D0 / ( cv_c(1:n_cry) * gamma_c(1:n_cry)   &
         * ( gamma_c(1:n_cry) - 1.D0 ) )

    bar_p_c(1:n_cry) = ( rho0_c(1:n_cry) * C0_c(1:n_cry) ** 2.d0 -              &
         gamma_c(1:n_cry) * p0_c(1:n_cry) ) / gamma_c(1:n_cry)


    IF (.NOT. (crystallization_model .EQ. 'Vitturi2010' ) ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong crystallization model chosen.'
       WRITE(*,*) 'Please choose between:'
       WRITE(*,*) ''
       WRITE(*,*) 'Vitturi2010'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    IF ( (crystallization_model .EQ. 'Vitturi2010' ) .AND.                      & 
         (.NOT.(n_cry .EQ. 1 ))  ) THEN

       WRITE(*,*) ''
       WRITE(*,*) 'Wrong number of crystal components inserted'
       WRITE(*,*) 'Using Vitturi2010, n_cry has to be 1'
       WRITE(*,*) ''

       CALL ABORT

    END IF

    ! ------- READ melt_parameters NAMELIST -------------------------------------
    READ(input_unit, melt_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist MELT_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    T0_m = C0_m **2.D0 / ( cv_m * gamma_m * ( gamma_m - 1.D0 ) )
    bar_p_m = ( rho0_m * C0_m ** 2.d0 - gamma_m * p0_m ) / gamma_m

    ! ------- READ viscosity_parameters NAMELIST --------------------------------
    READ(input_unit, viscosity_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist VISCOSITY_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    check_model = .FALSE.

    DO i=1,n_visc_melt_models

       IF ( TRIM(visc_melt_model) .EQ. TRIM(available_visc_melt_models(i)) ) THEN
          
          check_model = .TRUE.

       END IF

    END DO

    IF ( .NOT.check_model ) THEN

       WRITE(*,*) 'Wrong drag_funct_model chosen.'
       WRITE(*,*) 'Please choose between:'
       
       DO i=1,n_visc_melt_models
          
          WRITE(*,*) available_visc_melt_models(i)

       END DO

       STOP

    END IF


    IF ( visc_melt_model .EQ. 'Giordano_et_al2008' ) THEN
       
       tend1 = .FALSE.
       
       WRITE(*,*) 'search melt composition'
       
       melt_comp_search: DO
          
          READ(input_unit,*, END = 300 ) card
          
          IF( TRIM(card) == 'MELTCOMPOSITION' ) THEN
             
             EXIT melt_comp_search

          END IF
          
       END DO melt_comp_search
       
       READ(input_unit,*) 
       READ(input_unit,*) wt_init(1:12)
       
       IF ( verbose_level .GE. 1 ) WRITE(*,*) wt_init(1:12)

       GOTO 310
300    tend1 = .TRUE.
310    CONTINUE

       REWIND(input_unit)

    END IF

    check_model = .FALSE.

    DO i=1,n_theta_models

       IF ( TRIM(theta_model) .EQ. TRIM(available_theta_models(i)) ) THEN
          
          check_model = .TRUE.

       END IF

    END DO

    IF ( .NOT.check_model ) THEN

       WRITE(*,*) 'Wrong theta_model chosen.'
       WRITE(*,*) 'Please choose between:'
       
       DO i=1,n_theta_models
          
          WRITE(*,*) available_theta_models(i)

       END DO

       STOP

    END IF

    IF ( (theta_model .EQ. 'Vona_et_al2013_eq19' ) .OR.                         & 
         (theta_model .EQ. 'Vona_et_al2013_eq20' ) .OR.                         & 
         (theta_model .EQ. 'Vona_et_al2013_eq21' ) ) THEN

	WRITE(*,*) 'WARNING: bubbles_model not used.' 
        WRITE(*,*) 'Effect of bubbles is included by using Vona 2013 equations'	 

    END IF

    check_model = .FALSE.

    DO i=1,n_bubble_models

       IF ( TRIM(bubbles_model) .EQ. TRIM(available_bubble_models(i)) ) THEN
          
          check_model = .TRUE.

       END IF

    END DO

    IF ( .NOT.check_model ) THEN

       WRITE(*,*) 'Wrong bubbles_model chosen.'
       WRITE(*,*) 'Please choose between:'
       
       DO i=1,n_bubble_models
          
          WRITE(*,*) available_bubble_models(i)

       END DO

       STOP

    END IF


    ! ------- READ temperature_parameters NAMELIST ------------------------------
    READ(input_unit, temperature_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist TEMPERATURE_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    IF ( isothermal ) THEN

       T_in = fixed_temp

    END IF

    IF ( explosive_flag ) THEN
    
       ! ------- READ fragmentation_parameters NAMELIST ----------------------------
       READ(input_unit, fragmentation_parameters , IOSTAT = ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist FRAGMENTATION_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       
       IF (fragmentation_model .NE. 1 ) THEN
          
          WRITE(*,*) ''
          WRITE(*,*) 'Wrong fragmentation model chosen.'
          WRITE(*,*) 'Please choose between:'
          WRITE(*,*) ''
          WRITE(*,*) '1 (Volume_Fraction)'
          WRITE(*,*) ''
          
          CALL ABORT
          
       END IF

    END IF
       
    ! ------- READ external_water_parameters NAMELIST ---------------------------
    IF ( ext_water_flag ) THEN
       
       READ(input_unit, external_water_parameters , IOSTAT = ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist EXTERNAL_WATER_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       IF (ext_water_flag) THEN
          
          IF ( (.NOT. (aquifer_type .EQ. 'confined' ) ) .AND.                & 
               (.NOT. (aquifer_type .EQ. 'unconfined' ) ) ) THEN
             
             WRITE(*,*) ''
             WRITE(*,*) 'Wrong aquifer type chosen.'
             WRITE(*,*) 'Please choose between:'
             WRITE(*,*) ''
             WRITE(*,*) 'confined'
             WRITE(*,*) 'unconfined'
             WRITE(*,*) ''
             
             CALL ABORT
             
          END IF
          
       END IF

    END IF
       
    ! ------- READ source_parameters NAMELIST -----------------------------------
    READ(input_unit, source_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist SOURCE_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    IF ( ext_water_flag .OR. lateral_degassing_flag ) THEN

       ! ------- READ country_rock_parameters NAMELIST -----------------------------
       READ(input_unit, country_rock_parameters , IOSTAT = ios )
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist COUNTRY_ROCK_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       k_cr = 10 ** log10_k_cr

    END IF
       
    ! ------- READ relaxation_parameters NAMELIST -------------------------------
    ALLOCATE( log10_tau_d(n_gas) )
    ALLOCATE( log10_tau_c(n_cry) )

    READ(input_unit, relaxation_parameters , IOSTAT = ios )

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RELAXATION_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    
    drag_funct_coeff = 10.D0 ** log10_drag_funct_coeff

    IF ( p_relax_model .EQ. 'single' ) THEN
       
       log10_tau_p_coeff = 0.D0

    END IF

    tau_p_coeff = 10.D0 ** log10_tau_p_coeff


    tau_d(1:n_gas) = 10.D0 ** log10_tau_d(1:n_gas)
    tau_c(1:n_cry) = 10.D0 ** log10_tau_c(1:n_cry)


    check_model = .FALSE.

    DO i=1,n_drag_models

       IF ( TRIM(drag_funct_model) .EQ. TRIM(available_drag_models(i)) ) THEN
          
          check_model = .TRUE.

       END IF

    END DO

    IF ( .NOT.check_model ) THEN

       WRITE(*,*) 'Wrong drag_funct_model chosen.'
       WRITE(*,*) 'Please choose between:'
       
       DO i=1,n_drag_models
          
          WRITE(*,*) available_drag_models(i)

       END DO

       STOP

    END IF
       
    IF ( ( drag_funct_model .EQ. 'eval' ) .OR.                                  &
         ( drag_funct_model .EQ. 'Klug_and_Cashman' ) .OR.                      &
         ( drag_funct_model .EQ. 'drag' ) ) THEN
       
       READ(input_unit, bubbles_parameters , IOSTAT = ios )

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist BUBBLES_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
              
    ELSEIF ( ( drag_funct_model .EQ. 'darcy' ) .OR.                             & 
         ( drag_funct_model .EQ. 'forchheimer') .OR.                            & 
         ( drag_funct_model .EQ. 'forchheimer_wt')) THEN
       
       READ(input_unit, forchheimer_parameters , IOSTAT = ios )

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist FORCHHEIMER_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
       
          REWIND(input_unit)
          
       END IF
              
       bubble_number_density = 10.0D0 ** log10_bubble_number_density
       
    ELSEIF ( ( drag_funct_model .EQ. 'forchheimer_mod' ) .OR.                   &
         ( drag_funct_model .EQ. 'forchheimer_mod2' ) .OR.                      &
         ( drag_funct_model .EQ. 'forchheimer_mod3' ) ) THEN

       READ(input_unit, permeability_parameters , IOSTAT = ios )

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist PERMEABILITY_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF       
       
    END IF

    CLOSE( input_unit )

    bak_name = TRIM(run_name)//'.bak'

    OPEN(backup_unit,file=bak_name,status='unknown')

    WRITE(backup_unit, run_parameters )

    WRITE(backup_unit, geometry_parameters )

    WRITE(backup_unit, phases_parameters )

    IF ( method_of_moments_flag ) THEN

       WRITE(backup_unit, method_of_moments_parameters)

    END IF
    
    WRITE(backup_unit, steady_boundary_conditions )

    WRITE(backup_unit, exsolved_gas_parameters )

    WRITE(backup_unit, dissolved_gas_parameters )

    WRITE(backup_unit, crystals_parameters )

    WRITE(backup_unit, melt_parameters )

    WRITE(backup_unit, viscosity_parameters )

    WRITE(backup_unit, temperature_parameters )

    IF ( explosive_flag ) THEN
    
       WRITE(backup_unit, fragmentation_parameters )

    END IF
       
    IF ( ext_water_flag ) THEN
    
       WRITE(backup_unit, external_water_parameters )

    END IF

    WRITE(backup_unit, source_parameters )
       
    WRITE(backup_unit, country_rock_parameters )
       
    WRITE(backup_unit, relaxation_parameters )


    IF ( ( drag_funct_model .EQ. 'eval' ) .OR.                                  &
         ( drag_funct_model .EQ. 'Klug_and_Cashman' ) .OR.                      &
         ( drag_funct_model .EQ. 'drag' ) ) THEN
       
       WRITE(backup_unit, bubbles_parameters )
              
    ELSEIF ( ( drag_funct_model .EQ. 'darcy' ) .OR.                             & 
         ( drag_funct_model .EQ. 'forchheimer') .OR.                            & 
         ( drag_funct_model .EQ. 'forchheimer_wt')) THEN
       
       WRITE(backup_unit, forchheimer_parameters )
       
    ELSEIF ( ( drag_funct_model .EQ. 'forchheimer_mod' ) .OR.                   &
         ( drag_funct_model .EQ. 'forchheimer_mod2' ) .OR.                      &
         ( drag_funct_model .EQ. 'forchheimer_mod3' ) ) THEN

       WRITE(backup_unit, permeability_parameters  )
       
    END IF

    
    IF ( visc_melt_model .EQ. 'Giordano_et_al2008' ) THEN
       
       WRITE(backup_unit,*) '''MELTCOMPOSITION'''
       WRITE(backup_unit,*) 'SiO2  TiO2 Al2O3  FeO   MnO   MgO   CaO   Na2O  K2O   P2O5  H2O   F2O-1'     
       WRITE(backup_unit,107) wt_init(1:12)
       
107    FORMAT(12(f5.2,1x))
       
    END IF
    

    CLOSE(backup_unit)

  END SUBROUTINE read_param

  !******************************************************************************
  !> \brief Write the steady solution on the output unit
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !> \param[in]   zeta    vertical coordinate
  !> \param[in]   qp      array of physical variables
  !> \date 25/09/2012
  !******************************************************************************

  SUBROUTINE output_steady(zeta,qp,radius)

    ! External variables
    USE parameters, ONLY : n_vars , n_cry , n_gas

    USE parameters, ONLY : idx_p1 , idx_p2 , idx_u1 , idx_u2 , idx_T ,          &
         idx_alfa_first , idx_alfa_last , idx_beta_first , idx_beta_last

    USE geometry, ONLY : zeta_exit

    ! External subroutine
    USE constitutive, ONLY : sound_speeds
    USE constitutive, ONLY : frag_eff
    USE equations, ONLY : eval_qp2

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

       OPEN(steady_q_output_unit,FILE='test.out',STATUS='UNKNOWN')
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

    WRITE(steady_q_output_unit,1007,advance="no") zeta

    DO i = 1,n_vars

       WRITE(steady_q_output_unit,1007,advance="no") qp(i)

    END DO

    DO i = 1,1+n_cry+n_gas+n_gas+4

       WRITE(steady_q_output_unit,1007,advance="no") qp2(i)

    END DO
       
    WRITE(steady_q_output_unit,1007) radius

1007 FORMAT(1x,e15.8)

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

       mix_velocity = (1.0D0 - gasTotVolFrac) * rho1 / rho_mix * qp(idx_u1)

       DO j=1,n_gas

          mix_velocity = mix_velocity + qp(j) * rho2(j) / rho_mix * qp(idx_u2)

       END DO

       CLOSE(steady_p_output_unit)
       CLOSE(steady_q_output_unit)

       OPEN( dakota_unit,FILE='MAMMA.out',STATUS='UNKNOWN')

       pi = 4.D0 * ATAN(1.D0)

       WRITE(dakota_unit,*) 'Total Gas volume fraction', gasTotVolFrac
       WRITE(dakota_unit,*) 'Pressure 1',qp(idx_p1)
       WRITE(dakota_unit,*) 'Pressure 2',qp(idx_p2)
       WRITE(dakota_unit,*) 'Liquid/particles velocity',qp(idx_u1)
       WRITE(dakota_unit,*) 'Gas velocity',qp(idx_u2)
       WRITE(dakota_unit,*) 'Exit Temperature',qp(idx_T)
       WRITE(dakota_unit,*) 'Crystals volume fraction',SUM(qp(idx_beta_first:   &
            idx_beta_last))

       CALL sound_speeds(C_mix,mach) 
       CALL update_radius(zeta)

       WRITE(dakota_unit,*) 'Mach number',mach
       WRITE(dakota_unit,*) 'Exit mixture velocity', mix_velocity
       WRITE(dakota_unit,*) 'Mass flow rate', pi * radius * radius &
            * mix_velocity * rho_mix
       WRITE(dakota_unit,*) 'Fragmentation', frag_eff

       CLOSE( dakota_unit )

       OPEN( dakota_unit2,FILE='dakota.out',STATUS='UNKNOWN')

       pi = 4.D0 * ATAN(1.D0)

       WRITE(dakota_unit2,*) gasTotVolFrac
       WRITE(dakota_unit2,*) qp(idx_p1)
       WRITE(dakota_unit2,*) qp(idx_p2)
       WRITE(dakota_unit2,*) qp(idx_u1)
       WRITE(dakota_unit2,*) qp(idx_u2)
       WRITE(dakota_unit2,*) qp(idx_T)
       WRITE(dakota_unit2,*) mach
       WRITE(dakota_unit2,*) mix_velocity
       WRITE(dakota_unit2,*) pi * radius * radius * mix_velocity * rho_mix
       WRITE(dakota_unit2,*) frag_eff

       CLOSE( dakota_unit2 )

       exit_file = TRIM(run_name)//'_exit.std'

       OPEN( exit_unit,FILE=exit_file,STATUS='UNKNOWN')

       WRITE(exit_unit,*) 'Mass_flow_rate',mass_flow_rate
       WRITE(exit_unit,*) 'Gas volume fraction',1.D0-SUM(qp(idx_alfa_first:     &
            idx_alfa_last))
       WRITE(exit_unit,*) 'Pressure 1',qp(idx_p1)
       WRITE(exit_unit,*) 'Pressure 2',qp(idx_p2)
       WRITE(exit_unit,*) 'Crystals volume fraction',SUM(qp(idx_beta_first:     &
            idx_beta_last))
       WRITE(exit_unit,*) 'Liquid/particles velocity',qp(idx_u1)
       WRITE(exit_unit,*) 'Gas velocity',qp(idx_u2)
       WRITE(exit_unit,*) 'Mach number',mach

       CLOSE( exit_unit )

    END IF

  END SUBROUTINE output_steady

  !------------------------------------------------------------------------------
  !> \brief Numeric to String conversion
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !> \date 27/20/2009
  !> \param   k      integer to convert             (\b input)           
  !------------------------------------------------------------------------------

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

