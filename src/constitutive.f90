!********************************************************************************
!> \brief Constitutive equations
!> @author 
!> Mattia de' Michieli Vitturi
!> This module contains the constitutive equations for the closure of the
!> system. 
!> \date 04/02/2017
!********************************************************************************
MODULE constitutive

  USE geometry, ONLY : pi
  USE parameters, ONLY : verbose_level
  USE parameters, ONLY : n_eqns , n_vars, n_components
  USE parameters, ONLY : n_cry , n_gas , n_mom

  IMPLICIT none

  !--------- Constants for the equations of state -------------------------------
  REAL*8 :: cv_2         !< exsolved gas specific heat capacity at constant volume
  REAL*8 :: cv_m         !< melt specific heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_c(:)  !< crystals specific heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_d(:)  !< dissolved gas heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_g(:)  !< exsolved gas heat capacity at constant volume

  !~   REAL*8 :: C0_2        !< exsolved gas sound speed at atmospheric conditions
  REAL*8 :: C0_m             !< melt sound speed at atmospheric conditions
  REAL*8, ALLOCATABLE :: C0_c(:)   !< crystals sound speed at atm conditions
  REAL*8, ALLOCATABLE :: C0_d(:)   !< dissolved gas sound speed at atm conditions
  REAL*8, ALLOCATABLE :: C0_g(:)   !< exsolved gas sound speed at atm conditions

  REAL*8 :: gamma_2          !< exsolved gas adiabatic exponent 
  REAL*8 :: gamma_m          !< melt adiabatic exponent 
  REAL*8, ALLOCATABLE :: gamma_c(:)   !< crystals adiabatic exponent 
  REAL*8, ALLOCATABLE :: gamma_d(:)   !< dissolved gas adiabatic exponent 
  REAL*8, ALLOCATABLE :: gamma_g(:)   !< exsolved gas adiabatic exponent 

  REAL*8 :: rho0_2           !< exsolved gas reference density
  REAL*8 :: rho0_m           !< melt reference density
  REAL*8, ALLOCATABLE :: rho0_c(:)    !< crystals reference density
  REAL*8, ALLOCATABLE :: rho0_d(:)    !< dissolved gas reference density
  REAL*8, ALLOCATABLE :: rho0_g(:)    !< exsolved gas reference density

  REAL*8 :: T0_2            !< exsolved gas reference temperature
  REAL*8 :: T0_m            !< melt reference temperature
  REAL*8, ALLOCATABLE :: T0_c(:)     !< crystals gas reference temperature
  REAL*8, ALLOCATABLE :: T0_d(:)     !< dissolved gas reference temperature
  REAL*8, ALLOCATABLE :: T0_g(:)     !< exsolved gas reference temperature

  !~   REAL*8 :: p0_2       !< exsolved gas reference pressure
  REAL*8 :: p0_m            !< melt reference pressure
  REAL*8, ALLOCATABLE :: p0_c(:)     !< crystals reference pressure
  REAL*8, ALLOCATABLE :: p0_d(:)     !< dissolved gas reference pressure
  REAL*8, ALLOCATABLE :: p0_g(:)     !< exsolved gas reference pressure

  REAL*8 :: bar_e_2         !< exsolved gas formation energy
  REAL*8 :: bar_e_m         !< melt formation energy
  REAL*8, ALLOCATABLE :: bar_e_c(:)  !< crystals formation energy
  REAL*8, ALLOCATABLE :: bar_e_d(:)  !< dissolved gas formation energy
  REAL*8, ALLOCATABLE :: bar_e_g(:)  !< exsolved gas formation energy

  REAL*8 :: bar_p_m         !< melt cohesion pressure
  REAL*8, ALLOCATABLE :: bar_p_c(:)  !< crystals cohesion pressure
  REAL*8, ALLOCATABLE :: bar_p_d(:)  !< dissolved gas cohesion pressure
  REAL*8, ALLOCATABLE :: bar_p_g(:)  !< exsolved gas cohesion pressure

  !~   REAL*8 :: s0_2      !< exsolved gas reference entropy
  ! The reference entropies are used only for non-isothermal runs
  REAL*8 :: s0_m           !< melt reference entropy
  REAL*8, ALLOCATABLE :: s0_c(:)    !< crystals reference entropy
  REAL*8, ALLOCATABLE :: s0_d(:)    !< dissolved gas reference entropy
  REAL*8, ALLOCATABLE :: s0_g(:)    !< exsolved gas reference entropy

  !-----------------------------------------------------------------------------------------!
  REAL*8, ALLOCATABLE :: Pc_g(:)           !< critical gas pressure
  REAL*8, ALLOCATABLE :: Tc_g(:)           !< critical gas temperature
  REAL*8, ALLOCATABLE :: a_g(:)            !< parameter for the VDW EOS
  REAL*8, ALLOCATABLE :: b_g(:)            !< parameter for the VDW EOS  

  !> equation of state for gas\n
  !> - 'ideal'          => ideal gas law
  !> - 'VDW'            => Van der Waals equation of state
  !> .
  CHARACTER*20 :: gas_law
  !----------------------------------------------------------------------------------------!

  COMPLEX*16 :: cv_1       !< dis.gas+melt+crystals specific heat capacity at constant volume

  COMPLEX*16 :: e_mix      !< total internal energy

  COMPLEX*16 :: e_1        !< local specific internal energy of the melt-crystals phase
  COMPLEX*16 :: e_2        !< local specific internal energy of the exsolved gas
  COMPLEX*16 :: e_m        !< local specific internal energy of the melt
  COMPLEX*16, ALLOCATABLE :: e_c(:) !< local specific internal energy of the crystals
  COMPLEX*16, ALLOCATABLE :: e_d(:) !< local specific internal energy of the dissolved gas
  COMPLEX*16, ALLOCATABLE :: e_g(:) !< local specific internal energy of the exsolved gas

  COMPLEX*16 :: p_1        !< local pressure of the melt-crystals phase
  COMPLEX*16 :: p_2        !< local pressure of the exsolved gas

  ! Entropy is used only for non-isothermal runs
  COMPLEX*16 :: s_1        !< local specific entropy of the melt-crystals phase
  COMPLEX*16 :: s_2        !< local specific entropy of the exsolved gas
  COMPLEX*16 :: s_m        !< local specific entropy of the melt
  COMPLEX*16, ALLOCATABLE :: s_c(:) !< local specific entropy of the crystals
  COMPLEX*16, ALLOCATABLE :: s_d(:) !< local specific entropy of the dissolved gas
  COMPLEX*16, ALLOCATABLE :: s_g(:) !< local specific entropy of the exsolved gas

  COMPLEX*16 :: mu_1       !< free Gibbs energy of the melt-crystals phase
  COMPLEX*16 :: mu_2       !< free Gibbs energy of the exsolved gas
  COMPLEX*16 :: mu_m       !< free Gibbs energy of the melt
  COMPLEX*16, ALLOCATABLE :: mu_c(:)!< free Gibbs energy of the crystals
  COMPLEX*16, ALLOCATABLE :: mu_d(:)!< free Gibbs energy of the dissolved gas
  COMPLEX*16, ALLOCATABLE :: mu_g(:)!< free Gibbs energy of the exsolved gas

  COMPLEX*16 :: rho_1         !< dis_gas+melt+crystals phase local density
  COMPLEX*16 :: rho_2         !< exsolved gas local density
  COMPLEX*16 :: rho_m         !< melt local density
  COMPLEX*16, ALLOCATABLE :: rho_c(:)  !< crystals local density
  COMPLEX*16, ALLOCATABLE :: rho_d(:)  !< dissolved gas local density
  COMPLEX*16, ALLOCATABLE :: rho_g(:)  !< exsolved gas local density
  COMPLEX*16 :: rho_md        !< dis_gas+melt local density

  COMPLEX*16, ALLOCATABLE :: rhoB_c(:)     !< crystals bulk density
  COMPLEX*16 :: rhoB_m     !< melt bulk density

  COMPLEX*16 :: alfa_1     !< dis_gas+melt+crystals phase local volume fraction
  COMPLEX*16 :: alfa_2     !< total exsolved gas local volume fraction
  COMPLEX*16, ALLOCATABLE :: alfa_g(:)     !< exsolved gas phases local volume fraction

  COMPLEX*16 :: alfarho_2  !< bulk density of the exsolved gas

  COMPLEX*16 :: alfa_m_1                   !< melt volume fraction in phase 1
  COMPLEX*16, ALLOCATABLE :: alfa_d_1(:)   !< dissolved gas volume fractions in phase 1
  COMPLEX*16, ALLOCATABLE :: alfa_g_2(:)   !< exsolved gas volume fractions in phase 2

  COMPLEX*16, ALLOCATABLE :: beta(:)       !< crystal volume fraction in the melt-crystals phase
  COMPLEX*16, ALLOCATABLE :: beta_eq(:)    !< equil. cry. volume fraction in the melt-crystals phase

  COMPLEX*16, ALLOCATABLE :: mom_cry(:,:,:)  !< moments of the crystal referred to the melt-crystals phase

  REAL*8, ALLOCATABLE:: growth_mom(:,:,:)   !< moments of growth rate of crystals

  COMPLEX*16, ALLOCATABLE:: cry_shape_factor(:) !< shape factor of crystals

  REAL*8, ALLOCATABLE:: cry_shape_factor_in(:) !< shape factor of crystals (input)
  
  REAL*8, ALLOCATABLE:: T_m(:) !< liquidus temperature of crystals

  REAL*8, ALLOCATABLE :: T_u(:) !< temp of max growth rate of crystals
  
  REAL*8, ALLOCATABLE :: U_m(:) !< max growth rate of crystals

  REAL*8, ALLOCATABLE :: I_m(:) !< max nucleation rate of crystals

  REAL*8, ALLOCATABLE:: T_i(:) !< temp of max nucleation rate of crystals

  COMPLEX*16, ALLOCATABLE:: L0_cry(:,:) !< initial size of crystals (conduit bottom) 

  REAL*8, ALLOCATABLE :: L0_cry_in(:) !< initial size of phenocryst (conduit bottom)
  
  COMPLEX*16, ALLOCATABLE:: L_nucleus(:) !< size of new nucleus (input)

  REAL*8, ALLOCATABLE:: L_nucleus_in(:) !< size of new nucleus
  
  REAL*8, ALLOCATABLE :: cry_init_solid_solution(:,:) !< initial composition of phenocrysts
  
  REAL*8, ALLOCATABLE:: cry_current_solid_solution(:,:) !< composition of crystallizing minerals
  
  COMPLEX*16, ALLOCATABLE :: rhoB_components(:)   !< components bulk density
  
  COMPLEX*16, ALLOCATABLE :: sum_rhoB_components(:)   !< sum of components bulk density in the different crystals
  
  COMPLEX*16 :: u_1        !< melt-crystals phase local velocity
  COMPLEX*16 :: u_2        !< exsolved gas local velocity

  COMPLEX*16 :: x_1           !< melt-crystals phase local mass fraction
  COMPLEX*16 :: x_2           !< exsolved gas local mass fraction
  COMPLEX*16, ALLOCATABLE :: x_c(:)    !< crystals mass fraction (with respect to the mixture)
  COMPLEX*16 :: x_m           !< melt mass fraction (with respect to the mixture)
  COMPLEX*16, ALLOCATABLE :: x_d(:)    !< dissolved gas mass fraction (with respect to the mixture)
  COMPlEX*16, ALLOCATABLE :: x_g(:)    !< exsolved gas mass fraction (with respect to the mixture)

  COMPLEX*16, ALLOCATABLE :: x_d_md(:)     !< dissolved gas mass fraction in the melt+dis.gas phase
  COMPLEX*16, ALLOCATABLE :: x_d_md_eq(:)  !< equil. dis. gas mass fraction in the melt+dis.gas phase

  COMPLEX*16, ALLOCATABLE :: x_c_mc(:)  !< cry. mass fraction in the melt+cry phase

  COMPLEX*16, ALLOCATABLE :: x_c_1(:)      !< cristal mass fractions in phase 1
  COMPLEX*16 :: x_m_1             !< melt mass fraction in phase 1
  COMPLEX*16 :: x_md_1            !< melt+dis.gas mass fraction in phase 1
  COMPLEX*16, ALLOCATABLE :: x_d_1(:)      !< dissolved gas mass fractions in phase 1
  COMPLEX*16, ALLOCATABLE :: x_g_2(:)      !< exsolved gas mass fractions in phase 2

  COMPLEX*16 :: T          !< mixture local temperature
  COMPLEX*16 :: rho_mix    !< mixture local density
  COMPLEX*16 :: u_mix      !< mixture velocity
  COMPLEX*16 :: cv_mix     !< mixture specific heat capacity at constant volume
  COMPLEX*16 :: visc_mix   !< mixture viscosity

  COMPLEX*16 :: S          !< mixture entropy

  COMPLEX*16 :: C_1        !< first phase local sound speed
  COMPLEX*16 :: C_2        !< second phase local sound speed

  INTEGER :: n_drag_models

  CHARACTER (LEN=30), DIMENSION(20) :: available_drag_models
  
  !> drag function model\n
  !> - 'constant'        => drag_funct = 1
  !> - 'eval'            => drag_funct = f(q)
  !> - 'darcy'           => Darcy formulation from Eq. 16 Degruyter et al. 2012 
  !> - 'forchheimer'     => Forcheeimer from Eq. 16 Degruyter et al. 2012
  !> - 'drag'            => Stokes law for rising bubbles
  !> - 'single_velocity' => Single velocity model (u1=u2)
  !> .
  CHARACTER*30 :: drag_funct_model

  !> drag function for the relative velocity:\n
  COMPLEX*16 :: drag_funct

  !> coefficient for the drag function for the relative velocity:\n
  !> - drag_funct_coeff > 0   => positive drag (finite rate velocity relaxation)
  !> - drag_funct_coeff \f$\rightarrow +\infty \f$   => instantaneous relaxation
  !> - drag_funct_coeff = 0   => no velocity relaxation
  !> .
  REAL*8 :: drag_funct_coeff

  ! parameters of the Forchheimer model

  !> bubble number density (is referred to the liquid volume fraction)
  REAL*8 :: bubble_number_density
  REAL*8 :: log10_bubble_number_density

  !> tortuosity factor
  REAL*8 :: tortuosity_factor

  !> throat bubble ratio
  REAL*8 :: throat_bubble_ratio

  ! friction coefficient
  REAL*8 :: friction_coefficient

  ! drag coefficient
  REAL*8 :: C_D 

  ! ash particl size
  REAL*8 :: r_a

  !> pressure relaxation model\n
  !> - 'constant'    => tau_p = 1
  !> - 'eval'        => tau_p = f(q)
  !> .
  CHARACTER*20 :: p_relax_model

  !> pressure relaxation rate
  COMPLEX*16 :: tau_p

  !> pressure relaxation coefficient:\n
  !> - tau_p_coeff = 0     => instantaneous relaxation (sinlge pressure);
  !> - tau_p_coeff > 0     => finite rate relaxation;
  !> .
  REAL*8 :: tau_p_coeff

  !> crystallization parameter:\n
  !> - tau_c = 0     => instantaneous crystallization (equilibrium);
  !> - tau_c > 0     => finite rate crystallization;
  !> - tau_c = -1    => no crystallization.
  !> .
  REAL*8, ALLOCATABLE :: tau_c(:)

  !> chamber (equilibrium) crystal volume fraction
  REAL*8, ALLOCATABLE :: beta0(:)

  !> maximum crystal volume fraction
  REAL*8, ALLOCATABLE :: beta_max(:)

  !> Model for the equilibrium crystal volume fraction:\n
  !> - 'Vitturi2010'    => Eq. 4 of de' Michieli Vitturi et al. 2010;
  !> - 'None'           => Without crystallization
  !> .
  CHARACTER*20 :: crystallization_model

  !> index of fragmentation in the interval [0;1]
  REAL*8 :: frag_eff

  !> Parameter to choose the fragmentation model:\n
  !> - fragmentation_model = 1  => volume fraction
  !> - fragmentation_model = 2  => overpressure
  !> - fragmentation_model = 3  => strain rate
  !> .
  INTEGER :: fragmentation_model

  !> Threshold for the fragmentation
  REAL*8 :: frag_thr

  !> Fragmentation exponent
  REAL*8 :: tau_frag_exp

  !> fragmentation rate
  COMPLEX*16 :: tau_frag

  !> fragmentation coefficient:\n
  !> - tau_frag_coeff = 0     => instantaneous fragmentation;
  !> - tau_frag_coeff > 0     => finite rate fragmentation;
  !> .
  REAL*8 :: tau_frag_coeff

  !> exsolution parameter:\n
  !> - tau_d = 0     => instantaneous exsolution (equilibrium);
  !> - tau_d > 0     => finite rate exsolution;
  !> - tau_d = -1    => no exsolution.
  !> .
  REAL*8, ALLOCATABLE :: tau_d(:)

  !> String for exsolution model:\n
  !> - 'Henry' => Henry's law;
  !> - 'Polynomial' => Polinomial law;
  !> - 'Zhang' => Zhang model (Zhang, JVGR '99).
  !> .
  CHARACTER*20 :: exsol_model

  !> Solubility parameter for the Henry's and polynomial law
  REAL*8, ALLOCATABLE :: solub(:)

  REAL*8, ALLOCATABLE :: solub_exp(:)

  REAL*8, ALLOCATABLE :: x_ex_dis_in(:)

  !> gravitational acceleration
  REAL*8 :: grav

  !> country rock permeability
  REAL*8 :: log10_k_cr
  REAL*8 :: k_cr

  !> country rock density
  REAL*8 :: rho_cr

  !> gas viscosity
  REAL*8 :: visc_2

  !> melt viscosity
  COMPLEX*16 :: visc_melt

  !> melt+crystal viscosity
  COMPLEX*16 :: visc_1

  !> relative viscosity due to bubbles
  COMPLEX*16 :: visc_rel_bubbles

  INTEGER :: n_theta_models 

  CHARACTER (LEN=30), DIMENSION(20) :: available_theta_models
  
  !> Parameter to choose the model for the influence of crystal on the mixture:
  !> 'Lejeune_and_Richet1995'
  !> 'Dingwell1993'
  !> 'Melnik_and_Sparks1999'
  !> 'Costa2005'
  !> 'Melnik_and_Sparks2005'
  !> 'Vona_et_al2011'
  !> 'Vona_et_al2011_mod'
  !> 'Vona_et_al2013_eq19'
  !> 'Vona_et_al2013_eq20'
  !> 'Vona_et_al2013_eq21'
  !> 'Fixed_value'
  !> .  
  CHARACTER*30 :: theta_model

  !> Relative viscosity of the crystals
  COMPLEX*16 :: theta

  !> Fixed value for the relative viscosity of the crystals
  REAL*8 :: theta_fixed

  !> Coefficients for the relative viscosity models
  REAL*8 :: c1 , c2 , c3

  !> Lithostatic pressure
  REAL*8 :: p_lith

  !> Hydrostatic pressure
  REAL*8 :: p_hydro

  !> Elevation above the bottom for the evaluation of the lithostatic pressure
  REAL*8 :: zeta_lith

  INTEGER :: n_bubble_models

  CHARACTER (LEN=30), DIMENSION(20) :: available_bubble_models
  
  !> Parameter to choose the model for the influence of the bubbles on the mixture:\n 
  !> - 'Quane-Russel         -> For Campi-Flegrei we use 0.63 as reported in the 2004 Report (Task 2.2)
  !> - 'Eilers'              -> Mader et al. 2013
  !> - 'Sibree'              -> Mader et al. 2013
  !> - 'Taylor'              -> Mader et al. 2013
  !> - 'Mackenzie'           -> Mader et al. 2013
  !> - 'DucampRaj'           -> Mader et al. 2013
  !> - 'BagdassarovDingwell' -> Mader et al. 2013
  !> - 'Rahaman'             -> Mader et al. 2013
  !> .
  CHARACTER*20 :: bubbles_model

  !> Flag to choose the eruptive style:\n
  !> - explosive_flag = .TRUE.   => explosive eruption
  !> - explosive_flag = .FALSE.  => effusive eruption
  LOGICAL :: explosive_flag

  !> Magma permeability
  COMPLEX*16 :: permkc 
  COMPLEX*16 :: inv_permkc  

  REAL*8 :: alfa_switch
  REAL*8 :: a_2nd , b_2nd , perm0

  INTEGER :: n_visc_melt_models

  CHARACTER (LEN=30), DIMENSION(20) :: available_visc_melt_models
  
  !> Parameter to select the melt viscosity (bubbles and crystal-free) model:\n
  !> - 'Hess_and_Dingwell1996'
  !> - 'Whittington_et_al2008'
  !> - 'Romano_et_al2003t'
  !> - 'Romano_et_al2003p'
  !> - 'Giordano_et_al2008'
  !> - 'Giordano_et_al2009'
  !> - 'Di_Genova_et_al2013_eqn_3,5'
  !> - 'Di_Genova_et_al2013_eqn_4,5'
  !> .
  CHARACTER*30 :: visc_melt_model
  
  REAL*8, DIMENSION(12) :: wt_init

  !> Water density
  REAL*8 :: rho_w

  REAL*8 :: xa,xb,xc

CONTAINS

  SUBROUTINE initialize_models

    IMPLICIT NONE

    n_drag_models = 6
    
    available_drag_models(1) = 'constant'
    available_drag_models(2) = 'eval'
    available_drag_models(3) = 'darcy'
    available_drag_models(4) = 'forchheimer'
    available_drag_models(5) = 'drag'
    available_drag_models(6) = 'single_velocity'

    n_theta_models = 11

    available_theta_models(1) = 'Lejeune_and_Richet1995'
    available_theta_models(2) = 'Dingwell1993'
    available_theta_models(3) = 'Melnik_and_Sparks1999'
    available_theta_models(4) = 'Costa2005'
    available_theta_models(5) = 'Melnik_and_Sparks2005'
    available_theta_models(6) = 'Vona_et_al2011'
    available_theta_models(7) = 'Vona_et_al2011_mod'
    available_theta_models(8) = 'Vona_et_al2013_eq19'
    available_theta_models(9) = 'Vona_et_al2013_eq20'
    available_theta_models(10) = 'Vona_et_al2013_eq21'
    available_theta_models(11) = 'Fixed_value'

    n_bubble_models = 10

    available_bubble_models(1) = 'none'
    available_bubble_models(2) = 'Costa2007'
    available_bubble_models(3) = 'Quane-Russel'
    available_bubble_models(4) = 'Eilers'
    available_bubble_models(5) = 'Sibree'
    available_bubble_models(6) = 'Taylor'
    available_bubble_models(7) = 'Mackenzie'
    available_bubble_models(8) = 'DucampRaj'
    available_bubble_models(9) = 'BagdassarovDingwell'
    available_bubble_models(10) = 'Rahaman'

    n_visc_melt_models = 8

    available_visc_melt_models(1) = 'Hess_and_Dingwell1996'
    available_visc_melt_models(2) = 'Whittington_et_al2008'
    available_visc_melt_models(3) = 'Romano_et_al2003t'
    available_visc_melt_models(4) = 'Romano_et_al2003p'
    available_visc_melt_models(5) = 'Giordano_et_al2008'
    available_visc_melt_models(6) = 'Giordano_et_al2009'
    available_visc_melt_models(7) = 'Di_Genova_et_al2013_eqn_3,5'
    available_visc_melt_models(8) = 'Di_Genova_et_al2013_eqn_4,5'
    
  END SUBROUTINE initialize_models
  
  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE allocate_phases_parameters

    ALLOCATE( cv_c(n_cry) )
    ALLOCATE( C0_c(n_cry) )
    ALLOCATE( gamma_c(n_cry) )
    ALLOCATE( rho0_c(n_cry) )
    ALLOCATE( T0_c(n_cry) )
    ALLOCATE( p0_c(n_cry) )
    ALLOCATE( bar_e_c(n_cry) )
    ALLOCATE( bar_p_c(n_cry) )
    ALLOCATE( s0_c(n_cry) )
    ALLOCATE( e_c(n_cry) )
    ALLOCATE( s_c(n_cry) )
    ALLOCATE( mu_c(n_cry) )
    ALLOCATE( rho_c(n_cry) )
    ALLOCATE( rhoB_c(n_cry) )
    ALLOCATE( beta(n_cry) )
    ALLOCATE( beta_eq(n_cry) )
    ALLOCATE( x_c_mc(n_cry) )
    ALLOCATE( x_c(n_cry) )
    ALLOCATE( x_c_1(n_cry) )
    ALLOCATE( tau_c(n_cry) )
    ALLOCATE( beta0(n_cry) )
    ALLOCATE( beta_max(n_cry) )
    ALLOCATE( cv_d(n_gas) )
    ALLOCATE( C0_d(n_gas) )
    ALLOCATE( gamma_d(n_gas) )
    ALLOCATE( rho0_d(n_gas) )
    ALLOCATE( T0_d(n_gas) )
    ALLOCATE( p0_d(n_gas) )
    ALLOCATE( bar_e_d(n_gas) )
    ALLOCATE( bar_p_d(n_gas) )
    ALLOCATE( s0_d(n_gas) )
    ALLOCATE( e_d(n_gas) )
    ALLOCATE( s_d(n_gas) )
    ALLOCATE( mu_d(n_gas) )
    ALLOCATE( rho_d(n_gas) )
    ALLOCATE( alfa_d_1(n_gas) )

    ALLOCATE( x_d(n_gas) )
    ALLOCATE( x_d_md(n_gas) )
    ALLOCATE( x_d_md_eq(n_gas) )
    ALLOCATE( x_d_1(n_gas) )
    ALLOCATE( tau_d(n_gas) )
    ALLOCATE( solub(n_gas) )
    ALLOCATE( solub_exp(n_gas) )

    ALLOCATE( x_ex_dis_in(n_gas) )

    ALLOCATE( cv_g(n_gas) )
    ALLOCATE( C0_g(n_gas) )
    ALLOCATE( gamma_g(n_gas) )
    ALLOCATE( rho0_g(n_gas) )
    ALLOCATE( T0_g(n_gas) )
    ALLOCATE( p0_g(n_gas) )
    ALLOCATE( bar_e_g(n_gas) )
    ALLOCATE( bar_p_g(n_gas) )
    ALLOCATE( s0_g(n_gas) )
    ALLOCATE( e_g(n_gas) )
    ALLOCATE( s_g(n_gas) )
    ALLOCATE( mu_g(n_gas) )
    ALLOCATE( rho_g(n_gas) )
    ALLOCATE( x_g(n_gas) )
    ALLOCATE( x_g_2(n_gas) )
    ALLOCATE( alfa_g(n_gas) )
    ALLOCATE( alfa_g_2(n_gas) )

    ALLOCATE( Pc_g(n_gas) )
    ALLOCATE( Tc_g(n_gas) )
    ALLOCATE( a_g(n_gas) )
    ALLOCATE( b_g(n_gas) )

  END SUBROUTINE allocate_phases_parameters


  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Equation of state
  !
  !> This subroutine evaluates, from the mixture temperature T, the phase 
  !> densities \f$ \rho_i \f$ and the mass fractions, the values of 
  !> \f$ e_i, p_i, s_i, \mu_i \f$ for the two phases.
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE eos

    USE complexify 
    IMPLICIT none

    e_c(1:n_cry) = DCMPLX(bar_e_c(1:n_cry),0.D0) + DCMPLX(cv_c(1:n_cry),0.D0)   &
         * T + DCMPLX(bar_p_c(1:n_cry),0.D0) / rho_c(1:n_cry)

    e_m = DCMPLX(bar_e_m,0.D0) + DCMPLX(cv_m,0.D0) * T +                        &
         DCMPLX(bar_p_m,0.D0) / rho_m

    e_d(1:n_gas) = DCMPLX(bar_e_d(1:n_gas),0.D0) + DCMPLX(cv_d(1:n_gas),0.D0)   &
         * T + DCMPLX(bar_p_d(1:n_gas),0.D0) / rho_d(1:n_gas)

    e_1 = SUM( x_c_1(1:n_cry) * e_c(1:n_cry) ) + x_m_1 * e_m +                  &
         SUM( x_d_1(1:n_gas) * e_d(1:n_gas) )

    e_g(1:n_gas) = DCMPLX(bar_e_g(1:n_gas),0.D0) + DCMPLX(cv_g(1:n_gas),0.D0)   &
         * T - DCMPLX(a_g(1:n_gas),0.0) * rho_g(1:n_gas)

    e_2 = SUM( x_g_2(1:n_gas) * e_g(1:n_gas) )

    s_c(1:n_cry) = s0_c(1:n_cry) + cv_c(1:n_cry) * log( T / T0_c(1:n_cry) *     &
         ( rho0_c(1:n_cry) / rho_c(1:n_cry) )**( gamma_c(1:n_cry) - 1.D0 ) )

    s_m = s0_m + cv_m * log( T / T0_m * ( rho0_m / rho_m )**( gamma_m - 1.D0 ) )

    s_d(1:n_gas) = s0_d(1:n_gas) + cv_d(1:n_gas) * log( T / T0_d(1:n_gas) *     &
         ( rho0_d(1:n_gas) / rho_d(1:n_gas) )**( gamma_d(1:n_gas) - 1.D0 ) )

    s_1 = SUM( x_c_1(1:n_cry) * s_c(1:n_cry) ) + x_m_1 * s_m +                  &
         SUM( x_d_1(1:n_gas) * s_d(1:n_gas) )

    s_g(1:n_gas) = s0_g(1:n_gas) + cv_g(1:n_gas) * log( T / T0_g(1:n_gas) *     &
         ( ( rho0_g(1:n_gas) / rho_g(1:n_gas) ) * ( DCMPLX(1.D0,0.D0) -         &
         b_g(1:n_gas) * rho_g(1:n_gas) ) ) ** ( gamma_g(1:n_gas) - 1.D0 ) )

    s_2 = SUM( x_g_2(1:n_gas) * s_g(1:n_gas) )

    mu_m = e_m + p_1/rho_m - T * s_m

    mu_c(1:n_cry) = e_c(1:n_cry) + p_1/rho_c(1:n_cry) - T * s_c(1:n_cry)

    mu_d(1:n_gas) = e_d(1:n_gas) + p_1/rho_d(1:n_gas) - T * s_d(1:n_gas)

    mu_1 = SUM( x_c_1(1:n_cry) * mu_c(1:n_cry) ) + x_m_1 * mu_m +               &
         SUM( x_d_1(1:n_gas) * mu_d(1:n_gas) )

    mu_g(1:n_gas) = e_g(1:n_gas) + p_2/rho_g(1:n_gas) - T * s_g(1:n_gas)

    mu_2 = SUM( x_g_2(1:n_gas) * mu_g(1:n_gas) )

    S = x_1 * s_1 + x_2 * s_2

  END SUBROUTINE eos

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Local sound speeds
  !
  !> This subroutine evaluates the local sound speed of the different phases 
  !> as:\n
  !> \f$ C_i = \sqrt \frac{\partial p_i}{\partial \rho_i} \f$
  !> and the mixture sound speed as:\n
  !> \f$ C_{mix}= \sqrt 1/(K_{mix}\rho_{mix})
  !> where the compressibility K_mix is the volumetric average of the 
  !> compressibilities of the different phases.
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE sound_speeds(C_mix,mach) 

    USE complexify 
    IMPLICIT none

    REAL*8, INTENT(OUT) :: C_mix , mach

    !> Sound speeds of the phases
    COMPLEX*16 :: C2_c(1:n_cry) , C2_m , C2_d(1:n_gas) , C2_1 , C2_2
    COMPLEX*16 :: C2_g(1:n_gas)

    !> Compressibilites of the phases
    COMPLEX*16 :: K_1 , K_2 , K_mix
    COMPLEX*16 :: K_c(1:n_cry) , K_m , K_d(1:n_gas) , K_g(1:n_gas)

    !REAL*8 :: C_mix_Wallis

    C2_c(1:n_cry) = C0_c(1:n_cry)**2.D0 * ( rho_c(1:n_cry) / rho0_c(1:n_cry) )  &
         ** ( gamma_c(1:n_cry) - DCMPLX(1.0,0.D0) ) * CDEXP( ( s_c(1:n_cry) -   &
         s0_c(1:n_cry) ) / cv_c(1:n_cry) )

    C2_m = C0_m**2.D0 * ( rho_m / rho0_m ) **                                   &
         ( gamma_m - DCMPLX(1.0,0.D0) ) * CDEXP( ( s_m - s0_m ) / cv_m )

    C2_d(1:n_gas) = C0_d(1:n_gas)**2.D0 * ( rho_d(1:n_gas) / rho0_d(1:n_gas) )  &
         ** ( gamma_d(1:n_gas) - DCMPLX(1.0,0.D0) ) * CDEXP( ( s_d(1:n_gas) -   &
         s0_d(1:n_gas) ) / cv_d(1:n_gas) )

    K_c(1:n_cry) = 1.D0 / ( rho_c(1:n_cry) * C2_c(1:n_cry) )
    K_m = 1.D0 / ( rho_m * C2_m )
    K_d(1:n_gas) = 1.D0 / ( rho_d(1:n_gas) * C2_d(1:n_gas) )

    K_1 = alfa_m_1*K_m + SUM( beta(1:n_cry)*K_c(1:n_cry) )                      &
         + SUM( alfa_d_1(1:n_gas) * K_d(1:n_gas) )

    C2_1 = 1.D0 / ( K_1 * rho_1 )

    C_1 = CDSQRT( C2_1 )

    C2_g(1:n_gas) = cv_g(1:n_gas) * ( gamma_g(1:n_gas) - DCMPLX(1.0,0.D0) )     &
         * gamma_g(1:n_gas) * T / ( DCMPLX(1.0,0.D0) - b_g(1:n_gas) *           &
         rho_g(1:n_gas) )**2.0 - DCMPLX(2.0,0.D0) * a_g(1:n_gas) * rho_g(1:n_gas)

    K_g(1:n_gas) = 1.D0 / ( rho_g(1:n_gas) * C2_g(1:n_gas) )

    K_2 = SUM( alfa_g_2(1:n_gas) * K_g(1:n_gas) ) 

    C2_2 = 1.D0 / ( K_2 * rho_2 )

    C_2 = CDSQRT( C2_2 )

    K_mix = alfa_1*K_1 + alfa_2*K_2

    C_mix = DREAL( 1.D0 / CDSQRT( K_mix * rho_mix ) )

    ! Sound speed of the mixture (Wallis, 1969)
    ! C_mix = REAL( C_2 / SQRT(x_2) * ( x_2 + (1.d0 - x_2 ) * rho_2 / rho_1 ) )

    mach = DREAL( u_mix / C_mix )

    IF ( verbose_level .GE. 2 ) THEN

       WRITE(*,*) 'K_1,k_2',K_1,K_2
       WRITE(*,*) 'K_mix,rho_mix',K_mix,rho_mix
       WRITE(*,*) 'u_mix,C_mix',u_mix,C_mix

    END IF


  END SUBROUTINE sound_speeds

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Phases densities with MoM
  !
  !> This subroutine evaluates some densities needed of when the MoM is used.
  !> \f$ \rho_i = \frac{p_i + \bar{p}_i}{c_{v,i}(\gamma_i-1)T} \f$
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE eval_densities_mom

    IMPLICIT none

    rho_c(1:n_cry) = ( p_1 + bar_p_c(1:n_cry) ) / ( T * cv_c(1:n_cry) *         &
         ( gamma_c(1:n_cry) - DCMPLX(1.D0,0.D0) ) )
    
    rho_m = ( p_1 + bar_p_m ) / ( T * cv_m * ( gamma_m - DCMPLX(1.D0,0.D0) ) )
    rho_d(1:n_gas) = ( p_1 + bar_p_d(1:n_gas) ) / ( T * cv_d(1:n_gas) *         &
         ( gamma_d(1:n_gas) - DCMPLX(1.D0,0.D0) ) )

    rho_md = DCMPLX(1.D0,0.D0) / ( SUM( x_d_md(1:n_gas) / rho_d(1:n_gas) )      &
         + ( DCMPLX(1.D0,0.D0) - SUM( x_d_md(1:n_gas) ) ) / rho_m )
    
  END SUBROUTINE eval_densities_mom
  
  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Phases densities
  !
  !> This subroutine evaluates the densities of the two phases given the 
  !> pressures p1 and p2, the temperature T and the mass fraction $x_{d,md}$.
  !> \f$ \rho_i = \frac{p_i + \bar{p}_i}{c_{v,i}(\gamma_i-1)T} \f$
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE eval_densities

    IMPLICIT none

    COMPLEX*16 :: a,b,c
    COMPLEX*16 :: y1,y2,y3

    COMPLEX*16 :: coeff1 , coeff2 , coeff3 

    INTEGER :: i

    rho_c(1:n_cry) = ( p_1 + bar_p_c(1:n_cry) ) / ( T * cv_c(1:n_cry) *         &
         ( gamma_c(1:n_cry) - DCMPLX(1.D0,0.D0) ) )

    rho_m = ( p_1 + bar_p_m ) / ( T * cv_m * ( gamma_m - DCMPLX(1.D0,0.D0) ) )
    rho_d(1:n_gas) = ( p_1 + bar_p_d(1:n_gas) ) / ( T * cv_d(1:n_gas) *         &
         ( gamma_d(1:n_gas) - DCMPLX(1.D0,0.D0) ) )

    rho_md = DCMPLX(1.D0,0.D0) / ( SUM( x_d_md(1:n_gas) / rho_d(1:n_gas) )      &
         + ( DCMPLX(1.D0,0.D0) - SUM( x_d_md(1:n_gas) ) ) / rho_m )

    rho_1 = SUM( beta(1:n_cry) * rho_c(1:n_cry) ) + ( DCMPLX(1.D0,0.D0) -       &
         SUM( beta(1:n_cry) ) ) * rho_md

    IF ( gas_law .EQ. 'IDEAL' )  THEN

       rho_g(1:n_gas) = ( p_2 ) / ( T * cv_g(1:n_gas) *                         &
            ( gamma_g(1:n_gas) - DCMPLX(1.D0,0.D0) ) )

    ELSEIF (  gas_law .EQ. 'VDW' ) THEN

       DO i=1,n_gas

          a = a_g(i) * DCMPLX(-1.D0,0.D0)
          b = cv_g(i) * ( gamma_g(i) - DCMPLX(1.D0,0.D0) ) * T  + p_2 * b_g(i)
          c = - p_2

          coeff1 = a / ( a_g(i) * b_g(i) )
          coeff2 = b / ( a_g(i) * b_g(i) )
          coeff3 = c / ( a_g(i) * b_g(i) )

          CALL solve_cubic( coeff1 , coeff2 , coeff3 , y1 , y2 , y3 )

          rho_g(i) = y1

       END DO

    END IF

    rho_2 = SUM( alfa_g_2(1:n_gas) * rho_g(1:n_gas) ) 

  END SUBROUTINE eval_densities

  !******************************************************************************
  !> \brief Solution of a cubic
  !
  !> This subroutine solves a cubic polynomial in the form y^3+ay^2+by+c=0.
  !> When the coefficients are real, the solutions found are real.
  !> \param[in]     a     polynomial coefficient 
  !> \param[in]     b     polynomial coefficient 
  !> \param[in]     c     polynomial coefficient 
  !> \param[out]    y1    solution
  !> \param[out]    y2    solution
  !> \param[out]    y3    solution
  !> \date 20/06/2013
  !******************************************************************************

  SUBROUTINE solve_cubic(a,b,c,y1,y2,y3)

    USE complexify
    IMPLICIT NONE

    COMPLEX*16, INTENT(IN) :: a,b,c
    COMPLEX*16, INTENT(OUT) :: y1,y2,y3

    COMPLEX*16 :: Q,R,D_temp
    COMPLEX*16 :: phi , arg , u , v

    Q = ( 3.D0*b - a**2) / DCMPLX(9.D0,0.D0)
    R = ( 9.D0*a*b - 27.D0*c - 2.D0*a**3 ) / DCMPLX(54.D0,0.D0)

    D_temp = Q**3 + R**2

    IF ( D_temp .GE. 0 ) THEN

       arg = R + CDSQRT(D_temp)

       u = SIGN(1.D0,REAL(arg))*(abs(arg))**(1.D0/3.D0)

       arg = R - CDSQRT(D_temp)
       v = SIGN(1.D0,REAL(arg))*(abs(arg))**(1.D0/3.D0)

       y1 = u + v - a / DCMPLX(3.D0,0.D0)

    ELSE

       phi = ACOS(R/CDSQRT(-(Q**3)))

       y1 = 2.D0 * CDSQRT(-Q) * COS(phi/ DCMPLX(3.D0,0.D0) )                    &
            - a / DCMPLX(3.D0,0.D0)

       y2 = 2.D0 * CDSQRT(-Q) * COS((phi+2.D0*pi)/ DCMPLX(3.D0,0.D0) )          &
            - a / DCMPLX(3.D0,0.D0)

       y3 = 2.D0 * CDSQRT(-Q) * COS((phi-2.D0*pi)/ DCMPLX(3.D0,0.D0) )          &
            - a/ DCMPLX(3.D0,0.D0)

    END IF

  END SUBROUTINE solve_cubic


  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Equilibrium Dissolved gas
  !
  !> This subrotine evaluates the equilibrium dissolved gas mass fraction by 
  !> means of the Henry's law:\n
  !> \f$ x_{d,md}^{eq} = s\sqrt P\ (kg/m^3)\f$\n
  !> where \f$s\f$ is a constant. If the pressure
  !> is smaller than nucleation pressure then \f$ x_{d,md}^{eq} = x_{tot}\f$ (total 
  !> mass gas fraction at the inlet). 
  !> \param[in]     pres               liquid pressure
  !> \param[in]     temp               temperature
  !> \param[in]     xtot               total gas fractions in the mixture
  !> \param[in]     alfa_g_2           exsolved gas volume fractions in phase 2
  !> \param[in]     x_g                gas mass fractions
  !> \param[in]     rho_md             dis_gas+melt density
  !> \param[out]    x_d_md_eq          equilibrium dissolved gas
  !> \date 08/10/12       
  !******************************************************************************

  SUBROUTINE f_xdis_eq

    USE complexify 

    IMPLICIT NONE

    COMPLEX*16 :: aa , bb , cc , pp
    INTEGER :: j

    IF ( REAL(p_2) .LE. 0.D0 ) THEN

       x_d_md_eq = DCMPLX(0.D0, 0.D0)

    ELSE

       SELECT CASE ( exsol_model )

       CASE DEFAULT


       CASE ( 'Henry' )

          ! Henry's law

          x_d_md_eq(1:n_gas) = solub(1:n_gas) * (alfa_g_2(1:n_gas) * p_2 )      &
               ** solub_exp(1:n_gas)

       CASE ( 'Polynomial' )

          ! Polynomial law

          x_d_md_eq(1:n_gas) = solub(1:n_gas) * (alfa_g_2(1:n_gas) * p_2 )      &
            * (alfa_g_2(1:n_gas)*p_2) + solub_exp(1:n_gas)*(alfa_g_2(1:n_gas)*p_2)

          DO j=1,n_gas
  
             IF( (- solub_exp(j) / solub(j) / 2.D0) .LT. (alfa_g_2(j)*p_2) ) THEN
      
                x_d_md_eq(j) = ( - solub_exp(j) ** 2.D0 / solub(j) / 4.D0)
     
             ENDIF
  
          ENDDO

       CASE ( 'Zhang' )

          ! Zhang

          aa = 0.4874d0 - 608.0d0 / T + 489530.d0 / ( T * T )

          bb = -0.06062d0 + 135.6d0 / T - 69200.d0 / ( T * T )

          cc = 0.00253d0 - 4.154d0 / T + 1509.0d0 / ( T * T )

          pp = p_2 * 1.0D-6

          x_d_md_eq(1:n_gas) = 1.0D-2 * ( aa * CDSQRT(alfa_g_2(1:n_gas) * pp) +  &
             bb * (alfa_g_2(1:n_gas) * pp) + cc * (alfa_g_2(1:n_gas) * pp) ** 1.5D0 )
          
       END SELECT

    END IF

  END SUBROUTINE f_xdis_eq

  !*****************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Equilibrium Crystal content
  !
  !> This subrotine evaluates the equilibrium crystal volume fraction as a 
  !> function of the pressure:\n
  !> \f$ \beta = \min\left( beta_{max} , 0.45 + 0.55 (0.58815 (\frac{p}{1.D6})^{-0.5226} 
  !>                 \right) \f$,\n
  !> \param[in]     press            liquid pressure
  !> \param[in]     temp             magmatic mixture temperature
  !> \param[in]     x_d               dissolved gas mass fraction
  !> \param[in]     rho_c            density of crystal components
  !> \param[in]     rho_1            density of phase 1
  !> \param[out]    beta_eq          equilibrium crystal volume fraction
  !> \date 08/10/2012 
  !*****************************************************************************

  SUBROUTINE f_beta_eq

    USE complexify 
    IMPLICIT NONE

    COMPLEX*16 :: x_d_md_tot ,x_d_md_wt_tot 
    COMPLEX*16 :: p_1_bar, T_celsius
    ! COMPLEX*16 :: crystal_mass_fraction(1:n_cry)

    INTEGER :: j

    p_1_bar = p_1 / 1.0e5
    T_celsius = T - 273.D0

    x_d_md_tot = SUM( x_d_md(1:n_gas) )
    x_d_md_wt_tot = x_d_md_tot * 100.D0

    SELECT CASE ( crystallization_model )

    CASE DEFAULT

    CASE ( 'Vitturi2010' )

       DO j=1,n_cry
          
          !----------------------------------------------------------------------
          beta_eq(j)=beta0(j) + 0.55D0*( 0.58815D0*( p_1/1.D6 )**( -0.5226D0 ) )
          !----------------------------------------------------------------------
          
          beta_eq(j) = MAX( beta0(j) + 1D-15, beta_eq(j) )
          beta_eq(j) = MIN( beta_max(j) , beta_eq(j) )

       END DO

    CASE ( 'None' )

       DO j=1,n_cry
          
          !----------------------------------------------------------------------
          beta_eq(j) = beta0(j)
          !----------------------------------------------------------------------
          
          beta_eq(j) = MIN( beta_max(j) , beta_eq(j) )
          
       END DO

    END SELECT

  END SUBROUTINE f_beta_eq

  !******************************************************************************
  !> \brief Growth rate
  !
  !> This function evaluates the growth rate of crystals given the size. 
  !> \param[in]   i_cry crystal component index 
  !> \param[in]   L_in  crystal size
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION growth_rate(i_cry, L_in)
    
    IMPLICIT NONE

    COMPLEX*8 :: growth_rate

    REAL*8 :: r_growth_rate
    
    INTEGER, INTENT(IN) :: i_cry
    REAL*8, INTENT(IN) :: L_in
    
    r_growth_rate = U_m(i_cry) * ( T_m(i_cry) - T ) * T_u(i_cry) /                &
         ( ( T_m(i_cry) - T_u(i_cry) ) * T ) * DEXP( ( - ( T_u(i_cry) - DREAL(T)) &
         * T_m(i_cry) / ( ( T_m(i_cry) - T_u(i_cry) ) * DREAL(T) ) )  )

    growth_rate = CMPLX( MAX( r_growth_rate, 0.D0), 0.D0 )

  END FUNCTION growth_rate
  
  !******************************************************************************
  !> \brief Nucleation rate
  !
  !> This function evaluates the nucleation rate of crystals given the size. 
  !> \param[in]   i_cry crystal component index 
  !> \param[in]   L_in  crystal size
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION nucleation_rate(i_cry, L_in)
    
    IMPLICIT NONE

    COMPLEX*16 :: nucleation_rate

    REAL*8 :: r_nucleation_rate
    
    INTEGER, INTENT(IN) :: i_cry
    REAL*8, INTENT(IN) :: L_in
  
    r_nucleation_rate = I_m(i_cry) * DEXP( (( T_u(i_cry) / (T_m(i_cry) - T_u(i_cry)) ) &
        * ( ( T_m(i_cry) / T_i(i_cry) ) - ( T_m(i_cry) / DREAL(T) ) ) -             &
        ( ( T_m(i_cry) - T_i(i_cry) ) ** 3.0 ) / (T_m(i_cry) + 3.0 * T_i(i_cry) ) * &
        ( ( T_m(i_cry) / ( T_i(i_cry) * ( T_m(i_cry) - T_i(i_cry) ) ** 2 ) ) -      &
        ( T_m(i_cry) / ( DREAL(T) * ( T_m(i_cry) - DREAL(T) ) ** 2 ) ) ) ) )

    nucleation_rate= DCMPLX( MAX( r_nucleation_rate, 0.D0 ), 0.D0 ) 

  END FUNCTION nucleation_rate
 
  !******************************************************************************
  !> @author 
  !
  !> This subrotine updates kinetic parameters
  !> \date 13/03/12       
  !******************************************************************************
  
  SUBROUTINE update_kinetics ! It must be modified

    IMPLICIT NONE

    INTEGER :: i, j
  
    DO i=1,n_components

       DO j=1,n_cry

          IF(( i == 1 .AND. j == 1 ) ) THEN

             cry_current_solid_solution(i,j) = 1.0
  
          ELSEIF(( i .GT. 1 ) .AND. ( j .GT. 1 )) THEN

             cry_current_solid_solution(i,j) = 1.0 / (n_components - 1.0)

          ELSE
              
             cry_current_solid_solution(i,j) = 0.0

          ENDIF

       ENDDO

    END DO 

    DO i=1,n_cry

       T_m(i) = 1200.00

       T_u(i) = 1150.00

       T_i(i) = 1100.00

    ENDDO

  END SUBROUTINE update_kinetics
 
  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Lithostatic pressure
  !
  !> This subrotine evaluates the lithostatic pressure at depth zeta
  !> \date 13/03/12       
  !******************************************************************************

  SUBROUTINE lithostatic_pressure

    USE geometry, ONLY : zN

    USE complexify 
    IMPLICIT NONE

    p_lith = ( zN - zeta_lith ) * rho_cr * grav

  END SUBROUTINE lithostatic_pressure

  
  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Hydrostatic pressure
  !
  !> This subrotine evaluates the hydrostatic pressure at depth zeta
  !> \date 13/03/12       
  !******************************************************************************

  SUBROUTINE hydrostatic_pressure

    USE geometry

    USE complexify 
    IMPLICIT NONE

    p_hydro = ( zN - zeta_lith ) * rho_w * grav

  END SUBROUTINE hydrostatic_pressure


  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Pressure relaxation term
  !
  !> This subrotine evaluates the term for the relaxation of the pressures.
  !> \date 08/10/2012
  !> \param[out]   pressure_relaxation            relaxation term 
  !******************************************************************************

  SUBROUTINE press_relax_term( pressure_relaxation )

    USE complexify
    IMPLICIT NONE

    COMPLEX*16, INTENT(OUT) :: pressure_relaxation

    CALL mixture_viscosity

    SELECT CASE ( p_relax_model )

    CASE DEFAULT

       tau_p = DCMPLX( 1.D0 , 0.D0 )

    CASE ( 'eval' ) 

       IF ( ( EXPLOSIVE_FLAG ) .AND. ( frag_eff .GT. 0.D0 ) ) THEN

          tau_p = ( visc_mix ** ( 1.D0 - frag_eff ) * visc_2 ** frag_eff )      &
               / ( alfa_1 * alfa_2 * rho_mix )

       ELSE

          tau_p = visc_mix / ( alfa_1 * alfa_2 * rho_mix )

       END IF

    CASE ( 'eval2' ) 

       tau_p = ( 1.D0 - frag_eff ) * visc_mix /                                 &
            ( (alfa_1 ** alfa_1) * rho_mix ) + 0.001 * frag_eff

    CASE ( 'constant' )

       tau_p = DCMPLX( 1.D0 , 0.D0 )

    CASE ( 'single_pressure' )

       tau_p = DCMPLX( 1.D0 , 0.D0 )

    END SELECT

    tau_p = tau_p_coeff * tau_p

    pressure_relaxation = - ( p_2 - p_1 ) / tau_p

  END SUBROUTINE press_relax_term

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Velocity relaxation term
  !
  !> This subrotine evaluates the drag function for the relaxation of the 
  !> velocities
  !> \date 08/10/2012
  !> \param[out]   velocity_relaxation            relaxation term 
  !******************************************************************************

  SUBROUTINE vel_relax_term( velocity_relaxation )

    USE complexify
    IMPLICIT NONE

    COMPLEX*16, INTENT(OUT) :: velocity_relaxation

    COMPLEX*16 :: Rey 
    COMPLEX*16 :: diam
    COMPLEX*16 :: radius_bubble
    COMPLEX*16 :: effusive_drag , explosive_drag

    COMPLEX*16 :: k_1 , k_2

    COMPLEX*16 :: throat_radius
    
    REAL*8 :: frag_transition, t

    frag_transition = frag_thr + 0.05D0

    t = MIN(1.0D0,MAX(0.0D0, ( REAL(alfa_2) - frag_thr ) /                      &
         (frag_transition - frag_thr )))

    SELECT CASE ( drag_funct_model )

    CASE DEFAULT

       effusive_drag = DCMPLX(1.D0,0.D0)

       explosive_drag = effusive_drag

    CASE ( 'eval' )

       ! permeability model

       CALL f_permkc

       ! Please note the bnd is multiplied by alfa_1, because the value is
       ! referred to liquid only
       diam = ( alfa_2 / ( 4.0 / 3.0 * pi *  bubble_number_density              &
            * ( alfa_1 ) ) ) ** ( 1.D0 / 3.D0 )

       Rey = alfa_2 * rho_2 * diam * CDABS( u_2 - u_1 ) / visc_2

       effusive_drag = visc_2 * inv_permkc * ( 1.D0 + 1.D-2 * Rey / alfa_1 )

       explosive_drag = effusive_drag

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'Rey',Rey
          WRITE(*,*) 'visc_2',visc_2
          WRITE(*,*) 'rho_2',rho_2
          WRITE(*,*) 'inv_permkc',inv_permkc

       END IF

    CASE ( 'Klug_and_Cashman' )

       ! permeability model

       CALL f_permkc

       diam = ( alfa_2 / ( 4.0 / 3.0 * pi *  bubble_number_density              &
            * ( alfa_1 ) ) ) ** ( 1.D0 / 3.D0 )

       Rey = alfa_2 * rho_2 * diam * CDABS( u_2 - u_1 ) / visc_2

       effusive_drag = visc_2 * inv_permkc * ( 1.D0 + 1.D-2 * Rey / alfa_1 )

       explosive_drag = 3.D0 * C_D / ( 8.D0 * r_a ) * rho_2 * CDABS( u_2 - u_1 ) 

    CASE ( 'darcy' )

       ! Darcy formulation from Eq. 16 Degruyter et al. 2012 

       radius_bubble = ( alfa_2 / ( 4.0 / 3.0 * pi *  bubble_number_density     &
            * ( alfa_1 ) ) ) ** ( 1.D0 / 3.D0 )

       throat_radius = radius_bubble * throat_bubble_ratio

       k_1 = 0.125D0 * throat_radius ** 2.D0 * alfa_2 ** tortuosity_factor

       effusive_drag = visc_2 / k_1

       explosive_drag = 3.D0 * C_D / ( 8.D0 * r_a ) * rho_2 * CDABS( u_2 - u_1 ) 

    CASE ( 'forchheimer' )

       ! Forchheimer (Eq. 16 Degruyter et al. 2012)

       IF ( REAL(alfa_2) .GT. 0.0001 ) THEN

          radius_bubble = ( alfa_2 / ( 4.0 / 3.0 * pi *  bubble_number_density  &
               * ( alfa_1 ) ) ) ** ( 1.D0 / 3.D0 )

          throat_radius = radius_bubble * throat_bubble_ratio

          k_1 = 0.125D0 * throat_radius ** 2.D0 * alfa_2 ** tortuosity_factor

          k_2 = throat_radius / friction_coefficient * alfa_2 ** ( ( 1.D0 +     &
               3.D0 * tortuosity_factor ) / 2.D0 )

          effusive_drag = visc_2 / k_1 + rho_2 / k_2 * CDABS( u_2 - u_1 )

          explosive_drag = 3.D0 * C_D / ( 8.D0 * r_a ) * rho_2                  &
               * CDABS( u_2 - u_1 ) 

       ELSE

          effusive_drag = DCMPLX(1.0E20,0.0)

          explosive_drag = DCMPLX(0.0,0.0)

       END IF


    CASE ( 'forchheimer_wt' )

       ! Forchheimer (Eq. 16 Degruyter et al. 2012) without transit. zone (phi_t)

       IF ( REAL(alfa_2) .GT. 0.0001 ) THEN

          radius_bubble = ( alfa_2 / ( 4.0 / 3.0 * pi *  bubble_number_density  &
               * ( alfa_1 ) ) ) ** ( 1.D0 / 3.D0 )

          throat_radius = radius_bubble * throat_bubble_ratio

          k_1 = 0.125D0 * throat_radius ** 2.D0 * alfa_2 ** tortuosity_factor

          k_2 = throat_radius / friction_coefficient * alfa_2 ** ( ( 1.D0 +     &
               3.D0 * tortuosity_factor ) / 2.D0 )

          effusive_drag = visc_2 / k_1 + rho_2 / k_2 * CDABS( u_2 - u_1 )

          explosive_drag = 3.D0 * C_D / ( 8.D0 * r_a ) * rho_2                  &
               * CDABS( u_2 - u_1 ) 

       ELSE

          effusive_drag = DCMPLX(1.0E20,0.0)

          explosive_drag = DCMPLX(0.0,0.0)

       END IF

    CASE ('drag')

       radius_bubble = ( alfa_2 / ( 4.0 / 3.0 * pi * ( alfa_1 ) ) )             &
            ** ( 1.D0 / 3.D0 )

       Rey = 2.D0 * radius_bubble * CDABS( u_2 - u_1 ) / visc_1

       C_D = DREAL ( 24.0D0 / Rey * ( 1.0D0 + 1.0D0 / 8.0D0 * Rey**0.72) )

       effusive_drag = 3.D0 * C_D / ( 8.D0 * radius_bubble ) * rho_1            &
            * CDABS( u_2 - u_1 )

       explosive_drag = effusive_drag

    CASE ( 'constant' )

       effusive_drag = DCMPLX(1.D0,0.D0)

       explosive_drag = effusive_drag

    CASE ( 'single_velocity' )

       effusive_drag = DCMPLX(1.D0,0.D0)

       explosive_drag = effusive_drag

    END SELECT

    IF ( drag_funct_model .EQ. 'forchheimer' ) THEN

       drag_funct = effusive_drag ** ( 1.D0 - t ) *                             &
               ( explosive_drag + 1.D-10 ) ** t

    ELSE
    
       drag_funct = effusive_drag ** ( 1.D0 - frag_eff ) *                      &
          ( explosive_drag + 1.D-10 ) ** frag_eff

    END IF

    drag_funct = drag_funct_coeff * drag_funct

    velocity_relaxation = - drag_funct * ( u_1 - u_2 ) * ( rho_mix              &
         / ( rho_1 * rho_2 ) )

    IF ( drag_funct_model .EQ. 'eval' ) THEN

       velocity_relaxation = - drag_funct * ( u_1-u_2 ) / ( x_1 * x_2 * rho_mix )

    ELSE

       velocity_relaxation = - drag_funct * ( u_1 - u_2 ) * ( rho_mix           &
            / ( rho_1 * rho_2 ) )
       
    END IF

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*)  drag_funct , ( u_1 - u_2 ) , x_1 , x_2 , rho_mix 
       WRITE(*,*)  - drag_funct * ( u_1 - u_2 ) / ( x_1 * x_2 * rho_mix ) 
       WRITE(*,*) 'velocity_relaxation',velocity_relaxation

    END IF

  END SUBROUTINE vel_relax_term

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Magma permeability
  !
  !> This subroutine evaluates the permeability by means of the relationship
  !> from Klug and Cashman, 1996:\n
  !> \f$ perm_{kc}= perm_0*10.0^{-10.2(100*\alpha)^{\frac{0.014}{\alpha}}}  \f$\n
  !> for \f$\alpha\geq\overline{alpha}\f$ For \f$\alpha<\overline{alpha}\f$ a
  !> quadratic extrapolation is applied.
  !> \date 02/10/2012       
  !******************************************************************************

  SUBROUTINE f_permkc

    IMPLICIT NONE

    REAL*8 :: permkc_alfa_switch
    REAL*8 :: d_permkc_alfa_switch
    REAL*8 :: log_base

    alfa_switch = 0.03

    log_base = 1.D0 / DLOG10(DEXP(1.D0))

    permkc_alfa_switch = (perm0*10.d0**(-10.2d0*(100.d0*alfa_switch)            &
         **(0.014d0/alfa_switch)) )

    d_permkc_alfa_switch =  log_base * permkc_alfa_switch *                     &
         DLOG10(permkc_alfa_switch / perm0) * ( 0.014D0/alfa_switch**2) *       &
         ( 1 - DLOG(100*alfa_switch) )

    a_2nd = ( permkc_alfa_switch - d_permkc_alfa_switch * alfa_switch)          &
         / ( alfa_switch**2 - 2.D0 * alfa_switch**2.D0 )

    b_2nd = d_permkc_alfa_switch - 2.D0 * a_2nd*alfa_switch

    IF ( REAL(alfa_2) .LT. alfa_switch ) THEN

       permkc = a_2nd * alfa_2*alfa_2 + b_2nd * alfa_2

    ELSE

       permkc = perm0*10.d0**(-10.2d0*(100.d0*alfa_2) ** ( 0.014D0 / alfa_2 ) )

    ENDIF

    inv_permkc = 1.D0 / permkc

  END SUBROUTINE f_permkc


  !******************************************************************************
  !> \brief Mixture viscosity
  !
  !> This subrotine evaluates the viscosity of the mixture (melt+crystals+gas)
  !> below the fragmentation level
  !> \date 11/03/2013  
  !******************************************************************************

  SUBROUTINE mixture_viscosity

    USE complexify 
    IMPLICIT NONE

    CALL f_viscliq

    CALL f_bubbles

    visc_mix = visc_1 * visc_rel_bubbles

  END SUBROUTINE mixture_viscosity

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Magma viscosity
  !
  !> This subroutine evaluates the viscosity of the (bubble-free) magma using
  !> the viscosity of the melt phase  (bubbles and crystal-free magma) and the 
  !> relative viscosity \f$1\theta\f$ due to the presence of the crsytals.
  !> \date 20/02/2009
  !******************************************************************************

  SUBROUTINE f_viscliq

    USE complexify 
    IMPLICIT NONE

    ! Calculate viscosity of the melt (bubbles and crystal-free)
    CALL f_viscmelt

    ! Calculate relative visocisty due to crystals
    CALL f_theta

    ! Calculate viscosity of the melt+crystal
    visc_1 = visc_melt * theta

  END SUBROUTINE f_viscliq

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Melt viscosity
  !
  !> This subrotine evaluates the melt viscosity (bubbles and crystal-free) using 
  !> an empirical relationship of melt viscosity according to water concentration
  !> and temperature of melt (Equation 7 in Hess and Dingwell '96):\n
  !>
  !> where w is the concentration of H2O in weight percent.
  !> \date 13/11/08       
  !******************************************************************************

  SUBROUTINE f_viscmelt

    USE complexify 
    IMPLICIT NONE

    COMPLEX*16 :: w
    COMPLEX*16 :: visc1exp , visc2exp , visc3exp 

    COMPLEX*16, DIMENSION(12) :: wt

    COMPLEX*16, DIMENSION(12) :: norm_wt, xmf
    REAL*8, DIMENSION(12) :: mw
    REAL*8 :: A
    COMPLEX*16 :: gfw , B , C
    REAL*8, DIMENSION(10) :: bb
    COMPLEX*16, DIMENSION(10) :: bcf
    COMPLEX*16, DIMENSION(7) :: cc , ccf
    COMPLEX*16 :: siti , tial , fmm , nak 
    COMPLEX*16 :: b1 , b2 , b3 , b4 , b5 , b6 , b7 , b11 , b12 , b13 , c1 , c2 ,&
         c3 , c4 , c5 , c6 , c11 

    COMPLEX*16 :: x_d_md_tot

    INTEGER :: i

    x_d_md_tot = SUM( x_d_md(1:n_gas) )

    w = x_d_md_tot * 100.d0

    SELECT CASE ( visc_melt_model )

    CASE DEFAULT

    CASE ( 'Hess_and_Dingwell1996')

       !  melt viscosity of rhyolite reaches a maximum at 1.e12

       IF ( REAL(w) .LT. 0.03 ) THEN

          visc_melt = DCMPLX( 1.0D12 , 0.D0 )

       ELSE

          visc1exp = (-3.545D0 + 0.833D0 * CDLOG(w))

          visc2exp = (9601.0D0 - 2368.0D0 * CDLOG(w))

          visc3exp = visc2exp / ( T - ( 195.70D0 + 32.250D0 * CDLOG(w) ) )

          visc_melt = 10.0**(visc1exp+visc3exp)

       END IF

    CASE ( 'Whittington_et_al2008')

       !  melt viscosity of dacites

       IF ( REAL(w) .LT. 0.03 ) THEN

          visc_melt = DCMPLX( 1.0D12 , 0.D0 )

       ELSE

          visc1exp = ( -4.43D0 )

          visc2exp = (7618.13D0 - 17.25D0 * LOG10(w + 0.26))

          visc3exp = visc2exp / ( T - ( 406.1D0 - 292.6D0 * LOG10(w + 0.26) ) )

          visc_melt = 10.0**(visc1exp+visc3exp)

       END IF

    CASE ( 'Romano_et_al2003t')

       !  Campi-Flegrei - Romano et al. (2003) (AMS, trachytes)

       IF ( REAL(w) .LT. 0.03 ) THEN

          visc_melt = DCMPLX( 1.0D12 , 0.D0 )

       ELSE

          visc1exp = ( -3.5405D0 +0.14467D0 * CDLOG(w))

          visc2exp = ( 9618.9D0 - 498.79D0 * CDLOG(w))

          visc3exp = visc2exp / ( T - ( 191.78D0 -35.518D0 * CDLOG(w) ) )

          visc_melt = 10.0**(visc1exp+visc3exp)

       END IF

    CASE ( 'Romano_et_al2003p')

       !  Campi-Flegrei - Romano et al. (2003) (TDPH, phonolites)

       IF ( REAL(w) .LT. 0.03 ) THEN

          visc_melt = DCMPLX( 1.0D12 , 0.D0 )

       ELSE

          visc1exp = ( -5.8996D0 -0.2857D0 * CDLOG(w))

          visc2exp = ( 10775.0D0 - 394.83D0 * CDLOG(w))

          visc3exp = visc2exp / ( T - ( 148.71D0 -21.65D0 * CDLOG(w) ) )

          visc_melt = 10.0**(visc1exp+visc3exp)

       END IF

    CASE ( 'Giordano_et_al2008' )

       ! -------------------------- Giordano et al. (2008) ----------------------
       ! Call for dissolved water in te melt from the model and add it to the   
       ! composition of glass in wt

       DO i = 1,12

          wt(i) = DCMPLX( wt_init(i) , 0.D0 )

       END DO

       wt(11) = w

       ! Create the normalized wt. % distribution

       norm_wt = 100.D0 * ( wt / SUM(wt) )

       ! Change to molar fractions

       mw = [ 60.0843 , 79.8658 , 101.961276 , 71.8444 , 70.937449 , 40.3044 ,  &
            56.0774 , 61.97894 , 94.1960 , 141.9446 , 18.01528 , 18.9984 ]
       gfw = 100.D0 / SUM( norm_wt / mw )
       xmf = ( norm_wt / mw ) * gfw

       ! Model coefficients

       bb  = [159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43,     &
            -0.91, 17.62]
       cc  = [2.75, 15.72, 8.32, 10.2, -12.29, -99.54, 0.3]

       ! Load composition-based matrix for multiplication against 
       ! model-coefficients

       siti = xmf(1) + xmf(2)
       tial = xmf(2) + xmf(3)
       fmm  = xmf(4) + xmf(5) + xmf(6)
       nak  = xmf(8) + xmf(9)
       b1  = siti
       b2  = xmf(3)
       b3  = xmf(4) + xmf(5) + xmf(10)
       b4  = xmf(6)
       b5  = xmf(7)
       b6  = xmf(8) + xmf(11) + xmf(12)
       b7  = xmf(11) + xmf(12) + CDLOG( 1.D0 + xmf(11) )
       b11 = siti * fmm
       b12 = ( siti + xmf(3) + xmf(10) ) * ( nak + xmf(11) )
       b13 = xmf(3) * nak

       c1 = xmf(1)
       c2 = tial
       c3 = fmm
       c4 = xmf(7)
       c5 = nak
       c6 = CDLOG( 1.D0 + xmf(11) + xmf(12) )
       c11 = xmf(3) + fmm + xmf(7) - xmf(10)
       c11 = c11 * ( nak + xmf(11) + xmf(12) )

       bcf =   [b1, b2, b3, b4, b5, b6, b7, b11, b12, b13]
       ccf =   [c1, c2, c3, c4, c5, c6, c11]   

       ! Model main parameters

       A  = -4.55D0
       B  = SUM( bb * bcf )
       C  = SUM( cc * ccf )

       ! Calculates viscosity

       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt

    CASE ('Di_Genova_et_al2013_eqn_3,5')

       ! Viscosity for peralkaline rhyolites from Pantelleria Island 
       ! Di Genova et al. 2013 (Eqs. 3,5)
       
       ! Call for dissolved water in te melt from the model and add it to the   
       ! composition of glass in wt

       DO i = 1,12

          wt(i) = DCMPLX( wt_init(i) , 0.D0 )

       END DO


       wt(11) = w

       ! Create the normalized wt. % distribution

       norm_wt = 100.D0 * ( wt / SUM(wt) )

       ! Change to molar fractions

       mw = [ 60.0843 , 79.8658 , 101.961276 , 71.8444 , 70.937449 , 40.3044 ,  &
            56.0774 , 61.97894 , 94.1960 , 141.9446 , 18.01528 , 18.9984 ]
       gfw = 100.D0 / SUM( norm_wt / mw )
       xmf = ( norm_wt / mw ) * gfw       

       A = -4.55D0

       ! ------------------------------------------
       b1 = 4278.17D0
       b2 = 8.6D0
       B = b1 + b2 * xmf(11)

       c1 = 513.D0
       c2 = -245.3D0
       C = c1 + c2 * LOG10( 1.D0 + xmf(11) )

       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt

    CASE ('Di_Genova_et_al2013_eqn_4,5')

       ! Viscosity for peralkaline rhyolites from Pantelleria Island 
       ! Di Genova et al. 2013 (Eqs. 4,5)

       ! Call for dissolved water in te melt from the model and add it to the   
       ! composition of glass in wt

       DO i = 1,12

          wt(i) = DCMPLX( wt_init(i) , 0.D0 )

       END DO

       wt(11) = w

       ! Create the normalized wt. % distribution

       norm_wt = 100.D0 * ( wt / SUM(wt) )

       ! Change to molar fractions

       mw = [ 60.0843 , 79.8658 , 101.961276 , 71.8444 , 70.937449 , 40.3044 ,  &
            56.0774 , 61.97894 , 94.1960 , 141.9446 , 18.01528 , 18.9984 ]
       gfw = 100.D0 / SUM( norm_wt / mw )
       xmf = ( norm_wt / mw ) * gfw

       A = -4.55D0

       ! -------- New parametrization -------------------------------------------
       b3 = 10528.64D0
       b4 = -4672.21D0
       B = b3 + b4 * LOG10( 1.D0 + xmf(11) )

       c3 = 172.27D0
       c4 = 89.75D0
       C = c3 + c4 * LOG10( 1.D0 + xmf(11) )

       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt


    CASE ( 'Giordano_et_al2009' )

       ! Call for dissolved water in te melt from the model and add it to the   
       ! composition of glass in wt

       DO i = 1,12

          wt(i) = DCMPLX( wt_init(i) , 0.D0 )

       END DO


       wt(11) = w

       ! Create the normalized wt. % distribution

       norm_wt = 100.D0 * ( wt / SUM(wt) )

       ! Change to molar fractions

       mw = [ 60.0843 , 79.8658 , 101.961276 , 71.8444 , 70.937449 , 40.3044 ,  &
            56.0774 , 61.97894 , 94.1960 , 141.9446 , 18.01528 , 18.9984 ]
       gfw = 100.D0 / SUM( norm_wt / mw )
       xmf = ( norm_wt / mw ) * gfw

       A = -4.55D0

       b1 = 6101.0D0
       b2 = -63.66D0

       B = b1 + b2 * xmf(11)

       c1 = 567.0D0
       c2 = -160.3D0
       C = c1 + c2 * LOG10( 1.D0 +xmf(11) )

       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt

    END SELECT

  END SUBROUTINE f_viscmelt

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Crystal relative viscosity
  !
  !> This subrotine evaluates the relative viscosity \f$\theta\f$ due to the 
  !> influence of crystal on the mixture viscosity according to the model 
  !> specified in the input file with the parameter theta_model:\n
  !> - theta_model = 1 => \f$ \theta = \left( 1 - \frac{\beta}{0.7}\right)^{-3,4}\f$ 
  !>                      (Einstein-Roscoe with constants from Lejeune & Richter '95);
  !> - theta_model = 2 => \f$ \theta = \left( 1+0.75\frac{\beta/0.84}{1-\beta/0.84} \right)^2\f$
  !>                        (Dingwell et al '93);
  !> - theta_model = 3 => \f$ \theta = \theta_0 \f$ Melnik & Sparks '99
  !> - theta_model = 4 => \f$ \theta = \left[ 1 - c_1 erf\left(0.5 * \sqrt\pi (1+\frac{c_2}{(1-\beta)^{c_3}})
  !>                          \right)\right]^\frac{-2.5}{c_1}\f$  (Costa)
  !> - theta_model = 5 => \f$ \theta = c_1 \cdot 10 ^{arctan(c_2(\beta-c_3))+\frac{\pi}{2}} \f$
  !>                          (Melnik & Sparks '05).  
  !> .
  !> \date 02/10/2012       
  !******************************************************************************

  SUBROUTINE f_theta

    USE complexify 
    IMPLICIT NONE

    COMPLEX*16 :: arg_erf , t , errorf, var_phi

    REAL*8 :: omega , betas , theta0, beta_aux

    REAL*8 :: p , a1 , a2 , a3 , a4 , a5
    REAL*8 :: phi_star, csi, delta, Einstein_coeff

    SELECT CASE (theta_model)

    CASE DEFAULT

    CASE ( 'Lejeune_and_Richet1995')

       !---------- Einstein-Roscoe ( with constants from Lejeune & Richet '95)
       theta = ( 1.0D0 - ( SUM( beta(1:n_cry) ) / 0.7D0 ) ) ** ( -3.4D0 )

    CASE ('Dingwell1993')

       !---------- Dingwell & al '93
       
       beta_aux = MIN( REAL(SUM(beta(1:n_cry))), 0.8399999 )       
       
       theta = ( 1.0D0 + 0.75D0 * ( beta_aux / 0.84D0 ) / ( 1.D0 - &
            ( beta_aux / 0.84D0 ) ) ) **2.D0

    CASE ('Fixed_value')

       !-----------Constant theta
       theta = theta_fixed    

    CASE ('Costa2005')

       !-----------Costa
       c1 = 0.9995D0
       c2 = 0.4D0
       c3 = 1.D0

       p = 0.3275911D0
       a1 = 0.254829592D0 
       a2 = -0.284496736D0 
       a3 = 1.421413741D0
       a4 = -1.453152027D0
       a5 = 1.061405429D0

       arg_erf = 0.5D0 * DSQRT(pi) * SUM( beta(1:n_cry) ) * ( 1.D0 + c2 /       &
            ( 1.d0 - SUM( beta(1:n_cry) ) ) ** c3 )

       t = 1.D0 / ( 1.D0 + p * arg_erf )

       errorf = 1.D0 - ( a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5)&
            * CDEXP( - arg_erf ** 2 )

       theta = ( 1.D0 - c1 * errorf ) ** ( - 2.5D0 / c1 )

    CASE ('Melnik_and_Sparks2005')

       !-----------Melnik & Sparks 2005 Eq.16
       c2 = 8.6D0
       c3 = 0.69D0
       c1 = 1.D0 / (10.D0 ** ( ATAN( - c2 *  c3  ) + pi/2.D0 ) )

       theta = c1 * 10.D0 ** ( ATAN( c2 * ( SUM( beta(1:n_cry) ) - c3 ) )       &
            + pi/2.D0 ) 

    CASE ('Melnik_and_Sparks1999')

       omega = 20.6d0
       betas = 0.62d0
       theta0 = 1.D0 / (10.D0 ** ( ATAN( - omega *  betas  ) + pi/2.D0 ) )
       
       theta = theta0 * 10.D0 ** ( 0.5D0 * pi + ATAN( omega * ( SUM(            &
            beta(1:n_cry) ) - betas )))

    CASE ('Vona_et_al2011')

       !-----------Costa et al. 2009, Vona et al. 2011

       p = 0.3275911D0
       a1 = 0.254829592D0 
       a2 = -0.284496736D0 
       a3 = 1.421413741D0
       a4 = -1.453152027D0
       a5 = 1.061405429D0

       ! This is the value reported in Vona et al 2011, but according to Fig 9
       ! of that article, a best fitting is obtained with 0.29

       phi_star = 0.274D0
       csi = 0.0327D0
       delta = 13.D0 - 0.84D0
       Einstein_coeff = 2.8D0


       var_phi = SUM( beta(1:n_cry) * rho_c(1:n_cry) / rho_1 ) / phi_star

       arg_erf = 0.5D0 * DSQRT(pi) / (1.0D0 - csi) * var_phi * ( 1.D0 +         &
            var_phi ** 0.84D0 ) 

       t = 1.D0 / ( 1.D0 + p * arg_erf )

       errorf = 1.D0 - ( a1 * t + a2 * t**2.0D0 + a3 * t**3.0D0                 &
            + a4 * t**4.0D0 + a5 * t**5.0D0) * CDEXP( - arg_erf ** 2.0D0 )

       theta =  theta_fixed * (1.D0 + var_phi ** delta ) /                      &
            ( ( 1.D0 - (1.D0 - csi) * errorf ) ** (Einstein_coeff * phi_star) ) 


    CASE ('Vona_et_al2011_mod')

       !-----------Costa et al. 2009, Vona et al. 2011

       p = 0.3275911D0
       a1 = 0.254829592D0 
       a2 = -0.284496736D0 
       a3 = 1.421413741D0
       a4 = -1.453152027D0
       a5 = 1.061405429D0

       ! Modification of Vona et al. 2011 in order to better reproduce
       ! the viscosity at Stromboli.

       phi_star = 0.39D0
       csi = 0.03D0
       delta = 2.0D0 - 0.84D0
       Einstein_coeff = 2.8D0


       var_phi = SUM( beta(1:n_cry) * rho_c(1:n_cry) / rho_1 ) / phi_star

       arg_erf = 0.5D0 * DSQRT(pi) / (1.0D0 - csi) * var_phi * ( 1.D0 +         &
            var_phi ** 0.84D0 ) 

       t = 1.D0 / ( 1.D0 + p * arg_erf )

       errorf = 1.D0 - ( a1 * t + a2 * t**2.0D0 + a3 * t**3.0D0 + a4 * t**4.0D0 &
            + a5 * t**5.0D0) * CDEXP( - arg_erf ** 2.0D0 )

       theta =  theta_fixed * (1.D0 + var_phi ** delta )/( ( 1.D0 - (1.D0 - csi)&
            * errorf ) ** (Einstein_coeff * phi_star) ) 
    
    CASE ('Vona_et_al2013_eq19')

       theta = ( 1.D0 - SUM(beta(1:n_cry)) / (1.D0 - alfa_2 ) ) ** ( - 5.D0     &
            / 2.D0) * ( 1 - alfa_2 ) ** (- 1.D0)

    CASE ('Vona_et_al2013_eq20')

       theta = ( 1.D0 - SUM(beta(1:n_cry)) - alfa_2 )  ** ( - ( 5.D0 *          &
            SUM(beta(1:n_cry)) + 2.D0 * alfa_2 ) / ( 2.D0 * ( SUM(beta(1:n_cry))&
            + alfa_2 ) ) )

    CASE ('Vona_et_al2013_eq21')

       theta = ( 1.D0 - alfa_2 / (1.D0 - SUM(beta(1:n_cry)) ) ) ** ( -1.0 )     &
           * ( 1 - SUM(beta(1:n_cry)) ) ** (- 5.D0/2.D0)

    END SELECT

  END SUBROUTINE f_theta

  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Exsolved gas relative viscosity
  !
  !> This subrotine evaluates the relative viscosity due to the influence of 
  !> exsolved gas on the mixture viscosity according to the model 
  !> specified in the input file with the parameter bubbles_model:\n
  !> - bubbles_model = Costa 2007   => 
  !> - bubbles_model = Einstein     => 
  !> - bubbles_model = Quane-Russel => 
  !> - bubbles_model = Eilers       => 
  !> - bubbles_model = Sibree       =>
  !> .
  !> \date 02/10/2012       
  !******************************************************************************

  SUBROUTINE f_bubbles

    USE complexify 
    USE geometry, ONLY : radius
    IMPLICIT NONE

    REAL*8 :: Ca, gamma_Ca, alfa_aux

    IF(  (.NOT. (theta_model .EQ. 'Vona_et_al2013_eq19' ) ) .AND.               & 
         (.NOT. (theta_model .EQ. 'Vona_et_al2013_eq20' ) ) .AND.               & 
         (.NOT. (theta_model .EQ. 'Vona_et_al2013_eq21' ) ) ) THEN

       SELECT CASE ( bubbles_model )

       CASE DEFAULT

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0)

       CASE ( 'none' )

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0)

       CASE ( 'Costa2007' )

          gamma_Ca = 0.25D0

          Ca =  r_a * visc_melt * (8.00D0 * (u_mix * pi * radius**2)            &
               / (3.0D0 * pi * radius**3.0D0)) / gamma_Ca

          visc_rel_bubbles = (1.0D0/(1.0D0 + 25.0D0*Ca*Ca))                     & 
               * ( 1.0D0 / (1.0D0-alfa_2)                                       &
               + 25.0D0*Ca*Ca*((1.0D0-alfa_2)**(5.0D0/3.0D0)))

       CASE ( 'Quane-Russel' )

          ! For Campi-Flegrei we use 0.63 as reported in 2004 Report (Task2.2)

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * CDEXP( ( - 0.63 * alfa_2 )   &
               / ( 1.D0 - alfa_2 ) )

       CASE ( 'Eilers' )

          ! Eq. (17) Mader et al. 2013

          alfa_aux = MIN( REAL(alfa_2), 0.772818199 )

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * (1.0D0 + (1.25D0 * alfa_aux)   &
               / (1.0D0 - 1.29 * alfa_aux) ) ** 2.0D0

       CASE ( 'Sibree' )

          ! Eq. (18) Mader et al. 2013

          alfa_aux = MIN( REAL(alfa_2), 0.83330833 )

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * ( 1.0D0 / ( 1.0D0            &
               - (1.2 * alfa_aux)** 0.33333D0 ) )

       CASE ( 'Taylor' )

          ! Eq. (16) Mader et al. 2013

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * (1.0D0 + alfa_2)

       CASE ( 'Mackenzie' )

          ! Eq. (19) Mader et al. 2013

          alfa_aux = MIN( REAL(alfa_2), 0.599994 )

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * (1.0D0 - (5.D0 / 3.D0)       &
               * alfa_aux)

       CASE ( 'DucampRaj' )

          ! Eq. (21) Mader et al. 2013, using b = 3

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * CDEXP( -3.D0 * ( alfa_2      &
               / ( 1.D0 - alfa_2 ) ) )

       CASE ( 'BagdassarovDingwell' )

          ! Eq. (22) Mader et al. 2013

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * ( 1.D0 / (1.D0 + 22.4D0      &
               * alfa_2 ) )

       CASE ( 'Rahaman' )

          ! Eq. (20) Mader et al. 2013

          visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * CDEXP( - 11.2D0 * alfa_2 )

       END SELECT

    ELSE

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0)

    END IF

  END SUBROUTINE f_bubbles

  !*****************************************************************************
  !> @author 
  !> Giuseppe La Spina
  !> \brief Bottom exsolved gas
  !
  !> This subrotine evaluates the exsolved gas volumetric fraction given the  
  !> the total gas mass fraction and the dissolved gass mass fraction:\n
  !> \f$ \alpha_2= \frac{x_2(1-\beta)\rho_{md}}{(1-x_{tot})\rho_2
  !>                   + x_2(1-\beta)\rho_{md}}\f$,\n
  !> where \f$ x_2=x_{tot}-x_{d,md} \f$ is the mass fraction of the exsolved gas 
  !> with respect to the crystal-free magma.
  !> \date 13/03/12       
  !> \param[in]    xtot      total gas mass fraction
  !> \param[in]    xmax      dissolved gas mass fraction
  !> \param[in]    r_beta    crystals volume fraction
  !> \param[in]    r_rho_md  dis.gas+melt density
  !> \param[in]    r_rho_2   exsolved gas density
  !> \param[out]   r_alfa_2  exsolved gas volume fraction
  !> \date 02/12/12       
  !******************************************************************************

  SUBROUTINE f_alfa(xtot,xmax,r_beta,r_rho_md,r_rho_2,r_alfa_2)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: xtot(n_gas) , xmax(n_gas) 
    REAL*8, INTENT(IN) :: r_beta(n_cry)
    REAL*8, INTENT(IN) :: r_rho_md
    REAL*8, INTENT(IN) :: r_rho_2 
    REAL*8, INTENT(OUT) :: r_alfa_2(n_gas)


    REAL*8 :: r_x_g(n_gas)    !< exsolved gas mass fraction

    ! Mass fraction of the exsolved gas with respect to the crystal-free magma
    r_x_g(1:n_gas) = xtot(1:n_gas) - xmax(1:n_gas)

    r_alfa_2 = ( r_x_g * ( 1.D0 - SUM( r_beta(1:n_cry) ) ) * r_rho_md ) /       &
         ( ( 1.D0 - xtot) * r_rho_2 + r_x_g * ( 1.D0 - SUM( r_beta(1:n_cry) ) ) &
         * r_rho_md )

  END SUBROUTINE f_alfa

  !*****************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Bottom exsolved gas
  !
  !> This subrotine evaluates the exsolved gas volumetric fraction given the  
  !> the total gas mass fraction and the dissolved gass mass fraction:\n
  !> \f$ \alpha_2= \frac{x_2(1-\beta)\rho_{md}}{(1-x_{tot})\rho_2
  !>                   + x_2(1-\beta)\rho_{md}}\f$,\n
  !> where \f$ x_2=x_{tot}-x_{d,md} \f$ is the mass fraction of the exsolved gas 
  !> with respect to the crystal-free magma.
  !> \date 13/03/12
  !> \param[in]    r_p_2     gas pressure
  !> \param[in]    xtot      total gas mass fraction
  !> \param[in]    r_beta    crystals volume fraction
  !> \param[in]    r_rho_md  dis.gas+melt density
  !> \param[in]    r_rho_g   exsolved gas density
  !> \param[out]   r_alfa_g  exsolved gas volume fraction
  !> \date 02/12/12       
  !******************************************************************************

  SUBROUTINE f_alfa3(r_p_2,xtot,r_beta,r_rho_md,r_rho_g,r_alfa_g)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: xtot(n_gas) 
    REAL*8, INTENT(IN) :: r_beta(n_cry)
    REAL*8, INTENT(IN) :: r_p_2
    REAL*8, INTENT(IN) :: r_rho_md
    REAL*8, INTENT(IN) :: r_rho_g(n_gas) 
    REAL*8, INTENT(OUT) :: r_alfa_g(n_gas)


    !REAL*8 :: r_x_g(n_gas)    !< exsolved gas mass fraction
    REAL*8 :: r_x_d(n_gas)
    REAL*8 :: r_alfa_g_2_max(n_gas),best_alfa_g_2(n_gas)
    REAL*8 :: r_alfa_g_2(n_gas), r_alfa_g_2_new(n_gas)
    REAL*8 :: h,error,best_error
    REAL*8 :: A(n_gas,n_gas),num(n_gas), den(n_gas)
    REAL*8 :: r_rho_md_beta
    INTEGER :: iter,i,j
    INTEGER :: pivot(n_gas)
    INTEGER :: ok

    r_alfa_g = 1.0e-5;

    r_alfa_g_2_max = (xtot / solub)**(1.0 / solub_exp) / r_p_2;


    r_rho_md_beta = ( 1.D0 - SUM( r_beta(1:n_cry) ) ) * r_rho_md

    h = 1e-5;

    iter = 1
    error = 1.0
    r_alfa_g_2 = [0.0 , 1.0]

    best_error = 1.0
    best_alfa_g_2 = r_alfa_g_2

    DO WHILE ( r_alfa_g_2(1) .LT. r_alfa_g_2_max(1) - h)

       r_alfa_g_2 = [r_alfa_g_2(1) + h, 1.0 - (r_alfa_g_2(1) + h)]

       DO i = 1,n_gas

          r_x_d(i) = DREAL( solub(i) * (r_alfa_g_2(i) * p_2)**(solub_exp(i)) )

          DO j = 1,n_gas
             A(i,j) = - xtot(i) * r_rho_g(j) + r_rho_md_beta * &
                  ( xtot(i) - r_x_d(i) ) 
          END DO

          A(i,i) = A(i,i) + r_rho_g(i)

          r_alfa_g(i) = r_rho_md_beta * ( xtot(i) - r_x_d(i) )

       END DO

       call DGESV(n_gas, 1, A , n_gas, pivot, r_alfa_g , n_gas, ok)

       r_alfa_g_2_new = r_alfa_g / SUM( r_alfa_g )

       error = MAXVAL(ABS(r_alfa_g_2_new - r_alfa_g_2))

       DO i = 1,n_gas
          IF (r_alfa_g(i) .LT. 0.0) THEN
             error = 1.0
          END IF
       END DO

       DO i = 1,n_gas
          IF (r_x_d(i) .GT. xtot(i)) THEN
             error = 1.0
          END IF
       END DO

       IF ( error .LT. best_error) THEN
          best_error = error
          best_alfa_g_2 = r_alfa_g_2
       END IF

    END DO

    r_alfa_g_2 = best_alfa_g_2

    DO i = 1,n_gas

       r_x_d(i) = DREAL (solub(i) * (r_alfa_g_2(i) * p_2)**(solub_exp(i)) )

       DO j = 1,n_gas
          A(i,j) = - xtot(i) * r_rho_g(j) + r_rho_md_beta *                     &
               ( xtot(i) - r_x_d(i) ) 
       END DO

       A(i,i) = A(i,i) + r_rho_g(i)

       r_alfa_g(i) = r_rho_md_beta * ( xtot(i) - r_x_d(i) )

    END DO

    call DGESV(n_gas, 1, A , n_gas, pivot, r_alfa_g , n_gas, ok)

    num = r_alfa_g * r_rho_g + ( 1.0 - SUM( r_alfa_g ))* r_rho_md_beta * r_x_d
    den = SUM(r_alfa_g * r_rho_g) + ( 1.0 - SUM( r_alfa_g )) * r_rho_md_beta
    error = MAXVAL(ABS(xtot - num/den));

    IF( error .GT. 1E-010 ) THEN
       WRITE(*,*) 'Initial exsolved volatile content error!'
       WRITE(*,*) 'error='
       WRITE(*,*) error
       STOP
    END IF


  END SUBROUTINE f_alfa3

  !******************************************************************************
  !> \brief Additional moments computation
  !
  !> This subroutine compute the additional moments of crystal components using
  !> the quadrature formulas.
  !> \param[in]   xi     abscissas for the quadrature
  !> \param[out]  w      weights for the quadrature
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_additional_moments

    ! external variables
    USE parameters, ONLY : n_nodes
    
    IMPLICIT NONE

    REAL*8, DIMENSION(n_cry, n_nodes, 2) :: Li
    REAL*8, DIMENSION(n_cry, n_nodes, 2) :: w
    REAL*8, DIMENSION(n_cry,n_nodes) :: growth_rate_array

    INTEGER :: i , j , k

    DO i = 1, n_cry

       DO k = 1,2

          CALL wheeler_algorithm( REAL(mom_cry(i,:,k)), Li(i,:,k), w(i,:,k) )

       END DO

    END DO

    DO i = 1,n_cry

       DO j=1,n_nodes

          growth_rate_array(i,j) = growth_rate(i, Li(i,j,k))

      END DO

    END DO

    DO i=1,n_cry

       DO k = 1,2
       
          DO j=0,n_mom-1

             growth_mom(i,j,k) = SUM( growth_rate_array(i,:)*w(i,:,k) &
                 * Li(i,:,k)**j ) / mom_cry(i,j,k)

          END DO
             
       END DO

    END DO

    RETURN

  END SUBROUTINE eval_additional_moments

  !******************************************************************************
  !> \brief Wheeler algorithm
  !
  !> This subroutine compute quadrature approximation with the Wheeler algorithm.
  !> Given the moments m from 0 to 2*n_nodes-1 it finds the nodes xi and the  
  !> weights w. 
  !> Modified from Marchisio and Fox 2013
  !> \param[in]     m       moments
  !> \param[out]    xi      abscissas
  !> \param[out]    w       weights 
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE wheeler_algorithm(m,xi,w)

    USE parameters, ONLY : n_nodes

    IMPLICIT NONE

    REAL*8, DIMENSION(n_mom), INTENT(IN) :: m
    REAL*8, DIMENSION(n_nodes), INTENT(OUT) :: xi
    REAL*8, DIMENSION(n_nodes), INTENT(OUT) :: w

    REAL*8, DIMENSION(n_nodes,n_nodes) :: jacobi
    REAL*8, DIMENSION(n_nodes) :: D
    REAL*8, DIMENSION(n_nodes-1) :: E

    REAL*8, DIMENSION(n_nodes) :: a , b
    REAL*8, DIMENSION(n_nodes+1, 2*n_nodes) :: sigma_wheeler

    REAL*8, DIMENSION(n_nodes,n_nodes) :: evec

    INTEGER :: i , l , k, consdis

    REAL*8, DIMENSION(2*n_nodes-2) :: WORK
    INTEGER :: INFO
    CHARACTER*1 :: JOBZ
    INTEGER :: LDZ

    consdis = 0

    DO i= 2,n_mom-1

       IF( (m(i+1)/m(i) - m(i)/m(i-1)) .GT. 1E-5) THEN

          consdis = 1 

       ENDIF

    END DO

    IF( consdis .EQ. 0 ) THEN

       xi(1) = m(2)/m(1)
       w(1) = m(1)

       DO i= 2, n_nodes

          xi(i) = 0.0
          w(i) = 0.0

       ENDDO

    ELSE

       DO i=1,n_nodes+1

          DO l=1,2*n_nodes
          
             sigma_wheeler(i,l) = 0.0
   
          END DO
       
       END DO
       
       DO l=0,2*n_nodes-1
          
          sigma_wheeler(2,l+1) = m(l+1)
          
       END DO

       !
       ! compute coefficients for Jacobi matrix 
       !
       
       a(1) = m(2) / m(1)
       b(1) = 0.D0

       DO k=1,n_nodes-1
          
          DO l=k,2*n_nodes-k-1
             
             sigma_wheeler(k+2,l+1) = sigma_wheeler(k+1,l+2) - a(k) * sigma_wheeler(k+1,l+1) - b(k) *      &
                  sigma_wheeler(k,l+1)
             
             a(k+1) = -sigma_wheeler(k+1,k+1) / sigma_wheeler(k+1,k) + sigma_wheeler(k+2,k+2) /            &
                  sigma_wheeler(k+2,k+1)
             
             b(k+1) = sigma_wheeler(k+2,k+1) / sigma_wheeler(k+1,k)
             
          END DO
          
       END DO
  
       !
       ! compute Jacobi matrix
       !
       
       DO i=1,n_nodes
          
          jacobi(i,i) = a(i)
          
          D(i) = jacobi(i,i)
          
       END DO
       
       DO i=1,n_nodes-1
          
          jacobi(i,i+1) = -(ABS(b(i+1)))**0.5
          jacobi(i+1,i) = -(ABS(b(i+1)))**0.5
          
          E(i) = jacobi(i,i+1)
          
       END DO
       
       !
       ! compute eigenvalues and eigenvectors
       !
       
       JOBZ = 'V'    ! compute the eigenvectors
       
       LDZ = n_nodes
       
       CALL DSTEV(JOBZ, n_nodes, D, E, evec, LDZ, WORK, INFO)
       
       !
       ! return weights
       !
       
       DO i=1,n_nodes
          
          w(i) = evec(1,i)**2 * m(1)
          
       END DO
       
       !
       ! return abscissas
       !
       
       DO i=1,n_nodes
          
          xi(i) = D(i)
          
       END DO

    ENDIF

    RETURN

  END SUBROUTINE wheeler_algorithm

END MODULE constitutive
