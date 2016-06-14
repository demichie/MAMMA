!*******************************************************************************
!> \brief Constitutive equations
!*******************************************************************************
MODULE constitutive

  USE geometry, ONLY : pi
  USE parameters, ONLY : verbose_level
  USE parameters, ONLY : n_eqns , n_vars
  USE parameters, ONLY : n_cry , n_gas
  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  !--------- Constants for the equations of state -----------------------------
  REAL*8 :: cv_2         !< exsolved gas specific heat capacity at constant volume
  REAL*8 :: cv_m         !< melt specific heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_c(:)  !< crystals specific heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_d(:)  !< dissolved gas heat capacity at constant volume
  REAL*8, ALLOCATABLE :: cv_g(:)  !< exsolved gas heat capacity at constant volume

  !~   REAL*8 :: C0_2        !< exsolved gas sound speed at atmospheric conditions
  REAL*8 :: C0_m             !< melt sound speed at atmospheric conditions
  REAL*8, ALLOCATABLE :: C0_c(:)      !< crystals sound speed at atmospheric conditions
  REAL*8, ALLOCATABLE :: C0_d(:)      !< dissolved gas sound speed at atmospheric conditions
  REAL*8, ALLOCATABLE :: C0_g(:)      !< exsolved gas sound speed at atmospheric conditions

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

  REAL*8, ALLOCATABLE :: fit(:,:)


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

  !> bubble number density
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
  !> - 'alphaMelts'     => Equilibrium from fitting of alphaMelts results;
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
  !> - 'Zhang' => Zhang model (Zhang, JVGR '99).
  !> .
  CHARACTER*20 :: exsol_model

  !> Solubility parameter for the Henry's law
  REAL*8, ALLOCATABLE :: solub(:)

  REAL*8, ALLOCATABLE :: solub_exp(:)

  REAL*8, ALLOCATABLE :: x_ex_dis_in(:)

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

  !> gravitational acceleration
  REAL*8 :: grav

  !> country rock permeability
  REAL*8 :: k_cr

  !> contry rock density
  REAL*8 :: rho_cr

  !> gas viscosity
  REAL*8 :: visc_2

  !> melt viscosity
  COMPLEX*16 :: visc_melt

  !> melt+crystal viscosity
  COMPLEX*16 :: visc_1

  !> relative viscosity due to bubbles
  COMPLEX*16 :: visc_rel_bubbles

  !> Parameter to choose the model for the influence of crystal on the mixture 
  !> viscosity according:\n
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

  !> Parameter to choose the model for the influence of the bubbles on the mixture:\n 
  !> - 'Einstein' 
  !> - 'Quane-Russel  -> For Campi-Flegrei we use 0.63 as reported in the 2004 Report (Task 2.2)
  !> - 'Eilers'       -> Mader et al. 2013
  !> - 'Sibree'       -> Mader et al. 2013
  !> .
  CHARACTER*20 :: bubbles_model

  !> Flag to choose the eruptive style:\n
  !> - explosive = .TRUE.   => explosive eruption
  !> - explosive = .FALSE.  => effusive eruption
  LOGICAL :: explosive

  !> Magma permeability
  COMPLEX*16 :: permkc 
  COMPLEX*16 :: inv_permkc  

  REAL*8 :: alfa_switch
  REAL*8 :: a_2nd , b_2nd , perm0

  !> Parameter to select the melt viscosity (bubbles and crystal-free) model:\n
  !> - 'Hess_and_Dingwell1996'
  !> - 'Romano_et_al2003'
  !> - 'Giordano_et_al2008'
  !> - 'Giordano_et_al2009'
  !> - 'Di_Genova_et_al2013_eqn_3,5'
  !> - 'Di_Genova_et_al2013_eqn_4,5'
  !> .
  CHARACTER*30 :: visc_melt_model

  
  REAL*8, DIMENSION(12) :: wt_init


  !> Flag to activate the injection of external water:\n
  !> - ext_water = .TRUE.   => interaction with external water
  !> - ext_water = .FALSE.  => no interaction
  LOGICAL :: ext_water

  !> Flag for instantaneous vaporization:\n
  !> - inst_vaporization = .TRUE.   => instantaneous vaporization
  !> - inst_vaporization = .FALSE.  => injection as liquid
  LOGICAL :: inst_vaporization

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

  !> Water density
  REAL*8 :: rho_w

  REAL*8 :: xa,xb,xc

  ! number of coefficients for the alphaMelts fitting
  INTEGER :: n_coeffs

CONTAINS

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

    IF ( present(c_qp) ) THEN

       qp = c_qp

    ELSE

       qp = DCMPLX( r_qp , 0.D0 )

    END IF

    alfa_g(1:n_gas) = qp(1:n_gas)

    alfa_2 = SUM( alfa_g(1:n_gas) )

    alfa_g_2(1:n_gas) = alfa_g(1:n_gas) / alfa_2

    alfa_1 = 1.D0 - alfa_2
    p_1 = qp(n_gas+1)
    p_2 = qp(n_gas+2)
    u_1 = qp(n_gas+3)
    u_2 = qp(n_gas+4)
    T = qp(n_gas+5)
    beta(1:n_cry) = qp(n_gas+5+1:n_gas+5+n_cry)
    x_d_md(1:n_gas) = qp(n_gas+5+n_cry+1:n_gas+5+n_cry+n_gas)

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

    !WRITE(*,*) rho_1
    !WRITE(*,*) rho_2
    !WRITE(*,*) alfa_1
    !WRITE(*,*) alfa_2
    !WRITE(*,*) u_1
    !WRITE(*,*) u_2
    !WRITE(*,*) rho_mix
    !WRITE(*,*) u_mix
    !READ(*,*)

    rhoB_m = rho_mix * x_m
    rhoB_c(1:n_cry) = rho_mix * x_c(1:n_cry)

    cv_mix = x_m * cv_m + SUM( x_c(1:n_cry) * cv_c(1:n_cry) )                   &
         + SUM( x_d(1:n_gas) * cv_d(1:n_gas) )                                  &
         + SUM( x_g(1:n_gas) * cv_g(1:n_gas) )

  END SUBROUTINE phys_var_qp

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

    !mu_1 = e_1 + p_1/rho_1 - T * s_1
    !mu_2 = e_2 + p_2/rho_2 - T * s_2

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

    INTEGER :: idx 
    INTEGER :: i

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

    idx = 0

    flux(1:n_eqns) = DCMPLX(0.D0,0.D0)

    !---- Mixture Density -------------------------------------------------------
    idx = idx + 1

    flux(idx) = rho_mix * u_mix * radius**2

    !WRITE(*,*) 'Mixture density'
    !WRITE(*,*) flux(idx)

    !---- Volumetric Fraction First Phase ---------------------------------------
    idx = idx + 1

    IF ( p_relax_model .EQ. 'single' ) THEN

       flux(idx) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       flux(idx) = rho_mix * u_mix * alfa_1 * radius**2

    END IF

    !WRITE(*,*) 'Volume Fraction1'
    !WRITE(*,*) flux(idx)


    !---- Mass Fraction Exsolved Gas Phases -------------------------------------
    DO i=1,n_gas

       idx = idx + 1

       flux(idx) = alfa_g(i) * rho_g(i) * u_2 * radius**2

       !WRITE(*,*) 'Mass fraction exsolved gas'
       !WRITE(*,*) flux(idx)

    END DO

    !---- Mixture Momentum ------------------------------------------------------
    idx = idx + 1

    flux(idx) = ( alfa_1 * rho_1 * u_1**2 + alfa_2 * rho_2 * u_2 ** 2 +         &
         alfa_1 * p_1 + alfa_2 * p_2 ) * radius**2

    !WRITE(*,*) 'Mixture momentum'
    !WRITE(*,*) flux(idx)

    !---- Relative Velocity -----------------------------------------------------
    idx = idx + 1

    IF ( drag_funct_model .EQ. 'single_velocity' ) THEN

       flux(idx) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       flux(idx) = ( 0.5D0 * ( u_1 * u_1 - u_2 * u_2 )                          &
            + ( mu_1 - mu_2 ) ) * radius**2

    END IF

    !WRITE(*,*) u_1
    !WRITE(*,*) u_2
    !WRITE(*,*) e_1
    !WRITE(*,*) e_2
    !WRITE(*,*) p_1
    !WRITE(*,*) p_2
    !WRITE(*,*) rho_1
    !WRITE(*,*) rho_2
    !WRITE(*,*) s_1
    !WRITE(*,*) s_2
    !WRITE(*,*) T
    !READ(*,*)

    !WRITE(*,*) 'Relative velocity'
    !WRITE(*,*) flux(idx)

    !---- Total Energy ----------------------------------------------------------
    idx = idx + 1

    IF ( isothermal ) THEN

       flux(idx) = DCMPLX( 0.D0 , 0.D0 )

    ELSE

       flux(idx) = ( alfa_1 * rho_1 * u_1 * ( e_1 + p_1/rho_1                   &
            + 0.5D0 * u_1 * u_1 ) + alfa_2 * rho_2 * u_2 * ( e_2 + p_2/rho_2    &
            + 0.5D0 * u_2* u_2 ) - rho_mix * x_1 * x_2 * ( u_1 - u_2 )          &
            * ( s_1 - s_2 ) * T ) * radius**2
    END IF

    !WRITE(*,*) 'Mixture Energy'
    !WRITE(*,*) flux(idx)


    !----- Crystal Phases -------------------------------------------------------
    DO i=1,n_cry

       idx = idx + 1

       flux(idx) = ( alfa_1 * rho_c(i) * beta(i) * u_1 ) * radius**2

       !WRITE(*,*) 'Crystal volume fraction'
       !WRITE(*,*) flux(idx)

    END DO

    !----- Dissolved Gas Phases -------------------------------------------------
    DO i=1,n_gas

       idx = idx + 1

       flux(idx) = ( x_d_md(i) * alfa_1 *  ( rho_1 - SUM( beta(1:n_cry) *       &
            rho_c(1:n_cry) ) ) * u_1 ) * radius**2

       !WRITE(*,*) 'Dissolved gas mass fraction'
       !WRITE(*,*) flux(idx)

    END DO

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

    !--------------- Evaluate the relaxation terms --------------
    CALL eval_relaxation_terms( relaxation_term )

    !--------------- Evaluate the forces terms ------------------
    CALL eval_forces_terms( forces_term )

    !--------------- Evaluate the forces terms ------------------
    CALL eval_source_terms( source_term )

    IF ( verbose_level .GE. 3 ) THEN

       WRITE(*,*)'relaxation_term'
       WRITE(*,*)relaxation_term / ( radius )**2
    
       WRITE(*,*)'forces_term'
       WRITE(*,*)forces_term / ( radius )**2
    
       WRITE(*,*)'source_term'
       WRITE(*,*)source_term / ( radius )**2
    
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

    INTEGER :: i
    INTEGER :: idx

    idx = 0

    relaxation_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! relaxation term for the mixture density idx = 1
    idx = idx+1 

    relaxation_term(idx) = DCMPLX(0.D0,0.D0)


    ! relaxation term for the volume fraction equation idx = 2
    idx = idx+1 

    CALL press_relax_term( pressure_relaxation )

    relaxation_term(idx) = pressure_relaxation * radius**2


    ! relaxation term for exsolved gas idx = 2+1,2+n_gas

    CALL f_xdis_eq

    DO i=1,n_gas

       idx = idx + 1

       IF ( REAL( x_d_md(i) ) .LT. REAL( x_d_md_eq(i) ) ) THEN

          relaxation_term(idx) = DCMPLX(0.D0,0.D0)

       ELSE

          relaxation_term(idx) = ( x_d_md(i) - x_d_md_eq(i) ) * ( 1.D0 - alfa_2)&
               * ( rho_1 - SUM( beta(1:n_cry) * rho_c(1:n_cry) ) ) / tau_d(i)   &
               * radius **2

       END IF

    END DO

    ! relaxation term for the mixture momentum idx = 2+n_gas+1
    idx = idx+1 

    relaxation_term(idx) = DCMPLX(0.D0,0.D0)


    ! relaxation term for relative velocity idx = 2+n_gas+2
    idx = idx + 1

    CALL vel_relax_term( velocity_relaxation )

    relaxation_term(idx) = velocity_relaxation * radius**2

    ! relaxation term for the mixture energy idx =2+n_gas+3
    idx = idx+1 

    relaxation_term(idx) = DCMPLX(0.D0,0.D0)

    ! relaxation term for crystallization idx = 2+n_gas+3+1,2+n_gas+3+n_cry

    CALL f_beta_eq

    DO i=1,n_cry

       idx = idx + 1

       relaxation_term(idx) = - ( 1.d0 - alfa_2 ) * rho_c(i) * ( beta(i) -      &
            beta_eq(i) ) / tau_c(i) * radius ** 2

    END DO

    ! relaxation term for dissolved gas 
    ! idx = 2+n_gas+3+n_cry+1,2+n_gas+3+n_cry+n_gas

    CALL f_xdis_eq

    DO i=1,n_gas

       idx = idx + 1

       IF ( REAL( x_d_md(i) ) .LT. REAL( x_d_md_eq(i) ) ) THEN

          relaxation_term(idx) = DCMPLX(0.D0,0.D0)

       ELSE

          relaxation_term(idx) = - ( x_d_md(i) - x_d_md_eq(i) ) *               &
               ( 1.D0 - alfa_2 ) * ( rho_1 - SUM( beta(1:n_cry) * rho_c(1:n_cry)&
               ) ) / tau_d(i) * radius ** 2

       END IF

    END DO

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

    USE geometry, ONLY : radius
    USE complexify 
    IMPLICIT none

    COMPLEX*16, INTENT(OUT) :: force_term(n_eqns)

    COMPLEX*16 :: visc_force_1 , visc_force_2
    COMPLEX*16 :: visc_force_1_rel , visc_force_2_rel

    REAL*8 :: gas_wall_drag

    INTEGER :: i

    INTEGER :: idx

    REAL*8 :: frag_transition, t

    frag_transition = frag_thr + 0.05D0

    t = MIN(1.0D0,MAX(0.0D0, ( REAL(alfa_2) - frag_transition ) /               &
         (frag_thr - frag_transition)))

    idx = 0

    force_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! Term for the mixture density ----------------------------------------------

    idx = idx + 1

    force_term(idx) = DCMPLX(0.D0,0.D0)

    ! Term for the first phase volume fraction ----------------------------------

    idx = idx + 1

    force_term(idx) = DCMPLX(0.D0,0.D0)

    ! Terms for the exsolved gas phases -----------------------------------------

    DO i=1,n_gas

       idx = idx + 1

       force_term(idx) = DCMPLX(0.D0,0.D0)

    END DO

    ! Term for the mixture momentum equation ------------------------------------

    idx = idx + 1

    force_term(idx) = DCMPLX(0.0,0.0)

    force_term(idx) = - rho_mix * grav * radius**2.D0

    CALL mixture_viscosity

    visc_force_1 = - 8.D0 * visc_mix * u_1

    ! Turbulent gas-wall friction (Degruyter et al. 2012)

    gas_wall_drag = 0.03D0

    visc_force_2 = - gas_wall_drag / 4.D0 * radius * rho_2 * CDABS( u_2 ) * u_2 

    visc_force_1 = visc_force_1 * ( 1.D0 - frag_eff )
    visc_force_2 = visc_force_2 * frag_eff

    !visc_force_1 =  visc_force_1 * t
    !visc_force_2 =  visc_force_2 * ( 1.0 - t )


    force_term(idx) = force_term(idx) + visc_force_2 + visc_force_1 

    ! Term for the relative velocity equation ----------------------------------

    idx = idx + 1

    force_term(idx) = DCMPLX(0.0,0.0)

    IF ( drag_funct_model .EQ. 'single_velocity' ) THEN

       force_term(idx) = DCMPLX( 0.D0 , 0.D0 ) 

    ELSE

       visc_force_1_rel = visc_force_1 / ( ( 1.D0 - alfa_2 ) * rho_1 ) 

       visc_force_2_rel = - visc_force_2 / ( alfa_2 * rho_2 )

       force_term(idx) = force_term(idx) + visc_force_2_rel + visc_force_1_rel

    END IF

    ! Term for the mixture energy ----------------------------------------------

    idx = idx + 1

    IF ( isothermal ) THEN

       force_term(idx) = DCMPLX(0.d0,0.D0)

    ELSE

       force_term(idx) = DCMPLX(0.d0,0.D0)

       force_term(idx) = - rho_mix * u_mix * grav * radius**2

       force_term(idx) = force_term(idx) + 2.D0 * visc_force_2 * u_2            &
            + 2.D0 * visc_force_1 * u_1

    END IF

    ! Terms for the crystal phases ----------------------------------------------

    DO i= 1,n_cry

       idx = idx + 1

       force_term(idx) = DCMPLX(0.d0,0.D0)

    END DO

    ! Terms for the dissolved gas phases ----------------------------------------

    DO i = 1,n_gas

       idx = idx + 1

       force_term(idx) = DCMPLX(0.d0,0.D0)

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

    USE geometry, ONLY : radius
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
          
          q_lat = ( rho_2 * alfa_2 * k_cr * ( p_2 - p_lith ) ) /       &
               ( visc_2 * radius )
          
       ELSE
          
          q_lat = DCMPLX(0.D0,0.D0)
          
       END IF
       
    END IF


    water_mass_flux = DCMPLX(0.D0,0.D0)

    source_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    ! --- TOTAL MASS SOURCE TERM ------------------------------------------------

    idx = idx + 1

    IF ( lateral_degassing ) THEN

       source_term(idx) = - 2.D0 * q_lat * radius

    ELSE

       source_term(idx) = DCMPLX(0.D0,0.D0)

    END IF

    IF ( ext_water ) THEN

       IF ( ( zeta_lith .GT. min_z_influx ) .AND.                               &
            ( zeta_lith .LT. min_z_influx + delta_z_influx ) ) THEN

          IF ( total_water_influx .GT. 0.D0 ) THEN

             water_mass_flux = DCMPLX( total_water_influx/delta_z_influx , 0.D0)

          ELSE

             rho_w = 1000;

             CALL hydrostatic_pressure

             IF ( p_hydro .GE. REAL(p_1) ) THEN

                visc_w = 2.414D-5 * 10.D0 ** ( 247.8D0 / ( T_w - 140.D0 ) )

                water_mass_flux = ( 2.D0 * radius * 3.14D0) * rho_w * k_cr /    &
                     visc_w * ( p_hydro - p_1 ) / radius
                     
             ELSE

                water_mass_flux = DCMPLX( 0.D0 , 0.D0 )

             END IF

          END IF

       ELSE

          water_mass_flux = DCMPLX(0.D0,0.D0)

       END IF

       source_term(idx) = source_term(idx) + water_mass_flux

    END IF

    source_term(idx) = DCMPLX(0.D0,0.D0)

    ! --- FIRST PHASE VOLUME FRACTION -------------------------------------------
    idx = idx + 1

    source_term(idx) = DCMPLX(0.D0,0.D0)

    ! --- EXSOLVED GAS MASS FRACTIONS -------------------------------------------

    DO i = 1,n_gas

       idx = idx + 1

       IF ( lateral_degassing ) THEN

          source_term(idx) = - 2.D0 * q_lat * radius

       ELSE

          source_term(idx) = DCMPLX(0.D0,0.D0)

       END IF

    END DO

    ! --- Mixture Momentum ------------------------------------------------------
    idx = idx + 1

    IF ( lateral_degassing ) THEN

       source_term(idx) = - 2.D0 * q_lat * radius * u_2

    ELSE

       source_term(idx) = DCMPLX(0.D0,0.D0)

    END IF

    ! --- Relative Velocity -----------------------------------------------------
    idx = idx + 1

    source_term(idx) = DCMPLX(0.D0,0.D0)

    ! --- MIXTURE ENERGY --------------------------------------------------------
    idx = idx + 1

    IF ( isothermal ) THEN

       source_term(idx) = DCMPLX(0.D0,0.D0)

    ELSE

       IF ( lateral_degassing ) THEN

          source_term(idx) = - 2.D0 * q_lat * radius * alfa_2 * ( T * cv_2 +    &
               0.5D0 * u_2*u_2 )

       ELSE

          source_term(idx) = DCMPLX(0.D0,0.D0)

       END IF

       IF ( ext_water ) THEN

          IF ( ( zeta_lith .GT. min_z_influx ) .AND.                            &
               ( zeta_lith .LT. min_z_influx + delta_z_influx ) ) THEN

             IF ( inst_vaporization ) THEN

                heat_flux = - water_mass_flux * ( cv_d(1) * ( T_boiling - T_w ) &
                     + cv_2 * ( - T_boiling ) + lambda_w )

             ELSE

                heat_flux = water_mass_flux * cv_d(1) * T_w

             END IF

             source_term(idx) = source_term(idx) + heat_flux

          END IF

       END IF

    END IF

    ! --- CRYSTALS BULK DENSITY -------------------------------------------------

    DO i = 1,n_cry

       idx = idx + 1

       source_term(idx) = DCMPLX(0.D0,0.D0)

    END DO

    ! --- DISSOLVED GAS BULK DENSITY --------------------------------------------

    idx = idx + 1

    ! --- DISSOLVED WATER -------------------------------------------------------
    IF ( ext_water .AND. .NOT.inst_vaporization ) THEN

       source_term(idx) = water_mass_flux

    ELSE

       source_term(idx) = DCMPLX(0.D0,0.D0)

    END IF

    DO i = 2,n_gas

       idx = idx + 1

       source_term(idx) = DCMPLX(0.D0,0.D0)

    END DO


  END SUBROUTINE eval_source_terms


  !******************************************************************************
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \brief Explicit Forces terms
  !
  !> This subroutine evaluates the forces to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity)
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
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


    ! --- TOTAL MASS SOURCE TERM

    idx = idx + 1

    expl_forces_term(idx) = 0.D0

    ! --- FIRST PHASE VOLUME FRACTION
    idx = idx + 1

    expl_forces_term(idx) = 0.D0

    ! --- EXSOLVED GAS MASS FRACTIONS

    DO i = 1,n_gas

       idx = idx + 1

       expl_forces_term(idx) = 0.D0

    END DO

    ! --- Mixture Momentum
    idx = idx + 1

    expl_forces_term(idx) = DREAL( rho_mix * grav * radius ** 2 )

    ! --- Relative Velocity
    idx = idx + 1

    expl_forces_term(idx) = 0.D0

    ! --- MIXTURE ENERGY
    idx = idx + 1

    IF ( isothermal ) THEN

       expl_forces_term(idx) = 0.D0

    ELSE

       expl_forces_term(idx) = DREAL( rho_mix * u_mix * grav * radius **2 )

    END IF


    ! --- CRYSTALS BULK DENSITY

    DO i = 1,n_cry

       idx = idx + 1

       expl_forces_term(idx) = 0.D0

    END DO

    ! --- DISSOLVED GAS BULK DENSITY

    idx = idx + 1

    DO i = 1,n_gas

       idx = idx + 1

       expl_forces_term(idx) = 0.D0

    END DO

  END SUBROUTINE eval_explicit_forces

  !********************************************************************************
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
  !********************************************************************************

  SUBROUTINE f_xdis_eq

    USE complexify 

    IMPLICIT NONE


    COMPLEX*16 :: aa , bb , cc , pp
    COMPLEX*16 :: P_bar, log10Ds, Ds, Dcl
    REAL*8 :: beta_S(4)
    REAL*8 :: melt_volume

    COMPLEX*16 :: melt_mass
    COMPLEX*16 :: co2_exs_mass,h2o_exs_mass,S_exs_mass,Cl_exs_mass
    COMPLEX*16 :: S_dis_mass,Cl_dis_mass
    COMPLEX*16 :: gas_mass, Ratio_exs_dis
    COMPLEX*16 :: S_mass,Cl_mass


    IF ( REAL(p_2) .LE. 0.D0 ) THEN

       x_d_md_eq = DCMPLX(0.D0, 0.D0)

    ELSE

       SELECT CASE ( exsol_model )

       CASE DEFAULT


       CASE ( 'Henry' )

          ! Henry's law

          x_d_md_eq(1:n_gas) = solub(1:n_gas) * (alfa_g_2(1:n_gas) * p_2 )      &
               ** solub_exp(1:n_gas)

       CASE ( 'Zhang' )

          ! Zhang

          aa = 0.4874d0 - 608.0d0 / T + 489530.d0 / ( T * T )

          bb = -0.06062d0 + 135.6d0 / T - 69200.d0 / ( T * T )

          cc = 0.00253d0 - 4.154d0 / T + 1509.0d0 / ( T * T )

          pp = p_2 * 1.0D-6

          x_d_md_eq(1:n_gas) = 1.0D-2 * ( aa * CDSQRT(pp) + bb * pp + cc        &
               * pp ** 1.5D0 )

       CASE ( 'H20-CO2-S-CL' )

          ! Henry's law (s is global variable and is a constant for rhyolite)
          
          melt_volume=1.0 !%m3
          
          melt_mass=melt_volume*rho_md
          
          !water exsolved mass
          h2o_exs_mass=melt_mass*x_g(1)
          
          !co2 exsolved mass
          co2_exs_mass=melt_mass*x_g(2)
          
          !S exsolved mass
          S_exs_mass=melt_mass*x_g(3)
          
          !Cl exsolved mass
          Cl_exs_mass=melt_mass*x_g(4)
          
          !Total volatile mass
          gas_mass=co2_exs_mass+h2o_exs_mass+S_exs_mass+Cl_exs_mass 
          
          !Ratio_exs_dis is ratio of all gas mass to melt mass
          Ratio_exs_dis=gas_mass/melt_mass 
          
          !Pressure in bar
          P_bar = p_2 / 1e5;
          
          !fitting coefficient
          beta_S = [4.8095e+02, -4.5850e+02, 3.0240e-03, -2.7519e-07] 
          
          !Solubility for S and Cl            
          log10Ds = beta(1) * P_bar**(-0.0075) + beta(2) + beta(3) * P_bar      &
               + beta(4) * P_bar**(2.0)
          Ds = 10**(log10Ds) 
          Dcl = DCMPLX(2.D0, 0.D0)
          
          !Equilibrium dissolved water    
          x_d_md_eq(1) = solub(1) * (alfa_g_2(1) * p_2) ** solub_exp(1)
          !Equilibrium dissolved co2      
          x_d_md_eq(2) = solub(2) * (alfa_g_2(2) * p_2) ** solub_exp(2)

          !Total mass of S
          S_mass=melt_mass*x_ex_dis_in(3)
          
          !Dissolved mass of S 
          S_dis_mass=S_mass - (Ds*Ratio_exs_dis/(Ds*Ratio_exs_dis + 1.0))*S_mass
          
          !Equilibrium dissolved S
          x_d_md_eq(3)=S_dis_mass/melt_mass 
          
          !Total mass of S
          Cl_mass=melt_mass*x_ex_dis_in(4)
          
          !Dissolved mass of Cl 
          Cl_dis_mass=Cl_mass- (Dcl*Ratio_exs_dis/(Dcl* Ratio_exs_dis + 1.0))   &
               *Cl_mass
          
          !Equilibrium dissolved Cl
          x_d_md_eq(4)=Cl_dis_mass/melt_mass 
          
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
    COMPLEX*16 :: crystal_mass_fraction(1:n_cry)

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
       
    CASE ( 'alphaMelts') 

       crystal_mass_fraction(1:n_cry) = ( fit(1,1:n_cry) * p_1_bar * p_1_bar    &
            + fit(2,1:n_cry) * T_celsius * T_celsius                            &
            + fit(3,1:n_cry) * x_d_md_wt_tot * x_d_md_wt_tot                    &
            + fit(4,1:n_cry) * p_1_bar * T_celsius                              &
            + fit(5,1:n_cry) * T_celsius * x_d_md_wt_tot                        &
            + fit(6,1:n_cry) * x_d_md_wt_tot * p_1_bar                          &
            + fit(7,1:n_cry) * p_1_bar                                          &
            + fit(8,1:n_cry) * T_celsius                                        &
            + fit(9,1:n_cry) * x_d_md_wt_tot                                    &
            + fit(10,1:n_cry) ) / 100.D0
       
       DO j=1,n_cry
          
          beta_eq(j) = crystal_mass_fraction(j) * rho_1 / rho_c(j) 
          
          beta_eq(j) = MAX( beta0(j) + 1D-15, beta_eq(j) )
          beta_eq(j) = MIN( beta_max(j) , beta_eq(j) )

       END DO


    END SELECT

  END SUBROUTINE f_beta_eq


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

    p_hydro = (  min_z_influx + delta_z_influx - zeta_lith ) * rho_w * grav

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

    SELECT CASE (  p_relax_model )

    CASE DEFAULT

       tau_p = DCMPLX( 1.D0 , 0.D0 )

    CASE ( 'eval' ) 

       IF ( ( EXPLOSIVE ) .AND. ( frag_eff .GT. 0.D0 ) ) THEN

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

    !REAL*8 :: frag_transition, t

    SELECT CASE ( drag_funct_model )

    CASE DEFAULT

       effusive_drag = DCMPLX(1.D0,0.D0)

       explosive_drag = effusive_drag

    CASE ( 'eval' )

       ! permeability model

       CALL f_permkc

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

    drag_funct = effusive_drag ** ( 1.D0 - frag_eff ) *                         &
         explosive_drag  ** frag_eff

    drag_funct = drag_funct_coeff * drag_funct

    velocity_relaxation = - drag_funct * ( u_1 - u_2 ) * ( rho_mix              &
         / ( rho_1 * rho_2 ) )

    IF ( drag_funct_model .EQ. eval ) THEN

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

    CASE ( 'Romano_et_al2003')

       !  Campi-Flegrei - Romano et al. (2003)

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

       A = -4.55D0

       ! -------- Giordano et al. 2009 ------------------------------------------
       b1 = 4278.17D0
       b2 = 8.6D0
       B = b1 + b2 * w

       c1 = 513.D0
       c2 = -245.3D0
       C = c1 + c2 * CDLOG( 1.D0 + w )

       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt

    CASE ('Di_Genova_et_al2013_eqn_4,5')

       ! Viscosity for peralkaline rhyolites from Pantelleria Island 
       ! Di Genova et al. 2013 (Eqs. 4,5)

       A = -4.55D0

       ! -------- New parametrization -------------------------------------------
       b3 = 10528.64D0
       b4 = -4672.21D0
       B = b3 + b4 * CDLOG( 1.D0 + w )

       c3 = 172.27D0
       c4 = 89.75D0
       C = c3 + c4 * CDLOG( 1.D0 + w )


       visc_melt = A + B / ( T - C )
       visc_melt = 10.D0 ** visc_melt


    CASE ( 'Giordano_et_al2009' )

       A = -4.55D0

       b1 = 6101.0D0
       b2 = -63.66D0

       B = b1 + b2 * 3.5D0 * w
       !B = b1 + b2 * w


       c1 = 567.0D0
       c2 = -160.3D0
       C = c1 + c2 * CDLOG( 1.D0 + 3.5D0 * w  )
       !C = c1 + c2 * CDLOG( 1.D0 + w  )

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

    REAL*8 :: omega , betas , theta0

    REAL*8 :: p , a1 , a2 , a3 , a4 , a5
    REAL*8 :: phi_star, csi, delta, Einstein_coeff

    SELECT CASE (theta_model)

    CASE DEFAULT

    CASE ( 'Lejeune_and_Richet1995')

       !---------- Einstein-Roscoe ( with constants from Lejeune & Richet '95)
       theta = ( 1.0D0 - ( SUM( beta(1:n_cry) ) / 0.7D0 ) ) ** ( -3.4D0 )

    CASE ('Dingwell1993')

       !---------- Dingwell & al '93
       theta = ( 1.0D0 + 0.75D0 * ( ( SUM(beta(1:n_cry) ) / 0.84D0 ) / ( 1.D0 -     &
            ( SUM( beta(1:n_cry) ) / .84D0 ) ) ) ) **2.D0

    CASE ('Fixed_value')

       !-----------Constant theta
       theta = theta_fixed    ! Melnik & Sparks '99

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
       c1 = 1.4D0
       c2 = 8.6D0
       c3 = 0.69D0

       theta = c1 * 10.D0 ** ( ATAN( c2 * ( SUM( beta(1:n_cry) ) - c3 ) ) + pi ) 

    CASE ('Melnik_and_Sparks1999')

       omega = 20.6d0
       betas = 0.62d0
       theta0 = 3.5d0

       theta = theta0 * 10.D0 ** ( 0.5D0 * pi + ATAN( omega * ( SUM(           &
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

       arg_erf = 0.5D0 * DSQRT(pi) / (1.0D0 - csi) * var_phi * ( 1.D0 +        &
            var_phi ** 0.84D0 ) 

       t = 1.D0 / ( 1.D0 + p * arg_erf )

       errorf = 1.D0 - ( a1 * t + a2 * t**2.0D0 + a3 * t**3.0D0 + a4 * t**4.0D0 + a5 * t**5.0D0)&
            * CDEXP( - arg_erf ** 2.0D0 )

       theta =  theta_fixed * (1.D0 + var_phi ** delta )/( ( 1.D0 - (1.D0 - csi) * errorf ) ** &
            (Einstein_coeff * phi_star) ) 


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

    END SELECT

  END SUBROUTINE f_theta

  !*****************************************************************************
  !>
  !>
  !*****************************************************************************

  SUBROUTINE f_bubbles

    USE complexify 
    IMPLICIT NONE

    SELECT CASE ( bubbles_model )

    CASE DEFAULT

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0)

    CASE ( 'none' )

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0)

    CASE ( 'Einstein' )

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) / ( 1.d0 - alfa_2 )

    CASE ( 'Quane-Russel' )

       ! For Campi-Flegrei we use 0.63 as reported in the 2004 Report (Task 2.2)

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * CDEXP( ( - 0.63 * alfa_2 )      &
            / ( 1.D0 - alfa_2 ) )

    CASE ( 'Eilers' )

       ! Eq. (17) Mader et al. 2013

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * (1.0D0 + (1.25D0 * alfa_2)      &
            / (1.0D0 - 1.29 * alfa_2) ) ** 2.0D0

    CASE ( 'Sibree' )

       ! Eq. (18) Mader et al. 2013

       visc_rel_bubbles = DCMPLX(1.0D0,0.0D0) * (1.0D0 / (1.0D0 - (1.2 * alfa_2)&
            ** 0.33333D0 ))

    END SELECT


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

END MODULE constitutive

