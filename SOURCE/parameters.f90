!*********************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!*********************************************************************
MODULE parameters

  IMPLICIT NONE

  REAL*8 :: pi

  REAL*8 :: eps_newton        !< threshold for the convergence of the
                              !< Newton's method 
  REAL*8 :: max_dt            !< Largest time step allowed
  REAL*8 :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL*8 :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

  !> Flag to add the relaxation terms after the linear reconstruction:\n
  !> - T      => evaluate the relaxation terms
  !> - F      => reconstruction without the relaxation 
  !> .
  LOGICAL :: interfaces_relaxation

  INTEGER :: n_cry         !< Numbeer of crystal phases
  INTEGER :: n_gas         !< Numbeer of crystal phases

  INTEGER :: n_vars        !< Number of conservative variables
  INTEGER :: n_eqns        !< Number of equations


  INTEGER :: n_nh

  INTEGER :: n_RK     !< Runge-Kutta order

  !> Limiter for the slope in the linear reconstruction:\n
  !> - 'none'     => no limiter (constant value)
  !> - 'minmod'   => minmod sloe;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  CHARACTER(LEN=20) :: limiter

  !> Finite volume method:\n
  !> - 'LxF'       => lax-friedrichs scheme;
  !> - 'GFORCE '   => gforce scheme;
  !> - 'KT'        => Kurganov and Tadmor semidiscrete scheme;
  !> .
  CHARACTER(LEN=20) :: solver_scheme     

  REAL*8 :: theta             !< Van Leer limiter parameter
  REAL*8 :: t_start           !< initial time for the run
  REAL*8 :: t_end             !< end time for the run
  REAL*8 :: t_output          !< time of the next output
  REAL*8 :: dt_output         !< time interval for the output of the solution

  !> Parameter for numerical scheme:\n
  !> - alfa_impl = 1.0   => Euler Implicit
  !> - alfa_impl = 0.5   => Crank-Nicolson
  REAL*8, PARAMETER :: alfa_impl = 1.0D0

  !> Maximum iterations of the Newthon-Raphson solver
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL*8, PARAMETER :: eps_thr = 1.D-2

  REAL*8, PARAMETER :: tol_abs = 1.D-3
  REAL*8, PARAMETER :: tol_rel = 1.D-3

  INTEGER :: verbose_level

  !> Flag for the shooting technique:\n
  !> - T   => the code search for inlet velociy
  !> - F   => the velocity is provided as input
  !> .
  LOGICAL :: shooting

  !> Residual for the convergence of the shooting method. The solution is 
  !> accepted if one of these conditions is satisfied:
  !> - ( P_exit - P_out ) / P_out < eps_conv
  !> - ( Mach_exit - 1 ) < eps_conv
  !> - ( Mass_flow_rate_exit - Mass_flow_rate ) < eps_conv
  !> .
  REAL*8 :: eps_conv

  !> Flag for the preconditioning
  !> - T      => evaluate the preconditioning
  !> - F      => not evaluate
  !> .
  LOGICAL :: preconditioning

  !> Flag for the preconditioning on conservative or entropic variables
  !> - T      => entropic variables
  !> - F      => conservative variables
  !> .
  LOGICAL :: entropic_prec

  !> Flag for dual time stepping
  !> - T      => use the dual-time stepping
  !> - F      => not use the dual-time stepping
  !> .
  LOGICAL :: dual_time

  INTEGER :: max_iter_tau

  REAL*8 :: dtau_dt_ratio

  !> Flag for local time stepping
  !> - T      => use local time step
  !> - F      => use global time step
  !> .
  LOGICAL :: local_time_step

  LOGICAL :: residual_smoothing

  LOGICAL :: complex_step_speed

  LOGICAL :: exact_speed


END MODULE parameters
