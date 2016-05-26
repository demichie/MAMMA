!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!*********************************************************************
MODULE geometry

  IMPLICIT NONE

  !> Location of the centers of the control volume of the domain
  REAL*8, ALLOCATABLE :: z_comp(:)

  !> Location of the boundaries of the control volumes of the domain
  REAL*8, ALLOCATABLE :: z_stag(:)

  !> Radius at the boundaries of the control volumes of the domain
  REAL*8, ALLOCATABLE :: radius_stag(:)

  REAL*8 :: z0           !< Left (bottom) of the physical domain
  REAL*8 :: zN           !< Right (top) of the physical domain
  REAL*8 :: zeta_exit    !< Right (top) of the physical domain
  REAL*8 :: dz           !< Control volumes size

  REAL*8 :: pi

  REAL*8 :: radius       !< Effective radius

  !> geometry model\n
  !> - 'fixed'           => constant radius
  !> - 'linear'          => linear change in radius
  !> - 'trans1'          => two radiuses and an transition zone
  !> - 'trans2'          => a zone with high radius with a transition zone
  !> - 'external'        => radius profile read from external file
  !> .
  CHARACTER*30 :: radius_model

  REAL*8 :: radius_fixed !< Fixed value of the radius
  REAL*8 :: radius_min   !< Fixed value of the minimum radius (used in non cylindrical conduits)
  REAL*8 :: radius_max   !< Fixed value of the maximum radius (used in non cylindrical conduits)
  REAL*8 :: radius_z     !< Characteristic depth for radius models trans1 and trans2
  REAL*8 :: radius_z_sig !< Characteristic sigma for radius model trans1 and trans2

  INTEGER :: comp_cells  !< Number of control volumes in the computational domain
  INTEGER :: comp_interfaces !< Number of interfaces (comp_cells+1)

CONTAINS

  !*********************************************************************
  !> \brief Finite volume grid initialization
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !*********************************************************************

  SUBROUTINE init_grid

    IMPLICIT none

    INTEGER j      !> loop counter

    pi = 4.D0 * DATAN(1.D0)

    comp_interfaces = comp_cells+1

    ALLOCATE( z_comp(comp_cells) )
    ALLOCATE( z_stag(comp_interfaces) )
    ALLOCATE( radius_stag(comp_interfaces) )
    
    dz = ( zN - z0 ) / comp_cells

    z_stag(1) = z0
    z_comp(1) = z0 + 0.5 * dz
    
    DO j=1,comp_cells

       z_stag(j+1) = MIN( z_stag(j) + dz , zN )

       z_comp(j) = 0.5 * ( z_stag(j) + z_stag(j+1) )

    END DO

    SELECT CASE ( radius_model )
 
    CASE DEFAULT

       radius_stag(1:comp_interfaces) = radius_fixed

      
    CASE ('fixed' )

       radius_stag(1:comp_interfaces) = radius_fixed

    CASE ( 'linear' )

       DO j=1,comp_interfaces
          
          radius_stag(j) = radius_min + (radius_max - radius_min) * z_stag(j) / (zN - z0)
          
       END DO
       
    CASE ( 'trans1' )

       DO j=1,comp_interfaces
          
          IF( z_stag(j) > (radius_z + radius_z_sig) ) THEN
             
             radius_stag(j) = radius_max
             
          ELSEIF( z_stag(j) < (radius_z - radius_z_sig) ) THEN
             
             radius_stag(j) = radius_min
             
          ELSE

             radius_stag(j) = radius_min + (radius_max - radius_min) * &
                  ( z_stag(j) - radius_z + radius_z_sig) / (2.D0 * radius_z_sig)
             
          END IF
          
       END DO
 
    CASE ( 'trans2' )
      
       DO j=1,comp_interfaces
          
          IF( z_stag(j) > (radius_z + radius_z_sig) ) THEN
             
             radius_stag(j) = radius_min
             
          ELSEIF( z_stag(j) > (radius_z) ) THEN
             
             radius_stag(j) = radius_max - (radius_max - radius_min) * &
                  ( z_stag(j) - radius_z ) / ( radius_z_sig )
             
          ELSEIF( z_stag(j) > (radius_z - radius_z_sig) ) THEN
             
             radius_stag(j) = radius_min + (radius_max - radius_min) * &
                  ( z_stag(j) - radius_z + radius_z_sig) / ( radius_z_sig )
             
          ELSE
             
             radius_stag(j) = radius_min
             
          END IF
          
       END DO

    CASE ( 'external' )
       
       OPEN( UNIT=10, FILE='DataRadius.txt' )
       
       DO j=1,comp_interfaces
          
          READ(10,*) radius_stag(j)
          
       END DO
       
       CLOSE(10)
       
    END SELECT

  END SUBROUTINE init_grid

  SUBROUTINE update_radius(zeta)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: zeta
    REAL*8 :: coeff_interp
    INTEGER :: j , z_idx

    z_idx = 1

    DO j=1,comp_interfaces-1

       IF ( z_stag(j) < zeta ) z_idx = j

    END DO    

    ! zeta is between z_stag(j) and z_stag(j+1)
    
    coeff_interp = ( zeta - z_stag(z_idx) ) / ( z_stag(z_idx+1) - z_stag(z_idx) )  

    radius = coeff_interp * radius_stag(z_idx+1) + ( 1.D0 - coeff_interp ) *        &
         radius_stag(z_idx)

  END SUBROUTINE update_radius

END MODULE geometry
