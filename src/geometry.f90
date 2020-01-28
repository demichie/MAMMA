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

  !> Eccentricity at the boundaries of the control volumes of the domain
  REAL*8, ALLOCATABLE :: eccen_stag(:)

  REAL*8 :: z0           !< Left (bottom) of the physical domain
  REAL*8 :: zN           !< Right (top) of the physical domain
  REAL*8 :: zeta_exit    !< Right (top) of the physical domain
  REAL*8 :: dz           !< Control volumes size

  REAL*8 :: pi

  REAL*8 :: radius       !< Effective radius
  REAL*8 :: f_eccen_a    !< Eccentricity factor a (rate between ellipse perimeter and 2*pi*sqrt(Ra*Rb))
  REAL*8 :: f_eccen_b    !< Eccentricity factor b (sqrt(2*Ra*Ra*Rb*Rb/(Ra*Ra+Rb*Rb))/sqrt(Ra*Rb))
  REAL*8 :: d_qradius 	!< Derivative of square radius along zeta

  !> geometry model\n
  !> - 'fixed'           => constant equivalent radius, constant eccentricity
  !> - 'linear'          => linear change in equivalent radius, constant eccentricity
  !> - 'trans1'          => two equivalent radii and a transition zone, constant eccentricity
  !> - 'trans2'          => a zone with high equivalent radius, constant eccentricity 
  !> - 'trans3'          => constant equivalent radius in lower portion and linearly variable upper portion, constant eccentricity
  !> - 'external'        => equivalent radius profile read from external file, constant eccentricity
  !> - 'linear_eccen'	 => constant equivalent radius, linear change in eccentricity
  !> - 'trans1_eccen'    => constant equivalent radius, two eccentricities and a transitional zone
  !> - 'trans4_eccen'    => constant minimum semi-axis, two eccentricities and a transitional zone
  !> .
  CHARACTER*30 :: radius_model

  REAL*8 :: radius_fixed !< Fixed value of equivalent radius (used in fixed, linear_eccen and trans1_eccen)
  REAL*8 :: radius_min   !< Fixed value of the minimum equivalent radius
  REAL*8 :: radius_max   !< Fixed value of the maximum equivalent radius
  REAL*8 :: radius_z     !< Characteristic depth for models trans1, trans2 and trans3
  REAL*8 :: radius_z_sig !< Characteristic sigma for models trans1 and trans2
  REAL*8 :: eccen_fixed  !< Fixed eccentricity of the conduit
  REAL*8 :: eccen_base   !< Value of the base eccentricity for models trans1_eccen and trans4_eccen
  REAL*8 :: eccen_top    !< Value of the top eccentricity for models trans1_eccen and trans4_eccen
  REAL*8 :: eccen_z_base !< Characteristic depth for models trans1_eccen and trans4_eccen
  REAL*8 :: eccen_z_top  !< Characteristic depth for models trans1_eccen and trans4_eccen
  REAL*8 :: eccen_axis_b !< Characteristic depth for model trans4_eccen
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
    ALLOCATE( eccen_stag(comp_interfaces) )

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

       eccen_stag(1:comp_interfaces) = eccen_fixed

    CASE ('fixed' )

       radius_stag(1:comp_interfaces) = radius_fixed
       
       eccen_stag(1:comp_interfaces) = eccen_fixed

    CASE ( 'linear' )

       DO j=1,comp_interfaces
          
          radius_stag(j) = radius_min + (radius_max - radius_min) * z_stag(j) / (zN - z0)
          
       END DO

       eccen_stag(1:comp_interfaces) = eccen_fixed

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
 
       eccen_stag(1:comp_interfaces) = eccen_fixed
 
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

       eccen_stag(1:comp_interfaces) = eccen_fixed

    CASE ( 'trans3' )

       DO j=1,comp_interfaces
          
          IF( z_stag(j) < (radius_z) ) THEN
             
             radius_stag(j) = radius_min
             
          ELSE

             radius_stag(j) = radius_min + (radius_max - radius_min) * &
                  ( z_stag(j) - radius_z ) / (zN - radius_z)
             
          END IF
          
       END DO
       
       eccen_stag(1:comp_interfaces) = eccen_fixed       

    CASE ( 'external' )
       
       OPEN( UNIT=10, FILE='DataRadius.txt' )
       
       DO j=1,comp_interfaces
          
          READ(10,*) radius_stag(j)
          
       END DO
       
       CLOSE(10)

       eccen_stag(1:comp_interfaces) = eccen_fixed

    CASE ( 'linear_eccen' )

       radius_stag(1:comp_interfaces) = radius_fixed

       DO j=1,comp_interfaces
          
          eccen_stag(j) = eccen_base + (eccen_top - eccen_base) * z_stag(j) / (zN - z0)
          
       END DO

    CASE ( 'trans1_eccen' )

       radius_stag(1:comp_interfaces) = radius_fixed

       DO j=1,comp_interfaces
          
          IF( z_stag(j) > (eccen_z_top) ) THEN
             
             eccen_stag(j) = eccen_top
             
          ELSEIF( z_stag(j) < (eccen_z_base) ) THEN
             
             eccen_stag(j) = eccen_base
             
          ELSE

             eccen_stag(j) = eccen_base + (eccen_top - eccen_base) * 		&
                  ( z_stag(j) - eccen_z_base ) / (eccen_z_top - eccen_z_base)
             
          END IF
          
       END DO

    CASE ( 'trans4_eccen' )

       DO j=1,comp_interfaces
          
          IF( z_stag(j) > (eccen_z_top) ) THEN
             
             eccen_stag(j) = eccen_top
             
          ELSEIF( z_stag(j) < (eccen_z_base) ) THEN
             
             eccen_stag(j) = eccen_base
             
          ELSE

             eccen_stag(j) = eccen_base + (eccen_top - eccen_base) * &
                  ( z_stag(j) - eccen_z_base ) / (eccen_z_top - eccen_z_base)
             
          END IF

       	  radius_stag(j) = eccen_axis_b / ( ( 1 - eccen_stag(j) ** 2.D0 ) ** 0.25D0)
          
       END DO

    END SELECT

  END SUBROUTINE init_grid

  SUBROUTINE update_radius(zeta)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: zeta
    REAL*8 :: coeff_interp, eccen_aux1, eccen_aux2
    INTEGER :: j , z_idx

    z_idx = 1

    DO j=1,comp_interfaces-1

       IF ( z_stag(j) < zeta ) z_idx = j

    END DO    

    ! zeta is between z_stag(j) and z_stag(j+1)
    
    coeff_interp = ( zeta - z_stag(z_idx) ) / ( z_stag(z_idx+1) - z_stag(z_idx) )  

    radius = coeff_interp * radius_stag(z_idx+1) + ( 1.D0 - coeff_interp ) *        &
         radius_stag(z_idx)

    d_qradius = 2.D0 * radius * ( radius_stag(z_idx + 1) - radius_stag(z_idx) ) / ( z_stag(z_idx+1) - z_stag(z_idx) )  

    eccen_aux1 = coeff_interp * eccen_stag( z_idx + 1 ) + ( 1.D0 - coeff_interp ) *          &
         eccen_stag(z_idx)

    eccen_aux2 = (1.D0 - eccen_aux1 ** 2.D0) ** (0.5D0)

    f_eccen_a =  ( 3.D0 * (1.D0 + eccen_aux2 ) - ((3.D0 + eccen_aux2) *  &
		              (1.D0 + 3.D0*eccen_aux2 ))**(0.5D0) ) / (2.D0 * eccen_aux2**(0.5D0) )

    f_eccen_b =  SQRT( (2.D0 * eccen_aux2 ) / (1 + eccen_aux2 ** 2.D0) )

  END SUBROUTINE update_radius

END MODULE geometry
