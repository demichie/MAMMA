!********************************************************************************
!> \brief Melts fit module
!
!> This module contains the procedures for reading the fit of melts
!> 
!********************************************************************************

MODULE melts_fit_module

  USE parameters, ONLY : n_cry, n_components

  !> Fraction of composition explained by phases in the fit
  REAL*8 :: wt_tot_0

  !> Normalized fraction of phases in melts fit (i.e. without crystals)
  REAL*8, ALLOCATABLE :: wt_components_fit(:)

  !> Normalized fraction of phases at conduit bottom (with crystals)
  REAL*8, ALLOCATABLE :: wt_components_init(:)

  !> Proportion of oxides in each component
  REAL*8, ALLOCATABLE :: wt_oxide_components(:,:)

  !> Residual values of oxides
  REAL*8,  DIMENSION(12) :: wt_oxide_residual

  !> Type of geochemical system (1: AbAnDi)
  INTEGER :: type_system

  !> Pressure list
  REAL*8, ALLOCATABLE :: pressure_list(:)

  !> Pressure length
  INTEGER :: n_pressure

  !> Fraction list
  REAL*8, ALLOCATABLE :: fraction_list(:)

  !> Pressure length
  INTEGER :: n_fraction

  !> Coefficients of fitting
  REAL*8, ALLOCATABLE :: p00(:,:)
  REAL*8, ALLOCATABLE :: p10(:,:)
  REAL*8, ALLOCATABLE :: p20(:,:)
  REAL*8, ALLOCATABLE :: p01(:,:)
  REAL*8, ALLOCATABLE :: p02(:,:)
  REAL*8, ALLOCATABLE :: p11(:,:)
  REAL*8, ALLOCATABLE :: q00(:,:)
  REAL*8, ALLOCATABLE :: q10(:,:)
  REAL*8, ALLOCATABLE :: q20(:,:)
  REAL*8, ALLOCATABLE :: q01(:,:)
  REAL*8, ALLOCATABLE :: q02(:,:)
  REAL*8, ALLOCATABLE :: q11(:,:)
  REAL*8, ALLOCATABLE :: s00(:,:)
  REAL*8, ALLOCATABLE :: s10(:,:)
  REAL*8, ALLOCATABLE :: s20(:,:)
  REAL*8, ALLOCATABLE :: s01(:,:)
  REAL*8, ALLOCATABLE :: s02(:,:)
  REAL*8, ALLOCATABLE :: s11(:,:)
  REAL*8, ALLOCATABLE :: T_ref(:,:)

  !> Information about the relation between crystals and phases 
  COMPLEX*16, ALLOCATABLE :: rel_cry_components(:,:)	

CONTAINS

  !******************************************************************************
  !> @author 
  !> Alvaro Aravena
  !> \brief Read file of fitting
  !
  !******************************************************************************

  SUBROUTINE read_fit

    INTEGER :: i, j, k, l, bol_lines
  
    LOGICAL :: fexist

    INTEGER, PARAMETER :: fitting_unit = 18

    CHARACTER(LEN=100), DIMENSION(23) :: headings_string

    CHARACTER(LEN=100) :: lines

    REAL*8, DIMENSION(100) :: lecture_array

    REAL*8, DIMENSION(100,100) :: lecture_matrix

    CHARACTER(LEN=50) :: input_fitting      !< Name of the file with the fitting coefficients

    input_fitting = 'fitting_coefficient.inp'

    headings_string(1) = 'pressure_list'
    headings_string(2) = 'fraction_list'
    headings_string(3) = 'p00'
    headings_string(4) = 'p01'
    headings_string(5) = 'p02'
    headings_string(6) = 'p10'
    headings_string(7) = 'p20'
    headings_string(8) = 'p11'
    headings_string(9) = 'q00'
    headings_string(10) = 'q01'
    headings_string(11) = 'q02'
    headings_string(12) = 'q10'
    headings_string(13) = 'q20'
    headings_string(14) = 'q11'
    headings_string(15) = 's00'
    headings_string(16) = 's01'
    headings_string(17) = 's02'
    headings_string(18) = 's10'
    headings_string(19) = 's20'
    headings_string(20) = 's11'
    headings_string(21) = 'T_ref'
    headings_string(22) = 'type'

    INQUIRE (FILE=input_fitting,exist=fexist)
    
    IF ( .NOT. fexist ) THEN

       WRITE(*,*) 'Fitting file not found'
       STOP

    ENDIF

    DO i=1,2

       OPEN(fitting_unit,FILE=input_fitting,STATUS='old')

       DO j=1,100

          lecture_array(j) = 0.D0

       ENDDO

       DO j=1,500

          READ(fitting_unit, * ,IOSTAT=bol_lines), lines
 
          IF(bol_lines .EQ. 0) THEN

             IF(TRIM(lines) .EQ.  TRIM(headings_string(i))) THEN

                READ(fitting_unit, * , IOSTAT=bol_lines) (lecture_array(l), l=1,100)

                IF(i == 1) THEN

                   n_pressure = MINLOC(lecture_array, DIM=1) - 1

                   ALLOCATE(  pressure_list(n_pressure) )

                   pressure_list(1:n_pressure) = lecture_array(1:n_pressure)

                ELSE

                   n_fraction = MINLOC(lecture_array, DIM=1) - 1

                   ALLOCATE(fraction_list(n_fraction))

                   fraction_list(1:n_fraction) = lecture_array(1:n_fraction)

                ENDIF

                EXIT

             ENDIF

          ENDIF

       ENDDO

       CLOSE(fitting_unit)

    ENDDO

    ALLOCATE( p00(n_pressure, n_fraction), p01(n_pressure, n_fraction), p02(n_pressure, n_fraction))
    ALLOCATE( p11(n_pressure, n_fraction), p10(n_pressure, n_fraction), p20(n_pressure, n_fraction))
    ALLOCATE( q00(n_pressure, n_fraction), q01(n_pressure, n_fraction), q02(n_pressure, n_fraction))
    ALLOCATE( q11(n_pressure, n_fraction), q10(n_pressure, n_fraction), q20(n_pressure, n_fraction))
    ALLOCATE( s00(n_pressure, n_fraction), s01(n_pressure, n_fraction), s02(n_pressure, n_fraction))
    ALLOCATE( s11(n_pressure, n_fraction), s10(n_pressure, n_fraction), s20(n_pressure, n_fraction))
    ALLOCATE( T_ref(n_pressure, n_fraction) )

    DO i = 3,21

       OPEN(fitting_unit,FILE=input_fitting,STATUS='old')

       DO j=1,100

          lecture_array(j) = 0.D0

       ENDDO

       DO j=1,500

          READ(fitting_unit, * ,IOSTAT=bol_lines), lines
 
          IF(bol_lines .EQ. 0) THEN

             IF(TRIM(lines) .EQ.  TRIM(headings_string(i))) THEN

                DO k = 1,n_pressure

                   READ(fitting_unit, * , IOSTAT=bol_lines) (lecture_array(l), l=1,n_fraction)

                   lecture_matrix(k,1:n_fraction) = lecture_array(1:n_fraction)

                ENDDO

                IF( i .EQ. 3 ) THEN

                   p00 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 4 ) THEN

                   p01 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 5 ) THEN

                   p02 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 6 ) THEN

                   p10 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 7 ) THEN

                   p20 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 8 ) THEN

                   p11 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 9 ) THEN

                   q00 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 10 ) THEN

                   q01 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 11 ) THEN

                   q02 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 12 ) THEN

                   q10 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 13 ) THEN

                   q20 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 14 ) THEN

                   q11 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 15 ) THEN

                   s00 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 16 ) THEN

                   s01 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 17 ) THEN

                   s02 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 18 ) THEN

                   s10 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 19 ) THEN

                   s20 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 20 ) THEN

                   s11 = lecture_matrix(1:n_pressure,1:n_fraction)

                ELSEIF( i .EQ. 21 ) THEN

                   T_ref = lecture_matrix(1:n_pressure,1:n_fraction)

                ENDIF

                EXIT

             ENDIF

          ENDIF

       ENDDO

       CLOSE(fitting_unit)

    ENDDO

    i = 22

    OPEN(fitting_unit,FILE=input_fitting,STATUS='old')

    DO j=1,100

       lecture_array(j) = 0.D0

    ENDDO

    DO j=1,500

       READ(fitting_unit, * ,IOSTAT=bol_lines), lines
 
       IF(bol_lines .EQ. 0) THEN

          IF(TRIM(lines) .EQ.  TRIM(headings_string(i))) THEN

             READ(fitting_unit, * , IOSTAT=bol_lines), type_system

             EXIT

          ENDIF

       ENDIF

    ENDDO

    CLOSE(fitting_unit)

    IF(type_system .EQ. 1) THEN

       wt_oxide_components(1:12, 1) = [68.74, 0.0, 19.44, 0.0, 0.0, 0.0, 0.0, 11.82, 0.0, 0.0, 0.0, 0.0] / 100.0
       wt_oxide_components(1:12, 2) = [43.19, 0.0, 36.65, 0.0, 0.0, 0.0, 20.16, 0.0, 0.0, 0.0, 0.0, 0.0] / 100.0
       wt_oxide_components(1:12, 3) = [55.49, 0.0, 0.0, 0.0, 0.0, 18.61, 25.90, 0.0, 0.0, 0.0, 0.0, 0.0] / 100.0

       IF(( n_cry .NE. 2 ) .OR. (n_components .NE. 3 )) THEN

          WRITE(*,*) 'Number of components/crystals is not compatible with system AbAnDi.'
          STOP

       ENDIF

    ELSE

          WRITE(*,*) 'System type is not supported.'
          STOP       

    ENDIF

    DO i=1,n_components
    
       DO j=1,n_cry

          IF(type_system == 1) THEN

             IF(( i == 3 .AND. j == 2 ) .OR. (i .LT. 3 .AND. j .LT. 2)) THEN

                rel_cry_components(i,j) = DCMPLX( 1.D0, 0.D0 )	
 
             ELSE
              
                rel_cry_components(i,j) = DCMPLX( 0.D0, 0.D0 )	  

             ENDIF

          ENDIF

       ENDDO
 
    ENDDO

    RETURN

  END SUBROUTINE read_fit

END MODULE melts_fit_module
