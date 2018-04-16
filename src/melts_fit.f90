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

  !> Information about the relation between crystals and phases 
  COMPLEX*16, ALLOCATABLE :: rel_cry_components(:,:)	

CONTAINS

  !******************************************************************************
  !> \brief Read file of melts fit (per ora, ci sono valori fissi)
  !
  !******************************************************************************

  SUBROUTINE read_fit

    INTEGER :: i, j
  
    wt_tot_0 = 0.9

    DO i=1,n_components

       wt_components_fit(i) =  1.D0 / n_components 
    
       DO j=1,n_cry

          IF(( i == 1 .AND. j == 1 ) .OR. (i .GT. 1 .AND. j .GT. 1)) THEN

             rel_cry_components(i,j) = DCMPLX( 1.D0, 0.D0 )	
 
          ELSE
              
             rel_cry_components(i,j) = DCMPLX( 0.D0, 0.D0 )	  

          ENDIF

       ENDDO
 
    ENDDO

    RETURN

  END SUBROUTINE read_fit

END MODULE melts_fit_module
