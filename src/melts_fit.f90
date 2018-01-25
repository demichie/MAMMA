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


CONTAINS

  !******************************************************************************
  !> \brief Read file of melts fit (per ora, ci sono valori fissi)
  !
  !******************************************************************************

  SUBROUTINE read_fit

    INTEGER :: i
  
    wt_tot_0 = 0.9

    DO i=1,n_components

       wt_components_fit(i) = 1.0 / n_components

    END DO 

    RETURN

  END SUBROUTINE read_fit

END MODULE melts_fit_module
