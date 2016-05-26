!*********************************************************************
!> \mainpage   MAMMA - Finite volume solver
!> MAMMA is a FORTRAN90 code designed to solved a conservative model
!> for compressible two-phase flow by finite volume methods. The model
!> is one-dimensional with different phase velocities and pressures but 
!> a single temperature for the two phases.\n
!> The finite volume solver is based on a semidiscrete central scheme 
!> and it is not tied on the specific eigenstructure of the model.\n
!> Version 0.9:\n
!> - 1D transient flow
!> .
!> \authors Mattia de' Michieli Vitturi (*,**), Giuseppe La Spina (***)\n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: demichie@pi.ingv.it \n
!> (**) School of Earth and Space Exploration, Arizona State University \n
!>      Tempe, AZ, USA \n
!>      E-mail: mdemichi@asu.edu \n
!> (***) Dipartimento di Matematica "L.Tonelli", Universita' di Pisa \n
!>       Largo Bruno Pontecorvo, 5\n
!>       I-56127 Pisa, Italy\n
!>       E-mail: glaspina@gmail.com
!*********************************************************************

!> \brief Main Program 
PROGRAM MAMMA_main

  USE inpout, ONLY : init_param , read_param
  USE geometry, ONLY : init_grid
  USE steady_solver, ONLY : steady_shooting

  IMPLICIT NONE

  REAL*8 :: t1 , t2

  CALL cpu_time(t1)

  CALL init_param

  CALL read_param

  CALL init_grid
  
  CALL steady_shooting
  
  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'


END PROGRAM MAMMA_main

