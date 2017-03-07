!*********************************************************************
!> \mainpage   MAMMA - Magma Ascent Mathematical Modeling and Analysis
!> MAMMA is a FORTRAN90 code designed to solved a conservative model
!> for steady magma ascent in a volcanic conduit, described as a 
!> compressible multicomponent two-phase flow.\n 
!> The governing multiphase equations for two-phase compressible flow 
!> are derived from the theory of thermodynamically compatible systems 
!> (Romenski et al., 2010).\n 
!> The model is one-dimensional (with variable radius) , with 
!> different phase velocities and pressures but a single temperature 
!> for the two phases.\n
!> The solver uses a shooting technique: for a given magma-chamber 
!> pressure, the model adjusts the inlet flow rate until the 
!> prescribed boundary conditions at exit are reached.\n
!> Version 1.0:\n
!> - 1D steady flow for axysimmetric geometry with variable radius
!> .
!> \authors Mattia de' Michieli Vitturi (*), Giuseppe La Spina (**), 
!> Alvaro Aravena Ponce (***)\n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!> (**) School of Earth, Atmospheric and Environmental Sciences \n
!>       The University of Manchester, \n
!>       Oxford Road, Manchester, M13 9PL. \n
!>       E-mail: giuseppe.laspina@manchester.ac.uk \n
!> (***) Dipertimento di Scienze della Terra \n
!>       Universita di Firenze, \n
!>       Florence, Italy
!*********************************************************************

!> \brief Main Program 
PROGRAM MAMMA_main

  USE constitutive, ONLY : initialize_models
  USE inpout, ONLY : init_param , read_param
  USE geometry, ONLY : init_grid
  USE steady_solver, ONLY : steady_shooting

  IMPLICIT NONE

  REAL*8 :: t1 , t2

  CALL cpu_time(t1)

  CALL initialize_models
  
  CALL init_param

  CALL read_param

  CALL init_grid
  
  CALL steady_shooting
  
  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'


END PROGRAM MAMMA_main

