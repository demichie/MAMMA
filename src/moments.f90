!********************************************************************************
!> \brief Method of Moments module
!
!> This module contains the procedures and the variables for the method of 
!> moments.
!> \date 22/10/2013
!> @author 
!> Mattia de' Michieli Vitturi
!********************************************************************************

MODULE moments_module

  IMPLICIT NONE

  !> Number of nodes for the quadrature
  INTEGER :: n_nodes

  !> Maximum number of moments
  INTEGER, PARAMETER :: max_nmom = 20

  !> Number of moments
  INTEGER :: n_mom

CONTAINS

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

  SUBROUTINE wheeler_algorithm(m,distribution,xi,w)

    IMPLICIT NONE

    REAL*8, DIMENSION(n_mom), INTENT(IN) :: m
    CHARACTER(LEN=20), INTENT(IN) :: distribution
    REAL*8, DIMENSION(n_nodes), INTENT(OUT) :: xi
    REAL*8, DIMENSION(n_nodes), INTENT(OUT) :: w


    REAL*8, DIMENSION(n_nodes,n_nodes) :: jacobi
    REAL*8, DIMENSION(n_nodes) :: D
    REAL*8, DIMENSION(n_nodes-1) :: E

    REAL*8, DIMENSION(n_nodes) :: a , b
    REAL*8, DIMENSION(n_nodes+1, 2*n_nodes) :: sigma

    REAL*8, DIMENSION(n_nodes,n_nodes) :: evec

    INTEGER :: i , l , k

    REAL*8, DIMENSION(2*n_nodes-2) :: WORK
    INTEGER :: INFO
    CHARACTER*1 :: JOBZ
    INTEGER :: LDZ

    IF ( distribution .EQ. 'constant' ) THEN

       xi(1) = m(2)/m(1)
       w(1) = m(1)

    ELSE

       DO l=1,2*n_nodes-2
          
          sigma(1,l+1) = 0.0
          
       END DO
       
       DO l=0,2*n_nodes-1
          
          sigma(2,l+1) = m(l+1)
          
       END DO
       
       !
       ! compute coefficients for Jacobi matrix 
       !
       
       
       a(1) = m(2) / m(1)
       b(1) = 0.D0
       
       
       DO k=1,n_nodes-1
          
          DO l=k,2*n_nodes-k-1
             
             sigma(k+2,l+1) = sigma(k+1,l+2) - a(k) * sigma(k+1,l+1) - b(k) *      &
                  sigma(k,l+1)
             
             a(k+1) = -sigma(k+1,k+1) / sigma(k+1,k) + sigma(k+2,k+2) /            &
                  sigma(k+2,k+1)
             
             b(k+1) = sigma(k+2,k+1) / sigma(k+1,k)
             
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
       
    END IF

    RETURN

  END SUBROUTINE wheeler_algorithm

  !******************************************************************************
  !> \brief Moments correction
  !
  !> This subroutine compute uses the algorithm of McGraw for the correction of 
  !> an unrealizable moment set. The subroutine first analyzes the moment set and
  !> then, if it is unrealizable, identifies the moemnt that has to be changed 
  !> the least to make the set realizable. 
  !> Modified from Marchisio and Fox 2013
  !> \param[in,out]    m       moments
  !> \param[out]       iter    number of iterations required
  !> \author Mattia de' Michieli Vitturi
  !> \date 22/10/2013
  !******************************************************************************

  SUBROUTINE moments_correction(m,iter)

    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(INOUT) :: m
    INTEGER, INTENT(OUT) :: iter

    INTEGER :: iter_max

    REAL*8, DIMENSION(n_mom,n_mom) :: bk , bkk

    REAL*8, DIMENSION(n_mom,n_mom) :: d

    INTEGER :: i , j , k

    LOGICAL :: check_corr

    INTEGER :: k_star

    REAL*8, DIMENSION(n_mom) :: cos_quad_alfa

    REAL*8 :: lnck

    iter_max = 100

    !
    ! define the unitary correction matrix
    !
    DO i = 1,n_mom

       DO j = 1,n_mom

          bkk(i,j) = 0.0

       END DO

    END DO

    DO k = 1,n_mom

       DO i = 1,n_mom

          IF ( i.EQ.k ) THEN 

             bkk(i,1) = 1.D0

          ELSE

             bkk(i,1) = 0.0

          END IF

       END DO

       DO j=2,4

          DO i=1,n_mom-j+1

             bkk(i,j) = bkk(i+1,j-1) - bkk(i,j-1)

          END DO

       END DO

       DO i=1,n_mom

          bk(k,i) = bkk(i,4)

       END DO

    END DO

    check_corr = .TRUE.
    iter = 0

    ! start the correction loop

    d = 0.D0

    DO WHILE ( ( check_corr ) .AND. ( iter .LT. iter_max ) )

       check_corr = .FALSE.

       DO i = 1,n_mom

          d(i,1) = LOG(m(i))

       END DO

       DO j = 2,n_mom

          DO i = 1,n_mom-j+1

             d(i,j) = d(i+1,j-1) - d(i,j-1)

             IF ( ( j .EQ. 3 ) .AND. ( d(i,j) .LT. 0.D0 ) ) THEN 

                check_corr = .TRUE.

             END IF

          END DO
       END DO

       IF ( check_corr ) THEN

          iter = iter + 1
          k_star = 1

          DO k = 1,n_mom

             cos_quad_alfa(k) = ( ( DOT_PRODUCT( bk(k,:) , d(:,4) ) ) /         &
                  ( SQRT( DOT_PRODUCT( bk(k,:) , bk(k,:) ) ) *                  &
                  SQRT( DOT_PRODUCT( d(:,4) , d(:,4) ) ) ) ) **2.D0

             IF ( cos_quad_alfa(k) .GE. cos_quad_alfa(k_star) ) THEN

                k_star = k

             END IF

          END DO

          lnck = -( DOT_PRODUCT( bk(k_star,:) , d(:,4) ) ) /                    &
               DOT_PRODUCT( bk(k_star,:) , bk(k_star,:) )

          m(k_star) = DEXP(lnck) * m(k_star)

       END IF

    END DO

  END SUBROUTINE moments_correction

  !******************************************************************************
  !> \brief Moments correction Wright
  !
  !> This subroutine compute uses the algorithm of McGraw for the correction of 
  !> an unrealizable moment set. The subroutine first analyzes the moment set and
  !> then, if it is unrealizable, identifies the moemnt that has to be changed 
  !> the least to make the set realizable. 
  !> Modified from Marchisio and Fox 2013
  !> \param[in,out]    m       moments
  !> \param[out]       iter    number of iterations required
  !> \author Mattia de' Michieli Vitturi
  !> \date 22/10/2013
  !******************************************************************************

  SUBROUTINE moments_correction_wright(m)

    IMPLICIT NONE
    
    REAL*8, DIMENSION(:), INTENT(INOUT) :: m

    INTEGER :: i,j

    INTEGER :: k

    REAL*8 :: mu , sigma_var

    REAL*8, DIMENSION(n_mom) :: m1 , m2

    i = 2
    j = 3

    ! calculate the corresponding mean and variance 

    mu = (j/(i*j-i**2))*log(m(i+1)/m(1)) + (i/(i*j-j**2))*log(m(j+1)/m(1))
    
    sigma_var = ((2.0/(j**2))*log(m(j+1)/m(1))-(2.0/(i*j))*log(m(i+1)/m(1)))/(1.0-i/j)

    if (sigma_var .LT. 0) then

       sigma_var = 0.0
    end if

    ! calculate the new moments for the first distribution
    
    DO k=1,n_mom

       m1(k) = m(1)*dexp((k-1.0)*mu+((k-1.0)**2*sigma_var)/2.0)

    END DO


    ! select the two additional moments to be used for the first log-normal

    i=1
    j=3
	
    ! calculate the corresponding mean and variance 

    mu = (j/(i*j-i**2))*log(m(i+1)/m(1)) + (i/(i*j-j**2))*log(m(j+1)/m(1))

    sigma_var = ((2.0/(j**2))*log(m(j+1)/m(1))-(2.0/(i*j))*log(m(i+1)/m(1)))/(1.0-i/j)

    IF (sigma_var < 0) THEN

       sigma_var = 0.0

    END IF

    ! calculate the new moments for the second distribution

    DO k=1,n_mom
    
       m2(k) =  m(1)*exp((k-1.0)*mu+((k-1.0)**2*sigma_var)/2.0)

    END DO

    ! calculate the new set of moments from the arithmetic mean of the previous
    ! two sets of moments
    
    DO k=1,n_mom

       m(k) = (m1(k)+m2(k))/2.0

    END DO
   
  END SUBROUTINE moments_correction_wright



  !******************************************************************************
  !> \brief Beta function
  !
  !> This function evaluates the beta function B(z,w). This is the name used by 
  !> Legendre and Whittaker and Watson (1990) for the beta integral (also called 
  !> the Eulerian integral of the first kind).
  !> \param[in]   z    first parameter of the beta function
  !> \param[in]   w    second parameter of the beta function
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  REAL*8 FUNCTION beta_function(z,w)

    IMPLICIT NONE

    REAL*8, INTENT(IN):: w,z

    REAL*8:: gamln_z , gamln_w , gamln_zw

    gamln_z = gammln(z)
    gamln_w = gammln(w)
    gamln_zw = gammln(z+w)

    beta_function = DEXP( gamln_z + gamln_w - gamln_zw )

    RETURN

  END FUNCTION beta_function

  !******************************************************************************
  !> \brief Gamma function logarithm
  !
  !> This function returns the logarithm of the gamma function, gammaln(A) = 
  !> log(gamma(A)). Input must be nonnegative and real. Converted from Numerical
  !> Recipes.
  !> \param[in]   xx   argument of the function
  !> \date 22/10/2013
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  FUNCTION gammln(xx)

    IMPLICIT NONE

    REAL*8 :: gammln,xx

    INTEGER :: j

    REAL*8 :: ser , tmp , x , y , cof(14)

    cof(1) = 57.1562356658629235
    cof(2) = -59.5979603554754912
    cof(3) = 14.1360979747417471
    cof(4) = -0.491913816097620199
    cof(5) = 0.339946499848118887d-4
    cof(6) = 0.465236289270485756d-4
    cof(7) = -0.983744753048795646d-4
    cof(8) = 0.158088703224912494d-3
    cof(9) = -0.210264441724104883d-3
    cof(10) = 0.217439618115212643d-3
    cof(11) = -0.164318106536763890d-3
    cof(12) = 0.844182239838527433d-4
    cof(13) = -0.261908384015814087d-4
    cof(14) = .368991826595316234d-5

    IF (xx .LE. 0) THEN
       
       WRITE(*,*) "bad arg in gammln"
       gammln = 0.D0

    ELSE

       x = xx
       y = x
       
       tmp = x + 5.2421875D0
       tmp = ( x + 0.5D0 ) * log(tmp) - tmp
       
       ser = 0.999999999999997092
       
       DO j=1,14
          
          y = y+1.D0
          ser = ser + cof(j)/y
          
       END DO
       
       gammln = tmp + log(2.5066282746310005*ser/x)
       
    END IF

    RETURN

  END FUNCTION gammln

  !******************************************************************************
  !> \brief Binomial Coefficient
  !
  !> This subroutine evaluates the binomial coefficients c(n,m) and stores the 
  !> values in the matrix (array) c(n,m). The coefficient is used for the 
  !> computation of the moments of the normal distribution.
  !> From: http://www-users.york.ac.uk/~hcb1/fortranY1/fortran90/array15.f90
  !> \param[in]    nmax   max index of the binomial coefficent
  !> \param[out]   c      array of the binomial coefficients
  !> \date 03/05/2014
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************


  SUBROUTINE coefficient(nmax,c)
    IMPLICIT NONE

    ! This example evaluates the binomial coefficients c(n,m) and 
    ! stores the values in the matrix (array) c(n,m)		    		
    ! declarations
    INTEGER :: n,m
    INTEGER, INTENT(IN) :: nmax
    ! specify the dimensions of the two dimensional array 
    ! (matrix) c(n,m).
    INTEGER, DIMENSION(0:nmax,0:nmax), INTENT(OUT) :: c

    ! initialize array. Note how all elements are set to zero 
    ! by this one statement. Fortran90 has the facility to treat 
    ! arrays as single objects. 
    c = 0

    ! use recurrence relation c(n,m)=c(n-1,m-1) + c(n-1,m) 
    ! where m<n to generate elements

    c(0,0) = 1
    c(1,0) = 1
    c(1,1) = 1

    DO  n = 2,nmax

       c(n,0) = 1
       c(n,n) = 1

       DO  m = 1,n-1

          c(n,m) = c(n-1,m-1) + c(n-1,m)

       END DO

    END DO

    RETURN

  END SUBROUTINE coefficient

END MODULE moments_module
