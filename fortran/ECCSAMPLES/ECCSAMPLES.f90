MODULE ECCSAMPLESmod

use inversebetamod

 implicit none

 CONTAINS

 !=======================================================================
 SUBROUTINE ECCSAMPLES(transit,verbose,a,b,xe,xw,efinal,wfinal)

 ! === AUTHOR ===
 ! David M. Kipping, Columbia University
 ! d.kipping@columbia.edu

 ! === DESCRIPTION ===
 ! ECCSAMPLES transforms a sample drawn from a bivariate uniform distribution
 ! into a sample drawn from P(e,w|transit), where P(e,w|transit) is the a-priori
 ! joint probability distribution a planet's orbital eccentricity, e, and 
 ! argument of periastron, w, given that the planet is known to transit its
 ! host star. A fundamental assumption of this code is that the intinsic
 ! prior distribution of eccentricity, P(e), is described by a Beta distribution
 ! with shape parameters a & b. One may also draw directly from the intrinsic
 ! priors instead, by setting the logical "transit" to FALSE.

 ! === INPUTS ===
 ! transit = Logical, where TRUE => planet is known to transit its host star;
 !           FALSE => it is unknown whether the planets transits.
 ! a = First shape parameter of the Beta distribution describing P(e).
 !     Must be real and satisfy a>0.
 ! b = Second shape parameter of the Beta distribution describing P(e).
 !     Must be real and satisfy b>1.
 ! xe = Proxy parameter for eccentricity, e. Must be real and satisfy 0<xe<1.
 !      xe is assumed to be uniformly distributed.
 ! xw = Proxy parameter for argument of periastron, w. Must be real and satisfy 
 !      0<xw<1. xw is assumed to be uniformly distributed.
 ! verbose = Logical determining whether to provide verbose output.

 ! === OUTPUTS ===
 ! efinal = Transformed value of xe, an eccentricity variate drawn from the
 !          target distribution. Real and falls in the range 0<efinal<1.
 ! wfinal = Transformed value of xw, an argument of periastron variate drawn 
 !          from the target distribution. Real and falls in the range 
 !          0<wfinal<2pi.

 ! === REFERENCE ===
 ! Use of this code should cite:
 ! Kipping, D. M., 2014, MNRAS, 444, 2263-2269
 ! Kipping, D. M., 2013, MNRAS, 434, L51-55

 ! === DEPENDENCIES ===
 ! There several dependecies needed to correctly compile ECCSAMPLES:
 ! * f1 by Colavecchia & Gasaneo (2004) 
 !   [at http://cpc.cs.qub.ac.uk/summaries/ADSJ_v1_0.html]
 ! * SPECFUN by William Cody and Laura Stoltz 
 !   [at http://people.sc.fsu.edu/~jburkardt/f_src/specfun/specfun.html]
 ! * inversebeta, modified code David Kipping originally by John Burkardt
 !   [included with ECCSAMPLES]
 !
 ! f1 REFERENCES...
 !
 ! FD Colavecchia, G Gasaneo,
 ! f1: a code to compute Appell's F1 hypergeometric function,
 ! Computer Physics Communications,
 ! Volume 157, Number 1, 2004, pages 32-38.
 ! 
 ! RC Forrey,
 ! Computing the hypergeometric function,
 ! Journal of Computational Physics,
 ! Volume 137, 1997, pages 79-100.
 !
 ! SPECFUN REFERENCES...
 !
 ! S Zhang, J Jin,
 ! Computation of Special Functions,
 ! Wiley, 1996,
 ! ISBN: 0-471-11963-6,
 ! LC: QA351.C45
 !
 ! inversebeta REFERENCES...
 !
 ! Allan Macleod,
 ! Algorithm AS 245,
 ! A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
 ! Applied Statistics,
 ! Volume 38, Number 2, 1989, pages 397-402.
 !
 ! GW Cran, KJ Martin, GE Thomas,
 ! Remark AS R19 and Algorithm AS 109:
 ! A Remark on Algorithms AS 63: The Incomplete Beta Integral
 ! and AS 64: Inverse of the Incomplete Beta Integeral,
 ! Applied Statistics,
 ! Volume 26, Number 1, 1977, pages 111-114.
 !
 ! KL Majumder, GP Bhattacharjee,
 ! Algorithm AS 63:
 ! The incomplete Beta Integral,
 ! Applied Statistics,
 ! Volume 22, Number 3, 1973, pages 409-411.

 ! === COMPUTATION TIME ===
 ! Computation time is ~500us per sample using an Intel Core i7-3720QM
 ! processor with hyperthreading activated.

 implicit none

 INTEGER :: i, ifault
 LOGICAL, INTENT(IN) :: transit, verbose
 REAL(8), INTENT(IN) :: a, b, xe, xw
 REAL(8), INTENT(OUT) :: efinal, wfinal
 REAL(8) :: etrial, eprec, REval
 REAL(8) :: wtrial, wprec, RWval
 REAL(8) :: beta_log
 LOGICAL :: process
 REAL(8), PARAMETER :: etol = 1.0D-3 		   ! Tolerance in precision of e
 REAL(8), PARAMETER :: wtol = 6.283185307179586D-3 ! Tolerance in precision of w
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0

 ! Initialize with assumption of no errors
 process = .TRUE.

 ! Check uniform variates are valid
 IF( xe.GT.1.0D0 .OR. xe.LT.0.0D0 .OR. xw.GT.1.0D0 .OR. xw.LT.0.0D0 ) THEN
   IF( verbose ) THEN
     write(*,*) 'ERROR: xe and xw must be in the interval [0,1]'
   END IF
   process = .FALSE. ! Instruct code to abort
 END IF

 ! Check shape parameters are valid
 IF( a.LE.0.0D0 .OR. b.LE.1.0D0 ) THEN
   IF( verbose ) THEN
     write(*,*) 'ERROR: ECCSAMPLES requires a>0 and b>1'
   END IF
   process = .FALSE. ! Instruct code to abort
 END IF

 IF( process ) THEN
 
   ! State uniform random variates, if verbose
   IF ( verbose ) THEN
     write(*,*) 'xe = ',xe
     write(*,*) 'xw = ',xw
   END IF

   IF( transit ) THEN ! Draw samples from P(e,w|transit)
    
     IF( verbose ) THEN
       write(*,*) 'Planet assumed to be transiting'
     END IF

     ! Define an intial starting point for the e iteration
     etrial = xe*0.5D0	! Start Newton's method from e_1 = xe/2
     eprec = 1.0D0	! Precision of first trial set to unity
     i = 1
     IF ( verbose ) THEN
       write(*,*) 'i delta(e_i) e_i'
       write(*,*) 1,eprec,etrial
     END IF
     ! Perform Newton's method on e
     DO WHILE ( eprec .GT. etol )
       i = i + 1
       call eratio(a,b,etrial,xe,REval)
       eprec = DABS(REval)
       etrial = etrial - REval
       IF ( verbose ) THEN
         write(*,*) i,eprec,etrial
       END IF
     END DO
     ! Define final e value
     efinal = etrial
     IF ( verbose ) THEN
       write(*,*) 'e(final) = ',efinal
     END IF

     ! Define an intial starting point for the w iteration
     wtrial = xw*twopi	! Start Newton's method from w_1 = xw*2pi
     wprec = 1.0D0	! Precision of first trial set to unity
     i = 1
     IF ( verbose ) THEN
       write(*,*) 'i delta(w_i) w_i'
       write(*,*) 1,wprec,wtrial
     END IF
     ! Perform Newton's method on w
     DO WHILE ( wprec .GT. wtol )
       i = i + 1
       call wratio(efinal,wtrial,xw,RWval)
       wprec = DABS(RWval)
       wtrial = wtrial - RWval
       IF ( verbose ) THEN
         write(*,*) i,wprec,wtrial
       END IF
     END DO
     ! Define final w value
     wfinal = wtrial
     IF ( verbose ) THEN
       write(*,*) 'w(final) = ',wfinal
     END IF

   ELSE ! Draw samples from P(e,w)

     IF( verbose ) THEN
       write(*,*) 'It is assumed that it is unknown whether the planet transits'
     END IF

     ! e drawn from Beta distribution
     beta_log = alngam(a,ifault) + alngam(b,ifault) - alngam(a+b,ifault)
     efinal = xinbta (a,b,beta_log,xe,ifault )
     IF ( verbose .AND. ifault .NE. 0 ) THEN
       write(*,*) 'ERROR: ifault = ',ifault,'; problem in inversebeta.f'
     END IF
     IF ( verbose ) THEN
       write(*,*) 'e(final) = ',efinal
     END IF
   
     ! w drawn from uniform distribution
     wfinal = xw*twopi	! Start Newton's method from w_1 = xw*2pi
     IF ( verbose ) THEN
       write(*,*) 'w(final) = ',wfinal
     END IF

   END IF

 ELSE ! process = .FALSE.

   ! Invalid inputs, forcing improper outputs
   efinal = -1.0D0; wfinal = -1.0D0

 END IF

 END SUBROUTINE ECCSAMPLES
 !=======================================================================

 !=======================================================================
 SUBROUTINE eratio(a,b,e,x,Rval)

 ! Computes the ratio of [ f(e) ]/[ df(e)/de ], where f(e) = CDF(e) - x
 ! and CDF(e) is the cumulative density function of the probability P(e|b<1)
 ! This ratio is used to iteratively solve the equation f(e) = 0, for e, using
 ! Newton's method.

 ! The ratio, R, is shown in Kipping (2014) to be given by:
 ! R = (1/4) * (1-e^2) * (1-e)^(1-b) * e^(1-a) * T
 !
 ! where...
 ! T = ( e^(a+1) / (1+a) )*( 2Q + 3P - U ) + 4S
 ! P = 2F1[1+a,1-b,2+a;e]
 ! Q = 2F1[1+a,2-b,2+a;e]
 ! U = F1[1+a,-b,1,2+a,e,-e]
 ! S = B[a,1+b;e] - x * B[a,b-1] * 2F1[1,a,a+b-1;-1]
 !
 ! where...
 ! B[a,b] is the Beta function
 ! B[a,b;z] is the incomplete Beta function with argument z
 ! 2F1[a,b,c;z] is Gauss' hypergeometric function with argument z
 ! F1[a,b,b';c;x,y] is Appell's hypergeometric function with arguments x & y
 !

 implicit none

 REAL(8), INTENT(IN) :: a, b, e, x
 REAL(8) :: bt1, bt2, bix, hfre, hfim
 REAL(8) :: Sval, Pre, Pim, Qre, Qim, Uval, Tval
 COMPLEX(8) :: c1, c2, c3, c4, f1
 REAL(8), INTENT(OUT) :: Rval

 ! Get S
 call beta(a,b+1.0D0,bt1)
 call beta(a,b-1.0D0,bt2)
 call incob(a,b+1.0D0,e,bix)
 call hyp(-1.0D0,1.0D0,a,a+b-1.0D0,hfre,hfim)
 Sval = bix*bt1 - x*bt2*hfre

 ! Get U, P, Q
 c1 = 1.0D0 + a
 c2 = -b
 c3 = 1.0D0
 c4 = 2.0D0 + a
 Uval = f1(c1,c2,c3,c4,e,-e)
 call hyp(e,a+1.0D0,1.0D0-b,a+2.0D0,Pre,Pim)
 call hyp(e,a+1.0D0,2.0D0-b,a+2.0D0,Qre,Qim)

 ! Get Tval
 Tval = ( e**(1.0D0+a)/( 1.0D0+a ) )*( 3.0D0*Pre + 2.0D0*Qre - Uval )
 Tval = Tval + 4.0D0*Sval

 ! Get Rval (ratio)
 Rval = 0.25D0*(1.0D0-e**2)*(1.0D0-e)**(1.0D0-b)*e**(1.0D0-a)*Tval

 END SUBROUTINE eratio
 !=======================================================================

 !=======================================================================
 SUBROUTINE wratio(e,w,x,Rval)

 ! Computes the ratio of [ f(w|e) ]/[ df(w|e)/dw ], where f(w|e) = CDF(w|e) - x
 ! and CDF(w|e) is the cumulative density function of the probability P(w|e,b<1)
 ! This ratio is used to iterarively solve the equation f(w|e) = 0, for w, using
 ! Newton's method.

 ! The ratio, R, is shown in Kipping (2014) to be given by:
 ! R = [ w - 2*pi*x + e*(1 - cos(w) ) ] /
 !     [ 1 + e*sin(w) ]

 implicit none

 REAL(8), INTENT(IN) :: e, w, x
 REAL(8), INTENT(OUT) :: Rval
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0

 ! Get Rval (ratio)
 Rval = w - twopi*x + e*( 1.0D0 - DCOS(w) )
 Rval = Rval/( 1.0D0 + e*DSIN(w) )

 END SUBROUTINE wratio
 !=======================================================================

END MODULE ECCSAMPLESmod
