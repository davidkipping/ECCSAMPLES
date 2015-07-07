PROGRAM example_call

 use ECCSAMPLESmod

 implicit none

 INTEGER :: i
 INTEGER, PARAMETER :: n = 1E4 ! Number of samples to draw
 REAL(8) :: alpha, beta        ! Beta distribution shape parameters alpha & beta
 REAL(8) :: xe, xw             ! Random uniform variates
 REAL(8), DIMENSION(n) :: e, w
 LOGICAL, PARAMETER :: transit = .TRUE.
 LOGICAL, PARAMETER :: verbose = .FALSE.

 ! Define alpha and beta
 alpha = 0.867D0
 beta  = 3.030D0

 ! Generates samples
 DO i=1,n
   ! Get random uniform variates xe and xw
   call RANDOM_NUMBER(xe)
   call RANDOM_NUMBER(xw)
   call ECCSAMPLES(transit,verbose,alpha,beta,xe,xw,e(i),w(i))
 END DO

 ! Export samples
 OPEN(10,FILE='ECCSAMPLES.dat',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,n
   write(10,*) e(i),w(i)
 END DO
 CLOSE(10)

END PROGRAM example_call
