C ==============================================================
C === MODULE inversebeta ===

      MODULE inversebetamod
      implicit none
      CONTAINS

      function alngam ( xvalue, ifault )

c*********************************************************************72
c
cc ALNGAM computes the logarithm of the gamma function.
c
c  Modified:
c
c    30 March 1999
c
c  Author:
c
c    Allan Macleod
c    Modifications by John Burkardt
c
c  Reference:
c
c    Allan Macleod,
c    Algorithm AS 245,
c    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
c    Applied Statistics,
c    Volume 38, Number 2, 1989, pages 397-402.
c
c  Parameters:
c
c    Input, double precision XVALUE, the argument of the Gamma function.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    1, XVALUE is less than or equal to 0.
c    2, XVALUE is too big.
c
c    Output, double precision ALNGAM, the logarithm of the gamma function of X.
c
      implicit none

      double precision alngam
      double precision alr2pi
      parameter ( alr2pi = 0.918938533204673D+00 )
      integer ifault
      double precision r1(9)
      double precision r2(9)
      double precision r3(9)
      double precision r4(5)
      double precision x
      double precision x1
      double precision x2
      double precision xlge
      parameter ( xlge = 5.10D+06 )
      double precision xlgst
      parameter ( xlgst = 1.0D+30 )
      double precision xvalue
      double precision y

      data r1 /
     &  -2.66685511495D+00,
     &  -24.4387534237D+00,
     &  -21.9698958928D+00,
     &   11.1667541262D+00,
     &   3.13060547623D+00,
     &   0.607771387771D+00,
     &   11.9400905721D+00,
     &   31.4690115749D+00,
     &   15.2346874070D+00 /

      data r2 /
     &  -78.3359299449D+00,
     &  -142.046296688D+00,
     &   137.519416416D+00,
     &   78.6994924154D+00,
     &   4.16438922228D+00,
     &   47.0668766060D+00,
     &   313.399215894D+00,
     &   263.505074721D+00,
     &   43.3400022514D+00 /

      data r3 /
     &  -2.12159572323D+05,
     &   2.30661510616D+05,
     &   2.74647644705D+04,
     &  -4.02621119975D+04,
     &  -2.29660729780D+03,
     &  -1.16328495004D+05,
     &  -1.46025937511D+05,
     &  -2.42357409629D+04,
     &  -5.70691009324D+02 /

      data r4 / 
     &   0.279195317918525D+00, 
     &   0.4917317610505968D+00,
     &   0.0692910599291889D+00, 
     &   3.350343815022304D+00,
     &   6.012459259764103D+00 /

      x = xvalue
      alngam = 0.0D+00
c
c  Check the input.
c
      if ( xlgst .le. x ) then
        ifault = 2
        return
      end if

      if ( x .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0
c
c  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
c
      if ( x .lt. 1.5D+00 ) then

        if ( x .lt. 0.5D+00 ) then

          alngam = - dlog ( x )
          y = x + 1.0D+00
c
c  Test whether X < machine epsilon.
c
          if ( y .eq. 1.0D+00 ) then
            return
          end if

        else

          alngam = 0.0D+00
          y = x
          x = ( x - 0.5D+00 ) - 0.5D+00

        end if

        alngam = alngam + x * ((((
     &      r1(5)   * y 
     &    + r1(4) ) * y 
     &    + r1(3) ) * y
     &    + r1(2) ) * y 
     &    + r1(1) ) / ((((
     &                y 
     &    + r1(9) ) * y 
     &    + r1(8) ) * y
     &    + r1(7) ) * y 
     &    + r1(6) )

        return

      end if
c
c  Calculation for 1.5 <= X < 4.0.
c
      if ( x .lt. 4.0D+00 ) then

        y = ( x - 1.0D+00 ) - 1.0D+00

        alngam = y * ((((
     &      r2(5)   * x 
     &    + r2(4) ) * x 
     &    + r2(3) ) * x 
     &    + r2(2) ) * x
     &    + r2(1) ) / ((((
     &                x 
     &    + r2(9) ) * x 
     &    + r2(8) ) * x 
     &    + r2(7) ) * x
     &    + r2(6) )
c
c  Calculation for 4.0 <= X < 12.0.
c
      else if ( x .lt. 12.0D+00 ) then

        alngam = ((((
     &      r3(5)   * x 
     &    + r3(4) ) * x 
     &    + r3(3) ) * x 
     &    + r3(2) ) * x 
     &    + r3(1) ) / (((( 
     &                x 
     &    + r3(9) ) * x 
     &    + r3(8) ) * x 
     &    + r3(7) ) * x 
     &    + r3(6) )
c
c  Calculation for X >= 12.0.
c
      else

        y = dlog ( x )
        alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

        if ( x .le. xlge ) then

          x1 = 1.0D+00 / x
          x2 = x1 * x1

          alngam = alngam + x1 * ( ( 
     &           r4(3)   * 
     &      x2 + r4(2) ) * 
     &      x2 + r4(1) ) / ( ( 
     &      x2 + r4(5) ) * 
     &      x2 + r4(4) )

        end if

      end if

      return
      end function alngam

      function xinbta ( p, q, beta, alpha, ifault )

c*********************************************************************72
c
cc XINBTA computes inverse of the incomplete Beta function.
c
c  Modified:
c
c    10 January 2008
c
c  Author:
c
c    GW Cran, KJ Martin, GE Thomas
c    Modifications by John Burkardt
c
c  Reference:
c
c    GW Cran, KJ Martin, GE Thomas,
c    Remark AS R19 and Algorithm AS 109:
c    A Remark on Algorithms AS 63: The Incomplete Beta Integral
c    and AS 64: Inverse of the Incomplete Beta Integeral,
c    Applied Statistics,
c    Volume 26, Number 1, 1977, pages 111-114.
c
c  Parameters:
c
c    Input, double precision P, Q, the parameters of the incomplete
c    Beta function.
c
c    Input, double precision BETA, the logarithm of the value of
c    the complete Beta function.
c
c    Input, double precision ALPHA, the value of the incomplete Beta
c    function.  0 <= ALPHA <= 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    nonzero, an error occurred.
c
c    Output, double precision XINBTA, the argument of the incomplete
c    Beta function which produces the value ALPHA.
c
c  Local Parameters:
c
c    Local, double precision SAE, the most negative decimal exponent
c    which does not cause an underflow.
c
      implicit none

      double precision a
      double precision acu
      double precision adj
      double precision alpha
      double precision beta
      double precision betain
      double precision fpu
      double precision g
      double precision h
      integer iex
      integer ifault
      logical indx
      double precision p
      double precision pp
      double precision prev
      double precision q
      double precision qq
      double precision r
      double precision s
      double precision sae
      parameter ( sae = -37.0D+00 )
      double precision sq
      double precision t
      double precision tx
      double precision w
      double precision xin
      double precision xinbta
      double precision y
      double precision yprev

      fpu = 10.D+00**sae

      ifault = 0
      xinbta = alpha
c
c  Test for admissibility of parameters.
c
      if ( p .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( q .le. 0.0D+00 ) then
        ifault = 1
        return
      end if


      if ( alpha .lt. 0.0D+00 .or. 1.0D+00 .lt. alpha ) then
        ifault = 2
        return
      end if
c
c  If answer is easy to determine, return immediately.
c
      if ( alpha .eq. 0.0D+00 ) then
        xinbta = 0.0D+00
        return
      end if

      if ( alpha .eq. 1.0D+00 ) then
        xinbta = 1.0D+00
        return
      end if
c
c  Change tail if necessary.
c
      if ( 0.5D+00 .lt. alpha ) then
        a = 1.0D+00 - alpha
        pp = q
        qq = p
        indx = .true.
      else
        a = alpha
        pp = p
        qq = q
        indx = .false.
      end if
c
c  Calculate the initial approximation.
c
      r = dsqrt ( - dlog ( a * a ) )

      y = r - ( 2.30753D+00 + 0.27061D+00 * r ) 
     &  / ( 1.0D+00 + ( 0.99229D+00 + 0.04481D+00 * r ) * r )

      if ( 1.0D+00 .lt. p .and. 1.0D+00 .lt. q ) then

        r = ( y * y - 3.0D+00 ) / 6.0D+00
        s = 1.0D+00 / ( pp + pp - 1.0D+00 )
        t = 1.0D+00 / ( qq + qq - 1.0D+00 )
        h = 2.0D+00 / ( s + t )
        w = y * dsqrt ( h + r ) / h - ( t - s ) 
     &  * ( r + 5.0D+00 / 6.0D+00 - 2.0D+00 / ( 3.0D+00 * h ) )
        xinbta = pp / ( pp + qq * dexp ( w + w ) )

      else

        r = qq + qq
        t = 1.0D+00 / ( 9.0D+00 * qq )
        t = r * ( 1.0D+00 - t + y * dsqrt ( t ) )**3

        if ( t .le. 0.0D+00 ) then
          xinbta = 1.0D+00 
     &    - dexp ( ( dlog ( ( 1.0D+00 - a ) * qq ) + beta ) / qq )
        else

          t = ( 4.0D+00 * pp + r - 2.0D+00 ) / t
  
          if ( t .le. 1.0D+00 ) then
            xinbta = dexp ( ( dlog ( a * pp ) + beta ) / pp )
          else
            xinbta = 1.0D+00 - 2.0D+00 / ( t + 1.0D+00 )
          end if

        end if

      end if
c
c  Solve for X by a modified Newton-Raphson method,
c  using the function BETAIN.
c
      r = 1.0D+00 - pp
      t = 1.0D+00 - qq
      yprev = 0.0D+00
      sq = 1.0D+00
      prev = 1.0D+00

      if ( xinbta .lt. 0.0001D+00 ) then
        xinbta = 0.0001D+00
      end if

      if ( 0.9999D+00 .lt. xinbta ) then
        xinbta = 0.9999D+00
      end if

      iex = max ( - 5.0D+00 / pp**2 - 1.0D+00 / a**0.2D+00 - 13.0D+00, 
     &  sae )

      acu = 10.0D+00**iex

    7 continue

      y = betain ( xinbta, pp, qq, beta, ifault )

      if ( ifault .ne. 0 ) then
        ifault = 3
        return
      end if

      xin = xinbta
      y = ( y - a ) * dexp ( beta + r * dlog ( xin ) 
     &  + t * dlog ( 1.0D+00 - xin ) )

      if ( y * yprev .le. 0.0D+00 ) then
        prev = max ( sq, fpu )
      end if

      g = 1.0D+00

    9 continue

      adj = g * y
      sq = adj * adj

      if ( prev .le. sq ) then
        go to 10
      end if

      tx = xinbta - adj

      if ( 0.0D+00 .le. tx .and. tx .le. 1.0D+00 ) then
        go to 11
      end if

10    continue

      g = g / 3.0D+00
      go to 9

   11 continue

      if ( prev .le. acu ) then
        go to 12
      end if

      if ( y * y .le. acu ) then
        go to 12
      end if

      if ( tx .eq. 0.0D+00 .or. tx .eq. 1.0D+00 ) then
        go to 10
      end if

      if ( tx .eq. xinbta ) then
        go to 12
      end if

      xinbta = tx
      yprev = y
      go to 7

   12 continue

      if ( indx ) then 
        xinbta = 1.0D+00 - xinbta
      end if

      return
      end function xinbta

      END MODULE inversebetamod

      function betain ( x, p, q, beta, ifault )

c*********************************************************************72
c
cc BETAIN computes the incomplete Beta function ratio.
c
c  Modified:
c
c    06 January 2008
c
c  Author:
c
c    KL Majumder, GP Bhattacharjee
c    Modifications by John Burkardt
c
c  Reference:
c
c    KL Majumder, GP Bhattacharjee,
c    Algorithm AS 63:
c    The incomplete Beta Integral,
c    Applied Statistics,
c    Volume 22, Number 3, 1973, pages 409-411.
c
c  Parameters:
c
c    Input, double precision X, the argument, between 0 and 1.
c
c    Input, double precision P, Q, the parameters, which
c    must be positive.
c
c    Input, double precision BETA, the logarithm of the complete
c    beta function.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    nonzero, an error occurred.
c
c    Output, double precision BETAIN, the value of the incomplete
c    Beta function ratio.
c
      implicit none

      double precision acu
      parameter ( acu = 0.1D-14 )
      double precision ai
      double precision beta
      double precision betain
      double precision cx
      integer ifault
      logical indx
      integer ns
      double precision p
      double precision pp
      double precision psq
      double precision q
      double precision qq
      double precision rx
      double precision temp
      double precision term
      double precision x
      double precision xx

      betain = x
      ifault = 0
c
c  Check the input arguments.
c
      if ( p .le. 0.0D+00 .or. q .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( x .lt. 0.0D+00 .or. 1.0D+00 .lt. x ) then
        ifault = 2
        return
      end if
c
c  Special cases.
c
      if ( x .eq. 0.0D+00 .or. x .eq. 1.0D+00 ) then
        return
      end if
c
c  Change tail if necessary and determine S.
c
      psq = p + q
      cx = 1.0D+00 - x

      if ( p .lt. psq * x ) then
        xx = cx
        cx = x
        pp = q
        qq = p
        indx = .true.
      else
        xx = x
        pp = p
        qq = q
        indx = .false.
      end if

      term = 1.0D+00
      ai = 1.0D+00
      betain = 1.0D+00
      ns = int ( qq + cx * psq )
c
c  Use the Soper reduction formula.
c
      rx = xx / cx
      temp = qq - ai
      if ( ns .eq. 0 ) then
        rx = xx
      end if

10    continue

      term = term * temp * rx / ( pp + ai )
      betain = betain + term
      temp = dabs ( term )

      if ( temp .le. acu .and. temp .le. acu * betain ) then

        betain = betain * dexp ( pp * dlog ( xx ) 
     &  + ( qq - 1.0D+00 ) * dlog ( cx ) - beta ) / pp

        if ( indx ) then
          betain = 1.0D+00 - betain
        end if

        return

      end if

      ai = ai + 1.0D+00
      ns = ns - 1

      if ( ns .ge. 0 ) then
        temp = qq - ai
        if ( ns .eq. 0 ) then
          rx = xx
        end if
      else
        temp = psq
        psq = psq + 1.0D+00
      end if

      go to 10

      return
      end !function
