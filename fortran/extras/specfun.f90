
subroutine incob ( a, b, x, bix )

!*****************************************************************************80
!
!! INCOB computes the incomplete beta function Ix(a,b).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BIX, the function value.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bix
  real ( kind = 8 ) bt
  real ( kind = 8 ) dk(51)
  real ( kind = 8 ) fk(51)
  integer ( kind = 4 ) k
  real ( kind = 8 ) s0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) x

  s0 = ( a + 1.0D+00 ) / ( a + b + 2.0D+00 )
  call beta ( a, b, bt )

  if ( x <= s0 ) then

    do k = 1, 20
      dk(2*k) = k * ( b - k ) * x / &
        ( a + 2.0D+00 * k - 1.0D+00 ) / ( a + 2.0D+00 * k )
    end do

    do k = 0, 20
      dk(2*k+1) = - ( a + k ) * ( a + b + k ) * x &
        / ( a + 2.0D+00 * k ) / ( a + 2.0D+00 * k + 1.0D+00 )
    end do

    t1 = 0.0D+00
    do k = 20, 1, -1
      t1 = dk(k) / ( 1.0D+00 + t1 )
    end do
    ta = 1.0D+00 / ( 1.0D+00 + t1 )
    bix = x ** a * ( 1.0D+00 - x ) ** b / ( a * bt ) * ta

  else

    do k = 1, 20
      fk(2*k) = k * ( a - k ) * ( 1.0D+00 - x ) &
        / ( b + 2.0D+00 * k - 1.0D+00 ) / ( b + 2.0D+00 * k )
    end do

    do k = 0,20
      fk(2*k+1) = - ( b + k ) * ( a + b + k ) * ( 1.0D+00 - x ) &
        / ( b + 2.0D+00 * k ) / ( b + 2.0D+00 * k + 1.0D+00 )
    end do

    t2 = 0.0D+00
    do k = 20, 1, -1
      t2 = fk(k) / ( 1.0D+00 + t2 )
    end do
    tb = 1.0D+00 / ( 1.0D+00 + t2 )
    bix = 1.0D+00 - x ** a * ( 1.0D+00 - x ) ** b / ( b * bt ) * tb

  end if

  return
end

subroutine beta ( p, q, bt )

!*****************************************************************************80
!
!! BETA computes the Beta function B(p,q).
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, Q, the parameters.
!    0 < P, 0 < Q.
!
!    Output, real ( kind = 8 ) BT, the value of B(P,Q).
!
  implicit none

  real ( kind = 8 ) bt
  real ( kind = 8 ) gp
  real ( kind = 8 ) gpq
  real ( kind = 8 ) gq
  real ( kind = 8 ) p
  real ( kind = 8 ) ppq
  real ( kind = 8 ) q

  call gamma ( p, gp )
  call gamma ( q, gq )
  ppq = p + q
  call gamma ( ppq, gpq )
  bt = gp * gq / gpq

  return
end

subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
  implicit none

  real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, & 
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, & 
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, & 
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)
  real ( kind = 8 ) ga
  real ( kind = 8 ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = 8 ) )
      end do
      z = z - real ( m, kind = 8 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end
