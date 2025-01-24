!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
! Modified version from https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.f90
! * embedded everything into module
! * changed single 32 bit types to 64 bit

module Minpack
    use iso_fortran_env
    implicit none

    public :: lmder, lmder1, syminv

    private

    contains

subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

!*****************************************************************************80
!
!! CHKDER checks the gradients of M functions of N variables.
!
!  Discussion:
!
!    CHKDER checks the gradients of M nonlinear functions in N variables,
!    evaluated at a point X, for consistency with the functions themselves.
!
!    The user calls CHKDER twice, first with MODE = 1 and then with MODE = 2.
!
!    MODE = 1.
!      On input,
!        X contains the point of evaluation.
!      On output,
!        XP is set to a neighboring point.
!
!    Now the user must evaluate the function and gradients at X, and the
!    function at XP.  Then the function is called again:
!
!    MODE = 2.
!      On input,
!        FVEC contains the function values at X,
!        FJAC contains the function gradients at X.
!        FVECP contains the functions evaluated at XP.
!      On output,
!        ERR contains measures of correctness of the respective gradients.
!
!    The function does not perform reliably if cancellation or
!    rounding errors cause a severe loss of significance in the
!    evaluation of a function.  Therefore, none of the components
!    of X should be unusually small (in particular, zero) or any
!    other value which may cause loss of significance.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) M, is the number of functions.
!
!    Input, integer(int32) N, is the number of variables.
!
!    Input, real(real64) X(N), the point at which the jacobian is to be
!    evaluated.
!
!    Input, real(real64) FVEC(M), is used only when MODE = 2.
!    In that case, it should contain the function values at X.
!
!    Input, real(real64) FJAC(LDFJAC,N), an M by N array.  When MODE = 2,
!    FJAC(I,J) should contain the value of dF(I)/dX(J).
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Output, real(real64) XP(N), on output with MODE = 1, is a neighboring
!    point of X, at which the function is to be evaluated.
!
!    Input, real(real64) FVECP(M), on input with MODE = 2, is the function
!    value at XP.
!
!    Input, integer(int32) MODE, should be set to 1 on the first call, and
!    2 on the second.
!
!    Output, real(real64) ERR(M).  On output when MODE = 2, ERR contains
!    measures of correctness of the respective gradients.  If there is no
!    severe loss of significance, then if ERR(I):
!      = 1.0D+00, the I-th gradient is correct,
!      = 0.0D+00, the I-th gradient is incorrect.
!      > 0.5D+00, the I-th gradient is probably correct.
!      < 0.5D+00, the I-th gradient is probably incorrect.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) eps
  real(real64) epsf
  real(real64) epslog
  real(real64) epsmch
  real(real64) err(m)
  real(real64) fjac(ldfjac,n)
  real(real64) fvec(m)
  real(real64) fvecp(m)
  integer(int32) i
  integer(int32) j
  integer(int32) mode
  real(real64) temp
  real(real64) x(n)
  real(real64) xp(n)

  epsmch = epsilon ( epsmch )
  eps = sqrt ( epsmch )
!
!  MODE = 1.
!
  if ( mode == 1 ) then

    do j = 1, n
      temp = eps * abs ( x(j) )
      if ( temp == 0.0D+00 ) then
        temp = eps
      end if
      xp(j) = x(j) + temp
    end do
!
!  MODE = 2.
!
  else if ( mode == 2 ) then

    epsf = 100.0D+00 * epsmch
    epslog = log10 ( eps )

    err = 0.0D+00

    do j = 1, n
      temp = abs ( x(j) )
      if ( temp == 0.0D+00 ) then
        temp = 1.0D+00
      end if
      err(1:m) = err(1:m) + temp * fjac(1:m,j)
    end do

    do i = 1, m

      temp = 1.0D+00

      if ( fvec(i) /= 0.0D+00 .and. fvecp(i) /= 0.0D+00 .and. &
        abs ( fvecp(i)-fvec(i)) >= epsf * abs ( fvec(i) ) ) then
        temp = eps * abs ( (fvecp(i)-fvec(i)) / eps - err(i) ) &
          / ( abs ( fvec(i) ) + abs ( fvecp(i) ) )
      end if

      err(i) = 1.0D+00

      if ( epsmch < temp .and. temp < eps ) then
        err(i) = ( log10 ( temp ) - epslog ) / epslog
      end if

      if ( eps <= temp ) then
        err(i) = 0.0D+00
      end if

    end do

  end if

  return
end
subroutine dogleg ( n, r, lr, diag, qtb, delta, x )

!*****************************************************************************80
!
!! DOGLEG finds the minimizing combination of Gauss-Newton and gradient steps.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA, the
!    problem is to determine the convex combination X of the
!    Gauss-Newton and scaled gradient directions that minimizes
!    (A*X - B) in the least squares sense, subject to the
!    restriction that the euclidean norm of D*X be at most DELTA.
!
!    This function completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization of A.  That is, if A = Q*R, where Q has
!    orthogonal columns and R is an upper triangular matrix,
!    then DOGLEG expects the full upper triangle of R and
!    the first N components of Q'*B.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the matrix R.
!
!    Input, real(real64) R(LR), the upper triangular matrix R stored
!    by rows.
!
!    Input, integer(int32) LR, the size of the R array, which must be 
!    no less than (N*(N+1))/2.
!
!    Input, real(real64) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real(real64) QTB(N), the first N elements of the vector Q'* B.
!
!    Input, real(real64) DELTA, is a positive upper bound on the
!    euclidean norm of D*X(1:N).
!
!    Output, real(real64) X(N), the desired convex combination of the
!    Gauss-Newton direction and the scaled gradient direction.
!
  implicit none

  integer(int32) lr
  integer(int32) n

  real(real64) alpha
  real(real64) bnorm
  real(real64) delta
  real(real64) diag(n)
  real(real64) epsmch
  real(real64) gnorm
  integer(int32) i
  integer(int32) j
  integer(int32) jj
  integer(int32) k
  integer(int32) l
  real(real64) qnorm
  real(real64) qtb(n)
  real(real64) r(lr)
  real(real64) sgnorm
  real(real64) sum2
  real(real64) temp
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) x(n)

  epsmch = epsilon ( epsmch )
!
!  Calculate the Gauss-Newton direction.
!
  jj = ( n * ( n + 1 ) ) / 2 + 1

  do k = 1, n

     j = n - k + 1
     jj = jj - k
     l = jj + 1
     sum2 = 0.0D+00

     do i = j + 1, n
       sum2 = sum2 + r(l) * x(i)
       l = l + 1
     end do

     temp = r(jj)

     if ( temp == 0.0D+00 ) then

       l = j
       do i = 1, j
         temp = max ( temp, abs ( r(l)) )
         l = l + n - i
       end do

       if ( temp == 0.0D+00 ) then
         temp = epsmch
       else
         temp = epsmch * temp
       end if

     end if

     x(j) = ( qtb(j) - sum2 ) / temp

  end do
!
!  Test whether the Gauss-Newton direction is acceptable.
!
  wa1(1:n) = 0.0D+00
  wa2(1:n) = diag(1:n) * x(1:n)
  qnorm = enorm ( n, wa2 )

  if ( qnorm <= delta ) then
    return
  end if
!
!  The Gauss-Newton direction is not acceptable.
!  Calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
     temp = qtb(j)
     do i = j, n
       wa1(i) = wa1(i) + r(l) * temp
       l = l + 1
     end do
     wa1(j) = wa1(j) / diag(j)
  end do
!
!  Calculate the norm of the scaled gradient.
!  Test for the special case in which the scaled gradient is zero.
!
  gnorm = enorm ( n, wa1 )
  sgnorm = 0.0D+00
  alpha = delta / qnorm

  if ( gnorm /= 0.0D+00 ) then
!
!  Calculate the point along the scaled gradient which minimizes the quadratic.
!
    wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)

    l = 1
    do j = 1, n
      sum2 = 0.0D+00
      do i = j, n
        sum2 = sum2 + r(l) * wa1(i)
        l = l + 1
      end do
      wa2(j) = sum2
    end do

    temp = enorm ( n, wa2 )
    sgnorm = ( gnorm / temp ) / temp
!
!  Test whether the scaled gradient direction is acceptable.
!
    alpha = 0.0D+00
!
!  The scaled gradient direction is not acceptable.
!  Calculate the point along the dogleg at which the quadratic is minimized.
!
    if ( sgnorm < delta ) then

      bnorm = enorm ( n, qtb )
      temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta )
      temp = temp - ( delta / qnorm ) * ( sgnorm / delta) ** 2 &
        + sqrt ( ( temp - ( delta / qnorm ) ) ** 2 &
        + ( 1.0D+00 - ( delta / qnorm ) ** 2 ) &
        * ( 1.0D+00 - ( sgnorm / delta ) ** 2 ) )

      alpha = ( ( delta / qnorm ) * ( 1.0D+00 - ( sgnorm / delta ) ** 2 ) ) &
        / temp

    end if

  end if
!
!  Form appropriate convex combination of the Gauss-Newton
!  direction and the scaled gradient direction.
!
  temp = ( 1.0D+00 - alpha ) * min ( sgnorm, delta )

  x(1:n) = temp * wa1(1:n) + alpha * x(1:n)

  return
end
function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) N, is the length of the vector.
!
!    Input, real(real64) X(N), the vector whose norm is desired.
!
!    Output, real(real64) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer(int32) n
  real(real64) x(n)
  real(real64) enorm

  enorm = sqrt ( sum ( x(1:n) ** 2 ))

  return
end
function enorm2 ( n, x )

!*****************************************************************************80
!
!! ENORM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This routine was named ENORM.  It has been renamed "ENORM2",
!    and a simplified routine has been substituted.
!
!    The Euclidean norm is computed by accumulating the sum of
!    squares in three different sums.  The sums of squares for the
!    small and large components are scaled so that no overflows
!    occur.  Non-destructive underflows are permitted.  Underflows
!    and overflows do not occur in the computation of the unscaled
!    sum of squares for the intermediate components.
!
!    The definitions of small, intermediate and large components
!    depend on two constants, RDWARF and RGIANT.  The main
!    restrictions on these constants are that RDWARF^2 not
!    underflow and RGIANT^2 not overflow.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer(int32) N, is the length of the vector.
!
!    Input, real(real64) X(N), the vector whose norm is desired.
!
!    Output, real(real64) ENORM2, the Euclidean norm of the vector.
!
  implicit none

  integer(int32) n

  real(real64) agiant
  real(real64) enorm2
  integer(int32) i
  real(real64) rdwarf
  real(real64) rgiant
  real(real64) s1
  real(real64) s2
  real(real64) s3
  real(real64) x(n)
  real(real64) xabs
  real(real64) x1max
  real(real64) x3max

  rdwarf = sqrt ( tiny ( rdwarf ) )
  rgiant = sqrt ( huge ( rgiant ) )

  s1 = 0.0D+00
  s2 = 0.0D+00
  s3 = 0.0D+00
  x1max = 0.0D+00
  x3max = 0.0D+00
  agiant = rgiant / real ( n, real64 )

  do i = 1, n

    xabs = abs ( x(i) )

    if ( xabs <= rdwarf ) then

      if ( x3max < xabs ) then
        s3 = 1.0D+00 + s3 * ( x3max / xabs ) ** 2
        x3max = xabs
      else if ( xabs /= 0.0D+00 ) then
        s3 = s3 + ( xabs / x3max ) ** 2
      end if

    else if ( agiant <= xabs ) then

      if ( x1max < xabs ) then
        s1 = 1.0D+00 + s1 * ( x1max / xabs ) ** 2
        x1max = xabs
      else
        s1 = s1 + ( xabs / x1max ) ** 2
      end if

    else

      s2 = s2 + xabs ** 2

    end if

  end do
!
!  Calculation of norm.
!
  if ( s1 /= 0.0D+00 ) then

    enorm2 = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )

  else if ( s2 /= 0.0D+00 ) then

    if ( x3max <= s2 ) then
      enorm2 = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
    else
      enorm2 = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
    end if

  else

    enorm2 = x3max * sqrt ( s3 )

  end if

  return
end
subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

!*****************************************************************************80
!
!! FDJAC1 estimates an N by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This function computes a forward-difference approximation
!    to the N by N jacobian matrix associated with a specified
!    problem of N functions in N variables. If the jacobian has
!    a banded form, then function evaluations are saved by only
!    approximating the nonzero terms.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer(int32) n
!      real(real64) fvec(n)
!      integer(int32) iflag
!      real(real64) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) N, the number of functions and variables.
!
!    Input, real(real64) X(N), the point where the jacobian is evaluated.
!
!    Input, real(real64) FVEC(N), the functions evaluated at X.
!
!    Output, real(real64) FJAC(LDFJAC,N), the N by N approximate
!    jacobian matrix.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC, which
!    must not be less than N.
!
!    Output, integer(int32) IFLAG, is an error flag returned by FCN.  
!    If FCN returns a nonzero value of IFLAG, then this routine returns 
!    immediately to the calling program, with the value of IFLAG.
!
!    Input, integer(int32) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the
!    jacobian is not banded, set ML and MU to N-1.
!
!    Input, real(real64) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) n

  real(real64) eps
  real(real64) epsfcn
  real(real64) epsmch
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fvec(n)
  real(real64) h
  integer(int32) i
  integer(int32) iflag
  integer(int32) j
  integer(int32) k
  integer(int32) ml
  integer(int32) msum
  integer(int32) mu
  real(real64) temp
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) x(n)

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )
  msum = ml + mu + 1
!
!  Computation of dense approximate jacobian.
!
  if ( n <= msum ) then

     do j = 1, n

        temp = x(j)
        h = eps * abs ( temp )
        if ( h == 0.0D+00 ) then
          h = eps
        end if

        iflag = 1
        x(j) = temp + h
        call fcn ( n, x, wa1, iflag )

        if ( iflag < 0 ) then
          exit
        end if

        x(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h

     end do

  else
!
!  Computation of banded approximate jacobian.
!
     do k = 1, msum

        do j = k, n, msum
          wa2(j) = x(j)
          h = eps * abs ( wa2(j) )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
          x(j) = wa2(j) + h
        end do

        iflag = 1
        call fcn ( n, x, wa1, iflag )

        if ( iflag < 0 ) then
          exit
        end if

        do j = k, n, msum

           x(j) = wa2(j)

           h = eps * abs ( wa2(j) )
           if ( h == 0.0D+00 ) then
             h = eps
           end if

           fjac(1:n,j) = 0.0D+00

           do i = 1, n
             if ( j - mu <= i .and. i <= j + ml ) then
               fjac(i,j) = ( wa1(i) - fvec(i) ) / h
             end if
           end do

        end do

     end do

  end if

  return
end
subroutine fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )

!*****************************************************************************80
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This function computes a forward-difference approximation
!    to the M by N jacobian matrix associated with a specified
!    problem of M functions in N variables.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer(int32) n
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) M, is the number of functions.
!
!    Input, integer(int32) N, is the number of variables.  
!    N must not exceed M.
!
!    Input, real(real64) X(N), the point where the jacobian is evaluated.
!
!    Input, real(real64) FVEC(M), the functions evaluated at X.
!
!    Output, real(real64) FJAC(LDFJAC,N), the M by N approximate
!    jacobian matrix.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC, 
!    which must not be less than M.
!
!    Output, integer(int32) IFLAG, is an error flag returned by FCN.  
!    If FCN returns a nonzero value of IFLAG, then this routine returns 
!    immediately to the calling program, with the value of IFLAG.
!
!    Input, real(real64) EPSFCN, is used in determining a suitable
!    step length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) eps
  real(real64) epsfcn
  real(real64) epsmch
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fvec(m)
  real(real64) h
!  integer(int32) i
  integer(int32) iflag
  integer(int32) j
  real(real64) temp
  real(real64) wa(m)
  real(real64) x(n)

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )

  do j = 1, n

    temp = x(j)
    h = eps * abs ( temp )
    if ( h == 0.0D+00 ) then
      h = eps
    end if

    iflag = 1
    x(j) = temp + h
    call fcn ( m, n, x, wa, iflag )

    if ( iflag < 0 ) then
      exit
    end if

    x(j) = temp
    fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

  end do

  return
end
subroutine hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
  factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

!*****************************************************************************80
!
!! HYBRD seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  The user must provide a
!    subroutine which calculates the functions.  
!
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer(int32) n
!      real(real64) fvec(n)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) N, the number of functions and variables.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(N), the functions evaluated at the output X.
!
!    Input, real(real64) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, integer(int32) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, integer(int32) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the jacobian
!    is not banded, set ML and MU to at least n - 1.
!
!    Input, real(real64) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real(real64) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(int32) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(real64) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(int32) NPRINT, enables controlled printing of
!    iterates if it is positive.  In this case, FCN is called with IFLAG = 0 at
!    the beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG.
!    See the description of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, relative error between two consecutive iterates is at most XTOL.
!    2, number of calls to FCN has reached or exceeded MAXFEV.
!    3, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress, as measured by the improvement
!       from the last five jacobian evaluations.
!    5, iteration is not making good progress, as measured by the improvement
!       from the last ten iterations.
!
!    Output, integer(int32) NFEV, the number of calls to FCN.
!
!    Output, real(real64) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Output, real(real64) R(LR), the upper triangular matrix produced by
!    the QR factorization of the final approximate jacobian, stored rowwise.
!
!    Input, integer(int32) LR, the size of the R array, which must be no
!    less than (N*(N+1))/2.
!
!    Output, real(real64) QTF(N), contains the vector Q'*FVEC.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) lr
  integer(int32) n

  real(real64) actred
  real(real64) delta
  real(real64) diag(n)
  real(real64) epsfcn
  real(real64) epsmch
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fnorm
  real(real64) fnorm1
  real(real64) fvec(n)
  integer(int32) i
  integer(int32) iflag
  integer(int32) info
  integer(int32) iter
  integer(int32) iwa(1)
  integer(int32) j
  logical jeval
  integer(int32) l
  integer(int32) maxfev
  integer(int32) ml
  integer(int32) mode
  integer(int32) msum
  integer(int32) mu
  integer(int32) ncfail
  integer(int32) nslow1
  integer(int32) nslow2
  integer(int32) ncsuc
  integer(int32) nfev
  integer(int32) nprint
  logical pivot
  real(real64) pnorm
  real(real64) prered
  real(real64) qtf(n)
  real(real64) r(lr)
  real(real64) ratio
  logical sing
  real(real64) sum2
  real(real64) temp
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) wa3(n)
  real(real64) wa4(n)
  real(real64) x(n)
  real(real64) xnorm
  real(real64) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    return
  else if ( xtol < 0.0D+00 ) then
    return
  else if ( maxfev <= 0 ) then
    return
  else if ( ml < 0 ) then
    return
  else if ( mu < 0 ) then
    return
  else if ( factor <= 0.0D+00 ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
    return
  end if

  if ( mode == 2 ) then

    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do

  end if
!
!  Evaluate the function at the starting point
!  and calculate its norm.
!
  iflag = 1
  call fcn ( n, x, fvec, iflag )
  nfev = 1

  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( n, fvec )
!
!  Determine the number of calls to FCN needed to compute the jacobian matrix.
!
  msum = min ( ml + mu + 1, n )
!
!  Initialize iteration counter and monitors.
!
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
!
!  Beginning of the outer loop.
!
30 continue

    jeval = .true.
!
!  Calculate the jacobian matrix.
!
    iflag = 2
    call fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

    nfev = nfev + msum

    if ( iflag < 0 ) then
      go to 300
    end if
!
!  Compute the QR factorization of the jacobian.
!
    pivot = .false.
    call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
!
!  On the first iteration, if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then

        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do

      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm ( n, wa3 )
      delta = factor * xnorm
      if ( delta == 0.0D+00 ) then
        delta = factor
      end if

    end if
!
!  Form Q' * FVEC and store in QTF.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then
         temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
         qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
       end if

     end do
!
!  Copy the triangular factor of the QR factorization into R.
!
     sing = .false.

     do j = 1, n
        l = j
        do i = 1, j - 1
          r(l) = fjac(i,j)
          l = l + n - i
        end do
        r(l) = wa1(j)
        if ( wa1(j) == 0.0D+00 ) then
          sing = .true.
        end if
     end do
!
!  Accumulate the orthogonal factor in FJAC.
!
     call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
180    continue
!
!  If requested, call FCN to enable printing of iterates.
!
        if ( 0 < nprint ) then
          iflag = 0
          if ( mod ( iter - 1, nprint ) == 0 ) then
            call fcn ( n, x, fvec, iflag )
          end if
          if ( iflag < 0 ) then
            go to 300
          end if
        end if
!
!  Determine the direction P.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( n, wa2, wa4, iflag )
        nfev = nfev + 1

        if ( iflag < 0 ) then
          go to 300
        end if

        fnorm1 = enorm ( n, wa4 )
!
!  Compute the scaled actual reduction.
!
        actred = -1.0D+00
        if ( fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
        endif
!
!  Compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
          sum2 = 0.0D+00
          do j = i, n
            sum2 = sum2 + r(l) * wa1(j)
            l = l + 1
          end do
          wa3(i) = qtf(i) + sum2
        end do

        temp = enorm ( n, wa3 )
        prered = 0.0D+00
        if ( temp < fnorm ) then
          prered = 1.0D+00 - ( temp / fnorm ) ** 2
        end if
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        if ( 0.0D+00 < prered ) then
          ratio = actred / prered
        end if
!
!  Update the step bound.
!
        if ( ratio < 0.1D+00 ) then

          ncsuc = 0
          ncfail = ncfail + 1
          delta = 0.5D+00 * delta

        else

          ncfail = 0
          ncsuc = ncsuc + 1

          if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
            delta = max ( delta, pnorm / 0.5D+00 )
          end if

          if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
            delta = pnorm / 0.5D+00
          end if

        end if
!
!  Test for successful iteration.
!
!  Successful iteration.
!  Update X, FVEC, and their norms.
!
        if ( 0.0001D+00 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Determine the progress of the iteration.
!
        nslow1 = nslow1 + 1
        if ( 0.001D+00 <= actred ) then
          nslow1 = 0
        end if

        if ( jeval ) then
          nslow2 = nslow2 + 1
        end if

        if ( 0.1D+00 <= actred ) then
          nslow2 = 0
        end if
!
!  Test for convergence.
!
        if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
          info = 1
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Tests for termination and stringent tolerances.
!
        if ( maxfev <= nfev ) then
          info = 2
        end if

        if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
          info = 3
        end if

        if ( nslow2 == 5 ) then
          info = 4
        end if

        if ( nslow1 == 10 ) then
          info = 5
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Criterion for recalculating jacobian approximation
!  by forward differences.
!
        if ( ncfail == 2 ) then
          go to 290
        end if
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
        do j = 1, n
          sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
          wa2(j) = ( sum2 - wa3(j) ) / pnorm
          wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
          if ( 0.0001D+00 <= ratio ) then
            qtf(j) = sum2
          end if
        end do
!
!  Compute the QR factorization of the updated jacobian.
!
        call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
        call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
        call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
!
!  End of the inner loop.
!
        jeval = .false.
        go to 180

  290   continue
!
!  End of the outer loop.
!
     go to 30

  300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( n, x, fvec, iflag )
  end if

  return
end
subroutine hybrd1 ( fcn, n, x, fvec, tol, info )

!*****************************************************************************80
!
!! HYBRD1 seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRD.  The user must provide a
!    subroutine which calculates the functions.  
!
!    The jacobian is calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    19 August 2016
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer(int32) n
!      real(real64) fvec(n)
!      integer(int32) iflag
!      real(real64) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) N, the number of functions and variables.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(N), the functions evaluated at the output X.
!
!    Input, real(real64) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the
!    description of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN has reached or exceeded 200*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, the iteration is not making good progress.
!
  implicit none

!  integer(int32) lwa
  integer(int32) n

  real(real64) diag(n)
  real(real64) epsfcn
  real(real64) factor
  external fcn
  real(real64) fjac(n,n)
  real(real64) fvec(n)
  integer(int32) info
!  integer(int32) j
  integer(int32) ldfjac
  integer(int32) lr
  integer(int32) maxfev
  integer(int32) ml
  integer(int32) mode
  integer(int32) mu
  integer(int32) nfev
  integer(int32) nprint
  real(real64) qtf(n)
  real(real64) r((n*(n+1))/2)
  real(real64) tol
  real(real64) x(n)
  real(real64) xtol

  if ( n <= 0 ) then
    info = 0
    return
  end if

  if ( tol < 0.0D+00 ) then
    info = 0
    return
  end if

  xtol = tol
  maxfev = 200 * ( n + 1 )
  ml = n - 1
  mu = n - 1
  epsfcn = 0.0D+00
  diag(1:n) = 1.0D+00
  mode = 2
  factor = 100.0D+00
  nprint = 0
  info = 0
  nfev = 0
  fjac(1:n,1:n) = 0.0D+00
  ldfjac = n
  r(1:(n*(n+1))/2) = 0.0D+00
  lr = ( n * ( n + 1 ) ) / 2
  qtf(1:n) = 0.0D+00

  call hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
    factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end
subroutine hybrj ( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
  factor, nprint, info, nfev, njev, r, lr, qtf )

!*****************************************************************************80
!
!! HYBRJ seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRJ finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  The user must provide a
!    subroutine which calculates the functions and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!
!      subroutine fcn ( n, x, fvec, fjac, ldfjac, iflag )
!      integer(int32) ldfjac
!      integer(int32) n
!      real(real64) fjac(ldfjac,n)
!      real(real64) fvec(n)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) N, the number of functions and variables.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(N), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an N by N matrix, containing
!    the orthogonal matrix Q produced by the QR factorization
!    of the final approximate jacobian.
!
!    Input, integer(int32) LDFJAC, the leading dimension of the
!    array FJAC.  LDFJAC must be at least N.
!
!    Input, real(real64) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, integer(int32) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input/output, real(real64) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(int32) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(real64) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(int32) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG.  
!    See the description of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, relative error between two consecutive iterates is at most XTOL.
!    2, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    3, XTOL is too small.  No further improvement in
!       the approximate solution X is possible.
!    4, iteration is not making good progress, as measured by the
!       improvement from the last five jacobian evaluations.
!    5, iteration is not making good progress, as measured by the
!       improvement from the last ten iterations.
!
!    Output, integer(int32) NFEV, the number of calls to FCN 
!    with IFLAG = 1.
!
!    Output, integer(int32) NJEV, the number of calls to FCN 
!    with IFLAG = 2.
!
!    Output, real(real64) R(LR), the upper triangular matrix produced
!    by the QR factorization of the final approximate jacobian, stored rowwise.
!
!    Input, integer(int32) LR, the size of the R array, which must 
!    be no less than (N*(N+1))/2.
!
!    Output, real(real64) QTF(N), contains the vector Q'*FVEC.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) lr
  integer(int32) n

  real(real64) actred
  real(real64) delta
  real(real64) diag(n)
  real(real64) epsmch
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fnorm
  real(real64) fnorm1
  real(real64) fvec(n)
  integer(int32) i
  integer(int32) iflag
  integer(int32) info
  integer(int32) iter
  integer(int32) iwa(1)
  integer(int32) j
  logical jeval
  integer(int32) l
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) ncfail
  integer(int32) nslow1
  integer(int32) nslow2
  integer(int32) ncsuc
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  logical pivot
  real(real64) pnorm
  real(real64) prered
  real(real64) qtf(n)
  real(real64) r(lr)
  real(real64) ratio
  logical sing
  real(real64) sum2
  real(real64) temp
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) wa3(n)
  real(real64) wa4(n)
  real(real64) x(n)
  real(real64) xnorm
  real(real64) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    if ( iflag < 0 ) then
      info = iflag
    end if
    iflag = 0
    if ( 0 < nprint ) then
      call fcn ( n, x, fvec, fjac, ldfjac, iflag )
    end if
    return
  end if

  if ( ldfjac < n .or. &
       xtol < 0.0D+00 .or. &
       maxfev <= 0 .or. &
       factor <= 0.0D+00 .or. &
       lr < ( n * ( n + 1 ) ) / 2 ) then
    if ( iflag < 0 ) then
      info = iflag
    end if
    iflag = 0
    if ( 0 < nprint ) then
      call fcn ( n, x, fvec, fjac, ldfjac, iflag )
    end if
    return
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        if ( iflag < 0 ) then
          info = iflag
        end if
        iflag = 0
        if ( 0 < nprint ) then
          call fcn ( n, x, fvec, fjac, ldfjac, iflag )
        end if
        return
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( n, x, fvec, fjac, ldfjac, iflag )
  nfev = 1

  if ( iflag < 0 ) then
    if ( iflag < 0 ) then
      info = iflag
    end if
    iflag = 0
    if ( 0 < nprint ) then
      call fcn ( n, x, fvec, fjac, ldfjac, iflag )
    end if
    return
  end if

  fnorm = enorm ( n, fvec )
!
!  Initialize iteration counter and monitors.
!
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
!
!  Beginning of the outer loop.
!
  do

    jeval = .true.
!
!  Calculate the jacobian matrix.
!
    iflag = 2
    call fcn ( n, x, fvec, fjac, ldfjac, iflag )
    njev = njev + 1

    if ( iflag < 0 ) then
      info = iflag
      iflag = 0
      if ( 0 < nprint ) then
        call fcn ( n, x, fvec, fjac, ldfjac, iflag )
      end if
      return
    end if
!
!  Compute the QR factorization of the jacobian.
!
    pivot = .false.
    call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
!
!  On the first iteration, if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do
      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm ( n, wa3 )
      delta = factor * xnorm
      if ( delta == 0.0D+00 ) then
        delta = factor
      end if

    end if
!
!  Form Q'*FVEC and store in QTF.
!
    qtf(1:n) = fvec(1:n)

    do j = 1, n
      if ( fjac(j,j) /= 0.0D+00 ) then
        sum2 = 0.0D+00
        do i = j, n
          sum2 = sum2 + fjac(i,j) * qtf(i)
        end do
        temp = - sum2 / fjac(j,j)
        do i = j, n
          qtf(i) = qtf(i) + fjac(i,j) * temp
        end do
      end if
    end do
!
!  Copy the triangular factor of the QR factorization into R.
!
    sing = .false.

    do j = 1, n
      l = j
      do i = 1, j - 1
        r(l) = fjac(i,j)
        l = l + n - i
      end do
      r(l) = wa1(j)
      if ( wa1(j) == 0.0D+00 ) then
        sing = .true.
      end if
    end do
!
!  Accumulate the orthogonal factor in FJAC.
!
    call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
    if ( mode /= 2 ) then
      do j = 1, n
        diag(j) = max ( diag(j), wa2(j) )
      end do
    end if
!
!  Beginning of the inner loop.
!
    do
!
!  If requested, call FCN to enable printing of iterates.
!
      if ( 0 < nprint ) then

        iflag = 0
        if ( mod ( iter - 1, nprint ) == 0 ) then
          call fcn ( n, x, fvec, fjac, ldfjac, iflag )
        end if

        if ( iflag < 0 ) then
          info = iflag
          iflag = 0
          if ( 0 < nprint ) then
            call fcn ( n, x, fvec, fjac, ldfjac, iflag )
          end if
          return
        end if

      end if
!
!  Determine the direction P.
!
      call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
      wa1(1:n) = - wa1(1:n)
      wa2(1:n) = x(1:n) + wa1(1:n)
      wa3(1:n) = diag(1:n) * wa1(1:n)
      pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
      if ( iter == 1 ) then
        delta = min ( delta, pnorm )
      end if
!
!  Evaluate the function at X + P and calculate its norm.
!
      iflag = 1
      call fcn ( n, wa2, wa4, fjac, ldfjac, iflag )
      nfev = nfev + 1

      if ( iflag < 0 ) then
        info = iflag
        iflag = 0
        if ( 0 < nprint ) then
          call fcn ( n, x, fvec, fjac, ldfjac, iflag )
        end if
        return
      end if

      fnorm1 = enorm ( n, wa4 )
!
!  Compute the scaled actual reduction.
!
      actred = -1.0D+00
      if ( fnorm1 < fnorm ) then
        actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
      end if
!
!  Compute the scaled predicted reduction.
!
      l = 1
      do i = 1, n
        sum2 = 0.0D+00
        do j = i, n
          sum2 = sum2 + r(l) * wa1(j)
          l = l + 1
        end do
        wa3(i) = qtf(i) + sum2
      end do

      temp = enorm ( n, wa3 )
      prered = 0.0D+00
      if ( temp < fnorm ) then
        prered = 1.0D+00 - ( temp / fnorm ) ** 2
      end if
!
!  Compute the ratio of the actual to the predicted reduction.
!
      if ( 0.0D+00 < prered ) then
        ratio = actred / prered
      else
        ratio = 0.0D+00
      end if
!
!  Update the step bound.
!
      if ( ratio < 0.1D+00 ) then

        ncsuc = 0
        ncfail = ncfail + 1
        delta = 0.5D+00 * delta

      else

        ncfail = 0
        ncsuc = ncsuc + 1

        if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
          delta = max ( delta, pnorm / 0.5D+00 )
        end if

        if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
          delta = pnorm / 0.5D+00
        end if

      end if
!
!  Test for successful iteration.
!

!
!  Successful iteration.
!  Update X, FVEC, and their norms.
!
      if ( 0.0001D+00 <= ratio ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:n) = wa4(1:n)
        xnorm = enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
      end if
!
!  Determine the progress of the iteration.
!
      nslow1 = nslow1 + 1
      if ( 0.001D+00 <= actred ) then
        nslow1 = 0
      end if

      if ( jeval ) then
        nslow2 = nslow2 + 1
      end if

      if ( 0.1D+00 <= actred ) then
        nslow2 = 0
      end if
!
!  Test for convergence.
!
      if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
        info = 1
      end if

      if ( info /= 0 ) then
        iflag = 0
        if ( 0 < nprint ) then
          call fcn ( n, x, fvec, fjac, ldfjac, iflag )
        end if
        return
      end if
!
!  Tests for termination and stringent tolerances.
!
      if ( maxfev <= nfev ) then
        info = 2
      end if

      if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
        info = 3
      end if

      if ( nslow2 == 5 ) then
        info = 4
      end if

      if ( nslow1 == 10 ) then
        info = 5
      end if

      if ( info /= 0 ) then
        iflag = 0
        if ( 0 < nprint ) then
          call fcn ( n, x, fvec, fjac, ldfjac, iflag )
        end if
        return
      end if
!
!  Criterion for recalculating jacobian.
!
      if ( ncfail == 2 ) then
        exit
      end if
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
      do j = 1, n
        sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
        wa2(j) = ( sum2 - wa3(j) ) / pnorm
        wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
        if ( 0.0001D+00 <= ratio ) then
          qtf(j) = sum2
        end if
      end do
!
!  Compute the QR factorization of the updated jacobian.
!
      call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
      call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
      call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
!
!  End of the inner loop.
!
      jeval = .false.

    end do
!
!  End of the outer loop.
!
  end do

end
subroutine hybrj1 ( fcn, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! HYBRJ1 seeks a zero of N equations in N variables by Powell's method.
!
!  Discussion:
!
!    HYBRJ1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRJ.  The user
!    must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( n, x, fvec, fjac, ldfjac, iflag )
!      integer(int32) ldfjac
!      integer(int32) n
!      real(real64) fjac(ldfjac,n)
!      real(real64) fvec(n)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) N, the number of functions and variables.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(N), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer(int32) LDFJAC, the leading dimension of  FJAC.
!    LDFJAC must be at least N.
!
!    Input, real(real64) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at most
!    TOL.  TOL should be nonnegative.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) n

  real(real64) diag(n)
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fvec(n)
  integer(int32) info
!  integer(int32) j
  integer(int32) lr
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  real(real64) qtf(n)
  real(real64) r((n*(n+1))/2)
  real(real64) tol
  real(real64) x(n)
  real(real64) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  maxfev = 100 * ( n + 1 )
  xtol = tol
  mode = 2
  diag(1:n) = 1.0D+00
  factor = 100.0D+00
  nprint = 0
  lr = ( n * ( n + 1 ) ) / 2

  call hybrj ( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
    factor, nprint, info, nfev, njev, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end
subroutine lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

!*****************************************************************************80
!
!! LMDER minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer(int32) ldfjac
!      integer(int32) n
!      real(real64) fjac(ldfjac,n)
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) M, is the number of functions.
!
!    Input, integer(int32) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Input, real(real64) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real(real64) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, real(real64) GTOL.  Termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer(int32) MAXFEV.  Termination occurs when the number of
!    calls to FCN with IFLAG = 1 is at least MAXFEV by the end of an iteration.
!
!    Input/output, real(real64) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(int32) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(real64) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(int32) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of
!       squares are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer(int32) NFEV, the number of calls to FCN with
!    IFLAG = 1.
!
!    Output, integer(int32) NJEV, the number of calls to FCN with
!    IFLAG = 2.
!
!    Output, integer(int32) IPVT(N), defines a permutation matrix P
!    such that JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column
!    IPVT(J) of the identity matrix.
!
!    Output, real(real64) QTF(N), contains the first N elements of Q'*FVEC.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) actred
  real(real64) delta
  real(real64) diag(n)
  real(real64) dirder
  real(real64) epsmch
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fnorm
  real(real64) fnorm1
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gnorm
  real(real64) gtol
!  integer(int32) i
  integer(int32) iflag
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) iter
  integer(int32) j
  integer(int32) l
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  real(real64) par
  logical pivot
  real(real64) pnorm
  real(real64) prered
  real(real64) qtf(n)
  real(real64) ratio
  real(real64) sum2
  real(real64) temp
  real(real64) temp1
  real(real64) temp2
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) wa3(n)
  real(real64) wa4(m)
  real(real64) xnorm
  real(real64) x(n)
  real(real64) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    go to 300
  end if

  if ( m < n ) then
    go to 300
  end if

  if ( ldfjac < m &
    .or. ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 &
     .or. maxfev <= 0 .or. factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
  nfev = 1
  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
  do
!
!  Calculate the jacobian matrix.
!
    iflag = 2
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )

    njev = njev + 1

    if ( iflag < 0 ) then
      go to 300
    end if
!
!  If requested, call FCN to enable printing of iterates.
!
    if ( 0 < nprint ) then
      iflag = 0
      if ( mod ( iter - 1, nprint ) == 0 ) then
        call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
      end if
      if ( iflag < 0 ) then
        go to 300
      end if
    end if
!
!  Compute the QR factorization of the jacobian.
!
    pivot = .true.
    call qrfac ( m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2 )
!
!  On the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do
      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)

      xnorm = enorm ( n, wa3 )

      if ( xnorm == 0.0D+00 ) then
        delta = factor
      else
        delta = factor * xnorm
      end if

    end if
!
!  Form Q'*FVEC and store the first N components in QTF.
!
    wa4(1:m) = fvec(1:m)

    do j = 1, n

      if ( fjac(j,j) /= 0.0D+00 ) then
        sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
        temp = - sum2 / fjac(j,j)
        wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
      end if

      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)

    end do
!
!  Compute the norm of the scaled gradient.
!
    gnorm = 0.0D+00

    if ( fnorm /= 0.0D+00 ) then

      do j = 1, n
        l = ipvt(j)
        if ( wa2(l) /= 0.0D+00 ) then
          sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
          gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
        end if
      end do

    end if
!
!  Test for convergence of the gradient norm.
!
    if ( gnorm <= gtol ) then
      info = 4
      go to 300
    end if
!
!  Rescale if necessary.
!
    if ( mode /= 2 ) then
      do j = 1, n
        diag(j) = max ( diag(j), wa2(j) )
      end do
    end if
!
!  Beginning of the inner loop.
!
    do
!
!  Determine the Levenberg-Marquardt parameter.
!
      call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction p and x + p. calculate the norm of p.
!
      wa1(1:n) = - wa1(1:n)
      wa2(1:n) = x(1:n) + wa1(1:n)
      wa3(1:n) = diag(1:n) * wa1(1:n)

      pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
      if ( iter == 1 ) then
        delta = min ( delta, pnorm )
      end if
!
!  Evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      call fcn ( m, n, wa2, wa4, fjac, ldfjac, iflag )

      nfev = nfev + 1

      if ( iflag < 0 ) then
        go to 300
      end if

      fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
      if ( 0.1D+00 * fnorm1 < fnorm ) then
        actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
      else
        actred = - 1.0D+00
      end if
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
      do j = 1, n
        wa3(j) = 0.0D+00
        l = ipvt(j)
        temp = wa1(l)
        wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
      end do

      temp1 = enorm ( n, wa3 ) / fnorm
      temp2 = ( sqrt ( par ) * pnorm ) / fnorm
      prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
      dirder = - ( temp1 ** 2 + temp2 ** 2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
      if ( prered /= 0.0D+00 ) then
        ratio = actred / prered
      else
        ratio = 0.0D+00
      end if
!
!  Update the step bound.
!
      if ( ratio <= 0.25D+00 ) then

        if ( 0.0D+00 <= actred ) then
          temp = 0.5D+00
        end if

        if ( actred < 0.0D+00 ) then
          temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
        end if

        if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
          temp = 0.1D+00
        end if

        delta = temp * min ( delta, pnorm / 0.1D+00 )
        par = par / temp

      else

        if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
          delta = 2.0D+00 * pnorm
          par = 0.5D+00 * par
        end if

      end if
!
!  Successful iteration.
!
!  Update X, FVEC, and their norms.
!
      if ( 0.0001D+00 <= ratio ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:m) = wa4(1:m)
        xnorm = enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
      end if
!
!  Tests for convergence.
!
      if ( abs ( actred) <= ftol .and. &
        prered <= ftol .and. &
        0.5D+00 * ratio <= 1.0D+00 ) then
        info = 1
      end if

      if ( delta <= xtol * xnorm ) then
        info = 2
      end if

      if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) then
        info = 3
      end if

      if ( info /= 0 ) then
        go to 300
      end if
!
!  Tests for termination and stringent tolerances.
!
      if ( nfev >= maxfev ) then
        info = 5
      end if

      if ( abs ( actred ) <= epsmch .and. prered <= epsmch &
        .and. 0.5D+00 * ratio <= 1.0D+00 ) then
        info = 6
      end if

      if ( delta <= epsmch * xnorm ) then
        info = 7
      end if

      if ( gnorm <= epsmch ) then
        info = 8
      end if

      if ( info /= 0 ) then
        go to 300
      end if
!
!  End of the inner loop. repeat if iteration unsuccessful.
!
      if ( 0.0001D+00 <= ratio ) then
        exit
      end if

    end do
!
!  End of the outer loop.
!
  end do

  300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
  end if

  return
end
subroutine lmder1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! LMDER1 minimizes M functions in N variables by Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDER.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer(int32) ldfjac
!      integer(int32) n
!      real(real64) fjac(ldfjac,n)
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) M, the number of functions.
!
!    Input, integer(int32) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer(int32) LDFJAC, is the leading dimension of FJAC,
!    which must be no less than M.
!
!    Input, real(real64) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) diag(n)
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gtol
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  real(real64) qtf(n)
  real(real64) tol
  real(real64) x(n)
  real(real64) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < m ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 100 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  mode = 1
  nprint = 0

  call lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
subroutine lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
  diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

!*****************************************************************************80
!
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer(int32) m
!      integer(int32) n
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) M, the number of functions.
!
!    Input, integer(int32) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Input, real(real64) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real(real64) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  Therefore, XTOL
!    measures the relative error desired in the approximate solution.  XTOL
!    should be nonnegative.
!
!    Input, real(real64) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer(int32) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, real(real64) EPSFCN, is used in determining a suitable step 
!    length for the forward-difference approximation.  This approximation 
!    assumes that the relative errors in the functions are of the order of 
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real(real64) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(int32) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(real64) FACTOR, determines the initial step bound.  
!    This bound is set to the product of FACTOR and the euclidean norm of
!    DIAG*X if nonzero, or else to FACTOR itself.  In most cases, FACTOR should 
!    lie in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(int32) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN has reached or exceeded MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer(int32) NFEV, the number of calls to FCN.
!
!    Output, real(real64) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Output, integer(int32) IPVT(N), defines a permutation matrix P such
!    that JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
!    of the identity matrix.
!
!    Output, real(real64) QTF(N), the first N elements of Q'*FVEC.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) actred
  real(real64) delta
  real(real64) diag(n)
  real(real64) dirder
  real(real64) epsfcn
  real(real64) epsmch
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fnorm
  real(real64) fnorm1
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gnorm
  real(real64) gtol
  integer(int32) i
  integer(int32) iflag
  integer(int32) iter
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) j
  integer(int32) l
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) nprint
  real(real64) par
  logical pivot
  real(real64) pnorm
  real(real64) prered
  real(real64) qtf(n)
  real(real64) ratio
  real(real64) sum2
  real(real64) temp
  real(real64) temp1
  real(real64) temp2
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) wa3(n)
  real(real64) wa4(m)
  real(real64) x(n)
  real(real64) xnorm
  real(real64) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0

  if ( n <= 0 ) then
    go to 300
  else if ( m < n ) then
    go to 300
  else if ( ldfjac < m ) then
    go to 300
  else if ( ftol < 0.0D+00 ) then
    go to 300
  else if ( xtol < 0.0D+00 ) then
    go to 300
  else if ( gtol < 0.0D+00 ) then
    go to 300
  else if ( maxfev <= 0 ) then
    go to 300
  else if ( factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( m, n, x, fvec, iflag )
  nfev = 1

  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
30 continue
!
!  Calculate the jacobian matrix.
!
  iflag = 2
  call fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )
  nfev = nfev + n

  if ( iflag < 0 ) then
    go to 300
  end if
!
!  If requested, call FCN to enable printing of iterates.
!
  if ( 0 < nprint ) then
    iflag = 0
    if ( mod ( iter - 1, nprint ) == 0 ) then
      call fcn ( m, n, x, fvec, iflag )
    end if
    if ( iflag < 0 ) then
      go to 300
    end if
  end if
!
!  Compute the QR factorization of the jacobian.
!
  pivot = .true.
  call qrfac ( m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2 )
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if ( iter == 1 ) then

       if ( mode /= 2 ) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if ( wa2(j) == 0.0D+00 ) then
             diag(j) = 1.0D+00
           end if
         end do
       end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm ( n, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) then
         delta = factor
       end if
     end if
!
!  Form Q' * FVEC and store the first N components in QTF.
!
     wa4(1:m) = fvec(1:m)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then
         sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
         temp = - sum2 / fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       end if

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     end do
!
!  Compute the norm of the scaled gradient.
!
     gnorm = 0.0D+00

     if ( fnorm /= 0.0D+00 ) then

       do j = 1, n

         l = ipvt(j)

         if ( wa2(l) /= 0.0D+00 ) then
           sum2 = 0.0D+00
           do i = 1, j
             sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
           end do
           gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
         end if

       end do

     end if
!
!  Test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 300
     end if
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
200  continue
!
!  Determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( m, n, wa2, wa4, iflag )
        nfev = nfev + 1
        if ( iflag < 0 ) then
          go to 300
        end if
        fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
        else
          actred = -1.0D+00
        end if
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
!
        do j = 1, n
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( sqrt ( par ) * pnorm ) / fnorm
        prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
        dirder = - ( temp1 ** 2 + temp2 ** 2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        if ( prered /= 0.0D+00 ) then
          ratio = actred / prered
        end if
!
!  Update the step bound.
!
        if ( ratio <= 0.25D+00 ) then

           if ( actred >= 0.0D+00 ) then
             temp = 0.5D+00
           endif

           if ( actred < 0.0D+00 ) then
             temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
           end if

           if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
             temp = 0.1D+00
           end if

           delta = temp * min ( delta, pnorm / 0.1D+00  )
           par = par / temp

        else

           if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
             delta = 2.0D+00 * pnorm
             par = 0.5D+00 * par
           end if

        end if
!
!  Test for successful iteration.
!

!
!  Successful iteration. update X, FVEC, and their norms.
!
        if ( 0.0001D+00 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Tests for convergence.
!
        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 1
        end if

        if ( delta <= xtol * xnorm ) then
          info = 2
        end if

        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) info = 3

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Tests for termination and stringent tolerances.
!
        if ( maxfev <= nfev ) then
          info = 5
        end if

        if ( abs ( actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 6
        end if

        if ( delta <= epsmch * xnorm ) then
          info = 7
        end if

        if ( gnorm <= epsmch ) then
          info = 8
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
        if ( ratio < 0.0001D+00 ) then
          go to 200
        end if
!
!  End of the outer loop.
!
     go to 30

300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( m, n, x, fvec, iflag )
  end if

  return
end
subroutine lmdif1 ( fcn, m, n, x, fvec, tol, info )

!*****************************************************************************80
!
!! LMDIF1 minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDIF.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer(int32) n
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer(int32) M, the number of functions.
!
!    Input, integer(int32) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Input, real(real64) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN has reached or exceeded 200*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer(int32) m
  integer(int32) n

  real(real64) diag(n)
  real(real64) epsfcn
  real(real64) factor
  external fcn
  real(real64) fjac(m,n)
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gtol
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) ldfjac
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) nprint
  real(real64) qtf(n)
  real(real64) tol
  real(real64) x(n)
  real(real64) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 200 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  epsfcn = 0.0D+00
  mode = 1
  nprint = 0
  ldfjac = m

  call lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
    diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This function completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  
!
!    If, however, the limit of 10 iterations is reached, then the output 
!    PAR will contain the best value obtained so far.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of R.
!
!    Input/output, real(real64) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer(int32) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer(int32) IPVT(N), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the 
!    identity matrix.
!
!    Input, real(real64) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real(real64) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real(real64) DELTA, an upper bound on the euclidean norm
!    of D*X.  DELTA should be positive.
!
!    Input/output, real(real64) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real(real64) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real(real64) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer(int32) ldr
  integer(int32) n

  real(real64) delta
  real(real64) diag(n)
  real(real64) dwarf
  real(real64) dxnorm
  real(real64) gnorm
  real(real64) fp
!  integer(int32) i
  integer(int32) ipvt(n)
  integer(int32) iter
  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) nsing
  real(real64) par
  real(real64) parc
  real(real64) parl
  real(real64) paru
!  real(real64) qnorm
  real(real64) qtb(n)
  real(real64) r(ldr,n)
  real(real64) sdiag(n)
  real(real64) sum2
  real(real64) temp
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) x(n)
!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if
    if ( nsing < n ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    if ( iter == 0 ) then
      par = 0.0D+00
    end if
    return
  end if
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  if ( n <= nsing ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  end if
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  do j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta

  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
!
!  Beginning of an iteration.
!
  do
 
    iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
    if ( par == 0.0D+00 ) then
      par = max ( dwarf, 0.001D+00 * paru )
    end if

    wa1(1:n) = sqrt ( par ) * diag(1:n)

    call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm ( n, wa2 )
    temp = fp
    fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      exit
    end if
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      exit
    else if ( iter == 10 ) then
      exit
    end if
!
!  Compute the Newton correction.
!
    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      wa1(j) = wa1(j) / sdiag(j)
      temp = wa1(j)
      wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
    end do

    temp = enorm ( n, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
    if ( 0.0D+00 < fp ) then
      parl = max ( parl, par )
    else if ( fp < 0.0D+00 ) then
      paru = min ( paru, par )
    end if
!
!  Compute an improved estimate for PAR.
!
    par = max ( parl, par + parc )
!
!  End of an iteration.
!
  end do
!
!  Termination.
!
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end
subroutine lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

!*****************************************************************************80
!
!! LMSTR minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMSTR minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm
!    which uses minimal storage.
!
!    The user must provide a subroutine which calculates the functions and
!    the rows of the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the rows of the jacobian.
!    FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
!      integer(int32) m
!      integer(int32) n
!      real(real64) fjrow(n)
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
!    the jacobian at X, and return this vector in FJROW.
!    To terminate the algorithm, set the output value of IFLAG negative.
!
!    Input, integer(int32) M, the number of functions.
!
!    Input, integer(int32) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an N by N array.  The upper
!    triangle of FJAC contains an upper triangular matrix R such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    triangular part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Input, real(real64) FTOL.  Termination occurs when both the actual and
!    predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real(real64) XTOL.  Termination occurs when the relative error 
!    between two consecutive iterates is at most XTOL.  XTOL should be 
!    nonnegative.
!
!    Input, real(real64) GTOL. termination occurs when the cosine of the 
!    angle between FVEC and any column of the jacobian is at most GTOL in 
!    absolute value.  Therefore, GTOL measures the orthogonality desired 
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer(int32) MAXFEV.  Termination occurs when the number 
!    of calls to FCN with IFLAG = 1 is at least MAXFEV by the end of 
!    an iteration.
!
!    Input/output, real(real64) DIAG(N).  If MODE = 1, then DIAG is set 
!    internally.  If MODE = 2, then DIAG must contain positive entries that 
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(int32) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(real64) FACTOR, determines the initial step bound.  This 
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(int32) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the 
!    description of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer(int32) NFEV, the number of calls to FCN 
!    with IFLAG = 1.
!
!    Output, integer(int32) NJEV, the number of calls to FCN 
!    with IFLAG = 2.
!
!    Output, integer(int32) IPVT(N), defines a permutation matrix P such
!    that JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular.
!    Column J of P is column IPVT(J) of the identity matrix.
!
!    Output, real(real64) QTF(N), contains the first N elements of Q'*FVEC.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) actred
  real(real64) delta
  real(real64) diag(n)
  real(real64) dirder
  real(real64) epsmch
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) fnorm
  real(real64) fnorm1
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gnorm
  real(real64) gtol
  integer(int32) i
  integer(int32) iflag
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) iter
  integer(int32) j
  integer(int32) l
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  real(real64) par
  logical pivot
  real(real64) pnorm
  real(real64) prered
  real(real64) qtf(n)
  real(real64) ratio
  logical sing
  real(real64) sum2
  real(real64) temp
  real(real64) temp1
  real(real64) temp2
  real(real64) wa1(n)
  real(real64) wa2(n)
  real(real64) wa3(n)
  real(real64) wa4(m)
  real(real64) x(n)
  real(real64) xnorm
  real(real64) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    go to 340
  else if ( m < n ) then
    go to 340
  else if ( ldfjac < n ) then
    go to 340
  else if ( ftol < 0.0D+00 ) then
    go to 340
  else if ( xtol < 0.0D+00 ) then
    go to 340
  else if ( gtol < 0.0D+00 ) then
    go to 340
  else if ( maxfev <= 0 ) then
    go to 340
  else if ( factor <= 0.0D+00 ) then
    go to 340
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 340
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( m, n, x, fvec, wa3, iflag )
  nfev = 1

  if ( iflag < 0 ) then
    go to 340
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
   30 continue
!
!  If requested, call FCN to enable printing of iterates.
!
     if ( 0 < nprint ) then
       iflag = 0
       if ( mod ( iter-1, nprint ) == 0 ) then
         call fcn ( m, n, x, fvec, wa3, iflag )
       end if
       if ( iflag < 0 ) then
         go to 340
       end if
     end if
!
!  Compute the QR factorization of the jacobian matrix calculated one row
!  at a time, while simultaneously forming Q'* FVEC and storing
!  the first N components in QTF.
!
     qtf(1:n) = 0.0D+00
     fjac(1:n,1:n) = 0.0D+00
     iflag = 2

     do i = 1, m
       call fcn ( m, n, x, fvec, wa3, iflag )
       if ( iflag < 0 ) then
         go to 340
       end if
       temp = fvec(i)
       call rwupdt ( n, fjac, ldfjac, wa3, qtf, temp, wa1, wa2 )
       iflag = iflag + 1
     end do

     njev = njev + 1
!
!  If the jacobian is rank deficient, call QRFAC to
!  reorder its columns and update the components of QTF.
!
     sing = .false.
     do j = 1, n
       if ( fjac(j,j) == 0.0D+00 ) then
         sing = .true.
       end if
       ipvt(j) = j
       wa2(j) = enorm ( j, fjac(1,j) )
     end do

     if ( sing ) then

       pivot = .true.
       call qrfac ( n, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2 )

       do j = 1, n

         if ( fjac(j,j) /= 0.0D+00 ) then

           sum2 = dot_product ( qtf(j:n), fjac(j:n,j) )
           temp = - sum2 / fjac(j,j)
           qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp

         end if

         fjac(j,j) = wa1(j)

       end do

     end if
!
!  On the first iteration
!    if mode is 1,
!      scale according to the norms of the columns of the initial jacobian.
!    calculate the norm of the scaled X,
!    initialize the step bound delta.
!
     if ( iter == 1 ) then

       if ( mode /= 2 ) then

         diag(1:n) = wa2(1:n)
         do j = 1, n
           if ( wa2(j) == 0.0D+00 ) then
             diag(j) = 1.0D+00
           end if
         end do

       end if

       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm ( n, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) then
         delta = factor
       end if

     end if
!
!  Compute the norm of the scaled gradient.
!
     gnorm = 0.0D+00

     if ( fnorm /= 0.0D+00 ) then

       do j = 1, n
         l = ipvt(j)
         if ( wa2(l) /= 0.0D+00 ) then
           sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
           gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
         end if
       end do

     end if
!
!  Test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 340
     end if
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
240    continue
!
!  Determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)
        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( m, n, wa2, wa4, wa3, iflag )
        nfev = nfev + 1
        if ( iflag < 0 ) then
          go to 340
        end if
        fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
        else
          actred = -1.0D+00
        end if
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
        do j = 1, n
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( sqrt(par) * pnorm ) / fnorm
        prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
        dirder = - ( temp1 ** 2 + temp2 ** 2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
        if ( prered /= 0.0D+00 ) then
          ratio = actred / prered
        else
          ratio = 0.0D+00
        end if
!
!  Update the step bound.
!
        if ( ratio <= 0.25D+00 ) then

          if ( actred >= 0.0D+00 ) then
            temp = 0.5D+00
          else
            temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
          end if

          if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
            temp = 0.1D+00
          end if

          delta = temp * min ( delta, pnorm / 0.1D+00 )
          par = par / temp

        else

          if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
            delta = pnorm / 0.5D+00
            par = 0.5D+00 * par
          end if

        end if
!
!  Test for successful iteration.
!
        if ( ratio >= 0.0001D+00 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Tests for convergence, termination and stringent tolerances.
!
        if ( abs ( actred ) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 1
        end if

        if ( delta <= xtol * xnorm ) then
          info = 2
        end if

        if ( abs ( actred ) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) then
          info = 3
        end if

        if ( info /= 0 ) then
          go to 340
        end if

        if ( nfev >= maxfev) then
          info = 5
        end if

        if ( abs ( actred ) <= epsmch .and. prered <= epsmch &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 6
        end if

        if ( delta <= epsmch * xnorm ) then
          info = 7
        end if

        if ( gnorm <= epsmch ) then
          info = 8
        end if

        if ( info /= 0 ) then
          go to 340
        end if
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
        if ( ratio < 0.0001D+00 ) then
          go to 240
        end if
!
!  End of the outer loop.
!
     go to 30

  340 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( m, n, x, fvec, wa3, iflag )
  end if

  return
end
subroutine lmstr1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! LMSTR1 minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMSTR1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm
!    which uses minimal storage.
!
!    This is done by using the more general least-squares solver
!    LMSTR.  The user must provide a subroutine which calculates
!    the functions and the rows of the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    19 August 2016
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the rows of the jacobian.
!    FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
!      integer(int32) m
!      integer(int32) n
!      real(real64) fjrow(n)
!      real(real64) fvec(m)
!      integer(int32) iflag
!      real(real64) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
!    the jacobian at X, and return this vector in FJROW.
!    To terminate the algorithm, set the output value of IFLAG negative.
!
!    Input, integer(int32) M, the number of functions.
!
!    Input, integer(int32) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real(real64) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(real64) FVEC(M), the functions evaluated at the output X.
!
!    Output, real(real64) FJAC(LDFJAC,N), an N by N array.  The upper
!    triangle contains an upper triangular matrix R such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower triangular part of FJAC contains information generated
!    during the computation of R.
!
!    Input, integer(int32) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Input, real(real64) TOL. Termination occurs when the algorithm 
!    estimates either that the relative error in the sum of squares is at 
!    most TOL or that the relative error between X and the solution is at 
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer(int32) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer(int32) ldfjac
  integer(int32) m
  integer(int32) n

  real(real64) diag(n)
  real(real64) factor
  external fcn
  real(real64) fjac(ldfjac,n)
  real(real64) ftol
  real(real64) fvec(m)
  real(real64) gtol
  integer(int32) info
  integer(int32) ipvt(n)
  integer(int32) maxfev
  integer(int32) mode
  integer(int32) nfev
  integer(int32) njev
  integer(int32) nprint
  real(real64) qtf(n)
  real(real64) tol
  real(real64) x(n)
  real(real64) xtol

  if ( n <= 0 ) then
    info = 0
    return
  end if

  if ( m < n ) then
    info = 0
    return
  end if

  if ( ldfjac < n ) then
    info = 0
    return
  end if

  if ( tol < 0.0D+00 ) then
    info = 0
    return
  end if

  fvec(1:n) = 0.0D+00
  fjac(1:ldfjac,1:n) = 0.0D+00
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  maxfev = 100 * ( n + 1 )
  diag(1:n) = 0.0D+00
  mode = 1
  factor = 100.0D+00
  nprint = 0
  info = 0
  nfev = 0
  njev = 0
  ipvt(1:n) = 0
  qtf(1:n) = 0.0D+00

  call lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
subroutine qform ( m, n, q, ldq )

!*****************************************************************************80
!
!! QFORM produces the explicit QR factorization of a matrix.
!
!  Discussion:
!
!    The QR factorization of a matrix is usually accumulated in implicit
!    form, that is, as a series of orthogonal transformations of the
!    original matrix.  This routine carries out those transformations,
!    to explicitly exhibit the factorization constructed by QRFAC.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) M, is a positive integer input variable set
!    to the number of rows of A and the order of Q.
!
!    Input, integer(int32) N, is a positive integer input variable set
!    to the number of columns of A.
!
!    Input/output, real(real64) Q(LDQ,M).  Q is an M by M array.
!    On input the full lower trapezoid in the first min(M,N) columns of Q
!    contains the factored form.
!    On output, Q has been accumulated into a square matrix.
!
!    Input, integer(int32) LDQ, is a positive integer input variable 
!    not less than M which specifies the leading dimension of the array Q.
!
  implicit none

  integer(int32) ldq
  integer(int32) m
  integer(int32) n

  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) minmn
  real(real64) q(ldq,m)
  real(real64) temp
  real(real64) wa(m)

  minmn = min ( m, n )

  do j = 2, minmn
    q(1:j-1,j) = 0.0D+00
  end do
!
!  Initialize remaining columns to those of the identity matrix.
!
  q(1:m,n+1:m) = 0.0D+00

  do j = n + 1, m
    q(j,j) = 1.0D+00
  end do
!
!  Accumulate Q from its factored form.
!
  do l = 1, minmn

    k = minmn - l + 1

    wa(k:m) = q(k:m,k)

    q(k:m,k) = 0.0D+00
    q(k,k) = 1.0D+00

    if ( wa(k) /= 0.0D+00 ) then

      do j = k, m
        temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
        q(k:m,j) = q(k:m,j) - temp * wa(k:m)
      end do

    end if

  end do

  return
end
subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This function uses Householder transformations with optional column
!    pivoting to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  
!
!    The Householder transformation for column K, K = 1,2,...,min(M,N), 
!    is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  
!
!    The form of this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK routine.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) M, the number of rows of A.
!
!    Input, integer(int32) N, the number of columns of A.
!
!    Input/output, real(real64) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q, the non-trivial elements of
!    the U vectors described above.
!
!    Input, integer(int32) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer(int32) IPVT(LIPVT), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer(int32) LIPVT, the dimension of IPVT, which should 
!    be N if pivoting is used.
!
!    Output, real(real64) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real(real64) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer(int32) lda
  integer(int32) lipvt
  integer(int32) m
  integer(int32) n

  real(real64) a(lda,n)
  real(real64) acnorm(n)
  real(real64) ajnorm
  real(real64) epsmch
!  integer(int32) i
  integer(int32) i4_temp
  integer(int32) ipvt(lipvt)
  integer(int32) j
  integer(int32) k
  integer(int32) kmax
  integer(int32) minmn
  logical pivot
  real(real64) r8_temp(m)
  real(real64) rdiag(n)
  real(real64) temp
  real(real64) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(kmax) < rdiag(k) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j + 1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp ** 2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = - ajnorm

  end do

  return
end
subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This function completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of R.
!
!    Input/output, real(real64) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer(int32) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer(int32) IPVT(N), defines the permutation matrix P such 
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real(real64) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real(real64) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real(real64) X(N), the least squares solution.
!
!    Output, real(real64) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer(int32) ldr
  integer(int32) n

  real(real64) c
  real(real64) cotan
  real(real64) diag(n)
  integer(int32) i
  integer(int32) ipvt(n)
  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) nsing
  real(real64) qtb(n)
  real(real64) qtbpj
  real(real64) r(ldr,n)
  real(real64) s
  real(real64) sdiag(n)
  real(real64) sum2
  real(real64) t
  real(real64) temp
  real(real64) wa(n)
  real(real64) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t ** 2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k + 1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end
subroutine r1mpyq ( m, n, a, lda, v, w )

!*****************************************************************************80
!
!! R1MPYQ computes A*Q, where Q is the product of Householder transformations.
!
!  Discussion:
!
!    Given an M by N matrix A, this function computes A*Q where
!    Q is the product of 2*(N - 1) transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!    eliminate elements in the I-th and N-th planes, respectively.
!    Q itself is not given, rather the information to recover the
!    GV, GW rotations is supplied.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) M, the number of rows of A.
!
!    Input, integer(int32) N, the number of columns of A.
!
!    Input/output, real(real64) A(LDA,N), the M by N array.
!    On input, the matrix A to be postmultiplied by the orthogonal matrix Q.
!    On output, the value of A*Q.
!
!    Input, integer(int32) LDA, the leading dimension of A, which must not
!    be less than M.
!
!    Input, real(real64) V(N), W(N), contain the information necessary
!    to recover the Givens rotations GV and GW.
!
  implicit none

  integer(int32) lda
  integer(int32) m
  integer(int32) n

  real(real64) a(lda,n)
  real(real64) c
  integer(int32) i
  integer(int32) j
  real(real64) s
  real(real64) temp
  real(real64) v(n)
  real(real64) w(n)
!
!  Apply the first set of Givens rotations to A.
!
  do j = n - 1, 1, -1

     if ( 1.0D+00 < abs ( v(j) ) ) then
       c = 1.0D+00 / v(j)
       s = sqrt ( 1.0D+00 - c ** 2 )
     else
       s = v(j)
       c = sqrt ( 1.0D+00 - s ** 2 )
     end if

     do i = 1, m
        temp =   c * a(i,j) - s * a(i,n)
        a(i,n) = s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do
!
!  Apply the second set of Givens rotations to A.
!
  do j = 1, n - 1

     if ( abs ( w(j) ) > 1.0D+00 ) then
       c = 1.0D+00 / w(j)
       s = sqrt ( 1.0D+00 - c ** 2 )
     else
       s = w(j)
       c = sqrt ( 1.0D+00 - s ** 2 )
     end if

     do i = 1, m
        temp =     c * a(i,j) + s * a(i,n)
        a(i,n) = - s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do

  return
end
subroutine r1updt ( m, n, s, ls, u, v, w, sing )

!*****************************************************************************80
!
!! R1UPDT re-triangularizes a matrix after a rank one update.
!
!  Discussion:
!
!    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
!    N-vector V, the problem is to determine an orthogonal matrix Q such that
!
!      (S + U * V' ) * Q
!
!    is again lower trapezoidal.
!
!    This function determines Q as the product of 2 * (N - 1)
!    transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    where GV(I), GW(I) are Givens rotations in the (I,N) plane
!    which eliminate elements in the I-th and N-th planes,
!    respectively.  Q itself is not accumulated, rather the
!    information to recover the GV and GW rotations is returned.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) M, the number of rows of S.
!
!    Input, integer(int32) N, the number of columns of S.  
!    N must not exceed M.
!
!    Input/output, real(real64) S(LS).  On input, the lower trapezoidal
!    matrix S stored by columns.  On output S contains the lower trapezoidal
!    matrix produced as described above.
!
!    Input, integer(int32) LS, the length of the S array.  LS must be at
!    least (N*(2*M-N+1))/2.
!
!    Input, real(real64) U(M), the U vector.
!
!    Input/output, real(real64) V(N).  On input, V must contain the 
!    vector V.  On output V contains the information necessary to recover the
!    Givens rotations GV described above.
!
!    Output, real(real64) W(M), contains information necessary to
!    recover the Givens rotations GW described above.
!
!    Output, logical SING, is set to TRUE if any of the diagonal elements
!    of the output S are zero.  Otherwise SING is set FALSE.
!
  implicit none

  integer(int32) ls
  integer(int32) m
  integer(int32) n

  real(real64) cos
  real(real64) cotan
  real(real64) giant
  integer(int32) i
  integer(int32) j
  integer(int32) jj
  integer(int32) l
  real(real64) s(ls)
  real(real64) sin
  logical sing
  real(real64) tan
  real(real64) tau
  real(real64) temp
  real(real64) u(m)
  real(real64) v(n)
  real(real64) w(m)
!
!  GIANT is the largest magnitude.
!
  giant = huge ( giant )
!
!  Initialize the diagonal element pointer.
!
  jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
!
!  Move the nontrivial part of the last column of S into W.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
!  Rotate the vector V into a multiple of the N-th unit vector
!  in such a way that a spike is introduced into W.
!
  do j = n - 1, 1, -1

    jj = jj - ( m - j + 1 )
    w(j) = 0.0D+00

    if ( v(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of V.
!
      if ( abs ( v(n) ) < abs ( v(j) ) ) then
        cotan = v(n) / v(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
        cos = sin * cotan
        tau = 1.0D+00
        if ( abs ( cos ) * giant > 1.0D+00 ) then
          tau = 1.0D+00 / cos
        end if
      else
        tan = v(j) / v(n)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
        sin = cos * tan
        tau = sin
      end if
!
!  Apply the transformation to V and store the information
!  necessary to recover the Givens rotation.
!
      v(n) = sin * v(j) + cos * v(n)
      v(j) = tau
!
!  Apply the transformation to S and extend the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) - sin * w(i)
        w(i) = sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do

    end if

  end do
!
!  Add the spike from the rank 1 update to W.
!
   w(1:m) = w(1:m) + v(n) * u(1:m)
!
!  Eliminate the spike.
!
  sing = .false.

  do j = 1, n-1

    if ( w(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of the spike.
!
      if ( abs ( s(jj) ) < abs ( w(j) ) ) then

        cotan = s(jj) / w(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
        cos = sin * cotan

        if ( 1.0D+00 < abs ( cos ) * giant ) then
          tau = 1.0D+00 / cos
        else
          tau = 1.0D+00
        end if

      else

        tan = w(j) / s(jj)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
        sin = cos * tan
        tau = sin

      end if
!
!  Apply the transformation to S and reduce the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) + sin * w(i)
        w(i) = - sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do
!
!  Store the information necessary to recover the Givens rotation.
!
      w(j) = tau

    end if
!
!  Test for zero diagonal elements in the output S.
!
    if ( s(jj) == 0.0D+00 ) then
      sing = .true.
    end if

    jj = jj + ( m - j + 1 )

  end do
!
!  Move W back into the last column of the output S.
!
  l = jj
  do i = n, m
    s(l) = w(i)
    l = l + 1
  end do

  if ( s(jj) == 0.0D+00 ) then
    sing = .true.
  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real(real64) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer(int32) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real(real64) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer(int32), parameter :: i4_huge = 2147483647
  integer(int32) k
  real(real64) r8_uniform_01
  integer(int32) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, real64 ) * 4.656612875D-10

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(int32) M, the number of rows in A.
!
!    Input, integer(int32) N, the number of columns in A.
!
!    Input, real(real64) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer(int32) m
  integer(int32) n

  real(real64) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(int32) M, N, the number of rows and columns.
!
!    Input, real(real64) A(M,N), an M by N matrix to be printed.
!
!    Input, integer(int32) ILO, JLO, the first row and column to print.
!
!    Input, integer(int32) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer(int32), parameter :: incx = 5
  integer(int32) m
  integer(int32) n

  real(real64) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer(int32) i
  integer(int32) i2hi
  integer(int32) i2lo
  integer(int32) ihi
  integer(int32) ilo
  integer(int32) inc
  integer(int32) j
  integer(int32) j2
  integer(int32) j2hi
  integer(int32) j2lo
  integer(int32) jhi
  integer(int32) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), real64 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(int32) N, the number of components of the vector.
!
!    Input, real(real64) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer(int32) n

  real(real64) a(n)
  integer(int32) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine rwupdt ( n, r, ldr, w, b, alpha, c, s )

!*****************************************************************************80
!
!! RWUPDT computes the decomposition of triangular matrix augmented by one row.
!
!  Discussion:
!
!    Given an N by N upper triangular matrix R, this function
!    computes the QR decomposition of the matrix formed when a row
!    is added to R.  If the row is specified by the vector W, then
!    RWUPDT determines an orthogonal matrix Q such that when the
!    N+1 by N matrix composed of R augmented by W is premultiplied
!    by Q', the resulting matrix is upper trapezoidal.
!
!    The matrix Q' is the product of N transformations
!
!      G(N)*G(N-1)* ... *G(1)
!
!    where G(I) is a Givens rotation in the (I,N+1) plane which eliminates
!    elements in the (N+1)-st plane.  RWUPDT also computes the product
!    Q'*C where C is the (N+1)-vector (B,ALPHA).  
!
!    Q itself is not accumulated, rather the information to recover the G 
!    rotations is supplied.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of R.
!
!    Input/output, real(real64) R(LDR,N).  On input the upper triangular
!    part of R must contain the matrix to be updated.  On output R contains the
!    updated triangular matrix.
!
!    Input, integer(int32) LDR, the leading dimension of the array R.
!    LDR must not be less than N.
!
!    Input, real(real64) W(N), the row vector to be added to R.
!
!    Input/output, real(real64) B(N).  On input, the first N elements
!    of the vector C.  On output the first N elements of the vector Q'*C.
!
!    Input/output, real(real64) ALPHA.  On input, the (N+1)-st element
!    of the vector C.  On output the (N+1)-st element of the vector Q'*C.
!
!    Output, real(real64) C(N), S(N), the cosines and sines of the
!    transforming Givens rotations.
!
  implicit none

  integer(int32) ldr
  integer(int32) n

  real(real64) alpha
  real(real64) b(n)
  real(real64) c(n)
  real(real64) cotan
  integer(int32) i
  integer(int32) j
  real(real64) r(ldr,n)
  real(real64) rowj
  real(real64) s(n)
  real(real64) tan
  real(real64) temp
  real(real64) w(n)

  do j = 1, n

    rowj = w(j)
!
!  Apply the previous transformations to R(I,J), I=1,2,...,J-1, and to W(J).
!
    do i = 1, j - 1
      temp =   c(i) * r(i,j) + s(i) * rowj
      rowj = - s(i) * r(i,j) + c(i) * rowj
      r(i,j) = temp
    end do
!
!  Determine a Givens rotation which eliminates W(J).
!
    c(j) = 1.0D+00
    s(j) = 0.0D+00

    if ( rowj /= 0.0D+00 ) then

      if ( abs ( r(j,j) ) < abs ( rowj ) ) then
        cotan = r(j,j) / rowj
        s(j) = 0.5D+00 / sqrt ( 0.25d+00 + 0.25D+00 * cotan ** 2 )
        c(j) = s(j) * cotan
      else
        tan = rowj / r(j,j)
        c(j) = 0.5D+00 / sqrt ( 0.25d+00 + 0.25D+00 * tan ** 2 )
        s(j) = c(j) * tan
      end if
!
!  Apply the current transformation to R(J,J), B(J), and ALPHA.
!
      r(j,j) =  c(j) * r(j,j) + s(j) * rowj
      temp =    c(j) * b(j)   + s(j) * alpha
      alpha = - s(j) * b(j)   + c(j) * alpha
      b(j) = temp

    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer(int32) d
  integer(int32) h
  integer(int32) m
  integer(int32) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer(int32) n
  integer(int32) s
  integer(int32) values(8)
  integer(int32) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end


!!!!!! https://people.sc.fsu.edu/~jburkardt/f_src/asa007/asa007.f90

subroutine cholesky ( a, n, nn, u, nullty, ifault )

!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!    The missing initialization "II = 0" has been added to the code.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Michael Healy.
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix
!    stored by rows in lower triangular form as a one dimensional array,
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer ( kind = 4 ) NN, the dimension of A, (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, if NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(nn)
  real ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
!  real ( kind = 8 ) rsq
  real ( kind = 8 ) u(nn)
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end
subroutine syminv ( a, n, c, w, nullty, ifault )

!*****************************************************************************80
!
!! SYMINV computes the inverse of a symmetric matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 7:
!    Inversion of a Positive Semi-Definite Symmetric Matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 198-199.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored
!    by rows in lower triangular form as a one dimensional array, in the 
!    sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
!    inverse if A is singular, stored using the same storage scheme employed
!    for A.  The program is written in such a way that A and U can share 
!    storage.
!
!    Workspace, real ( kind = 8 ) W(N).
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY 
!    is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, N < 1.
!    2, A is not positive semi-definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) c((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mdiag
  integer ( kind = 4 ) ndiag
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x

  ifault = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  nrow = n
!
!  Compute the Cholesky factorization of A.
!  The result is stored in C.
!
  nn = ( n * ( n + 1 ) ) / 2

  call cholesky ( a, n, nn, c, nullty, ifault )

  if ( ifault /= 0 ) then
    return
  end if
!
!  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
!  of C, row by row starting with the last row.
!  IROW = the row number,
!  NDIAG = location of last element in the row.
!
  irow = nrow
  ndiag = nn

  do
!
!  Special case, zero diagonal element.
!
    if ( c(ndiag) == 0.0D+00 ) then

      l = ndiag
      do j = irow, nrow
        c(l) = 0.0D+00
        l = l + j
      end do

    else

      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do

      icol = nrow
      jcol = nn
      mdiag = nn

      do

        l = jcol

        if ( icol == irow ) then
          x = 1.0D+00 / w(irow)
        else
          x = 0.0D+00
        end if

        k = nrow

        do while ( irow < k )

          x = x - w(k) * c(l)
          k = k - 1
          l = l - 1

          if ( mdiag < l ) then
            l = l - k + 1
          end if

        end do

        c(l) = x / w(irow)

        if ( icol <= irow ) then
          exit
        end if

        mdiag = mdiag - icol
        icol = icol - 1
        jcol = jcol - 1

      end do

    end if

    ndiag = ndiag - irow
    irow = irow - 1

    if ( irow <= 0 ) then
      exit
    end if

  end do

  return
end

end module
