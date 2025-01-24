

MODULE utils
  use globals
  use ieps
  implicit none

  ! logical :: print_enabled = .true.
  ! logical :: warnings_enabled = .true.

CONTAINS
  
  FUNCTION  get_condensed_m(z) result(m)
    ! returns condensed m where the ones not needed are filled with 0 (returns same size as z)
    type(inum), intent(in) :: z(:)
    integer :: m(size(z)), pos, i 
    m = 1
    pos = 1
    do i = 1,size(z)
      if(abs(z(i)) < zero) then
        if(i == size(z)) then
          pos = pos + 1
        else 
          m(pos) = m(pos) + 1
        end if
      else 
        pos = pos + 1
      end if
    end do
    m(pos:) = 0
  END FUNCTION get_condensed_m

  FUNCTION get_condensed_z(m, z_in) result(z_out)
    ! returns condensed z vector
    integer :: m(:), i, pos
    type(inum) :: z_in(:), z_out(size(m))
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(i) = z_in(pos)
    end do
  END FUNCTION get_condensed_z

  FUNCTION  get_flattened_z(m,z_in) result(z_out)
    ! returns flattened version of z based on m and z
    integer :: m(:), i, pos
    type(inum) :: z_in(:), z_out(sum(m))
    z_out = izero
    pos = 0
    do i=1,size(m)
      pos = pos + m(i)
      z_out(pos) = z_in(i)
    end do
  END FUNCTION get_flattened_z

  FUNCTION find_amount_trailing_zeros(z) result(res)
    type(inum) :: z(:)
    integer :: res, i
    res = 0
    do i = size(z), 1, -1
      if( abs(z(i)) < zero ) then
        res = res + 1
      else
        exit
      end if
    end do
  END FUNCTION find_amount_trailing_zeros

  FUNCTION find_marker(v) result(res)
    type(inum) :: v(:)
    integer res
    do res=1,size(v)
      if(v(res)%i0 == marker%i0) then
        return
      endif
    enddo
  END FUNCTION find_marker

  FUNCTION find_first_zero(v) result(res)
    ! returns index of first zero, or -1 if there is no zero
    integer :: v(:), i, res
    res = -1
    do i = 1,size(v)
      if(v(i) == 0) then
        res = i
        return
      end if
    end do
  END FUNCTION find_first_zero

  FUNCTION min_index(v)
    ! returns the index of the smallest element in v
    real(kind=prec) :: v(:), minimum
    integer :: min_index, i
    min_index = 1
    minimum = 1e15
    do i = 1,size(v)
      if(v(i) < minimum .and. v(i) > zero) then
        minimum = v(i)
        min_index = i
      end if
    end do
  END FUNCTION min_index

  FUNCTION zeroes(n) result(res)
    integer :: n
    type(inum) :: res(n)
    res = izero
  END FUNCTION zeroes

  RECURSIVE FUNCTION factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res, i
    res = 1
    do i=n,1,-1
      res = res*i
    enddo
  END FUNCTION factorial

  ! Adapted from https://rosettacode.org/wiki/Evaluate_binomial_coefficients#Fortran
  ! published under GNU Free Documentation License 1.2
  !
  ! This uses approximately (1.55 * n - 2.5) bit integers. This means
  ! that we can go up to ~ 83 for 128bit and ~ 42 on 64bit compilers.
  ! This is actually fine because the Bernoulli numbers this is used
  ! for are already O(10^19) and O(10^60) resp.
  FUNCTION binom(n, r)
    integer, intent(in) :: n, r
    real(kind=prec) :: binom

#ifdef __GFORTRAN__
    integer(16) :: num, den
    integer, parameter :: nmax = 82
#else
    integer(8) :: num, den
    integer, parameter :: nmax = 42
#endif

    integer i, k
    integer, parameter :: primes(20) = (/2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71/)

    if(n > nmax) then
      print*, "Can't calculate this large", n, nmax
    endif
    num = 1
    den = 1
    do i=0,r-1
      num = num*(n-i)
      den = den*(i+1)
      if (i > 0) then
        ! Divide out common prime factors
        do k=1,size(primes)
          if (mod(i,primes(k)) == 0) then
            num = num/primes(k)
            den = den/primes(k)
          end if
        end do
      end if
    end do
    binom = real(num/den, kind=prec)
  END FUNCTION binom
  ! subroutine print(s1,s2,s3,s4,s5)
  !   character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
  !   if(print_enabled) then
  !     print*, s1, s2, s3, s4, s5
  !   end if
  ! end subroutine print

  ! subroutine warn(s1,s2,s3,s4,s5)
  !   character(len = *), intent(in), optional :: s1, s2, s3, s4, s5
  !   if(warnings_enabled) then
  !     print*, 'Warning: ', s1, s2, s3, s4, s5
  !   end if
  ! end subroutine warn

END MODULE utils
