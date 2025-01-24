! newton.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Straightforward implementation of the Newton-Raphson method
!     for finding roots of an equation
!

! modified to use double precision for MCFM

module newton_raphson
    use types

    implicit none

contains

subroutine find_root( f, xinit, tol, maxiter, result, success )

    interface
        function f(x)
            use types
            real(dp), intent(in) :: x
            real(dp) :: f
        end function f
    end interface

    real(dp), intent(in)     :: xinit
    real(dp), intent(in)     :: tol
    integer, intent(in)  :: maxiter
    real(dp), intent(out)    :: result
    logical, intent(out) :: success

    real(dp)                 :: eps = 1.0e-4
    real(dp)                 :: fx1
    real(dp)                 :: fx2
    real(dp)                 :: fprime
    real(dp)                 :: x
    real(dp)                 :: xnew
    integer              :: i

    result  = 0.0
    success = .false.

    x = xinit
    do i = 1,max(1,maxiter)
        fx1    = f(x)
        fx2    = f(x+eps)
        !write(*,*) i, fx1, fx2, eps
        fprime = (fx2 - fx1) / eps

        xnew   = x - fx1 / fprime

        if ( abs(xnew-x) <= tol ) then
            success = .true.
            result  = xnew
            exit
        endif

        x = xnew
        !write(*,*) i, x
     enddo

end subroutine find_root

end module

