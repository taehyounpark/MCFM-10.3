
MODULE globals
  implicit none

#ifndef KINDREAL
#define KINDREAL selected_real_kind(15,32)
#endif
#ifndef KINDINT
#define KINDINT 4
#endif

  integer, parameter :: prec = KINDREAL
  integer, parameter :: ikin = KINDINT
  real(kind=prec), parameter :: zero = 10._prec**(-precision(1._prec))     ! values smaller than this count as zero
  real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec

  ! The following parameters control the accuracy of the evaluation
  real(kind=prec), protected :: MPLdelta = zero          ! if the MPL sum changes less then del it is truncated.
  real(kind=prec), protected ::  Lidelta = zero          ! like MPLdelta but for polylogs
  real(kind=prec), protected :: HoelderCircle = 1.1_prec ! when to apply Hoelder convolution?
  integer, parameter :: PolyLogCacheSize(2) = (/ 5, 100 /)
        ! = (/ mmax, n /). At most n polylogs with weight mmax will be cached

  complex(kind=prec), parameter :: i_ = (0.,1._prec)
  integer :: verb = 0

CONTAINS 
  
#ifdef DEBUG
  SUBROUTINE parse_cmd_args
    integer :: i
    character(len=32) :: arg
    i = 0
    do
      call get_command_argument(i, arg)
      if (len_trim(arg) == 0) exit

      ! parse verbosity
      if(trim(arg) == '-verb') then
        call get_command_argument(i+1,arg)
        read(arg,*) verb               ! str to int
      end if

      i = i+1
    end do
  END SUBROUTINE parse_cmd_args
#endif

  SUBROUTINE SET_OPTIONS(mpldel, lidel, hcircle)
    real(kind=prec), optional :: hcircle, mpldel, lidel
    if (present(MPLdel)) MPLdelta = mpldel
    if (present( Lidel)) LiDelta  =  lidel
    if (present(hcircle)) HoelderCircle = hcircle
  END SUBROUTINE

END MODULE globals
