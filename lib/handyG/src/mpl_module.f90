 
MODULE mpl_module
  use globals
  use ieps
  implicit none

#ifndef KINDREAL
#define KINDREAL 8
#endif
#if KINDREAL==16
  real(kind=prec), parameter :: underflowalert = 1.e-600_prec
  real(kind=prec), parameter ::  overflowalert = 1.e+600_prec
#else
  real(kind=prec), parameter :: underflowalert = 1.e-250_prec
  real(kind=prec), parameter ::  overflowalert = 1.e+250_prec
#endif

#ifdef MPL_CACHE
  ! Max weight, no. cache
  integer,parameter :: cachesize(2) = (/ 4, 2500 /)
  type el
    integer :: m(cachesize(1))
    complex(kind=prec) :: c(cachesize(1))
    complex(kind=prec) :: ans
  end type el

  type(el) cache(cachesize(1),abs(cachesize(2)))
  integer cacheind(cachesize(1))
#endif
CONTAINS 

  FUNCTION MPL_converges(m,x)
    ! checks if the MPL converges 
    complex(kind=prec) :: x(:)
    integer :: m(:)
    logical :: MPL_converges
    MPL_converges = .false.
    if(abs(product(x)) < 1) then
      if(m(1) /= 1 .or. abs(x(1) - 1) < zero) then
        MPL_converges = .true.
      end if
    end if
  END FUNCTION MPL_converges

#ifdef MPL_CACHE
  FUNCTION CHECK_CACHE(m, x, res)
    integer m(:)
    complex(kind=prec) :: x(:), res
    logical check_cache
    integer k, j, i
    j = size(x)
    if (j<=4) then
      do k=1,cacheind(j)
        if (all(cache(j,k)%m(:j) == m)) then
          do i=1,j
            if( abs(cache(j,k)%c(i)-x(i)).gt.zero ) goto 123
          enddo
          res = cache(j,k) % ans
          check_cache = .true.
          return
        endif
123     continue
      enddo
    endif
    check_cache = .false.
  END FUNCTION
#endif



  FUNCTION MPL(m, x) result(res)
    integer :: m(:)
    complex(kind=prec) :: x(:)
    complex(kind=prec) :: res
    complex(kind=prec) :: t(size(x)), cpow
    integer(kind=ikin) :: q, j, k, ipow
#ifdef MPL_CACHE
    if (check_cache(m,x,res)) return
#endif

    j = size(x)


#ifdef DEBUG
    if(verb >= 70) print*, 'called MPL(',m,',',x,')'
#endif
    if(size(m) /= size(x)) then
      print*, 'Error: m and x must have the same length'
    end if
    res=0
    q=0
    t=0
    do while(.true.)
      res = t(1)
      q = q+1
      if (q < 0) exit

      cpow = x(j)**q
      ipow = q**m(j)

      if(ipow .eq. 0) exit
      if (abs(cpow).lt.underflowalert) exit
      if (abs(cpow).gt. overflowalert) exit

      t(j) = t(j) + cpow/ipow
      do k=1,j-1
        ipow = (k+q)**m(j-k)
        cpow = x(j-k)**(k+q)
        if(ipow .eq. 0) exit
        if (abs(cpow).lt.underflowalert) exit
        if (abs(cpow).gt. overflowalert) exit
        t(j-k) = t(j-k) + t(j-k+1) * cpow/ipow
      enddo

      if (mod(q,2_ikin) .eq. 1) then
        if (abs(t(1)-res).lt.MPLdelta/100.) exit
      endif
    enddo
    res = t(1)


#ifdef MPL_CACHE
    if (j<=cachesize(1)) then
      if (cacheind(j)+1 > abs(cachesize(2))) then
        if (cachesize(2) < 0) print*,"MPL cache for depth=",j," is full. Try increasing cachesize(2)"
        return
      endif
      cacheind(j) = cacheind(j) + 1
      cache(j,cacheind(j)) = el( reshape(m, (/cachesize(1)/), pad=[      0            ]), &
                                 reshape(x, (/cachesize(1)/), pad=[cmplx(0.,kind=prec)]), &
                                 res )
    endif
#endif
  END FUNCTION MPL

END MODULE mpl_module
