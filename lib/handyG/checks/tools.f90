MODULE TTOOLS
  use globals
  use mpl_module
  use gpl_module
  real(kind=prec) :: tol = 8.0e-7
  logical :: tests_successful = .true.
contains

  subroutine iprint(imsg, typ)
    character(len=*) imsg
    character(len=200) msg
    integer :: typ
    character(len=5),parameter :: red   = char(27)//'[31m'
    character(len=5),parameter :: green = char(27)//'[32m'
    character(len=5),parameter :: orange= char(27)//'[33m'
    character(len=4),parameter :: norm  = char(27)//'[0m'
    character       ,parameter :: cr    = achar(13)
    integer, save :: prevlen, prevtype
    character(len=200), save :: prevmsg

    if ( (prevtype == -1).and.(typ .ne. -1) ) then
      msg = prevmsg(1:prevlen) // trim(imsg)
    else
      msg = trim(imsg)
    endif

    select case(typ)
      case(0)
        print*,green //'[PASS]'//norm//' '//trim(msg)
      case(1)
        print*, red  //'[FAIL]'//norm//' '//trim(msg)
      case(2)
        print*, red  //'[FATL]'//norm//' '//trim(msg)
        stop 1
      case(3)
        print*,orange//'[WARN]'//norm//' '//trim(msg)
      case(4)
        print*,orange//'[INFO]'//norm//' '//trim(msg)
      case(-1)
        write(*,'(a)',advance='no')' [    ]'//' '//trim(msg)//cr
    end select
    prevtype = typ
    prevlen  = len_trim(msg)
    prevmsg  = msg
  end subroutine



  function readint(prev,i)
    integer i, readint, st
    character(len=32) :: arg
    character(len=*) :: prev
    i=i+1
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) call iprint("Argument "//prev//" requires a number",2)
    read(arg,*,iostat=st) readint
    if (st .ne. 0) call iprint("For argument "//prev//": "//trim(arg)//" is not a number",2)
  end function


  subroutine check(res, ref, ans, ttol)
    complex(kind=prec) :: res, ref
    real(kind=prec) :: delta, mytol
    real(kind=prec), optional :: ttol
    character(len=40) :: msg
    logical, optional :: ans

    mytol = tol
    if (present(ttol)) mytol = ttol


    delta = abs(res-ref)
    if(delta < mytol) then
      write(msg, 900) delta
      call iprint(trim(msg), 0)
      if (present(ans)) ans = .true.
    else
      write(msg, 900) delta
      call iprint(trim(msg), 1)
      tests_successful = .false.
      if (present(ans)) ans = .false.
    end if

900 format(" with delta = ",ES10.3)
  end subroutine check

  subroutine test_one_MPL(m,x,ref, test_id)
    integer :: m(:)
    complex(kind=prec) :: x(:), ref, res
    character(len=*) :: test_id

    call iprint('  testing MPL '//test_id//' ...',-1)
    res = MPL(m,x)
    call check(res,ref)
  end subroutine test_one_MPL


  subroutine test_one_condensed(m,z,y,k,ref,test_id)
    integer :: m(:), k
    complex(kind=prec) :: z(:), y, res, ref
    character(len=*) :: test_id

    call iprint('  testing GPL '//test_id//' ...',-1)
    res = G_condensed(m,toinum(z),inum(y,di0),k)
    call check(res,ref)
  end subroutine test_one_condensed

  subroutine test_one_flat(z,ref,test_id, ans)
    complex(kind=prec) :: z(:), res, ref
    character(len=*) :: test_id
    logical, optional :: ans

    call iprint('  testing GPL '//test_id//' ...',-1)
    res = G(z)
    if (present(ans)) then
      call check(res,ref,ans)
    else
      call check(res,ref)
    endif
  end subroutine test_one_flat

END MODULE TTOOLS
