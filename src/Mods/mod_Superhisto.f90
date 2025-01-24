!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module Superhisto
    use omp_lib
    use types
    implicit none

    private

    type sh_histogram
        character(len=:), allocatable :: title
        real(dp) :: xmin, xmax, delx
        integer :: nbins = 0
        real(dp), allocatable :: xx(:)
        real(dp), allocatable :: xxsq(:)

        ! for real emission temporary binning of all
        ! dipole contributions
        real(dp), allocatable :: tmp(:)

        ! for plotting additional histogram columns
        real(dp), allocatable :: extras(:,:)

        logical :: customBinning
        real(dp), allocatable :: xarr(:)

        contains
        
        procedure :: init_custom => shinit_custom
        procedure :: init => shinit
        procedure :: initialized => shinitialized
        procedure :: book => shbook
        procedure :: tmpbook => shtmpbook
        procedure :: tmpcommit => shtmpcommit
        procedure :: reset => shreset
        procedure :: tmpreset => shtmpreset
        procedure :: add => shadd
        procedure :: iterprocess => shiterprocess
        procedure :: printme => shprint
        procedure :: serialize => serializeHistogram
        procedure :: deserialize => deserializeHistogram
    end type sh_histogram

    public :: sh_histogram
    public :: shbook
    public :: shreset
    public :: shinit
    public :: shinitialized
    public :: shadd
    public :: shtmpbook
    public :: shtmpcommit
    public :: shtmpreset

    ! for txt type writes
    public :: shwrite
    public :: shwritepdf

    ! for top type writes
    public :: shwritetop

    contains

    subroutine serializeHistogram(histo, unit)
        implicit none

        class(sh_histogram), intent(in) :: histo
        integer, intent(in) :: unit

        write (unit) histo%nbins
        if (histo%nbins > 0) then
            write (unit) len(histo%title)
            write (unit) histo%title
            write (unit) histo%xmin, histo%xmax, histo%delx
            write (unit) histo%customBinning
            if (histo%customBinning) then
                write (unit) size(histo%xarr)
                write (unit) histo%xarr
            endif
            write (unit) histo%xx
            write (unit) histo%xxsq
        endif

    end subroutine

    subroutine deserializeHistogram(histo, unit)
        implicit none
        integer, intent(in) :: unit
        class(sh_histogram), intent(out) :: histo

        integer :: tlen, nbins
        character(:), allocatable :: title
        real(dp) :: xmin, xmax, delx
        integer :: xarrsize

        read (unit) nbins
        if (nbins > 0) then
            read (unit) tlen
            allocate(character(tlen) :: title)
            read (unit) title
            read (unit) xmin, xmax, delx

            read (unit) histo%customBinning
            if (histo%customBinning) then
                read (unit) xarrsize
                allocate(histo%xarr(xarrsize))
                read (unit) histo%xarr
                call histo%init_custom(title,histo%xarr)
            else
                call histo%init(title,xmin,xmax,delx)
            endif

            read (unit) histo%xx
            read (unit) histo%xxsq
        endif

    end subroutine

    subroutine shprint(this)
        implicit none
        class(sh_histogram), intent(in) :: this

        integer :: k

        write (*,*) this%title
        write (*,*) "underflow ", this%xx(0), this%xxsq(0)
        if (this%customBinning) then
            do k=1,this%nbins
                write(*,*) this%xarr(k), this%xarr(k+1), this%xx(k), this%xxsq(k)
            enddo
        else
            do k=1,this%nbins
                write(*,*) this%xmin+(k-1)*this%delx, this%xx(k), this%xxsq(k)
            enddo
        endif
        write (*,*) "overflow ", this%xx(this%nbins+1), this%xxsq(this%nbins+1)
        write (*,*) "sum", sum(this%xx), sqrt(sum(this%xxsq**2))
        write (*,*) ""

    end subroutine

    subroutine shiterprocess(this, ha, ncall)
        implicit none
        class(sh_histogram), intent(inout) :: this
        class(sh_histogram), intent(in) :: ha
        real(dp), intent(in) :: ncall
        real(dp), allocatable :: sigsq(:)

        integer :: k

        ! this will be 1-indexed through the automatic allocation
        sigsq = ((ncall*ha%xxsq - ha%xx**2)/(ncall-1))

        do k=0,ha%nbins+1
            if (ha%xx(k) /= 0._dp) then
                this%xx(k) = this%xx(k) + ha%xx(k) / sigsq(k+1)
                this%xxsq(k) = this%xxsq(k) + 1._dp / sigsq(k+1)
            endif
        enddo

    end subroutine

    subroutine shadd(this, histo)
        implicit none

        class(sh_histogram), intent(inout) :: this
        class(sh_histogram), intent(in) :: histo

        ! some basic sanity check
        if (this%nbins /= histo%nbins) then 
           write(6,*) 'Abort in mod_Superhisto'
           stop
        endif
        this%xx = this%xx + histo%xx
        this%xxsq = this%xxsq + histo%xxsq

    end subroutine

    function shinitialized(histo)
        implicit none
        logical :: shinitialized
        class(sh_histogram), intent(in) :: histo

        if ( histo%nbins > 0 ) then
            shinitialized = .true.
        else
            shinitialized = .false.
        endif

    end function

    subroutine shinit_custom(newsh, title, xarr)
        implicit none
        class(sh_histogram), intent(inout) :: newsh
        character(len=*), intent(in) :: title
        real(dp), intent(in) :: xarr(:)

        integer :: n 

        newsh%title = title
        newsh%xmin = xarr(1)
        n = size(xarr)
        newsh%xmax = xarr(n)
        newsh%delx = 0._dp
        newsh%customBinning = .true.
        newsh%xarr = xarr

        newsh%nbins = n-1

        if (.not. allocated(newsh%xx)) then
            allocate(newsh%xx(0:n), source=0._dp)
            allocate(newsh%xxsq(0:n), source=0._dp)
            allocate(newsh%tmp(0:n), source=0._dp)
        else
            newsh%xx = 0._dp
            newsh%xxsq = 0._dp
            newsh%tmp = 0._dp
        endif


    end subroutine shinit_custom

    subroutine shinit(newsh, title, xmin, xmax, delx)
        implicit none
        class(sh_histogram), intent(inout) :: newsh
        character(len=*), intent(in) :: title
        real(dp), intent(in) :: xmin, xmax, delx
        integer :: nbins

        newsh%title = title
        newsh%xmin = xmin
        newsh%xmax = xmax
        newsh%delx = delx
        newsh%customBinning = .false.

        nbins = ceiling((xmax-xmin)/delx)
        newsh%nbins = nbins

        !write (*,*) "initializing histo ", title, xmin, xmax, delx, nbins, &
                !"for thread", omp_get_thread_num()

        if (allocated(newsh%xx) .eqv. .false.) then
            allocate(newsh%xx(0:nbins+1), source=0._dp)
            allocate(newsh%xxsq(0:nbins+1), source=0._dp)
            allocate(newsh%tmp(0:nbins+1), source=0._dp)
        else
            newsh%xx = 0._dp
            newsh%xxsq = 0._dp
            newsh%tmp = 0._dp
        endif

    end subroutine

    subroutine shbook(hist, x, wgt)
        use ieee_arithmetic
        use omp_lib
        implicit none
        class(sh_histogram), intent(inout) :: hist
        real(dp), intent(in) :: x, wgt

        integer :: ibin,j

        if (ieee_is_nan(x)) then
            write (*,*) "NaN in shtmpbook, please fix observable ", hist%title
            return
        endif

        if (hist%customBinning) then
            if (x < hist%xmin) then
                ibin = 0
            else if (x > hist%xmax) then
                ibin = hist%nbins + 1
            else if (x == hist%xmin) then
                ibin = 1
            else
                loopbins: do j=1,hist%nbins
                    if (x <= hist%xarr(j+1)) then
                        ibin = j
                        exit loopbins
                    endif
                enddo loopbins
            endif
        else
            if (x < hist%xmin) then
                ibin = 0
            else if (x > hist%xmax) then
                ibin = hist%nbins + 1
            else if (x == hist%xmin) then
                ibin = 1
            else
                ibin = ceiling((x - hist%xmin)/hist%delx)
            endif
        endif

        hist%xx(ibin) = hist%xx(ibin) + wgt
        hist%xxsq(ibin) = hist%xxsq(ibin) + wgt**2

    end subroutine

    subroutine shtmpbook(hist, x, wgt)
        use ieee_arithmetic
        implicit none
        class(sh_histogram), intent(inout) :: hist
        real(dp), intent(in) :: x, wgt

        integer :: ibin,j

        if (ieee_is_nan(x)) then
            write (*,*) "NaN in shtmpbook, please fix observable ", hist%title
            return
        endif

        if (hist%customBinning) then
            if (x < hist%xmin) then
                ibin = 0
            else if (x > hist%xmax) then
                ibin = hist%nbins + 1
            else if (x == hist%xmin) then
                ibin = 1
            else
                loopbins: do j=1,hist%nbins
                    if (x <= hist%xarr(j+1)) then
                        ibin = j
                        exit loopbins
                    endif
                enddo loopbins
            endif
        else
            if (x < hist%xmin) then
                ibin = 0
            else if (x > hist%xmax) then
                ibin = hist%nbins + 1
            else if (x == hist%xmin) then
                ibin = 1
            else
                ibin = ceiling((x - hist%xmin)/hist%delx)
            endif
        endif

        hist%tmp(ibin) = hist%tmp(ibin) + wgt
    end subroutine

    subroutine shtmpcommit(this)
        implicit none
        class(sh_histogram), intent(inout) :: this

        this%xx = this%xx + this%tmp
        this%xxsq = this%xxsq + this%tmp**2

        this%tmp = 0._dp
    end subroutine

    subroutine shreset(hist)
        implicit none
        class(sh_histogram), intent(inout) :: hist

        if (hist%nbins > 0) then
            hist%xx = 0._dp
            hist%xxsq = 0._dp
            hist%tmp = 0._dp
        endif
    end subroutine

    subroutine shtmpreset(hist)
        implicit none
        class(sh_histogram), intent(inout) :: hist

        if (hist%nbins > 0) then
            hist%tmp = 0._dp
        endif
    end subroutine

    subroutine shwrite(this, iounit)
        implicit none
        type(sh_histogram), intent(in) :: this
        integer, intent(in) :: iounit
        character :: tab = char(9)
        integer :: k

        write (iounit,'(A,A)') "# ", this%title
        write (iounit,'(A,G15.8,A,G15.8)') "# underflow"//tab, this%xx(0), tab, this%xxsq(0)
        write (iounit,'(A,G15.8,A,G15.8)') "# overflow"//tab, this%xx(this%nbins+1), tab, this%xxsq(this%nbins+1)
        write (iounit,'(A,G15.8,A,G15.8)') "# sum"//tab, sum(this%xx), tab, sqrt(sum(this%xxsq**2))

        if (this%customBinning) then
            if (allocated(this%extras)) then
                write (iounit,'(A)') "# xmin"//tab//"xmax"//tab//"taucutuncert"//tab// &
                            "fitchisquare"//tab//"fitresult"//tab//"fiterror"
                do k=1,this%nbins
                    write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8)') &
                             this%xarr(k), &
                        tab, this%xarr(k+1), &
                        tab, this%extras(1,k), &
                        tab, this%extras(2,k), &
                        tab, this%xx(k), &
                        tab, this%xxsq(k)
                enddo
            else
                write (iounit,'(A)') "# xmin"//tab//"xmax"//tab//"cross"//tab// &
                            "numerror"
                do k=1,this%nbins
                    write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8)') this%xarr(k), tab, &
                        this%xarr(k+1), tab, this%xx(k), tab, this%xxsq(k)
                enddo
            endif
        else
            if (allocated(this%extras)) then
                write (iounit,'(A)') "# xmin"//tab//"xmax"//tab//"taucutuncert"//tab// &
                            "fitchisquare"//tab//"fitresult"//tab//"fiterror"
                do k=1,this%nbins
                    write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8)') &
                             this%xmin+(k-1)*this%delx, &
                        tab, this%xmin+k*this%delx, &
                        tab, this%extras(1,k), &
                        tab, this%extras(2,k), &
                        tab, this%xx(k), &
                        tab, this%xxsq(k)
                enddo
            else
                write (iounit,'(A)') "# xmin"//tab//"xmax"//tab//"cross"//tab// &
                            "numerror"
                do k=1,this%nbins
                    write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8)') this%xmin+(k-1)*this%delx, tab, &
                        this%xmin+k*this%delx, tab, this%xx(k), tab, this%xxsq(k)
                enddo
            endif
        endif

    end subroutine

    subroutine shwritepdf(histCentral,histSymm,histPlus,histMinus,iounit)
        implicit none
        type(sh_histogram), intent(in) :: histCentral, histSymm, histPlus, histMinus
        integer, intent(in) :: iounit
        character :: tab = char(9)
        integer :: k

        write (iounit,'(A,A)') "# ", histCentral%title
        write (iounit,'(A,G15.8,A,G15.8)') "# underflow"//tab, histCentral%xx(0), tab, histCentral%xxsq(0)
        write (iounit,'(A,G15.8,A,G15.8)') "# overflow"//tab, histCentral%xx(histCentral%nbins+1), &
            tab, histCentral%xxsq(histCentral%nbins+1)
        write (iounit,'(A,G15.8,A,G15.8)') "# sum"//tab, sum(histCentral%xx), tab, sqrt(sum(histCentral%xxsq**2))
        write (iounit,'(A)') "# xmin"//tab//"xmax"//tab//"cross"//tab// &
                    "numerror"//tab//"pdfsymm"//tab//"pdfplus"//tab//"pdfminus"
        if (histCentral%customBinning) then
            do k=1,histCentral%nbins
                write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A)') &
                    histCentral%xarr(k), tab, &
                    histCentral%xarr(k+1), tab, &
                    histCentral%xx(k), tab, histCentral%xxsq(k), tab, &
                    histSymm%xx(k), tab, &
                    histPlus%xx(k), tab, &
                    histMinus%xx(k), tab
            enddo
        else
            do k=1,histCentral%nbins
                write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8,A)') &
                    histCentral%xmin+(k-1)*histCentral%delx, tab, &
                    histCentral%xmin+k*histCentral%delx, tab, &
                    histCentral%xx(k), tab, histCentral%xxsq(k), tab, &
                    histSymm%xx(k), tab, &
                    histPlus%xx(k), tab, &
                    histMinus%xx(k), tab
            enddo
        endif

    end subroutine

    subroutine shwritetop(this, iounit)
        implicit none
        type(sh_histogram), intent(in) :: this
        integer, intent(in) :: iounit
        character :: tab = char(9)
        integer :: k

        write (iounit,'(A,A)') "# ", this%title
        write (iounit,'(A,G15.8,A,G15.8)') "# underflow"//tab, this%xx(0), tab, this%xxsq(0)
        write (iounit,'(A,G15.8,A,G15.8)') "# overflow"//tab, this%xx(this%nbins+1), tab, this%xxsq(this%nbins+1)
        write (iounit,'(A,G15.8,A,G15.8)') "# sum"//tab, sum(this%xx), tab, sqrt(sum(this%xxsq**2))

        if (allocated(this%extras)) then
            write (iounit,'(A)') "# x"//tab//"taucutuncert"//tab// &
                        "fitchisquare"//tab//"fitresult"//tab//"fiterror"
            do k=1,this%nbins
                write (iounit,'(G15.8,A,G15.8,A,G15.8,A,G15.8,A,G15.8)') &
                         this%xmin+(k-0.5)*this%delx, &
                    tab, this%extras(1,k), &
                    tab, this%extras(2,k), &
                    tab, this%xx(k), &
                    tab, this%xxsq(k)
            enddo
        else
            write (iounit,'(A)') "# x"//tab//tab//"cross"//tab//"numerror"
            do k=1,this%nbins
                write (iounit,'(G15.8,A,G15.8,A,G15.8)') this%xmin+(k-0.5)*this%delx, tab, &
                    this%xx(k), tab, this%xxsq(k)
            enddo
        endif


    end subroutine


end module
