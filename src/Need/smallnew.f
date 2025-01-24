!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module SafetyCuts
      contains

      function passed_smallnew(p,npart, cutoff_in)
c JC, Sep. 2018: applies a cut on dimensionless variables
c       --- cut on Ei/100 GeV [typical hard scale]
c       --- cut on sij/Ei/Ej, i.e. dot product of directions
      use types
      implicit none
      logical:: passed_smallnew
      include 'mxpart.f'
      include 'cutoff.f'
      include 'taucut.f'
      include 'debug.f'

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: npart
      real(dp), intent(in), optional :: cutoff_in

      integer :: j,k
      real(dp) :: dot,checkmin,Eref
c     integer, save :: ipp
c!$omp threadprivate(ipp)

      if (dosingcheck) then
        passed_smallnew = .true.
        return
      endif

      checkmin=1.e9_dp    ! a big number
      Eref=1.e2_dp

c loops over particles in positions ipp and above since
c they are (usually) the partons that cause problems
c      do j=ipp,npart+2
c Changed to start at 3 for additional safety
      do j=3,npart+2
         if (p(j,4) == 0._dp) then
           passed_smallnew = .false.
           return
         endif
c        write(6,*) 'Ej',j,abs(p(j,4)/1.e2_dp)
c        write(6,*) 'j1',j,abs(dot(p,j,1)/p(j,4)/p(1,4))
c        write(6,*) 'j2',j,abs(dot(p,j,2)/p(j,4)/p(2,4))
        checkmin=min(checkmin,abs(p(j,4)/Eref))
        checkmin=min(checkmin,abs(dot(p,j,1)/p(j,4)/p(1,4)))
        checkmin=min(checkmin,abs(dot(p,j,2)/p(j,4)/p(2,4)))
        if (j == npart+2) cycle
        do k=j+1,npart+2
           if (p(k,4) == 0._dp) then
             passed_smallnew = .false.
             return
           endif
c           write(6,*) 'jk',j,k,abs(dot(p,j,k)/p(j,4)/p(k,4))
          checkmin=min(checkmin,abs(dot(p,j,k)/p(j,4)/p(k,4)))
        enddo
      enddo
c      stop

c for processes such as X+jet, ensure that pt(X) is not abnormally small
c (leading to Gram determinant issues in virtual)
c      if ((ipp == 5) .and. (npart > 2)) then
c        checkmin=min(checkmin,pttwo(3,4,p)/Eref)
c      endif

      if (present(cutoff_in)) then
          if (checkmin < cutoff_in) then
              passed_smallnew = .false.
              return
          else
              passed_smallnew = .true.
              return
          endif
      endif

c DISABLED the following check since it (at least) causes problems with
c simple processes such as W and Z production at NNLO
c for the virtual contribution, always cut below 10^-7 for stability,
c      if (kpart == kvirt) then
c        if (checkmin < 1.e-7_dp) then
c            passed_smallnew = .false.
c            return
c        endif
c      endif

c make cut
      if (checkmin < cutoff) then
          passed_smallnew = .false.
          return
      endif
c for non-SCET calculations do not need such a stringent cutoff
      if (.not. usescet) then
        if (checkmin < 1.e-6_dp) then
            passed_smallnew = .false.
            return
        endif
      endif

      passed_smallnew = .true.
      return

      end function passed_smallnew

      end module SafetyCuts
