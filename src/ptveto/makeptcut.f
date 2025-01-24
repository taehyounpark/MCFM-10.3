!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine makeptcut(pparton,pjet,isub,passed,nd)
          use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
      include 'kprocess.f'
      include 'hdecaymode.f'
      include 'jetlabel.f'
      logical:: passed
      integer:: j,isub
      real(dp):: pparton(mxpart,4),pjet(mxpart,4),tau,tauc,pt
      integer, intent(in) :: nd
      integer, save :: ipp
!$omp threadprivate(ipp)

      logical :: bin
      common/bin/bin

c--- determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp)  /=  'pp') .and. (plabel(ipp)  /=  ''))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in makeptcut.f'
          stop
        endif
c        write(6,*) 'found ipp=',ipp
      endif

      scetreweight(:) = 0._dp

      passed=.false.
      includeTaucutgrid(nd) = .true.

c--- if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart  ==  ipp-3+ntau) .and. (isub  ==  0))  then
        passed=.true.
        return
      endif

c calculate lead-jet pt
      tau=0._dp
      do j=ipp,ipp+jets-1
      tau=max(tau,pt(j,pjet))
      enddo

      if (dynamictau) then
        tauc=getdynamictau(pjet)
      else
        tauc=taucut
      endif

      if (bin .and. doMultitaucut) then
          ! no sampled value for taus will pass cuts, if smaller than the smallest taucut
          if (tau  < smallestTaucut*(tauc/taucut)) then
              scetreweight(:) = 0._dp
              includeTaucutgrid(nd) = .false.
              return
          endif

          ! otherwise compute "weights" for other taucuts
          do j=1,size(tcutarray)
              if (tau < tcutarray(j)*(tauc/taucut)) then
                  scetreweight(j) = 0._dp
              else
                  scetreweight(j) = 1._dp
              endif
          enddo

          ! and postpone this cut for later in the *int routines
          ! and in the plotting
          if (tau < tauc) then
              includeTaucutgrid(nd) = .false.
          endif
      else
          if (tau < tauc) return
      endif

c--- at this point cut has been passed
      passed=.true.
c      if (failedndp) then
c         write(6,*) 'WARNING, phase space point has passed, '
c         write(6,*) 'both tau cut AND has failed ndp check'
c      endif

      return
      end


