!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qtlumxmsq(p,xx,z1,z2,QB,order,central,hard_routine,xmsq)
        use LHAPDF
      use SCET
      use MCFMStorage, only: currentPart
      use Integration
c Wrapper routine for new implementation of QT subtraction
      implicit none
!      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'first.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'scet_const.f'
      include 'taucut.f'
      include 'nflav.f'
      include 'kprocess.f'
      include 'beamtype.f'
      include 'ipsgen.f'
      integer:: m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),xx(2),
     & tSb1(0:2),tSb2(0:4),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5,0:1),tBb1(-5:5,0:1),
     & tBa2(-5:5,0:2),tBb2(-5:5,0:2),
     & z1,z2,QB(2),
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & qtgetxmsq

      logical, intent(in) :: central
      real(dp) :: origtaucut
      external born_msq, hard_routine

      if (first) then
        first=.false.
        call Gammafill(nflav)
      endif

c to include gg->h
      if (kcase==kggfus0) then
        call qtsoft(0,order,tSb1,tSb2)
      else
        call qtsoft(1,order,tSb1,tSb2)
      endif

      if ((currentpart == nnloResVetoBelow) .and. (ipsgen == 1)) then
! for the pt-veto case, this is the matching correction and msq2 is not
! included in the fixed-order expansion of the resummed result
        call hard_routine(p,1,msq0,msq1,msq2)
        msq2(:,:)=0._dp
      else
! regular case
      call hard_routine(p,order,msq0,msq1,msq2)
      endif

      if (order >= 0) then
        call fdist(ih1,xx(1),facscale,tBa0,1)
        call fdist(ih2,xx(2),facscale,tBb0,2)
      endif
      if (order >= 1) then
        call qtbeam1(ih1,z1,xx(1),tBa1,1)
        call qtbeam1(ih2,z2,xx(2),tBb1,2)
      endif
      if (order >= 2) then
        call qtbeam2(ih1,z1,xx(1),tBa2,1)
        call qtbeam2(ih2,z2,xx(2),tBb2,2)
      endif

      xmsq = qtgetxmsq(p,order,tSb1,tSb2,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = qtgetxmsq(p,order,tSb1,tSb2,
     &                tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


      function qtgetxmsq(p,order,tSb1,tSb2,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
      use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'scet_const.f'
      include 'taucut.f'
      include 'kprocess.f'
      real(dp):: qtgetxmsq
      integer:: j,k,order
      real(dp):: p(mxpart,4),
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & tSb1(-1:1),tSb2(-1:3),hard(2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5,0:1),tBb1(-5:5,0:1),
     & tBa2(-5:5,0:2),tBb2(-5:5,0:2),
     & bit,qtassemble
      real(dp):: tauc

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

      qtgetxmsq=zip
      do j=-nf,nf
      do k=-nf,nf

c WARNING: this is meant to catch gg contributions that only enter for
c the first time at NNLO, i.e. everything so far except for gg->h process
      if ((j == 0) .and. (k == 0) .and. (order == 2) .and. (kcase /= kggfus0)) then
        qtgetxmsq = qtgetxmsq + tBa0(j)*tBb0(k)*ason4pi**2*msq2(j,k)
        cycle
      endif

      if (msq0(j,k) == zip) cycle

      hard(1)=msq1(j,k)/msq0(j,k)
      hard(2)=msq2(j,k)/msq0(j,k)

      if (onlypowcorr) then
        bit=zip
      else
        bit=qtassemble(order,tauc,
     &   tBa0(j),tBb0(k),tBa1(j,:),tBb1(k,:),
     &   tBa2(j,:),tBb2(k,:),tSb1,tSb2,hard)
      endif

c!!!---- no power corrections for now!
c!!      if ((incpowcorr) .or. (onlypowcorr)) then
c!!        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
c!!      endif

      qtgetxmsq=qtgetxmsq+bit*msq0(j,k)

      enddo
      enddo

      return
      end
