!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine GLYlumxmsq(p,xx,z1,z2,QB,order,central,hard_routine,xmsq)
      use LHAPDF
      use SCET
      use MCFMStorage, only: currentPart
      use Integration
c Wrapper routine for GLY implementation of QT subtraction
      implicit none
!      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'first.f'
      include 'facscale.f'
      include 'taucut.f'
      include 'nflav.f'
      include 'beamtype.f'
      include 'ipsgen.f'
      integer:: m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),
     & xx(2),F1(0:1,0:2),F2(0:1,0:4),hard(2),qsq,
     & Ba0(-5:5),Bb0(-5:5),
     & Ba1(-5:5,0:2),Bb1(-5:5,0:2),
     & Ba2(-5:5,0:4),Bb2(-5:5,0:4),
     & z1,z2,QB(2),
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & GLYgetxmsq

      logical, intent(in) :: central
      real(dp) :: origtaucut
      external born_msq, hard_routine

      if (first) then
        first=.false.
        call Gammafill(nflav)
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

      qsq=two*(p(1,4)*p(2,4)-p(1,1)*p(2,1)
     &        -p(1,2)*p(2,2)-p(1,3)*p(2,3))

      if (order >= 0) then
        call fdist(ih1,xx(1),facscale,Ba0,1)
        call fdist(ih2,xx(2),facscale,Bb0,2)
      endif
      if (order >= 1) then
        call GLYF1(qsq,F1)
        call GLYbeam1(ih1,z1,xx(1),Ba1,1)
        call GLYbeam1(ih2,z2,xx(2),Bb1,2)
      endif
      if (order >= 2) then
        call GLYF2(qsq,F2)
        call GLYbeam2(ih1,z1,xx(1),Ba2,1)
        call GLYbeam2(ih2,z2,xx(2),Bb2,2)
      endif

      xmsq = GLYgetxmsq(p,order,F1,F2,hard,
     &              Ba0,Bb0,Ba1,Bb1,Ba2,Bb2,msq0,msq1,msq2)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = GLYgetxmsq(p,order,F1,F2,hard,
     &              Ba0,Bb0,Ba1,Bb1,Ba2,Bb2,msq0,msq1,msq2)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


      function GLYgetxmsq(p,order,F1,F2,hard,
     &              Ba0,Bb0,Ba1,Bb1,Ba2,Bb2,msq0,msq1,msq2)
      use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'taucut.f'
      include 'qcdcouple.f'
      include 'transitionlabels.f'
      include 'kprocess.f'
      real(dp):: GLYgetxmsq
      integer:: j,k,order,pp
      real(dp):: p(mxpart,4),
     & F1(0:1,0:2),F2(0:1,0:4),hard(2),
     & Ba0(-5:5),Bb0(-5:5),
     & Ba1(-5:5,0:2),Bb1(-5:5,0:2),
     & Ba2(-5:5,0:4),Bb2(-5:5,0:4),
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & bit,GLYassemble
      real(dp):: tauc

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

      GLYgetxmsq=zip
      do j=-nf,nf
      do k=-nf,nf

c WARNING: this is meant to catch gg contributions that only enter for
c the first time at NNLO, i.e. everything so far except for gg->h process
      if ((j == 0) .and. (k == 0) .and. (order == 2) .and. (kcase /= kggfus0)) then
        GLYgetxmsq = GLYgetxmsq + Ba0(j)*Bb0(k)*ason4pi**2*msq2(j,k)
        cycle
      endif

      if (msq0(j,k) == zip) cycle

      hard(1)=msq1(j,k)/msq0(j,k)
      hard(2)=msq2(j,k)/msq0(j,k)

      if ((j == 0) .and. (k == 0)) then
        pp=gg
      else
        pp=qq
      endif

      if (onlypowcorr) then
        bit=zip
      else
        bit=GLYassemble(order,tauc,pp,
     &   Ba0(j),Bb0(k),Ba1(j,:),Bb1(k,:),
     &   Ba2(j,:),Bb2(k,:),F1,F2,hard)
      endif

c!!!---- no power corrections for now!
c!!      if ((incpowcorr) .or. (onlypowcorr)) then
c!!        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
c!!      endif

      GLYgetxmsq=GLYgetxmsq+bit*msq0(j,k)

      enddo
      enddo

      return
      end
