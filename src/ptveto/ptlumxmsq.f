!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine ptlumxmsq(p,xx,z1,z2,QB,order,central,hard_routine,xmsq)
      use LHAPDF
      use SCET
      use MCFMStorage, only: currentPart
      use Integration
c Wrapper routine for new implementation of pt subtraction
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
      include 'Rcut.f'
      include 'Lnu.f'
      include 'ipsgen.f'
      integer:: m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),xx(2),
     & soft(0:2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5),
     & z1,z2,QB(2),Lperp,
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & ptgetxmsq,tauc

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

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

!      do m=1,3
!      Lnu=10**m

c to include gg->h
      if (kcase==kggfus0) then
        call softveto(0,order,tauc,Rcut,soft)
      else
        call softveto(1,order,tauc,Rcut,soft)
      endif
      
      Lperp=log(facscale/tauc)

      if (order >= 0) then
        call fdist(ih1,xx(1),facscale,tBa0,1)
        call fdist(ih2,xx(2),facscale,tBb0,2)
      endif
      if (order >= 1) then
        call ptbeam1(ih1,z1,xx(1),Lperp,tBa1,1)
        call ptbeam1(ih2,z2,xx(2),Lperp,tBb1,2)
      endif
      if (order >= 2) then
        call ptbeam2(ih1,z1,xx(1),Lperp,tBa2,1)
        call ptbeam2(ih2,z2,xx(2),Lperp,tBb2,2)
      endif

      xmsq = ptgetxmsq(p,order,soft,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
     
!      write(6,*) 'Lnu,xmsq',Lnu,xmsq
!      enddo
!      pause

! This section is a repetition of the above for now;
! could be done more elegantly
      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  if (dynamictau) then
                    tauc=getdynamictau(p)
                  else
                    tauc=taucut
                  endif

                  if (kcase==kggfus0) then
                    call softveto(0,order,tauc,Rcut,soft)
                  else
                    call softveto(1,order,tauc,Rcut,soft)
                  endif

                  Lperp=log(facscale/tauc)

                  if (order >= 1) then
                    call ptbeam1(ih1,z1,xx(1),Lperp,tBa1,1)
                    call ptbeam1(ih2,z2,xx(2),Lperp,tBb1,2)
                  endif
                  if (order >= 2) then
                    call ptbeam2(ih1,z1,xx(1),Lperp,tBa2,1)
                    call ptbeam2(ih2,z2,xx(2),Lperp,tBb2,2)
                  endif

                  scetreweight(m) = ptgetxmsq(p,order,soft,
     &              tBa0,tBb0,tBa1,tBb1,tBa2,tBb2,msq0,msq1,msq2)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


      function ptgetxmsq(p,order,soft,
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
      include 'kpart.f'
      real(dp):: ptgetxmsq
      integer:: j,k,order
      real(dp):: p(mxpart,4),
     & msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),
     & soft(0:2),hard(2),
     & tBa0(-5:5),tBb0(-5:5),
     & tBa1(-5:5),tBb1(-5:5),
     & tBa2(-5:5),tBb2(-5:5),
     & bit,ptassemble

      ptgetxmsq=zip
      do j=-nf,nf
      do k=-nf,nf

c WARNING: this is meant to catch gg contributions that only enter for
c the first time at NNLO, i.e. everything so far except for gg->h process
      if ((j == 0) .and. (k == 0) .and. (order == 2) .and. (kcase /= kggfus0)) then
        ptgetxmsq = ptgetxmsq + tBa0(j)*tBb0(k)*ason4pi**2*msq2(j,k)
        cycle
      endif

      if (msq0(j,k) == zip) cycle

      hard(1)=msq1(j,k)/msq0(j,k)
      hard(2)=msq2(j,k)/msq0(j,k)

      if (onlypowcorr) then
        bit=zip
      else
        if (coeffonly) then
          bit=zip
        else
          bit=tBa0(j)*tBb0(k)
        endif
        if ( (order == 1) .or.
     &      ((order == 2) .and. (coeffonly .eqv. .false.)) ) then
          bit=bit+ason4pi*(
     &       (hard(1)+soft(1))*tBa0(j)*tBb0(k)
     &       +tBa0(j)*tBb1(k)+tBa1(j)*tBb0(k))
        endif
        if (order > 1) then
          bit=bit+ason4pi**2*(
     &       (hard(2)+soft(2))*tBa0(j)*tBb0(k)
     &      +(hard(1)+soft(1))*(tBa1(j)*tBb0(k)+tBa0(j)*tBb1(k))
     &      +hard(1)*soft(1)*tBa0(j)*tBb0(k)
     &      +tBa0(j)*tBb2(k)+tBa2(j)*tBb0(k)
     &      +tBa1(j)*tBb1(k))
        endif
      endif

c!!!---- no power corrections for now!
c!!      if ((incpowcorr) .or. (onlypowcorr)) then
c!!        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
c!!      endif

      ptgetxmsq=ptgetxmsq+bit*msq0(j,k)

      enddo
      enddo

      return
      end
