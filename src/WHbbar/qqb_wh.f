!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wh(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: p(mxpart,4)
      real(dp):: s,prop,fac,qqbWH,qbqWH,s56
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,msqhtautau,msqhbb

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      msq(:,:)=zip

c---calculate the 2 W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
c--- Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          s56=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(s56)
      elseif (hdecaymode == 'bqba') then
          s56=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(s56)
      else
        write(6,*) 'Unimplemented process in qqb_wh'
        stop
      endif
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

c--- adjust for fixed H->bb BR if necessary
      if ((FixBrHbb) .and. (hdecaymode == 'bqba')) then
        fac=fac*GamHbb/GamHbb0
      endif

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)


      do j=-nf,nf
      do k=-nf,nf
        if ((j > 0) .and. (k < 0)) msq(j,k)=Vsq(j,k)*qqbWH
        if ((j < 0) .and. (k > 0)) msq(j,k)=Vsq(j,k)*qbqWH
      enddo
      enddo

      return
      end

