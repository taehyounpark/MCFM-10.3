!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_whas2(p,msq)
      implicit none
      include 'types.f'
c---
c--- O(as^2) contributions that do not factorize over LO matrix element
c---
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
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'hdecaymode.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'taucut.f'
      integer:: j,k
      real(dp):: p(mxpart,4),s,fac,hdecay,prop,sH,
     & msqhtautau,msqhbb,msqhgamgam,msq(-nf:nf,-nf:nf)
      real(dp):: qqb_wh_HtopEFT,qqbWHtoploop,qbqWHtoploop

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---calculate the two W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
c--- adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
            fac=fac*GamHbb/GamHbb0
          endif
      elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)
      else
          write(6,*) 'Unimplemented decay mode in qqb_whas2.f'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      qqbWHtoploop=aveqq*qqb_wh_HtopEFT(1,2,3,4,p)*hdecay
      qbqWHtoploop=aveqq*qqb_wh_HtopEFT(2,1,3,4,p)*hdecay
      if ((FixBrHbb).and.(hdecaymode=='bqba')) then
         qqbWHtoploop=qqbWHtoploop*GamHbb/GamHbb0
         qbqWHtoploop=qbqWHtoploop*GamHbb/GamHbb0
      endif

      msq(:,:)=zip
      do j=-nf,nf
      do k=-nf,nf
      if (j*k >= 0) cycle ! skip gluons, qq, aa

      if ((j > 0) .and. (k < 0)) then
           msq(j,k)=Vsq(j,k)*qqbWHtoploop
      elseif ((j < 0) .and. (k > 0)) then
           msq(j,k)=Vsq(j,k)*qbqWHtoploop

      endif

      enddo
      enddo

      return
      end


