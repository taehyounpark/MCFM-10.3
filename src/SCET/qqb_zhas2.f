!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_zhas2(p,msq)
      implicit none
      include 'types.f'
c---
c--- O(as^2) contributions that do not factorize over LO matrix element
c---
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z
c                           |    |
c                           |     ->fermion(p3)+antifermion(p4)
c                           |
c                            ---> b(p5)+b(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'hdecaymode.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'taucut.f'
      integer:: j,k
      real(dp):: p(mxpart,4),s,fac,hdecay,prop,
     & msqgamgam,bit,msq(-nf:nf,-nf:nf),ggZH,
     & qqb_ZH_VItop,qqb_ZH_VIItop,sH

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &     -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
c     Deal with Higgs decay
      if (hdecaymode == 'tlta') then
         sH=s(5,6)+2._dp*mtau**2
         call htautaudecay(p,5,6,hdecay)
      elseif (hdecaymode == 'bqba') then
         sH=s(5,6)+2*mb**2
         call hbbdecay(p,5,6,hdecay)
c---  adjust for fixed H->bb BR if necessary
         if (FixBrHbb) then
            fac=fac*GamHbb/GamHbb0
         endif
      elseif (hdecaymode == 'gaga') then
         sH=s(5,6)
         hdecay=msqgamgam(hmass)
      elseif (hdecaymode == 'wpwm') then
         sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
         call hwwdecay(p,5,6,7,8,hdecay)
      else
         write(6,*) 'Unimplemented process in gg_hgg_v',hdecaymode
         stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      call gg_zh(p,ggZH)

      do j=-nf,nf
      k=-j
      bit=zip

      if (j > 0) then
        msq(j,k)=aveqq*hdecay*(qqb_ZH_VItop(1,2,3,4,p,abs(j))
     &                        +qqb_ZH_VIItop(1,2,3,4,p,abs(j)))
      elseif (j < 0) then
        msq(j,k)=aveqq*hdecay*(qqb_ZH_VItop(2,1,3,4,p,abs(k))
     &                        +qqb_ZH_VIItop(2,1,3,4,p,abs(k)))
      else
        msq(j,k)=ggZH
      endif

      enddo

      if (FixBrHbb .and. (hdecaymode=='bqba')) then
         msq(j,k)=msq(j,k)*GamHbb/GamHbb0
      endif

      return
      end

