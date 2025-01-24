!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine dkqqb_wh_v(p,msq)
      implicit none
      include 'types.f'
c-----Virtual corrections for Higgs decay to b-bbar
c-----Matrix element squared averaged over initial colors and spins
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
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'hbbparams.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & s,prop,fac,qqbWH,qbqWH,s56,hdecay,decayv
c----statement function
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
c----end statement function

      call qqb_wh(p,msq)
      call hbbdecay_v(p,5,6,decayv)

      msq(:,:)=msq(:,:)*decayv

c--- adjust for fixed H->bb BR if necessary
      if (FixBrHbb) then
        msq(:,:)=msq(:,:)*GamHbb0/GamHbb1
      endif

      return
      end

