!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dkqqb_zh_v(p,msq)
      implicit none
      include 'types.f'
c-----Virtual corrections for Higgs decay to b-bbar
c-----Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z
c                           |    |
c                           |    --> e^-(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c----Formula taken from Braaten and Leveille, PR D22, 715, 1980

      include 'nf.f'
      include 'mxpart.f'
      include 'hbbparams.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),decayv

      call qqb_zh(p,msq)
      call hbbdecay_v(p,5,6,decayv)

      msq(:,:)=msq(:,:)*decayv

c--- adjust for fixed H->bb BR if necessary
      if (FixBrHbb) then
         msq(:,:)=msq(:,:)*GamHbb0/GamHbb1
      endif

      return

      end
