!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_2jet_swap(pin,msq)
      implicit none
      include 'types.f'
c--- this is just a wrapper routine to qqb_dirgam_g,
c--- that interchanges p3 and p4 before the call
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: pin(mxpart,4),msq(-nf:nf,-nf:nf),
     & pswap(mxpart,4)

      pswap(:,:)=pin(:,:)
      pswap(3,:)=pin(4,:)
      pswap(4,:)=pin(3,:)

      call qqb_2jet(pswap,msq)

      return
      end

