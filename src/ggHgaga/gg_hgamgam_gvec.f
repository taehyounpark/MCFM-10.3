!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hgamgam_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
c  in is the label of the momentum contracted with n
      integer:: in
      real(dp):: msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      real(dp):: n(4),nDn,p(mxpart,4)

      msq(:,:)=0._dp

      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      call gg_hgamgam(p,msqt)

      msq(0,0)=-0.5_dp*nDn*msqt(0,0)

      return
      end


