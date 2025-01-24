!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine donothing_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: in
      real(dp), intent(in) :: p(mxpart,4), n(4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
      msq(:,:) = 0._dp
      end


