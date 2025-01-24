!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine storeptilde(nd,p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'ptilde.f'

      integer, intent(in) :: nd
      real(dp), intent(in) :: p(mxpart,4)

      ptilde(nd,:,:) = p(:,:)

      end

