!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine getptilde(nd,p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'ptilde.f'
      integer:: nd,i,j
      real(dp), intent(out) :: p(mxpart,4)

      do j=1,4
        do i=1,npart+2
        p(i,j)=ptilde(nd,i,j)
        enddo
      enddo

      return
      end

