!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine getptildejet(nd,pjet)
          use types
      implicit none
      include 'mxpart.f'
      include 'npart.f'
      include 'ptilde.f'
      integer, intent(in) :: nd
      real(dp), intent(out) :: pjet(mxpart,4)

      pjet(1:npart+2,1:4) = ptildejet(1:npart+2,1:4,nd)
      pjet(npart+3,1:4) = 0._dp   ! ensure next entry is zero

      return
      end

