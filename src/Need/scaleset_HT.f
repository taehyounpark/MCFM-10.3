!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_HT(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- partonic HT (i.e. scalar sum of particle pt's, whether or not they pass cuts)
      include 'mxpart.f'
      include 'npart.f'
      integer:: j
      real(dp):: p(mxpart,4),mu0,pt

      mu0=0._dp
      do j=3,npart+2
      mu0=mu0+pt(j,p)
      enddo

      return
      end

