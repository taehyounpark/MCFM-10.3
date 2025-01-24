!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_ptphoton(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  pt(photon)
      include 'mxpart.f'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0,pt

      if    ((kcase==kWgamma) .or.
     &       (kcase==kZgamma)) then
        mu0=pt(5,p)
      elseif((kcase==kdirgam) .or. (kcase==kgamjet)
     &   .or.(kcase==kgam_2j) .or. (kcase==kgam_3j)) then
        mu0=pt(3,p)
      else
        write(6,*) 'dynamicscale pt(photon)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

