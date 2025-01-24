!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_Msqpt345sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt345^2), where M is the mass of the particle (345)

      include 'mxpart.f'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==ktt_bbl) .or.
     &   (kcase==ktt_bbu) .or.
     &   (kcase==ktt_bbh)) then
        mu0=(p(3,4)+p(4,4)+p(5,4))**2-(p(3,3)+p(4,3)+p(5,3))**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt345^2)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

