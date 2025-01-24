!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_Msqpt3456sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt3456^2), where M is the mass of the particle (3456)

      include 'mxpart.f'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kWWqqbr) .or.
     &   (kcase==kWZbbar) .or.
     &   (kcase==kZZlept) .or.
     &   (kcase==kVVlept) .or.
     &   (kcase==kWW_jet) .or.
     &   (kcase==kWZ_jet) .or.
     &   (kcase==kZZ_jet)) then
        mu0=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2-(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt3456^2)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

