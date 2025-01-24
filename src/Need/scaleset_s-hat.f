!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_shat(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      include 'mxpart.f'
      include 'kpart.f'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if ((kpart==klord .or. kcase==kWgaj_a .or. kcase==kWga_ew .or.
     &     kcase==kW_only)) then
        mu0=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     &     -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale s-hat not supported beyond LO.'
        stop
      endif

      return
      end

