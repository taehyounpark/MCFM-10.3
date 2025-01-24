!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_HTprime(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt34^2)+pt5+pt6, where M is the mass of the particle (34)
      include 'mxpart.f'
      include 'kprocess.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,pt34,pttwo,pt

      if((kcase==khttjet) .or.
     &   (kcase==khttscl) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==kh2jmas) .or.
     &   (kcase==kh2jscl) .or.
     &   (kcase==kgagajj) .or.
     &   (kcase==kggfus2) .or.
     &   (kcase==kggfus0)) then

        pt34=pttwo(3,4,p)
        mu0=sqrt(mass3**2+pt34**2)

        if (p(5,4) > 1.d-8) then
          mu0=mu0+pt(5,p)
        endif

        if (p(6,4) > 1.d-8) then
          mu0=mu0+pt(6,p)
        endif

        if (p(7,4) > 1.d-8) then
          mu0=mu0+pt(7,p)
        endif

      else
        write(6,*) 'dynamicscale sqrt(HTprime)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

