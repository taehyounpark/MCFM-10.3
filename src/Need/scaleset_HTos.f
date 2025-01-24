!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_HTos(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt34^2)+sqrt(M^2+pt56^2)+pt7+pt8,
c---  where M34 and M56 are the masses of the particles (34) and (56)

      include 'mxpart.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,ptsum,pttwo,pt

      mu0=0._dp

      if (n3 == 1) then
c (34) is a massive particle)
        ptsum=pttwo(3,4,p)
        mu0=mu0+sqrt(mass3**2+ptsum**2)
      else
        mu0=mu0+pt(3,p)+pt(4,p)
      endif

      if (n2 == 1) then
c (56) is a massive particle)
        ptsum=pttwo(5,6,p)
        mu0=mu0+sqrt(mass2**2+ptsum**2)
      else
        if (p(5,4) > 1.e-8_dp) mu0=mu0+pt(5,p)
        if (p(6,4) > 1.e-8_dp) mu0=mu0+pt(6,p)
      endif

c no more resonances expected
      if (p(7,4) > 1.e-8_dp)  mu0=mu0+pt(7,p)
      if (p(8,4) > 1.e-8_dp)  mu0=mu0+pt(8,p)

      return
      end

