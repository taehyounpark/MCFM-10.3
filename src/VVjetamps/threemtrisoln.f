!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine threemtrisoln(S1,S2,P1DP2,
     & gam,al01,al02,al11,al12,al21,al22,D)
c--- given two momenta, computes and returns all relevant quantities
c---  involving "flat" momenta for three-mass triangle solution
c--- note that contents of za, zb are ameliorated to
c--- also contain spinor products corresponding to "flat" momenta
      implicit none
      include 'types.f'
      real(dp):: S1,S2,P1DP2,S1S2,
     & gam(2),D(2),al01(2),al02(2),al11(2),al12(2),al21(2),al22(2)
      integer i

      S1S2=S1*S2

c--- solutions for gamma, with some protection for numerical accuracy
      gam(1)=-P1DP2-sign(sqrt(P1DP2**2-S1S2),P1DP2)
      gam(2)=S1S2/gam(1)

c--- denominator
      do i=1,2
      D(i)=gam(i)**2-S1*S2
      al01(i)=S1*(gam(i)-S2)/D(i)
      al02(i)=S2*(gam(i)-S1)/D(i)
      al11(i)=al01(i)-S1/gam(i)
      al12(i)=al02(i)-1d0
      al21(i)=al01(i)-1d0
      al22(i)=al02(i)-S2/gam(i)
      enddo
      return
      end

