!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function fmtfull(s12,s34,s56)
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      include 'scalarselect.f'
      complex(dp):: fmtfull
      real(dp):: mtsq,s12,s34,s56,Del12,Del34,Del56,Del3
      complex(dp):: bdiff
      include 'cplx.h'

c     Implementation of BDKW Nucl Phys. 513, 3 (1998)
c     Eqn. (B.15) upgraded so that is gives the full mass
c     dependence.
      Del12=s12-s34-s56
      Del34=s34-s12-s56
      Del56=s56-s12-s34
      Del3=s12*Del12+s34*Del34+s56*Del56
      mtsq=mt**2

      fmtfull=
     & ((three*s12*s56*s34*Del34/Del3**2-(s12*s56-mtsq*del34)/Del3))
     & *(-loopI3(s12,s34,s56,mtsq,mtsq,mtsq,musq,0))
     & +((three*s56*Del56/Del3**2-half/Del3)*s12)
     & *bdiff(s34,s12,mtsq)
     & +((three*s12*Del12/Del3**2-half/Del3)*s56)
     & *bdiff(s34,s56,mtsq)
     & -cplx2(half*Del34/Del3,zip)
      end

