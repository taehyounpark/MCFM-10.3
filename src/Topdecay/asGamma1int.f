!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function asGamma1int(omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1int

c--   Author R.K. Ellis April 2012
c--   Integrand for NLO width with W-offshell
      include 'constants.f'
      real(dp):: omsq,mt,xi,ga,besq,asGamma1
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      asGamma1int=ga*xi/pi
     & /((1._dp-xi*omsq)**2+ga**2)*asGamma1(mt,besq,omsq)
      return
      end
