!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Gamma0int(omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0int

c--   Author R.K. Ellis April 2012
c--   Integrand for width with W-offshell
      include 'constants.f'
      real(dp):: omsq,mt,besq,xi,ga,Gamma0
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      Gamma0int=(ga*xi/pi)/((1._dp-xi*omsq)**2+ga**2)
     & *Gamma0(mt,besq,omsq)
      return
      end

