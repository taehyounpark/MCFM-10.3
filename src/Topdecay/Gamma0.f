!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Gamma0(mt,besq,omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0
c--   Author: John M. Campbell and R.K. Ellis, January 2012
c--   Taken from formula (2) of
c--   Fermilab-PUB-12-078-T

      include 'constants.f'
      include 'ewcouple.f'
      real(dp):: mt,omsq,besq,Gammainfty,f,P3b
      Gammainfty=GF*mt**3/(8._dp*rt2*pi)
      P3b=0.5_dp*sqrt(1._dp+omsq**2+besq**2-2._dp*(omsq+besq+omsq*besq))
      f=(1._dp-besq)**2+omsq*(1._dp+besq)-2._dp*omsq**2
      Gamma0=Gammainfty*2._dp*P3b*f
      return
      end

