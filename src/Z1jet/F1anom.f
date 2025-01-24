!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function F1anom(s12,s45,mt2,musq)
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scalarselect.f'
      complex(dp)::F1anom
      complex(dp)::lnrat
      real(dp)::s12,s45,mt2,musq

      if (mt2  ==  zero) then
      F1anom=one/two/(s45-s12)*(1.d0+s45/(s45-s12)*lnrat(-s12,-s45))
      return
      else
      F1anom=one/two/(s45-s12)
     & *(one+two*mt2*loopI3(s12,0.d0,s45,mt2,mt2,mt2,musq,0)
     & +(s45/(s45-s12))
     & *(loopI2(s45,mt2,mt2,musq,0)-loopI2(s12,mt2,mt2,musq,0)))
      endif

      return
      end

