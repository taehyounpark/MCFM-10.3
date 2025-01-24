!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function F3anom(s12,s45,mt2,musq)
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scalarselect.f'
      complex(dp)::F3anom
      complex(dp)::lnrat
      real(dp)::s12,s45,mt2,musq

      if (mt2  ==  zero) then
      F3anom=one/two/(s12-s45)*(1.d0-s12/(s12-s45)*lnrat(-s12,-s45))
      return
      else
      F3anom=one/two/(s12-s45)
     & *(one+two*mt2*loopI3(s12,0.d0,s45,mt2,mt2,mt2,musq,0)
     & -(s12/(s12-s45))*(loopI2(s45,mt2,mt2,musq,0)-loopI2(s12,mt2,mt2,musq,0)))
      endif

      return
      end

