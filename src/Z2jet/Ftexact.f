!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Ftexact(s23,mt2)
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scalarselect.f'
      complex(dp)::Ftexact,bb1
      real(dp)::s23,mt2,musq,s

c--- musq is irrelevant, set it to some value
      musq=mt2

      bb1 = loopI2(s23,mt2,mt2,musq,0)
      if(s23/mt2 < 0.25_dp) then
         s = s23/mt2
         bb1 = log(musq/mt2)+s/6._dp*(1._dp+s/10._dp*(1._dp+s/7._dp*(1._dp+s/6._dp)))
      endif
      Ftexact=
     & -(one+six*mt2*loopI3(s23,0._dp,0._dp,mt2,mt2,mt2,musq,0)
     &   +12._dp*mt2/s23*(bb1-loopI2(0._dp,mt2,mt2,musq,0)))

c expansion of entire contribution
      if(s23/mt2 < 0.25_dp) then
        s = s23/mt2
        Ftexact = s/20._dp*(1._dp+2._dp/21._dp*s*(1._dp+s/8._dp))
      endif

      return
      end
