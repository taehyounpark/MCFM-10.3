!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Ftriangle(x)
      implicit none
      include 'types.f'
      complex(dp):: Ftriangle
      include 'constants.f'
      real(dp):: x,y
      complex(dp):: arg
      Ftriangle=czip

      y=1._dp-4._dp*x
      if (y > 0._dp) then
      arg=cmplx((1._dp+sqrt(y))/(1._dp-sqrt(y)),kind=dp)
      Ftriangle=+chalf*(log(arg)-impi)**2
      elseif (y <= 0._dp) then
      Ftriangle=-ctwo*cmplx(asin(half/sqrt(x)),kind=dp)**2
      endif
      return
      end

