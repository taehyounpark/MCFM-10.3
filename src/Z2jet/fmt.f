!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function fmt(s12,s34,s56)
      implicit none
      include 'types.f'
      complex(dp):: fmt
      include 'masses.f'
      real(dp):: s12,s34,s56
c The full BDK expression keeps terms up to 1/mt^4 (II.15) but here we
c we will only use up to 1/mt^2, as done elsewhere
c      fmt=((1._dp+(2._dp*s34+s12+s56)/15._dp/mt**2)/(24._dp*mt**2))
      fmt=((1._dp)/(24._dp*mt**2))
      return
      end

