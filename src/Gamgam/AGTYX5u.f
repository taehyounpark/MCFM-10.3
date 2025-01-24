!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function AGTYX5u(s,t,Lu)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.16
      include 'types.f'
      real(dp):: AGTYX5u
      real(dp)::s,t,Lu
      AGTYX5u=32._dp/9._dp*Lu**2*(t**2+s**2)/(s*t)
      return
      end
