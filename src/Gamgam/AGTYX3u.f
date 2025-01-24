!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AGTYX3u(s,t,Lx,Ly,Lu)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.14
      include 'types.f'
      real(dp):: AGTYX3u
      include 'constants.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX3u= (1._dp/ 18*Lx**2+ (-2._dp/ 9*Ly+2._dp/ 9*Lu )*Lx
     & -4._dp/9*Ly*Lu+1._dp/18*pisq+2._dp/9*Lu**2+2._dp/9*Ly**2)
     & *(t**2+s**2)/(s*t)
      return
      end
