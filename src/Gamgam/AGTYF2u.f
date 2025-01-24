!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AGTYF2u(s,t,Lu)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. A.22
      include 'types.f'
      include 'constants.f'
      real(dp):: AGTYF2u,s,t,Lu

      AGTYF2u=(4._dp*pisq/27._dp+32._dp*Lu**2/9._dp-160._dp*Lu/27._dp)
     & *(t**2 + s**2)/(s*t)
      return
      end
