!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AGTYX4u(s,t,Lx,Ly,Lu)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.15
      include 'types.f'
      real(dp):: AGTYX4u
      include 'constants.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX4u=32._dp/3._dp*Lu
     & *(Lx+Lx**2-2._dp*LY-2._dp*Ly*Lx+2._dp*Ly**2+pisq)
     & +8._dp/3._dp*Lu
     & *(3._dp*pisq-6._dp*Ly-6._dp*Ly*Lx-14._dp
     & +3._dp*Lx**2+6._dp*Ly**2+3._dp*Lx)*(t**2+s**2)/(s*t)
     & -8._dp/3._dp*Lu
     & *(-pisq+2._dp*Ly*Lx-Lx**2+3._dp*Lx)*(t**2-s**2)/(s*t)

      return
      end

