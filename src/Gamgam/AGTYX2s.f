!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function AGTYX2s(t,u,Lx,Ly,Ls)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.5
      include 'types.f'
      real(dp):: AGTYX2s,X2sx,t,u,Lx,Ly,Ls
      AGTYX2s=+X2sx(t,u,Lx,Ly,Ls)+X2sx(u,t,Ly,Lx,Ls)
      return
      end

      function X2sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: X2sx
      real(dp)::t,u,Lx,Ly,Ls

      X2sx=-1._dp/9._dp*(2*Ls+Ly+Lx)
     & *(11*Ls+3*Lx**2+6*Lx*Ly+10*Ly+10*Lx-3*pisq)*(t/u)
     & + (-2._dp/3._dp*Lx*Ly-2._dp/3._dp*Ly**2-2._dp/
     & 3*Lx**2*Ly-4._dp/3._dp*Ly**2*Ls-4._dp/
     & 3*Lx*Ls-2._dp/3._dp*Ly**3)
c +   t <---> u
      return
      end
