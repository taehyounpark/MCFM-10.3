!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function AGTYX1s(t,u,Lx,Ly,Ls)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.4
      include 'types.f'
      real(dp):: AGTYX1s,X1sx,t,u,Lx,Ly,Ls
      AGTYX1s=X1sx(t,u,Lx,Ly,Ls)+X1sx(u,t,Ly,Lx,Ls)
      return
      end

      function X1sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: X1sx
      real(dp)::t,u,Lx,Ly,Ls

      X1sx=2._dp/3._dp*(3*Ly+Ly**2-7._dp+2*Lx**2)*(2*Ls+Ly+Lx )*(t/u)
     & +(4._dp/3._dp*Lx**2*Ly+8._dp/3._dp*Lx*Ls+4._dp/3*Lx*Ly
     & +4._dp/3._dp*Ly**3+4._dp/3._dp*Ly**2+8._dp/3._dp*Ly**2*Ls)
c + t <---> u
      return
      end
