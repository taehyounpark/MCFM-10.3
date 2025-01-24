!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function AGTYX4s(t,u,Lx,Ly,Ls)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.7
      include 'types.f'
      real(dp):: AGTYX4s,X4sx,t,u,Lx,Ly,Ls
      AGTYX4s=+X4sx(t,u,Lx,Ly,Ls)+X4sx(u,t,Ly,Lx,Ls)
      return
      end

      function X4sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: X4sx
      real(dp)::t,u,Lx,Ly,Ls

      X4sx= (-32._dp/3._dp*Ly*pisq+32._dp/
     & 3._dp*Ls*Lx**2+16._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-16._dp*pisq+16._dp/3._dp*Ls*Ly**2-112._dp/
     & 3._dp*Ls)*t/u
     & +(32._dp/3._dp*Ls*Lx**2+32._dp/3._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-32._dp/3._dp*pisq)
c      + t <---> u

      return
      end
