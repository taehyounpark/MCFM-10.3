!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function AGTYX5s(t,u,Ls)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. B.8
      include 'types.f'
      real(dp):: AGTYX5s,X5sx,t,u,Ls
      AGTYX5s=X5sx(t,u,Ls)+X5sx(u,t,Ls)
      return
      end

      function X5sx(t,u,Ls)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: X5sx
      real(dp)::t,u,Ls

      X5sx= (32._dp/9._dp*Ls**2+32._dp/9._dp*pisq)*t/u
c       + t <---> u
      return
      end
