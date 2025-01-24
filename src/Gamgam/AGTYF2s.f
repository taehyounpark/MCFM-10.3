!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function AGTYF2s(t,u,Ls)
      implicit none
c     Results taken from hep-ph/0201274
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c     Eq. A.11
      include 'types.f'
      real(dp):: AGTYF2s,F2sx,t,u,Ls

      AGTYF2s=F2sx(t,u,Ls)
     &       +F2sx(u,t,Ls)

      return
      end


      function F2sx(t,u,Ls)
      implicit none
      include 'types.f'
      real(dp):: F2sx
      include 'constants.f'
      real(dp):: t,u,Ls

      F2sx=(((-160*Ls)/27. + (32*Ls**2)/9. - (92*pisq)/27.)*t)/u

c     + t <---> u
      return
      end
