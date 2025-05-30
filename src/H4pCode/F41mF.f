!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function F41mF(psq,s,t)
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      complex(dp):: F41mF
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 21)

c      include 'scale.f'
      real(dp):: s,t,psq
c      real(dp):: den
c      complex(dp):: loopI4,loopI3
      complex(dp):: Lsm1

c      den=s*t
c--- MODIFIED: added factor of 1/2 and removed finite pieces from poles
c      F41mF=
c     & +den*loopI4(0._dp,0._dp,0._dp,psq,s,t,0._dp,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -   s*loopI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,0)
c     & -   t*loopI3(0._dp,0._dp,t,0._dp,0._dp,0._dp,musq,0)
c     & + psq*loopI3(0._dp,0._dp,psq,0._dp,0._dp,0._dp,musq,0)

c--- NOTE: checked on 8/30/09 that Lsm1 == (expression above)
      F41mF=Lsm1(-s,-psq,-t,-psq)

      return
      end

