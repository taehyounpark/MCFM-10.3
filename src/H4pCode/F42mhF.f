!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function F42mhF(psq,qsq,s,t)
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      complex(dp):: F42mhF

c      include 'scale.f'
c      include 'epinv.f'
      real(dp):: s,t,psq,qsq
c      real(dp):: den
c      complex(dp):: loopI4,loopI3
      complex(dp):: Lsm1_2mht
c      integer:: ep

c      den=s*t
c--- note: added a factor of 1/2 here
c      F42mhF=
c     & +den*loopI4(0._dp,0._dp,psq,qsq,s,t,0._dp,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  s*loopI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  t*loopI3(0._dp,0._dp,t,0._dp,0._dp,0._dp,musq,0)
c     & +psq*loopI3(0._dp,0._dp,psq,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & +qsq*loopI3(0._dp,0._dp,qsq,0._dp,0._dp,0._dp,musq,0)/2._dp

c--- NOTE: checked on 8/30/09 that Lsm1_2mht == (expression above)
      F42mhF=Lsm1_2mht(s,t,psq,qsq)

      return
      end

