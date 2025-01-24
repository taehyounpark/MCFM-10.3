!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function F31m(s)
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      complex(dp):: F31m

      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      real(dp):: s
c      complex(dp):: loopI3
      complex(dp):: lnrat
c      integer:: ep
c      F31m=czip
c      do ep=-2,0
c      F31m=F31m+s*epinv**(-ep)*loopI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep)
c      enddo

c--- NOTE: checked on 8/31/09 that this agrees with the expression above
      F31m=epinv*epinv2-epinv*lnrat(-s,musq)+0.5_dp*lnrat(-s,musq)**2

      return
      end

