!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function F33m(p1sq,p2sq,p3sq)
      implicit none
      include 'types.f'
      complex(dp):: F33m
      include 'scalarselect.f'
c      include 'scale.f'
      real(dp):: p1sq,p2sq,p3sq
c      complex(dp):: loopI3
      complex(dp):: I3m

c--- NOTE: checked on 8/30/09 that loopI3 == -I3m
c---       and F33m is defined to be (-1)*(scalar integral)
c      F33m=-loopI3(p1sq,p2sq,p3sq,0._dp,0._dp,0._dp,musq,0)
      F33m=I3m(p1sq,p2sq,p3sq)

      return
      end

