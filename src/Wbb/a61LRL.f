!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function a61LRL(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61LRL

c---  Corresponds to all outgoing
c     q(j1,-)+Q(j3,+)+e(j6,-)+q~(j4)+Q~(j2)+e~(j5)
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a61
      real(dp):: prop
      prop=s(5,6)/sqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
      a61LRL=a61('pp',j1,j2,j3,j4,j5,j6,zb,za)*prop
      return
      end

