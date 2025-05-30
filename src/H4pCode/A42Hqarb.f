!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- NB: phi-dagger amplitudes related directly to phi amplitudes
c---     using parity and charge conjugation
      function A42Hqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A42Hqarbmppm

      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A42phiqarbmppm

      A42Hqarbmppm=A42phiqarbmppm(j1,j2,j3,j4,za,zb)
     &            +A42phiqarbmppm(j2,j1,j4,j3,zb,za)

      return
      end

      function A42Hqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A42Hqarbmpmp

      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A42phiqarbmpmp

      A42Hqarbmpmp=A42phiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            +A42phiqarbmpmp(j2,j1,j4,j3,zb,za)

      return
      end
