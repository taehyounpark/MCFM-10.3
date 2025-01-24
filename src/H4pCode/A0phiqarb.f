!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function A0phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiqarbmppm

c----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phiqarbmppm=-za(j1,j4)**2/(za(j1,j2)*za(j3,j4))
      return
      end


      function A0phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiqarbmpmp

c----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiqarbmppm
      A0phiqarbmpmp=-A0phiqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end

      function A0phidqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidqarbmppm

c----Expresssion of Eq. (A.4) of arXiv:0906008
c----with factor of i removed
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      A0phidqarbmppm=-zb(j2,j3)**2/(zb(j1,j2)*zb(j3,j4))
      return
      end

      function A0phidqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidqarbmpmp

c----Expresssion of Eq. (A.4) of arXiv:0906008
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: A0phidqarbmppm
      integer:: j1,j2,j3,j4
      A0phidqarbmpmp=-A0phidqarbmppm(j1,j2,j4,j3,za,zb)
      return
      end
