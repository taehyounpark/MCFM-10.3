!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYA0qedpppm(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA0qedpppm

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- arXiv:1101.5947 [hep-ph], Eq. (A8)
c---- These are the twiddle functions arXiv:1101.5947 [hep-ph], Eq. (91)
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      integer:: e1,p2,p3,e4
      BSYA0qedpppm=-mt*zab(p3,q1,p2)
     & *(za(e1,e4)*zab(p3,q1,p2)-zb(p2,p3)*za(p3,e1)*za(p3,e4))
     & /(zab(p2,q1,p2)*zab(p3,q1,p3))

c--- NB:flip sign to agree with BSY code
      BSYA0qedpppm=-BSYA0qedpppm

      return
      end
