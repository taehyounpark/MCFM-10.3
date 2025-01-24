!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA0qqppmp

c-----Authors: John Campbell and Keith Ellis, March 2012
c---- arXiv:1101.5947 [hep-ph], Eq. (A9), fully Badger compliant
c---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'mxpart.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      integer:: e1,p2,p3,e4
      complex(dp):: zabe4q4p2

      zabe4q4p2=-zab(e4,q1,p2)-za(e4,p3)*zb(p3,p2)

      BSYA0qqppmp=mt*(za(e1,p3)*zabe4q4p2+za(e4,p3)*zab(e1,q1,p2))
     & /(s(p2,p3))

      return
      end
