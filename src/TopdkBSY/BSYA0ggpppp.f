!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYA0ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA0ggpppp

c-----Authors: John Campbell and Keith Ellis, March 2012
c---- arXiv:1101.5947 [hep-ph], Eq. (A5), fully Badger-compliant
c---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      integer:: e1,p2,p3,e4
      BSYA0ggpppp=mt**3*za(e1,e4)*zb(p2,p3)/(za(p2,p3)*zab(p2,q1,p2))

      return
      end
