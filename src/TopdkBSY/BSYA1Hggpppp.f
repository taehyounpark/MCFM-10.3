!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYA1Hggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1Hggpppp

c-----Authors: John Campbell and Keith Ellis, March 2012
c---- arXiv:1101.5947 [hep-ph], Eq. (100), fully Badger-compliant
c---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      include 'massiveintegrals.f'
      real(dp):: s23,mt2
      integer:: e1,p2,p3,e4

      s23=s(p2,p3)
      mt2=mt**2
      BSYA1Hggpppp=
     & -2d0*mt*(za(e1,e4)*zab(p2,q1,p2)-za(p2,e1)*za(p3,e4)*zb(p2,p3))
     & /(za(p2,p3)**3*zb(p2,p3))
     & *(s23*mt2*I3m2x3x41+2d0*mt2*F2m23+s23/6d0)

      return
      end

