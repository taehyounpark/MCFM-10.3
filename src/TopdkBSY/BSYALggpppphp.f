!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYALggpppphp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYALggpppphp

c-----Authors: John Campbell and Keith Ellis, March 2012
c---- arXiv:1101.5947 [hep-ph], Eq. (92), fully Badger-compliant
c---- higher-point (3- and 4-point) contributions only
c---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'massiveintegrals.f'
      real(dp):: mt3
      integer:: e1,p2,p3,e4,j

c-----setup variable controlling integrals to be used,
c-----depending on whether p2=2 or 3
      j=p2-1

      mt3=mt**3
      BSYALggpppphp=
     & +I41x2x3x4(j)*mt3*za(e1,e4)*zb(p2,p3)**2

      return
      end
