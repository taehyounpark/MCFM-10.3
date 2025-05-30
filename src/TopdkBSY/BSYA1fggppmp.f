!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BSYA1fggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYA1fggppmp

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- arXiv:1101.5947 [hep-ph], (limit of) Eq. (101)
c---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      integer:: e1,p2,p3,e4

      BSYA1fggppmp=czip
      return
      end

