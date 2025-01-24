!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function A1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A1fggpppp

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- taken from arXiv:1101.5947 [hep-ph], Eq. (100)
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: BSYA1fggpppp
      integer:: e1,p2,p3,e4

      A1fggpppp=BSYA1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
