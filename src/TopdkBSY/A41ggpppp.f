!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function A41ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A41ggpppp

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- taken from arXiv:1101.5947 [hep-ph], Eq. (4)
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: ALggpppp,ARggpppp,
     & BSYA1fggpppp,BSYA1Hggpppp
      integer:: e1,p2,p3,e4
      real(dp):: nlf,nhf
      parameter(nlf=5d0,nhf=1d0)

      A41ggpppp=ALggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     & -Ncinv**2*ARggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYA1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)*nlf/xn
     & -BSYA1Hggpppp(e1,p2,p3,e4,za,zb,zab,zba)*nhf/xn

      return
      end

