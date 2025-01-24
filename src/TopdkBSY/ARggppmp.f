!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function ARggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: ARggppmp

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      real(dp):: mt2,xlog,beta
      complex(dp):: BSYA0ggppmp,BSYARggppmp,VR
      integer:: e1,p2,p3,e4

      mt2=mt**2
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
      VR=cplx1(0.5d0*epinv)
     & -cplx1(epinv)*(1d0-2d0*mt2/s(p2,p3))*xlog
      ARggppmp=VR*BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &           +BSYARggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
