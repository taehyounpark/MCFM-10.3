!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function ALslcggpppm(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: ALslcggpppm

c-----Authors: John Campbell and Keith Ellis, November 2011
c---- arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2,xlog,beta
      complex(dp):: BSYA0qedpppm,BSYALslggpppm,Vslc,lnrat
      integer:: e1,p2,p3,e4

      mt2=mt**2
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
      Vslc=cplx1(epinv)*
     & (+s(1,p2)/s(p2,p3)*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     &  +s(1,p3)/s(p2,p3)*(lnrat(musq,-s(1,p3))+lnrat(mt2,-s(1,p3)))
     & +lnrat(musq,-s(p2,p3))+(1d0-2d0*mt2/s(p2,p3))*xlog)

      ALslcggpppm=Vslc*BSYA0qedpppm(e1,p2,p3,e4,za,zb,zab,zba)
     &                +BSYALslggpppm(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
