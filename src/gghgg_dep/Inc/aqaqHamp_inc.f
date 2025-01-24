!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c     Implementation of Eq.~(11.6) from arXiv:2002.04018 v2
      use fillformfactor_generic
      implicit none
      include 'Inc/zprods_decl.f'

      integer::p1,p2,p3,p4
      real(dp)::mtsq
      complex(dp)::FL,FT,zab2,zba2,amp(2,2)
      real(dp)::ssum,sprod
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)

      sprod=s(p1,p2)*s(p3,p4)
      ssum=s(p2,p4)+s(p2,p3)+s(p1,p4)+s(p1,p3)
      call fillformfactor(p1,p2,p3,p4,mtsq,FL,FT)
      amp(1,1)=2._dp*FL*za(p2,p4)*zb(p1,p3)
     & +FT/sprod*(zab2(p2,p4,p3,p1)*zab2(p4,p2,p1,p3)
     & +za(p2,p4)*zb(p1,p3)*ssum)
      amp(2,2)=2._dp*FL*za(p1,p3)*zb(p2,p4)
     & +FT/sprod*(zba2(p2,p4,p3,p1)*zba2(p4,p2,p1,p3)
     & +za(p1,p3)*zb(p2,p4)*ssum)
      amp(1,2)=2._dp*FL*za(p2,p3)*zb(p1,p4)
     & +FT/sprod*(zab2(p2,p4,p3,p1)*zba2(p4,p2,p1,p3)
     & +za(p2,p3)*zb(p1,p4)*ssum)
      amp(2,1)=2._dp*FL*za(p1,p4)*zb(p2,p3)
     & +FT/sprod*(zab2(p4,p2,p1,p3)*zba2(p2,p4,p3,p1)
     & +za(p1,p4)*zb(p2,p3)*ssum)
      return
