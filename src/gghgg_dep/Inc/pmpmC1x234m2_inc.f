!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.16) from arXiv:2002.04018 v2


      include 'Inc/zprods_decl.f'
      complex(dp)::pmpmC1x234m2_res
      integer p1,p2,p3,p4
      real(dp)::s12p13p14
      complex(dp)::zab2,zab1243,zab1342,zab1234

      zab2(p1,p2,p3,p4)=(za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4))

      s12p13p14=s(p1,p2)+s(p1,p3)+s(p1,p4)
      zab1243=zab2(p1,p2,p4,p3)
      zab1342=zab2(p1,p3,p4,p2)
      zab1234=zab2(p1,p2,p3,p4)

      pmpmC1x234m2_res=
     & +4._dp*s12p13p14/(zab1342*zab1234)
     & *((za(p2,p4)**2/(za(p2,p3)*za(p3,p4))
     &   -zb(p1,p3)**2/(zb(p1,p2)*zb(p1,p4)))
     & +zab1243/(zab1342)*(
     &  za(p1,p4)*za(p2,p4)/(za(p1,p2)*za(p3,p4))
     & +zb(p1,p3)*zb(p2,p3)/(zb(p1,p2)*zb(p3,p4)))
     & +zab1243/(zab1234)*(
     &  za(p1,p2)*za(p2,p4)/(za(p1,p4)*za(p2,p3))
     & +zb(p1,p3)*zb(p3,p4)/(zb(p1,p4)*zb(p2,p3))))
     & -(8._dp*zb(p1,p3)**4)
     & /(zb(p1,p2)*zb(p1,p4)*zb(p2,p3)*zb(p3,p4)*s12p13p14)

      return

