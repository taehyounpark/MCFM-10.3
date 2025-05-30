!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
c     This is a calculation for the process 1_g^+ 2_g^- 3_g^+ 4_g^- +H
c     Implementation of Eq.~(6.15) from arXiv:2002.04018 v2
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s24,s34,s234,s12ps13ps14
      complex(dp):: pmpmC1x234m0_res,zab2,zab1243,zab1342,zab1234

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s234=s23+s24+s34
      s12ps13ps14=s12+s13+s14

      zab1243=zab2(p1,p2,p4,p3)
      zab1342=zab2(p1,p3,p4,p2)
      zab1234=zab2(p1,p2,p3,p4)
      pmpmC1x234m0_res=s12ps13ps14
     & *(-2._dp*s234*zab1243**2
     & *(zb(p2,p3)**2*zab1234**2+zb(p3,p4)**2*zab1342**2)
     & /(zab1234**3*zab1342**3*zb(p2,p3)*zb(p3,p4))

     & +zb(p1,p3)**2*zb(p3,p4)
     & /(zab1234*zb(p1,p4)*zb(p2,p3)*zb(p2,p4))
     & +zb(p1,p3)**2*zb(p3,p2)
     & /(zab1342*zb(p1,p2)*zb(p4,p3)*zb(p4,p2))

     & +za(p2,p4)**2*za(p1,p4)
     & /(zab1342*za(p1,p3)*za(p1,p2)*za(p3,p4))
     & +za(p2,p4)**2*za(p1,p2)
     & /(zab1234*za(p1,p3)*za(p1,p4)*za(p3,p2)))

     & -2._dp*zb(p1,p3)**4/(zb(p1,p2)*zb(p1,p4)*zb(p3,p2)*zb(p3,p4))
      return
