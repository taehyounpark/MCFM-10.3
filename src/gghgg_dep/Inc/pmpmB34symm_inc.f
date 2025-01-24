!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s12,s13,s14,s23,s34,s24,s123,s124,s134,s234,s1234,Delta
      complex(dp):: pmpmB34symm_res,zab2

      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s123=s12+s13+s23
      s124=s12+s14+s24
      s134=s13+s14+s34
      s234=s23+s24+s34
      s1234=s12+s13+s14+s23+s24+s34

      Delta=(s1234-s12-s34)**2-4._dp*s12*s34

      pmpmB34symm_res =
     &  4._dp*zb(p1,p3)*za(p1,p4)**2*s134
     & /(za(p1,p3)*zab2(p1,p3,p4,p1)*zab2(p1,p3,p4,p2)**2)
     & +4._dp*zb(p1,p4)**2*za(p3,p4)*zab2(p1,p3,p4,p1)
     & *zab2(p4,p1,p2,p3)*(two*s123+s124)
     &  /(zb(p1,p2)*zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)**2*Delta)
     &+12._dp*za(p1,p2)*zb(p1,p4)**2*za(p3,p4)*zab2(p4,p1,p2,p3)*s134
     &  /(zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)**2*Delta)
     & +4._dp*zab2(p4,p1,p2,p3)*zab2(p4,p1,p3,p4)
     & /(zab2(p1,p3,p4,p2)**2*zab2(p3,p1,p2,p4))
     & +8._dp*zb(p1,p3)*za(p3,p4)*zab2(p1,p2,p3,p4)*zab2(p4,p1,p2,p3)
     & *(s234-s134)/(zab2(p1,p3,p4,p2)**2*zab2(p3,p1,p2,p4)*Delta)
     &-2._dp*zb(p1,p4)*za(p2,p3)*zab2(p4,p1,p2,p3)
     & /(zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)**2)
     &+4._dp*zb(p1,p4)*za(p2,p3)*zab2(p1,p3,p4,p1)*zab2(p4,p1,p2,p3)
     & *(s123-s124)/(zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)**2*Delta)
     &+3._dp*zab2(p2,p3,p4,p1)*zab2(p4,p1,p2,p3)*(s123-s124)
     & *(s234-s134)*(s134+s234)
     & /(zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)*Delta**2)
     & +zab2(p2,p3,p4,p1)*zab2(p4,p1,p2,p3)
     &  *(s13-5._dp*s14-5._dp*s23+s24-14._dp*s34)
     &  /(zab2(p1,p3,p4,p2)*zab2(p3,p1,p2,p4)*Delta)

      return


