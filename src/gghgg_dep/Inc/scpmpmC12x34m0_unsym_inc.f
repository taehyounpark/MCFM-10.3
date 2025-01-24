!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      real(dp):: s1234,s12,s13,s14,s23,s24,s34,Delta,Pi_mom,zab3123
      complex(dp):: scpmpmC12x34m0_unsym_res,zab2,zab4123,zab3142,
     & zab1234,zab3124,zab1342,zab2341

c  Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      Pi_mom(p1,p2,p3,p4)=s(p1,p3)+s(p2,p3)-s(p1,p4)-s(p2,p4)

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s1234 = s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4)+s(p3,p4)
      zab3123 = s13+s23

      zab1234 = zab2(p1,p2,p3,p4)
      zab3124 = zab2(p3,p1,p2,p4)
      zab3142 = zab2(p3,p1,p4,p2)
      zab1342 = zab2(p1,p3,p4,p2)
      zab2341 = zab2(p2,p3,p4,p1)
      zab4123 = zab2(p4,p1,p2,p3)
      Delta = (s1234-s12-s34)**2-4._dp*s12*s34

      scpmpmC12x34m0_unsym_res=
     & 2._dp*(za(p2,p3)**3*zb(p3,p4)*zab3123
     & *(zab3123*zb(p2,p3)-zb(p1,p2)*za(p1,p4)*zb(p3,p4)))
     & /(za(p1,p2)*zab1342*zab3124**3)

     & +2._dp*(za(p2,p3)**2*zb(p3,p4)*zab4123*(-2*s23-s24))
     & /(za(p1,p2)*zab1342*zab3124**2)

     & +2._dp*(zb(p1,p2)*za(p2,p3)**2*zb(p3,p4)
     & *(2*s12*(s23-s14-s34)+2*s13*s23
     & +2*s23**2+s14*s34-s23*s34
     & +2*zb(p1,p2)*za(p1,p3)*za(p2,p4)*zb(p3,p4)
     & ))/(zab1342*zab3124**3)

     & +2._dp*(zb(p1,p4)**2
     & *(za(p1,p4)*za(p2,p4)*(2*(s13-s24)-3*(s34+s14)-4*(s12+s23))
     & -2*za(p1,p3)*zb(p2,p3)*za(p2,p4)**2
     & +3*zb(p1,p3)*za(p1,p4)**2*za(p2,p3)))
     & /(zab1342*zab3124**2)

     & +s14**2*(s12*(6*s13-2*s14+2*s23+2*s24)-s14**2+s23**2)
     & /(zab1342**2*zab3124**2)

     &-8._dp*(s12*s13*s14*s23)
     & /(zab1342**2*zab3124**2)

     &+4._dp*(s14*s1234*zab2341*zab4123)/(zab1342*zab3124*Delta)
     &+4._dp*(za(p1,p2)*zb(p1,p3)*zab2341*zab4123*zab3142)
     & /(zab1342*zab3124*Delta)

     & +(zab1234*zab2341*zab3142*zab4123
     & *(Pi_mom(p4,p3,p2,p1)*Pi_mom(p1,p2,p3,p4)+Delta))
     & /(zab1342**2*zab3124**2*Delta)

     & -3._dp*(s1234*zab2341*zab4123
     & *Pi_mom(p4,p3,p2,p1)*Pi_mom(p1,p2,p3,p4)*(s13+s14+s23+s24))
     & /(2._dp*zab1342*zab3124*Delta**2)

     & +5._dp*(s1234*zab2341*zab4123*(s13+s14+s23+s24))
     & /(2._dp*zab1342*zab3124*Delta)

     & -4._dp*(zab2341*zab4123)/(zab1342*zab3124)
      return
