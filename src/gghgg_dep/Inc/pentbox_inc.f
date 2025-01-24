!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp):: C(5),detSx16,s12,s13,s14,s23,s24,s34,s234,s123,s1234,
     & Del4

      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)

      s234=s23+s24+s34
      s123=s12+s13+s23
      s1234=s12+s13+s14+s23+s24+s34

      Del4=(s12*s34-s13*s24-s14*s23)**2-4._dp*s13*s14*s23*s24
      detSx16=s12*s23*s34*(s14*s23-(s12+s13)*(s24+s34))+mtsq*Del4

      C(1)=s23*s34
     & *(s12*(s234-s23)+s123*(s234-s34)-s23*(s1234-s34))
      C(2)=s34*(s1234*s23*(s123-2*s12)
     & +s123*(s34*(s123-s23)+s12*(s234+s23)-s234*s123))
      C(3)=(s1234*s23-s234*s123)
     & *(s12*s234+s34*s123-s234*s123+s23*(s1234-s12-s34))
      C(4)=s12*(s1234*s23*(s234-2*s34)
     & +s234*(s12*(s234-s23)+s34*(s123+s23)-s234*s123))
      C(5)=s12*s23
     & *(s34*(s123-s23)+s234*(s123-s12)-s23*(s1234-s12))

      C(:)=-0.5_dp*C(:)/detSx16

      return
