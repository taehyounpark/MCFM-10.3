!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      zbaba222(p1,p2,p3,p4,p5,p6,p7,p9)=
     & +zb(p1,p2)*zaba22(p2,p4,p5,p6,p7,p9)
     & +zb(p1,p3)*zaba22(p3,p4,p5,p6,p7,p9)
      Gammat(p7,p9)=
     & zbaba222(p7,p1,p2,p3,p4,p5,p6,p9)-zbaba222(p7,p5,p6,p3,p4,p1,p2,p9)
      Delta3(p1,p2,p3,p4)=(s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4))**2
     & -4*s(p1,p2)*s(p3,p4)
      ss3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
