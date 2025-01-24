!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      complex(dp)::zab2,zaa22,zbb22
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaa22(p1,p2,p3,p4,p5,p6)=zab2(p1,p2,p3,p4)*za(p4,p6)
     &                        +zab2(p1,p2,p3,p5)*za(p5,p6)
      zbb22(p1,p2,p3,p4,p5,p6)=zab2(p4,p2,p3,p1)*zb(p4,p6)
     &                        +zab2(p5,p2,p3,p1)*zb(p5,p6)
