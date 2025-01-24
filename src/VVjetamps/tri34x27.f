!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tri34x27(p1,p2,p3,p4,p5,p6,p7,za,zb,tot,tot_sl)
      implicit none
      include'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      real(dp):: s156
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: tot,tot_sl
      complex(dp)::zab2,zaa22,zbb22
c     statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaa22(p1,p2,p3,p4,p5,p6)=zab2(p1,p2,p3,p4)*za(p4,p6)
     &                        +zab2(p1,p2,p3,p5)*za(p5,p6)
      zbb22(p1,p2,p3,p4,p5,p6)=zab2(p4,p2,p3,p1)*zb(p4,p6)
     &                        +zab2(p5,p2,p3,p1)*zb(p5,p6)
c     end statement functions


      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)

      tot=-0.5_dp*za(p1,p5)**2*za(p3,p4)*zab2(p1,p2,p7,p4)
     & *zbb22(p4,p2,p7,p3,p4,p2)
     & /(za(p2,p7)*za(p5,p6)*zab2(p7,p3,p4,p2)
     & *zaa22(p1,p5,p6,p2,p7,p1)*s156)
     & -0.5_dp*za(p1,p5)**2*zab2(p1,p2,p7,p4)*zbb22(p4,p2,p7,p3,p4,p2)
     & /(za(p2,p7)*zb(p3,p4)*za(p5,p6)*zab2(p7,p3,p4,p2)
     & *zaa22(p1,p5,p6,p2,p7,p1))
     & -0.5_dp*za(p1,p3)*(za(p1,p5)*zb(p2,p7))**2
     & *zaa22(p3,p2,p7,p3,p4,p7)
     & /(za(p3,p4)*za(p5,p6)*zab2(p7,p3,p4,p2)
     & *zaa22(p1,p5,p6,p2,p7,p1)*s156)
     & +(za(p1,p5)**2*zb(p2,p7)*za(p3,p7)*(zb(p2,p4)*zab2(p1,p3,p4,p7)
     & +0.5_dp*zb(p2,p7)*zab2(p1,p5,p6,p4)))
     & /(za(p5,p6)*zab2(p7,p3,p4,p2)*zaa22(p1,p5,p6,p2,p7,p1)*s156)
     & +za(p1,p5)**2*zb(p2,p4)*zb(p2,p7)
     & *(-0.5_dp*(s(p1,p2)+s(p2,p5)+s(p2,p6))*za(p1,p3)
     & -za(p1,p2)*zb(p2,p4)*za(p3,p4))
     & /(za(p5,p6)*zab2(p7,p3,p4,p2)*zaa22(p1,p5,p6,p2,p7,p1)*s156)
     & -1.5_dp*za(p1,p3)*za(p1,p5)**2*zb(p2,p7)*zb(p4,p7)
     & /(za(p5,p6)*zaa22(p1,p5,p6,p2,p7,p1)*s156)

      tot_sl=tot

      return
      end

