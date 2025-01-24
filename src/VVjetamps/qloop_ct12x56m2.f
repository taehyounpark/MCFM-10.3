!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_ct12x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.28)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_ct12x56m2
      real(dp)::s12,s56,s34,s347,delta
      complex(dp)::zab2,zba2,zaba22,zbab22,za7x5634x7,zb7x5634x7,
     & zb6x5612x7,zb6x1256x7,zb2x3456x7

c     statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
c     end statement functions

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s347=s(p3,p4)+s(p3,p7)+s(p4,p7)
      delta=(s347-s12-s56)**2-4*s12*s56
      za7x5634x7=zaba22(p7,p5,p6,p3,p4,p7)
      zb7x5634x7=zbab22(p7,p5,p6,p3,p4,p7)
      zb6x5612x7=zbab22(p6,p5,p6,p1,p2,p7)
      zb6x1256x7=zbab22(p6,p1,p2,p5,p6,p7)
      zb2x3456x7=zbab22(p2,p3,p4,p5,p6,p7)

      qloop_ct12x56m2=-za(p5,p3)*zb(p6,p4)
     & *(-za(p1,p7)*zab2(p7,p3,p4,p2)*delta
     & +zab2(p1,p5,p6,p2)*za7x5634x7
     & *(s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6)))
     & /(s12*za7x5634x7**2)

     & -za(p5,p3)*zb(p4,p7)*(zb2x3456x7*za(p7,p1)
     & -zb(p2,p7)*zaba22(p7,p5,p6,p3,p4,p1))
     & *(zb6x5612x7-zb6x1256x7)
     & /(s12*za7x5634x7*zb7x5634x7)

     & +zb(p2,p6)*za(p5,p3)
     & *(-2*zb(p2,p4)*(s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6))
     & -2*zb(p3,p4)*zab2(p3,p5,p6,p2)+2*zb(p2,p4)*za(p5,p6)*zb(p5,p6))
     & /(zb(p1,p2)*za7x5634x7)

     & +za(p1,p5)*za(p5,p3)*zb(p5,p6)
     & *(2*zab2(p1,p2,p7,p4)-4*za(p1,p2)*zb(p2,p4))
     & /(za(p1,p2)*za7x5634x7)
      qloop_ct12x56m2=qloop_ct12x56m2/(s56*s34)

      return
      end
