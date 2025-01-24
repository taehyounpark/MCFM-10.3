!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_c12x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.23)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp)::qloop_c12x56m2,qloop_c12x56m2_symbit,
     & qloop_c12x56m2_asybit

      qloop_c12x56m2=
     &          +qloop_c12x56m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &          +qloop_c12x56m2_symbit(p5,p6,p3,p4,p1,p2,p7,za,zb)
     &          +qloop_c12x56m2_asybit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &          -qloop_c12x56m2_asybit(p5,p6,p3,p4,p1,p2,p7,za,zb)

      return
      end


      function qloop_c12x56m2_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.26)
c symmetric part under 12 <-> 56
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c12x56m2_symbit
      real(dp)::s3,delta
      complex(dp)::zab2,zba2,zbab22,zaba22,
     & za7x1256x1,za7x1256x7,zb7x1256x7,zab1562,zab5126,
     & za3x5612x3,zb2x5612x7

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      delta=(s3(p3,p4,p7)-s(p1,p2)-s(p5,p6))**2-4*s(p1,p2)*s(p5,p6)
      za7x1256x1=zaba22(p7,p1,p2,p5,p6,p1)
      za7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      zb7x1256x7=zbab22(p7,p1,p2,p5,p6,p7)
      za3x5612x3=zaba22(p3,p5,p6,p1,p2,p3)
      zb2x5612x7=zbab22(p2,p5,p6,p1,p2,p7)
      zab1562=zab2(p1,p5,p6,p2)
      zab5126=zab2(p5,p1,p2,p6)
      qloop_c12x56m2_symbit=4*delta
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *(zb(p1,p2)*zb(p4,p7)*(za(p1,p3)*za(p5,p7)-za(p1,p7)*za(p3,p5))
     & +zb(p3,p4)*za(p3,p5)*zab2(p3,p6,p7,p2)
     & +zb(p3,p4)*za(p3,p5)*za(p3,p7)*zb(p7,p2))
     & /(zb7x1256x7*za7x1256x7**2)

     & +8*zb(p3,p4)*za(p3,p6)*zab5126
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *((s3(p5,p6,p1)-s3(p5,p6,p2))
     & *(za(p3,p5)*zb(p5,p2)-za(p3,p7)*zb(p7,p2))
     & -2*zab1562*(za(p3,p5)*zb(p5,p1)-za(p3,p7)*zb(p7,p1)))
     & /(zb7x1256x7*za7x1256x7**2)

     & +4*zb(p3,p4)*za(p3,p5)*(s3(p1,p2,p5)-s3(p1,p2,p6))
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *(2*zab2(p3,p6,p7,p1)*zab1562
     & -zab2(p3,p6,p7,p2)*(s3(p5,p6,p1)-s3(p5,p6,p2)))
     & /(zb7x1256x7*za7x1256x7**2)

     & -4*zb(p3,p4)*za(p1,p3)*zb(p2,p6)*za(p3,p7)*za(p5,p7)*delta
     & /(za7x1256x7**2)

     & -12*zb(p4,p7)*zab2(p3,p1,p2,p7)*zab1562*zab5126
     & *(s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6))
     & /(zb7x1256x7*za7x1256x7)

     & -4*za(p1,p3)*zb(p2,p7)*zb(p4,p7)*zab5126
     & *(s3(p5,p6,p1)-s3(p5,p6,p2))
     & *(s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p2,p3)
     & +s(p2,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6))
     & /(zb7x1256x7*za7x1256x7)

      qloop_c12x56m2_symbit=qloop_c12x56m2_symbit
     & +8*zb(p4,p7)*zab1562*zab5126
     & *(-zab2(p3,p1,p2,p7)*s(p1,p2)
     & +za(p1,p2)*zb(p1,p7)*zb(p2,p7)*za(p3,p7)
     & -za(p1,p3)*zb(p1,p7)*(s(p3,p7)+s(p4,p7)))
     & /(zb7x1256x7*za7x1256x7)

     & -2*za(p1,p3)*zb(p4,p3)*zb(p2,p7)*za(p3,p5)*zb(p6,p7)*delta
     & /(zb7x1256x7*za7x1256x7)

     & +4*za(p1,p3)*zb(p3,p4)*zab5126
     & *(zb(p2,p7)*(s3(p5,p6,p1)-s3(p5,p6,p2))
     & *(-zab2(p3,p4,p6,p7)+za(p3,p5)*zb(p5,p7))
     & +2*zb(p1,p7)*zab1562*zab2(p3,p4,p6,p7))
     & /(zb7x1256x7*za7x1256x7)

     & -2*za(p1,p3)*zb(p3,p4)*zb(p2,p7)*za(p3,p5)*zb(p6,p7)
     & *(s3(p1,p2,p5)-s3(p1,p2,p6))
     & *(s3(p5,p6,p1)-s3(p5,p6,p2))
     & /(zb7x1256x7*za7x1256x7)

      qloop_c12x56m2_symbit=qloop_c12x56m2_symbit
     & +16*zb(p2,p7)*zb(p4,p7)*zab5126
     & *(0.5_dp*za(p1,p2)*zb(p1,p4)*za(p3,p4)*zab1562
     & +0.5_dp*za(p1,p3)*za(p3,p4)*zb(p3,p4)
     & *(s3(p5,p6,p1)-s3(p5,p6,p2))-0.25_dp*za(p1,p2)*zab2(p3,p5,p6,p2)
     & *(s3(p5,p6,p1)-s3(p5,p6,p2)))
     & /(zb7x1256x7*za7x1256x7)

     & -8*za(p3,p7)**2*zb(p3,p4)*delta
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *(zb2x5612x7*za(p7,p5)-zb(p2,p7)*zab2(p7,p1,p2,p6)*za(p6,p5))
     & /(zb7x1256x7*za7x1256x7**3)

     & -4*zb(p4,p7)**2*za(p3,p4)*delta
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *(zb2x5612x7*za(p7,p5)-zb(p2,p7)*zab2(p7,p1,p2,p6)*za(p6,p5))
     & /(zb7x1256x7**2*za7x1256x7**2)

     & +4*za3x5612x3*zb(p3,p4)
     & *(2*zab1562*zab5126*(s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6))
     & +za(p1,p5)*zb(p2,p6)*delta)
     & /(za7x1256x7*delta)

      qloop_c12x56m2_symbit=qloop_c12x56m2_symbit
     & /(16*s(p1,p2)*s(p3,p4)*s(p5,p6))
      return
      end


      function qloop_c12x56m2_asybit(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.27)
c antisymmetric part under 12 <-> 56
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      complex(dp)::qloop_c12x56m2_asybit,zab2,zba2,zbab22,zaba22,
     & za7x1256x7,za7x5612x1,zb7x1256x7,za7x1256x1,
     & zb7x5612x4,zb7x1256x4,za7x5612x3,za7x1256x3
      integer::p1,p2,p3,p4,p5,p6,p7
      real(dp)::s3,delta

c      include 'texbit.f'
c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      delta=(s3(p3,p4,p7)-s(p1,p2)-s(p5,p6))**2-4*s(p1,p2)*s(p5,p6)
      za7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      zb7x1256x7=zbab22(p7,p1,p2,p5,p6,p7)
      za7x5612x1=zaba22(p7,p5,p6,p1,p2,p1)
      za7x1256x1=zaba22(p7,p1,p2,p5,p6,p1)
      zb7x5612x4=zbab22(p7,p5,p6,p1,p2,p4)
      zb7x1256x4=zbab22(p7,p1,p2,p5,p6,p4)
      za7x5612x3=zaba22(p7,p5,p6,p1,p2,p3)
      za7x1256x3=zaba22(p7,p1,p2,p5,p6,p3)
      qloop_c12x56m2_asybit =za(p3,p7)*delta
     & *(-2*za(p1,p5)*zb(p2,p6)*zab2(p7,p1,p2,p4)
     & -2*za(p1,p3)*zb(p2,p6)*zb(p3,p4)*za(p5,p7)
     & -4*zb(p1,p2)*za(p1,p5)*za(p1,p7)*zb(p4,p6))
     & /(za7x1256x7**2)

     & +za(p1,p3)*zb(p2,p4)*zab2(p5,p1,p2,p6)
     & *(-4*(s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6))
     & +8*(s3(p3,p4,p7)-s(p1,p2)))
     & /(za7x1256x7)

     & -za(p1,p5)*zb(p2,p7)*zb(p6,p7)
     & *(zb7x5612x4-zb7x1256x4)
     & *(za7x5612x3-za7x1256x3)
     & /(za7x1256x7*zb7x1256x7)

     & +2*za(p1,p5)*zb(p2,p7)*zb(p6,p7)*zab2(p3,p1,p2,p4)*delta
     & /(za7x1256x7*zb7x1256x7)

     & -4*zb(p2,p4)*za(p3,p5)*zb(p6,p7)
     & *(za7x5612x1-za7x1256x1)/(za7x1256x7)

      qloop_c12x56m2_asybit=qloop_c12x56m2_asybit
     & /(16*s(p1,p2)*s(p3,p4)*s(p5,p6))
      return
      end



