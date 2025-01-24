!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_c7x34m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.34,4.36)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp)::qloop_c7x34m2,qloop_c7x34m2_sym,qloop_c7x34m2_antisym

      qloop_c7x34m2=(+qloop_c7x34m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &               +qloop_c7x34m2_sym(p5,p6,p3,p4,p1,p2,p7,za,zb)
     &               +qloop_c7x34m2_antisym(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &               -qloop_c7x34m2_antisym(p5,p6,p3,p4,p1,p2,p7,za,zb))
      return
      end


      function qloop_c7x34m2_sym(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.35)
c symmetric part under 12 <-> 56
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c7x34m2_sym
      real(dp):: s3,s12x56,s347,s567,s134
      complex(dp)::zab2,zba2,zbab22,zaba22,
     & za7x1256x7,za7x1256x1,za1x3426x5,zb7x1256x7,zb2x3415x6

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions


      za7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      za7x1256x1=zaba22(p7,p1,p2,p5,p6,p1)
      zb7x1256x7=zbab22(p7,p1,p2,p5,p6,p7)
      za1x3426x5=zaba22(p1,p3,p4,p2,p6,p5)
      zb2x3415x6=zbab22(p2,p3,p4,p1,p5,p6)
      s12x56=s(p1,p5)+s(p1,p6)+s(p2,p5)+s(p2,p6)
      s347=s3(p3,p4,p7)
      s567=s3(p5,p6,p7)
      s134=s3(p1,p3,p4)

      qloop_c7x34m2_sym=4*s(p3,p4)*za(p1,p2)*zb(p5,p6)
     & *zb(p2,p7)**2*za(p3,p7)*zb(p4,p7)*za(p5,p7)*zab2(p5,p3,p4,p7)
     & /(zab2(p7,p3,p4,p7)*za7x1256x7*zb7x1256x7)

     & -4*za(p1,p7)*za(p3,p7)**2*za(p5,p7)*zb(p1,p2)*zb(p3,p4)
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *s12x56/za7x1256x7**3

     & -2*zb(p4,p7)**2*zb(p6,p7)*zab2(p1,p3,p4,p7)*zb(p1,p2)*za(p3,p4)
     & *za(p5,p6)
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & /(za7x1256x7*zb7x1256x7**2)

     & -za(p3,p7)*zb(p4,p7)
     & *(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)-zb(p6,p7)*za7x1256x1)
     & *s12x56
     & *(2*za1x3426x5+za(p1,p5)*(s347+s134+s(p1,p3)+s(p1,p4)))*zb(p1,p2)
     & /(za7x1256x7**2*zb7x1256x7)

     & -za(p3,p7)*za(p5,p7)*s12x56
     & *(za(p1,p3)*zb(p2,p6)*(s347+s(p2,p3)+s(p2,p4))
     & +za(p1,p3)*zb2x3415x6
     & -2*zb(p2,p7)*za(p3,p7)*zab2(p1,p2,p5,p6))*zb(p4,p3)
     & /(za7x1256x7**2)

     & -2*za(p3,p7)**2*(zb(p6,p5)*zab2(p5,p1,p2,p7)*za(p7,p1)
     & -zb(p6,p7)*za7x1256x1)
     & *(zb(p1,p2)*za(p1,p5)-zb(p2,p6)*za(p5,p6))*zb(p4,p3)
     & /(za7x1256x7**2)

     & +zb(p2,p7)*za(p3,p7)*zb(p4,p7)*zab2(p5,p3,p4,p7)
     & *(-za(p1,p5)*zb(p5,p6)*s12x56
     & +2*za(p1,p2)*zb(p2,p6)*s567)/(za7x1256x7*zb7x1256x7)

     & +zb(p2,p7)*zb(p4,p7)*zab2(p5,p3,p4,p7)
     & *(2*zb(p2,p3)*za(p3,p5)*za(p3,p7)
     & +4*zb(p2,p4)*za(p3,p7)*za(p4,p5)
     & -2*zb(p2,p4)*za(p3,p5)*za(p4,p7))*za(p1,p2)*zb(p5,p6)
     & /(za7x1256x7*zb7x1256x7)

     & -za(p1,p3)*za(p3,p7)
     & *(zb(p2,p6)*zab2(p5,p3,p4,p7)
     & +2*zb(p2,p7)*zab2(p5,p1,p2,p6))*zb(p4,p3)
     & /(za7x1256x7)

     & +za(p1,p2)*zb(p2,p6)*za(p3,p7)
     & *(2*zb(p2,p4)*zab2(p5,p3,p4,p7)+zb(p1,p2)*za(p1,p5)*zb(p4,p7))
     & /(za7x1256x7)

      qloop_c7x34m2_sym=qloop_c7x34m2_sym/(4*s(p1,p2)*s(p3,p4)*s(p5,p6))

      return
      end


      function qloop_c7x34m2_antisym(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.36)
c antisymmetric part under 12 <-> 56
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_c7x34m2_antisym
      real(dp):: s3,s127,s567
      complex(dp)::zab2,zba2,zaba22,za7x1256x7

c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      za7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      s127=s3(p1,p2,p7)
      s567=s3(p5,p6,p7)

      qloop_c7x34m2_antisym =
     & -za(p1,p5)*zb(p2,p4)*za(p3,p7)*zab2(p7,p3,p4,p7)
     & *(zab2(p7,p3,p4,p6)*s567-zb(p6,p7)*za7x1256x7)
     & /za7x1256x7**2
     & +za(p3,p7)*za(p5,p7)*zab2(p7,p3,p4,p7)*zb(p3,p4)
     & *(za(p1,p3)*zb(p2,p6)*s127+zb(p2,p4)*za(p3,p4)*zab2(p1,p3,p4,p6))
     & /za7x1256x7**2

      qloop_c7x34m2_antisym=
     & qloop_c7x34m2_antisym/(4*s(p1,p2)*s(p3,p4)*s(p5,p6))

      return
      end
