!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_ct34x56m2(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.24)
c     symmetric part under 3<-->5,4 <-> 6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_ct34x56m2
      real(dp)::s34,s56,s127,Delta
      complex(dp)::zba2,zaba22,za7x1256x7,za1x5627x1
c     statement functions
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions
      za7x1256x7=zaba22(p7,p1,p2,p5,p6,p7)
      za1x5627x1=zaba22(p1,p5,p6,p2,p7,p1)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      Delta=(s127-s34-s56)**2-4*s34*s56
      qloop_ct34x56m2=za(p3,p5)*zb(p4,p6)
     & *(za(p1,p7)**2*Delta-2*za1x5627x1*za7x1256x7)
     & /(za(p1,p2)*s34*s56*za7x1256x7**2)

      return
      end
