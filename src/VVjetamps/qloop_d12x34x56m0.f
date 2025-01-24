!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_d12x34x56m0(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.7)

      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d12x34x56m0,qloop_d12x34x56m0_symbit
      qloop_d12x34x56m0=
     &  qloop_d12x34x56m0_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     & +qloop_d12x34x56m0_symbit(p5,p6,p3,p4,p1,p2,p7,za,zb)
      return
      end

      function qloop_d12x34x56m0_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_d12x34x56m0_symbit,zab2,zaba22,
     & za7x5634x1,za7x5612x7
      real(dp)::s127,s567,s12,s56
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & zab2(p1,p2,p3,p4)*za(p4,p6)+zab2(p1,p2,p3,p5)*za(p5,p6)
      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s12=s(p1,p2)
      s56=s(p5,p6)

      za7x5634x1=zaba22(p7,p5,p6,p3,p4,p1)
      za7x5612x7=zaba22(p7,p5,p6,p1,p2,p7)

      qloop_d12x34x56m0_symbit=
     & -(za(p5,p7)*zab2(p7,p1,p2,p4)*za7x5634x1)**2
     & *(s127*s567-s12*s56)
     & /(4*zb(p3,p4)*za(p1,p2)*za(p5,p6)*za7x5612x7**4)
      return
      end
