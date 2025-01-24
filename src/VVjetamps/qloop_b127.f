!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qloop_b127(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
c     This is an implementation of deLaurentis, Campbell and Ellis
c     arXiv:2203.17170, Eq(4.41)
c     this function returns
c     part of qloop_b127 symmetric under (3,4)<-->(5,6)
c     i.e. +(qloop_b127(1,2,3,4,5,6,7)
c           +qloop_b127(1,2,5,6,3,4,7))
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: qloop_b127,qloop_b127_symbit
      qloop_b127=
     &  +qloop_b127_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
     &  +qloop_b127_symbit(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function qloop_b127_symbit(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7
      real(dp)::s3,s127,s356,s456,s17p27,s34,s56
      complex(dp)::qloop_b127_symbit,zab2,zba2,zbab22,zaba22,
     & zbab77,zaba77,zaba11,zab3564,
     & zab7124,zab5346,Delta
c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      zbab22(p1,p2,p3,p4,p5,p6)=
     & zb(p1,p2)*zab2(p2,p4,p5,p6)+zb(p1,p3)*zab2(p3,p4,p5,p6)
      zaba22(p1,p2,p3,p4,p5,p6)=
     & za(p1,p2)*zba2(p2,p4,p5,p6)+za(p1,p3)*zba2(p3,p4,p5,p6)
c     end statement functions

      s456=s3(p4,p5,p6)
      s356=s3(p3,p5,p6)
      s127=s3(p1,p2,p7);
      s17p27=s(p1,p7)+s(p2,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      Delta=(s127-s34-s56)**2-4*s34*s56
      zaba11=zaba22(p1,p5,p6,p2,p7,p1)
      zab3564=zab2(p3,p5,p6,p4)
      zab7124=zab2(p7,p1,p2,p4)
      zab5346=zab2(p5,p3,p4,p6)
      zbab77= zbab22(p7,p3,p4,p5,p6,p7)
      zaba77= zaba22(p7,p3,p4,p5,p6,p7)

      qloop_b127_symbit=
     &  -(za(p1,p7)*za(p5,p7)*zab7124)**2*s127*zbab77
     & /(za(p1,p2)*zb(p3,p4)*za(p5,p6)*s17p27**2*zaba77**3)

     & +zaba11*(zb(p4,p6)*za(p5,p6)+zb(p3,p4)*za(p3,p5))**2
     & /(2*delta*za(p1,p2)*zb(p3,p4)*za(p5,p6)*zaba77)

     & -(za(p1,p3)*za(p1,p7)*(s456-s356)*zab5346*zab7124)
     & /(2*delta*za(p1,p2)*zaba77**2)
     & +(zab2(p1,p2,p7,p4)*za(p1,p7)*(s127-s56)*zab5346*zab7124)
     & /(2*delta*za(p1,p2)*zb(p3,p4)*zaba77**2)

     & -(zaba22(p1,p3,p5,p1,p2,p7)*za(p1,p7)*zab3564*zab5346)
     & /(2*delta*za(p1,p2)*zaba77**2)

     & -(za(p1,p5)*za(p1,p7)*za(p5,p7)*s127*(s127-s34)*zab3564)
     & /(2*delta*za(p1,p2)*za(p5,p6)*zaba77**2)

     & +za(p1,p7)*za(p5,p7)*zab7124
     & *(zb(p2,p7)*zb(p4,p6)*za(p5,p6)-zb(p1,p2)*za(p1,p5)*zb(p4,p7)
     & +zb(p2,p7)*zb(p3,p4)*za(p3,p5))
     & /(2*zb(p3,p4)*za(p5,p6)*s17p27*zaba77**2)

     & +(3*za(p1,p7)*zb(p2,p7)*zb(p4,p7)*za(p5,p7)**2*s127*zab7124)
     & /(2*zb(p3,p4)*za(p5,p6)*s17p27**2*zaba77**2)

     & +(za(p1,p5)*za(p1,p7)*zb(p4,p7)*za(p5,p7)*zab7124)
     & /(2*za(p1,p2)*zb(p3,p4)*za(p5,p6)*zaba77**2)

     & -3*(za(p1,p7)**2*delta-4*zaba11*zaba77)*s127*zab3564*zab5346
     & /(4*delta**2*za(p1,p2)*zaba77**2)

      return
      end
