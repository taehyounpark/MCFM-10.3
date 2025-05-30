!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function TWbbpZab(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWbbpZab
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
c     Author: R.K. Ellis Feb, 2013
c     written by program WbbZdiags.frm
c     These are the Abelian diagrams with both the W and Z coming
c     off light line, only three diagrams with Z emitted before W
c     Calculation is performed for LH light-line (perforce because of W)
c     Calculation is performed for RH bbbar-line
c     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s34,s56,s134,s234,s567,s568,s5678
      complex(dp):: d9p,d10p
c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
c     end statement functions
      s34=s(p3,p4)
      s56=s(p5,p6)
      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s567=s3(p5,p6,p7)
      s568=s3(p5,p6,p8)
      s5678=s4(p5,p6,p7,p8)
      d9p= + s34**(-1)*s56**(-1)*s134**(-1)*s567**(-1) * ( za(p1,p3)*
     &    za(p2,p8)*za(p5,p6)*zb(p1,p4)*zb(p1,p6)*zb(p6,p7)*s5678**(-1)
     &     + za(p1,p3)*za(p2,p8)*za(p5,p7)*zb(p1,p4)*zb(p1,p7)*zb(p6,p7
     &    )*s5678**(-1) - za(p2,p8)*za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p4
     &    ,p6)*zb(p6,p7)*s5678**(-1) - za(p2,p8)*za(p3,p4)*za(p5,p7)*
     &    zb(p1,p4)*zb(p4,p7)*zb(p6,p7)*s5678**(-1) )
      d9p = d9p + s34**(-1)*s56**(-1)*s134**(-1)*s568**(-1) * ( za(p1,
     &    p3)*za(p2,p5)*za(p5,p8)*zb(p1,p4)*zb(p1,p7)*zb(p5,p6)*
     &    s5678**(-1) - za(p1,p3)*za(p2,p8)*za(p5,p8)*zb(p1,p4)*zb(p1,
     &    p7)*zb(p6,p8)*s5678**(-1) - za(p2,p5)*za(p3,p4)*za(p5,p8)*zb(
     &    p1,p4)*zb(p4,p7)*zb(p5,p6)*s5678**(-1) + za(p2,p8)*za(p3,p4)*
     &    za(p5,p8)*zb(p1,p4)*zb(p4,p7)*zb(p6,p8)*s5678**(-1) )

      d10p= + s34**(-1)*s56**(-1)*s234**(-1)*s567**(-1) * (  - za(p2,p3
     &    )*za(p2,p8)*za(p5,p6)*zb(p1,p6)*zb(p2,p4)*zb(p6,p7)*
     &    s5678**(-1) - za(p2,p3)*za(p2,p8)*za(p5,p7)*zb(p1,p7)*zb(p2,
     &    p4)*zb(p6,p7)*s5678**(-1) - za(p2,p3)*za(p3,p8)*za(p5,p6)*zb(
     &    p1,p6)*zb(p3,p4)*zb(p6,p7)*s5678**(-1) - za(p2,p3)*za(p3,p8)*
     &    za(p5,p7)*zb(p1,p7)*zb(p3,p4)*zb(p6,p7)*s5678**(-1) )
      d10p = d10p + s34**(-1)*s56**(-1)*s234**(-1)*s568**(-1) * (  -
     &    za(p2,p3)*za(p2,p5)*za(p5,p8)*zb(p1,p7)*zb(p2,p4)*zb(p5,p6)*
     &    s5678**(-1) + za(p2,p3)*za(p2,p8)*za(p5,p8)*zb(p1,p7)*zb(p2,
     &    p4)*zb(p6,p8)*s5678**(-1) - za(p2,p3)*za(p3,p5)*za(p5,p8)*zb(
     &    p1,p7)*zb(p3,p4)*zb(p5,p6)*s5678**(-1) + za(p2,p3)*za(p3,p8)*
     &    za(p5,p8)*zb(p1,p7)*zb(p3,p4)*zb(p6,p8)*s5678**(-1) )

      TWbbpZab=d9p+d10p
      return
      end
