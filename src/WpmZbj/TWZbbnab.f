!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function TWZbbnab(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWZbbnab
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
c     Author: R.K. Ellis Feb, 2013
c     written by program WZbbdiags.frm
c     These are the Nabelian diagrams with both the W and Z coming
c     off light line, only three diagrams with Z emitted before W
c     Calculation is performed for LH light-line (perforce because of W)
c     Calculation is performed for LH bbbar-line
c     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s34,s56,s78,s178,s278,s3456
      complex(dp):: d7,d8
c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
c     end statement functions
      s34=s(p3,p4)
      s56=s(p5,p6)
      s78=s(p7,p8)
      s278=s3(p2,p7,p8)
      s178=s3(p1,p7,p8)
      s3456=s4(p3,p4,p5,p6)
      d7= + s3456**(-1)*s34**(-1)*s56**(-1)*s78**(-1)*s178**(-1) * (
     &     - za(p1,p7)*za(p2,p3)*za(p3,p5)*zb(p1,p3)*zb(p1,p8)*zb(p4,p6
     &    ) + za(p1,p7)*za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zb(p3,
     &    p6) + za(p1,p7)*za(p2,p3)*za(p4,p5)*zb(p1,p4)*zb(p1,p8)*zb(p4
     &    ,p6) - za(p1,p7)*za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zb(
     &    p4,p6) - za(p1,p7)*za(p2,p5)*za(p3,p5)*zb(p1,p6)*zb(p1,p8)*
     &    zb(p4,p5) - za(p1,p7)*za(p2,p5)*za(p3,p6)*zb(p1,p6)*zb(p1,p8)
     &    *zb(p4,p6) + za(p2,p3)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(p3,p6
     &    )*zb(p4,p8) - za(p2,p3)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(p3,
     &    p8)*zb(p4,p6) + za(p2,p3)*za(p4,p5)*za(p7,p8)*zb(p1,p8)*zb(p4
     &    ,p6)*zb(p4,p8) - za(p2,p4)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(
     &    p4,p6)*zb(p4,p8) - za(p2,p5)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*
     &    zb(p4,p5)*zb(p6,p8) - za(p2,p5)*za(p3,p6)*za(p7,p8)*zb(p1,p8)
     &    *zb(p4,p6)*zb(p6,p8) )

      d8= + s3456**(-1)*s34**(-1)*s56**(-1)*s78**(-1)*s278**(-1) * (
     &    za(p2,p3)*za(p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p2,p8)*zb(p4,p6)
     &     - za(p2,p3)*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*zb(p3,p6
     &    ) - za(p2,p3)*za(p2,p7)*za(p4,p5)*zb(p1,p4)*zb(p2,p8)*zb(p4,
     &    p6) + za(p2,p4)*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*zb(p4
     &    ,p6) + za(p2,p5)*za(p2,p7)*za(p3,p5)*zb(p1,p6)*zb(p2,p8)*zb(
     &    p4,p5) + za(p2,p5)*za(p2,p7)*za(p3,p6)*zb(p1,p6)*zb(p2,p8)*
     &    zb(p4,p6) - za(p2,p7)*za(p3,p5)*za(p3,p7)*zb(p1,p3)*zb(p4,p6)
     &    *zb(p7,p8) + za(p2,p7)*za(p3,p5)*za(p3,p7)*zb(p1,p4)*zb(p3,p6
     &    )*zb(p7,p8) - za(p2,p7)*za(p3,p5)*za(p4,p7)*zb(p1,p4)*zb(p4,
     &    p6)*zb(p7,p8) - za(p2,p7)*za(p3,p5)*za(p5,p7)*zb(p1,p6)*zb(p4
     &    ,p5)*zb(p7,p8) - za(p2,p7)*za(p3,p6)*za(p5,p7)*zb(p1,p6)*zb(
     &    p4,p6)*zb(p7,p8) + za(p2,p7)*za(p3,p7)*za(p4,p5)*zb(p1,p4)*
     &    zb(p4,p6)*zb(p7,p8) )

      TWZbbnab=d7+d8
      return
      end
