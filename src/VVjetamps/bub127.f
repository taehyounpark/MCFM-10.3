!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub127(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p2,p5,p6,p4,p3,p1,p7,Qt,v1,v2,v3
      real(dp):: p(mxpart,4),pin(mxpart,4),s56,s34,s167,p53DP
      real(dp)::a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,zaQb,zaPQa,zbPQb,zaPQPb,res,bub127
      parameter(Qt=8)
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=-(zab(v1,p2,v2)+zab(v1,p1,v2)+zab(v1,p7,v2))
c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)-s167/s56*(zab(v1,p5,v2)+zab(v1,p6,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=-(
     & za(v1,p2)*zaQb(v2,p2)+za(v1,p1)*zaQb(v2,p1)+za(v1,p7)*zaQb(v2,p7)
     & )

c--- [v1|PQ|v2]
      zbPQb(v1,v2)=-(
     & zb(v1,p2)*zaQb(p2,v2)+zb(v1,p1)*zaQb(p1,v2)+zb(v1,p7)*zaQb(p7,v2)
     & )

c--- <v1|PQP|v2]
      zaPQPb(v1,v2)=-(
     & -za(v1,p2)*zbPQb(v2,p2)-za(v1,p1)*zbPQb(v2,p1)
     & -za(v1,p7)*zbPQb(v2,p7)
     & )

c--- copy momentum array
      p(:,:)=pin(:,:)

      s56=s(p5,p6)
      s34=s(p4,p3)
      s167=s(p2,p1)+s(p2,p7)+s(p1,p7)
      p53DP=0.5_dp
     & *(s(p5,p2)+s(p5,p1)+s(p5,p7)+s(p6,p2)+s(p6,p1)+s(p6,p7))

      a=s56
      b=2*p53DP
      c=s167
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qt,:)=-lamm*(p(p2,:)+p(p1,:)+p(p7,:))-(p(p5,:)+p(p6,:))*s167/s56
c      Qsq=s167+s167**2/s56-2d0*p53Dp*s167/s56

      if (swapz .eqv. .false.) then
        call spinoru(8,p,za,zb)
      else
        call spinoru(8,p,zb,za)
      endif

      Res=s167**3/s56*(
     & zb(Qt,p4)**2*za(p6,Qt)*za(p3,Qt)*za(p1,Qt)*zaPb(Qt,p5)
     & /2d0/(zaPb(Qt,Qt)**3*za(p2,Qt)*zaPb(Qt,p4)))
     &   -s167**3/s56*zb(Qt,p4)*za(p1,Qt)
     & /(zaPb(Qt,Qt)**3*za(p2,Qt)*zaPQa(p1,Qt))
     & *(zaQb(Qt,p5)*zaPQa(p3,Qt)*(
     &   za(p2,p6)*za(p1,Qt)/za(p2,Qt)+za(p6,Qt)*zaPQa(p1,p1)
     &   /zaPQa(p1,Qt)+za(p6,Qt)*zaPb(p1,Qt)/zaPb(Qt,Qt))
     &  +za(p6,Qt)*za(p1,Qt)/zaPb(Qt,Qt)*(
     &   zaQb(Qt,p5)*zaPQPb(p3,Qt)-zbPQb(Qt,p5)*zaPQa(p3,Qt)))
     &   -s167/s56*(
     & zb(Qt,p5)*zaPb(p6,Qt)*zaPb(Qt,p4)**2*zaPb(p3,Qt)*zaPb(p1,Qt)
     & /2d0/(zb(Qt,p4)*zaPb(Qt,Qt)**3*zaPb(p2,Qt)))
     &   +s167**2/s56*zaPb(Qt,p4)*zaPb(p1,Qt)
     & /(zaPb(Qt,Qt)**3*zaPb(p2,Qt)*zaPQPb(p1,Qt))
     & *(zbPQb(Qt,p5)*zaPQPb(p3,Qt)*(za(p2,p6)*zaPb(p1,Qt)/zaPb(p2,Qt)
     &   -za(p1,Qt)*zaPb(p6,Qt)/zaPb(Qt,Qt)
     &   +zaPb(p6,Qt)*zaPQa(p1,p1)/zaPQPb(p1,Qt))
     &  +zaPb(p6,Qt)*zaPb(p1,Qt)/zaPb(Qt,Qt)*(
     &   -za(p3,Qt)*zbPQb(Qt,p5)*zaPQa(p1,Qt)/za(p1,Qt)
     &   +zaQb(Qt,p5)*zaPQPb(p3,Qt)))
     &   -s167**2/s56*zb(p4,p2)*za(p2,p6)/zaPb(p2,p2)
     & *(s167*zb(p4,p2)*za(p2,p3)*za(p2,p1)*zaPb(p2,p5)
     &   /2d0/(zaPb(p2,p2)*zaPb(p2,p4)*zaPQa(p2,p2))
     &  +s167*za(p2,p1)**2*zaQb(p2,p5)*zaPQa(p3,p2)
     &   /(zaPQa(p2,p2)**2*zaPQa(p1,p2))
     &  -za(p3,p1)*zaPb(p1,p5)/zaPQa(p1,p2))
     &   +s167/s56*zb(p4,p5)*zaPb(p6,p4)*zaPb(p3,p4)*zaPb(p1,p4)
     & /2d0/(zbPQb(p4,p4)*zaPb(p2,p4))
     &   -1d0/s56
     & *zb(p4,p5)*za(p5,p6)*s167**2*zb(Qt,p5)*za(p1,Qt)*za(p3,Qt)
     & /(zaPb(Qt,Qt)**2*za(p2,Qt))
     &   +1d0/s56
     & *zb(p4,p5)*za(p5,p6)*s167*zaPb(Qt,p5)*zaPb(p1,Qt)*zaPb(p3,Qt)
     & /(zaPb(Qt,Qt)**2*zaPb(p2,Qt))
     &   -1d0/s56
     & *zb(p4,p5)*za(p5,p6)*s167**2*zb(p5,p2)*za(p3,p2)*za(p2,p1)
     & /(zaPQa(p2,p2)*zaPb(p2,p2))

c--- multiply by overall factor (including W propagators)
      res=res/(s56*s34)
     &   *za(p2,p1)/(za(p2,p7)*za(p1,p7))

c--- overall sign by hand
      bub127=-res

      return
      end

