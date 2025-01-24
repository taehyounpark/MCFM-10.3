!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub234(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p2,p4,p3,p5,p6,p1,p7,Qtil,Qtilsl,v1,v2,v3
      real(dp):: p(mxpart,4),pin(mxpart,4),sa,sb,twopbDa,t,
     & a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,zaQb,zaPQa,zbPQb,zaPQPb,
     & zaQPa,zbQPb,res,bub234
      integer Qta1,Qta2,Qta3,Qtb1,Qtb2
      parameter(Qtil=8,Qtilsl=9)
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)+zab(v1,Qta3,v2)
c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)+sa/sb*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2)
     &+za(v1,Qta3)*zaQb(v2,Qta3))

      zaQPa(v1,v2)=-zaPQa(v2,v1)

c--- [v1|PQ|v2]
      zbPQb(v1,v2)=(
     & zb(v1,Qta1)*zaQb(Qta1,v2)+zb(v1,Qta2)*zaQb(Qta2,v2)
     & +zb(v1,Qta3)*zaQb(Qta3,v2))

      zbQPb(v1,v2)=-zbPQb(v2,v1)

c--- <v1|PQP|v2]
      zaPQPb(v1,v2)=(
     & -za(v1,Qta1)*zbPQb(v2,Qta1)-za(v1,Qta2)*zbPQb(v2,Qta2)
     & -za(v1,Qta3)*zbPQb(v2,Qta3)
     & )

c      s(p2,p4)=(p(p2,4)+p(p4,4))**2-(p(p2,1)+p(p4,1))**2
c     &        -(p(p2,2)+p(p4,2))**2 -(p(p2,3)+p(p4,3))**2
      t(v1,v2,v3)=s(v1,v2)+s(v1,v3)+s(v2,v3)

c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P and Q
c---  P=p(Qta1)+p(Qta2)+p(Qta3), Q=P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
      Qta1=p2
      Qta2=p4
      Qta3=p3

      Qtb1=p1
      Qtb2=p7

      sb=s(Qtb1,Qtb2)
      sa=+s(Qta1,Qta2)+s(Qta1,Qta3)+s(Qta2,Qta3)
      twopbDa=s(Qta1,Qtb1)+s(Qta2,Qtb1)+s(Qta3,Qtb1)
     &       +s(Qta1,Qtb2)+s(Qta2,Qtb2)+s(Qta3,Qtb2)

      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamm
     &           +(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamp
     &           +(p(Qtb1,:)+p(Qtb2,:))*sa/sb

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      res=
     &(za(p2,p1)*za(p6,p1)*zaPb(p2,p5)*zaPb(p3,p2)**2*zb(p3,p4))/
     & (2._dp*sa*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)**2) +
     &(za(p2,p1)*za(p6,p1)*zaPb(p2,p5)*zaPb(p3,p2)
     &*zaPb(p3,p7)*zb(p3,p4))/
     & (sa*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPb(p2,p7)) +
     &(3*za(p6,p1)*zaPb(p3,p7)**2*zaPb(p7,p5)*zb(p3,p4))/
     & (2._dp*sa*za(p2,p7)*zaPb(p7,p7)**2) +
     &(za(p3,p7)*za(p6,p1)*zaPb(p3,p7)*zaPb(p1,p7)
     &*zaPb(p7,p5)*zb(p3,p4))/
     & (sa*za(p2,p7)*za(p1,p7)*zaPb(p7,p7)**2) -
     &(za(p3,p7)*za(p6,p1)*zaPb(p2,p5)*zaPb(p3,p7)
     &*zaPb(p1,p7)*zb(p3,p4))/
     & (sa*za(p2,p7)*za(p1,p7)*zaPb(p2,p7)*zaPb(p7,p7)) +
     &(sa*za(p6,Qtil)*za(p1,Qtil)**2*zaPb(p3,Qtil)**2*zaPb(Qtil,p5)*
     &   zb(p3,p4))/
     & (2._dp*sb*za(p2,Qtil)*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**3) +
     &(sa*za(p2,p6)*za(p2,p1)**2*zaPb(p2,p5)*zaPb(p3,p2)**2*zb(p3,p4))/
     & (2._dp*sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)**2*zaPQa(p2,p2)) +
     &(sa*za(p1,Qtil)**2*zaPb(p3,Qtil)*zaPb(p6,Qtil)*zaPb(Qtil,p5)*
     &   zaPQa(p3,Qtil)*zb(p3,p4))/
     & (sb*za(p2,Qtil)*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**4) -
     &(sa*za(p6,Qtil)*za(p1,Qtil)*zaPb(p3,Qtil)*zaPb(Qtil,p5)*
     &   zaPQa(p3,Qtil)*zb(p3,p4))/
     & (sb*za(p2,Qtil)*za(p7,Qtil)**2*zaPb(Qtil,Qtil)**3) +
     &(sa*za(p2,p1)*za(p6,Qtil)*za(p1,Qtil)*zaPb(p3,Qtil)*
     &   zaPb(Qtil,p5)*zaPQa(p3,Qtil)*zb(p3,p4))/
     & (sb*za(p2,Qtil)**2*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**3) +
     &(sa*za(p6,p7)*za(p1,p7)*zaPb(p3,p7)*zaPb(p7,p5)*zaPQa(p3,p7)*
     &   zb(p3,p4))/(sb*za(p2,p7)*zaPb(p7,p7)*zaPQa(p7,p7)**2) +
     &(sa*za(p6,p7)*za(p1,p7)*zaPb(p3,p7)**2*zaPb(p7,p5)*zb(p3,p4))/
     & (2._dp*sb*za(p2,p7)*zaPb(p7,p7)**2*zaPQa(p7,p7)) +
     &(sa*za(p6,Qtil)*za(p1,Qtil)**2*zaPb(p3,Qtil)*zaPb(Qtil,p5)*
     &   zaPQPb(p3,Qtil)*zb(p3,p4))/
     & (sb*za(p2,Qtil)*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**4) +
     &(sa*za(p3,Qtil)*zaPb(p6,Qtil)*zaPb(p1,Qtil)**2*zaPb(Qtil,p5)*
     &   zaPQPb(p3,Qtil)*zb(p3,p4))/
     & (sb*za(p1,p7)*zaPb(p2,Qtil)*zaPb(p7,Qtil)*zaPb(Qtil,Qtil)**4) -
     &(sa*za(p2,p6)*za(p2,p1)**2*zaPb(p2,p5)*zaPb(p3,p2)*zaQPa(p2,p3)*
     &zb(p3,p4))/(sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPQa(p2,p2)**2) -
     &(za(p2,p6)*za(p6,p1)*zaPb(p2,p7)*zaPb(p3,p2)**2
     &*zb(p3,p4)*zb(p6,p5))/
     & (2._dp*sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p2,p2)**2) -
     &(za(p2,p6)*za(p6,p1)*zaPb(p3,p2)*zaPb(p3,p7)*zb(p3,p4)*zb(p6,p5))/
     & (sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p2,p2)) -
     &(3*za(p6,p1)*za(p6,p7)*zaPb(p3,p7)**2*zb(p3,p4)*zb(p6,p5))/
     & (2._dp*sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p7,p7)) -
     &(za(p2,p1)*za(p6,p1)*zaPb(p2,p7)*zaPb(p3,p2)**2
     &*zb(p3,p4)*zb(p1,p5))/
     & (2._dp*sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p2,p2)**2) -
     &(za(p2,p1)*za(p6,p1)*zaPb(p3,p2)*zaPb(p3,p7)*zb(p3,p4)*zb(p1,p5))/
     & (sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p2,p2)) -
     &(3*za(p6,p1)*za(p1,p7)*zaPb(p3,p7)**2*zb(p3,p4)*zb(p1,p5))/
     & (2._dp*sa*t(p5,p6,p1)*za(p2,p7)*zaPb(p7,p7)) +
     &(sa*za(p2,p3)*za(p2,p1)**2*za(p6,p1)*zaPb(p2,p5)*zaPb(p3,p2)*
     &   zb(p3,p4)*zb(p7,p1))/
     & (sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPb(p2,p7)*zaPQa(p2,p2)) +
     &(sa*za(p3,p7)*za(p6,p1)*za(p1,p7)*zaPb(p3,p7)*zaPb(p7,p5)*
     &zb(p3,p4)*zb(p7,p1))/(sb*za(p2,p7)*zaPb(p7,p7)**2*zaPQa(p7,p7)) +
     &(za(p6,p1)*za(p1,Qtil)**2*zaPb(p3,Qtil)
     &*zaPb(Qtil,p5)*zaPQa(p3,Qtil)*zb(p3,p4)*zb(p7,p1))/
     & (sb*za(p2,Qtil)*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     &   zaQb(Qtil,p7)) + (sa**2*za(p3,Qtil)**2*zaPb(p6,Qtil)*
     &   zaPb(p1,Qtil)**2*zb(p3,p4)*zb(Qtil,p5))/
     & (2._dp*sb*za(p1,p7)*zaPb(p2,Qtil)*zaPb(p7,Qtil)
     & *zaPb(Qtil,Qtil)**3)-
     &(sa**2*za(p6,Qtil)*za(p1,Qtil)**2*zaPb(p3,Qtil)*zaPQa(p3,Qtil)*
     &   zb(p3,p4)*zb(Qtil,p5))/
     & (sb*za(p2,Qtil)*za(p1,p7)*za(p7,Qtil)*zaPb(Qtil,Qtil)**4) -
     &(sa**2*za(p3,Qtil)*zaPb(p6,Qtil)*zaPb(p1,Qtil)**2*
     &   zaPQa(p3,Qtil)*zb(p3,p4)*zb(Qtil,p5))/
     & (sb*za(p1,p7)*zaPb(p2,Qtil)*zaPb(p7,Qtil)*zaPb(Qtil,Qtil)**4) -
     &(sa**2*za(p3,Qtil)*za(p6,Qtil)*zaPb(p1,Qtil)**2*zaPQPb(p3,Qtil)*
     &   zb(p3,p4)*zb(Qtil,p5))/
     & (sb*za(p1,p7)*zaPb(p2,Qtil)*zaPb(p7,Qtil)*zaPb(Qtil,Qtil)**4) -
     &(sa**2*za(p3,Qtil)*zaPb(p6,Qtil)*zaPb(p1,Qtil)*zaPQPb(p3,Qtil)*
     &   zb(p3,p4)*zb(Qtil,p5))/
     & (sb*zaPb(p2,Qtil)*zaPb(p7,Qtil)**2*zaPb(Qtil,Qtil)**3) +
     &(sa**2*za(p2,p1)*za(p3,Qtil)*zaPb(p6,Qtil)*zaPb(p1,Qtil)*
     &   zaPQPb(p3,Qtil)*zb(p3,p4)*zb(Qtil,p5))/
     & (sb*za(p1,p7)*zaPb(p2,Qtil)**2*zaPb(p7,Qtil)*zaPb(Qtil,Qtil)**3)+
     &(sa*za(p3,Qtil)*za(p6,p1)*zaPb(p1,Qtil)**2*zaPQPb(p3,Qtil)*
     &   zb(p3,p4)*zb(p7,p1)*zb(Qtil,p5))/
     & (sb*za(p1,p7)*zaPb(p2,Qtil)*zaPb(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     &   zbPQb(Qtil,p7)) + (sa*za(p3,p7)*za(p6,p1)*zaPb(p3,p7)*
     &   zaPb(p1,p7)**2*zb(p3,p4)*zb(p7,p5)*zb(p7,p1))/
     & (sb*za(p1,p7)*zaPb(p2,p7)*zaPb(p7,p7)**2*zbQPb(p7,p7))

c--- multiply by overall factor (including W propagators)
      res=res/(za(p4,p3)*zb(p3,p4)*za(p5,p6)*zb(p6,p5))

c--- overall sign by hand
      bub234=-res

      return
      end

