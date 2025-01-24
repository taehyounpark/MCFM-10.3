!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub17(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p2,p5,p6,p4,p3,p1,p7,v1,v2,v3,
     & Qta1,Qta2,Qtb1,Qtb2
      integer, parameter:: Qtil=8,Qtilsl=9
      real(dp):: p(mxpart,4),pin(mxpart,4),a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,zaQb,zaPQa,zbPQb,zaPQPb,
     & zaQPa,zbQPb,res,bub17
      real(dp):: sa,sb,twopbDa,t
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)
c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)+sa/sb*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2))

      zaQPa(v1,v2)=-zaPQa(v2,v1)

c--- [v1|PQ|v2]
      zbPQb(v1,v2)=(
     & zb(v1,Qta1)*zaQb(Qta1,v2)+zb(v1,Qta2)*zaQb(Qta2,v2))

      zbQPb(v1,v2)=-zbPQb(v2,v1)

c--- <v1|PQP|v2]
      zaPQPb(v1,v2)=(
     & -za(v1,Qta1)*zbPQb(v2,Qta1)-za(v1,Qta2)*zbPQb(v2,Qta2)
     & )

      t(v1,v2,v3)=s(v1,v2)+s(v1,v3)+s(v2,v3)

c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P and Q
c---  P=p(Qta1)+p(Qta2)+p(Qta3), Q=P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
      Qta1=p1
      Qta2=p7

      Qtb1=p4
      Qtb2=p3
      sb=s(Qtb1,Qtb2)
      sa=s(Qta1,Qta2)
      twopbDa=s(Qtb1,Qta1)+s(Qtb2,Qta1)+s(Qtb1,Qta2)+s(Qtb2,Qta2)
      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:))*lamm+(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamp+(p(Qtb1,:)+p(Qtb2,:))*sa/sb

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      res=
     &-(sa*za(p2,p6)*za(p2,p1)**2*zaPb(p2,p5)*zaPb(p2,p4)
     &*zaPb(p3,p2)**2)/
     &(2._dp*sb*za(p2,p3)*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)**2
     &*zaPQa(p2,p2)) -
     &(sa*za(p6,p3)*za(p3,p1)**2*zaPb(p3,p5)*zaPb(p3,p4))/
     &(2._dp*sb*za(p2,p3)*za(p3,p7)*za(p1,p7)*zaPQa(p3,p3)) -
     &(sa*za(p2,p1)**3*zaPb(p2,p4)*zaPb(p3,p2)*zaPQa(p2,p6)*
     &  zaPQPb(p2,p5))/
     &(sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPQa(p2,p2)**2
     &*zaPQa(p2,p1)) -
     &(sa*za(p2,p1)**2*za(p6,p1)*zaPb(p2,p4)*zaPb(p3,p2)*zb(p5,p2))/
     &(sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPQa(p2,p1)) +
     &(sa*za(p2,p1)**3*zaPb(p2,p4)*zaPb(p3,p2)*zaPQa(p2,p6)*zb(p5,p2))/
     &(sb*za(p2,p7)*za(p1,p7)*zaPb(p2,p2)*zaPQa(p2,p2)*zaPQa(p2,p1)) -
     &(sa*za(p6,p3)*za(p3,p1)**2*zaPb(p3,p4)*zb(p5,p2))/
     &(2._dp*sb*za(p3,p7)*za(p1,p7)*zaPQa(p3,p3)) -
     &(sa*za(p6,p3)*za(p3,p1)*zb(p3,p4)*zb(p1,p5))/
     &(sb*za(p2,p7)*zaQb(p7,p1)) +
     &(sa*za(p6,p7)*za(p3,p1)**2*zaPb(p7,p5)*zb(p1,p4))/
     &(2._dp*sb*za(p2,p7)*za(p3,p7)*za(p1,p7)*zaQb(p7,p1)) +
     &(sa*za(p2,p1)*za(p6,p1)*za(p3,p1)*zb(p5,p2)*zb(p1,p4))/
     &(sb*za(p2,p7)*za(p1,p7)*zaQb(p1,p1)) -
     &(sa*za(p3,p1)*za(p1,p7)*zaQb(p2,p1)*zaQb(p6,p1)*zb(p5,p2)*
     &zb(p1,p4))/(sb*za(p2,p7)*zaQb(p1,p1)*zaQb(p7,p1)**2) +
     &(sa*za(p6,p7)*za(p3,p1)**2*zb(p5,p2)*zb(p1,p4))/
     &(2._dp*sb*za(p3,p7)*za(p1,p7)*zaQb(p7,p1)) -
     &(sa*za(p2,p1)*za(p6,p1)*zaPb(p2,p4)*zaPb(p3,p2)*zb(p7,p5))/
     &(sb*za(p2,p7)*zaPb(p2,p2)*zaPQa(p2,p1)) +
     &(sa*za(p6,p1)*za(p3,p1)*zb(p1,p4)*zb(p7,p5))/
     &(sb*za(p2,p7)*zaQb(p1,p1))
      res=res -
     &(sa**2*za(p2,p6)*za(p3,p1)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p2))/
     &(sb*t(p2,p5,p6)*za(p1,p7)**2*zaQb(p1,p1)*zb(p7,p1)) +
     &(sa**2*za(p2,p6)*zaQb(p3,p1)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p2))/
     &(sb*t(p2,p5,p6)*zaQb(p1,p1)*zaQb(p7,p1)**2*zb(p7,p1)) +
     &(sa**2*za(p2,p6)*za(p3,p7)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p2)**2)/
     &(2._dp*sb*t(p2,p5,p6)*za(p1,p7)**2
     & *zaQb(p7,p1)*zb(p1,p2)*zb(p7,p1))-
     &(sa**2*za(p5,p6)*za(p3,p1)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p5))/
     &(sb*t(p2,p5,p6)*za(p1,p7)**2*zaQb(p1,p1)*zb(p7,p1)) +
     &(sa**2*za(p5,p6)*zaQb(p3,p1)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p5))/
     &(sb*t(p2,p5,p6)*zaQb(p1,p1)*zaQb(p7,p1)**2*zb(p7,p1)) +
     &(sa**2*za(p5,p6)*za(p3,p7)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p5)**2)/
     &(2._dp*sb*t(p2,p5,p6)*za(p1,p7)**2
     &*zaQb(p7,p1)*zb(p1,p5)*zb(p7,p1))-
     &(za(p2,p1)**2*za(p6,p3)*zaPb(p2,p5)*zaPb(p3,p2)
     &*zb(p3,p4)*zb(p7,p1))/
     &(sb*za(p2,p7)*zaPb(p2,p2)*zaPQa(p2,p2))
      res=res +
     &(za(p6,p3)*za(p1,Qtil)**3*zaPb(p3,Qtil)*zaPQPb(Qtil,p5)*zb(p3,p4)*
     &zb(p7,p1))/
     &(sb*za(p2,Qtil)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*zaQPa(p1,Qtil)) +
     &(za(p6,p3)*za(p3,p1)*zb(p5,p2)*zb(p3,p4)*zb(p7,p1))/
     &(sb*zaQb(p7,p1)) + (za(p6,p3)*za(p1,Qtil)**3*zaPb(p3,Qtil)*
     &zaQPa(p2,Qtil)*zb(p5,p2)*zb(p3,p4)*zb(p7,p1))/
     &(sb*za(p2,Qtil)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*zaQPa(p1,Qtil)) +
     &(sa**2*za(p2,p6)*zaPb(p3,p2)**2*zb(p5,p2)*zb(p3,p4)*zb(p7,p2)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)*zaPb(Qtil,p2)*zb(p1,p2)*
     &zb(p7,p1)*zb(Qtil,p2)) + (sa**2*za(p5,p6)*zaPb(p3,p5)**2*
     &zb(p5,p2)*zb(p3,p4)*zb(p7,p5)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)*zaPb(Qtil,p5)*zb(p1,p5)*
     &zb(p7,p1)*zb(Qtil,p5)) - (sa**2*za(p2,p6)*zaPb(Qtil,p2)*
     &zaQb(p3,Qtil)**2*zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**3*zaQb(p1,Qtil)*zb(Qtil,p1)**2)
      res=res
     &-(sa**2*za(p5,p6)*zaPb(Qtil,p5)*zaQb(p3,Qtil)**2*zb(p5,p2)*
     &zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**3
     &*zaQb(p1,Qtil)*zb(Qtil,p1)**2) +
     &(sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p2,Qtil)
     &*zaQb(p6,Qtil)*zb(p5,p2)*
     &  zb(p7,p1)*zb(Qtil,p4)*zb(Qtil,p7)**2)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**3*zaQb(p1,Qtil)
     &*zb(Qtil,p1)**2)
     &- (sa**2*za(p2,p6)*za(p1,Qtil)*zaPb(Qtil,p2)*zaQb(p3,Qtil)**2*
     &zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1)) -
     &(sa**2*za(p5,p6)*za(p1,Qtil)*zaPb(Qtil,p5)*zaQb(p3,Qtil)**2*
     &zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1)) -
     &(sa**2*za(p2,p6)*zaPb(Qtil,p2)*zaQb(p3,Qtil)**2*zaQb(p1,p7)*
     &zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**3*zaQb(p1,Qtil)**2*zb(p7,p1)*
     &zb(Qtil,p1)) - (sa**2*za(p5,p6)
     &*zaPb(Qtil,p5)*zaQb(p3,Qtil)**2*
     &zaQb(p1,p7)*zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**3*zaQb(p1,Qtil)**2*zb(p7,p1)*
     &zb(Qtil,p1)) + (sa**2*za(p2,p6)*zaPb(p3,Qtil)**2*
     &zaPb(Qtil,p2)**2*zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)*zaPb(Qtil,Qtil)**3*zb(p7,p1)*
     &zb(Qtil,p2)*zb(Qtil,p1)) +
     &(sa**2*za(p5,p6)*zaPb(p3,Qtil)**2*zaPb(Qtil,p5)**2*zb(p5,p2)*
     &zb(p3,p4)*zb(Qtil,p7)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)*zaPb(Qtil,Qtil)**3*zb(p7,p1)*
     &zb(Qtil,p5)*zb(Qtil,p1)) +
     &(sa**2*za(p3,Qtil)**2*zaPb(p6,Qtil)*zb(p5,p2)*zb(Qtil,p4)*
     &zb(Qtil,p7)**2)/
     &(2._dp*sb*zaPb(p3,Qtil)*zaPb(Qtil,Qtil)**3*zb(Qtil,p1)) +
     &(sa**2*za(p3,Qtil)*za(p1,p7)*za(p1,Qtil)*zaQb(p2,Qtil)*
     &zaQb(p6,Qtil)*zb(p5,p2)*zb(p7,p1)*zb(Qtil,p4)*zb(Qtil,p7)**2)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1)) -
     &(sa**2*za(p2,p1)*za(p3,Qtil)*za(p1,p7)
     &*zaQb(p2,Qtil)*zaQb(p6,Qtil)*
     &zb(p5,p2)*zb(p7,p1)*zb(Qtil,p4)*zb(Qtil,p7)**2)/
     &(sb*zaPb(p2,Qtil)**2*zaPb(Qtil,Qtil)**3*zaQb(p1,Qtil)
     &*zb(Qtil,p1))
     &- (sa**3*za(p3,Qtil)**2*zaPb(p6,Qtil)*zb(Qtil,p5)*zb(Qtil,p4)*
     &  zb(Qtil,p7)**2)/
     &(2._dp*sb*zaPb(p2,Qtil)*zaPb(p3,Qtil)
     &*zaPb(Qtil,Qtil)**3*zb(Qtil,p1))
     &+ (sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p2,Qtil)*zaQb(p6,Qtil)*
     &  zaQb(p1,p4)*zb(p5,p2)*zb(Qtil,p7)**3)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**3
     &*zaQb(p1,Qtil)**2*zb(Qtil,p1))
     &- (sa*za(p6,p3)*za(p3,Qtil)*za(p1,p7)*zaQb(p2,Qtil)*zb(p5,p2)*
     &  zb(p3,p4)*zb(Qtil,p7)**3)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**2
     &*zaQb(p1,Qtil)*zb(Qtil,p1)) +
     &(2*sa**2*za(p2,p6)*zaPb(Qtil,p2)*zaQb(p3,Qtil)*zaQPa(p3,Qtil)*
     &zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**3)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(p7,p1)*
     &zb(Qtil,p1)) + (2*sa**2*za(p5,p6)*zaPb(Qtil,p5)*zaQb(p3,Qtil)*
     &zaQPa(p3,Qtil)*zb(p5,p2)*zb(p3,p4)*zb(Qtil,p7)**3)/
     &(sb*t(p2,p5,p6)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(p7,p1)*
     &zb(Qtil,p1)) - (sa**2*za(p3,Qtil)*za(p1,p7)*zaPQPb(Qtil,p5)*
     &zaQb(p6,Qtil)*zb(Qtil,p4)*zb(Qtil,p7)**3)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1))-
     &(sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p6,Qtil)*zaQPa(p2,Qtil)*
     &zb(p5,p2)*zb(Qtil,p4)*zb(Qtil,p7)**3)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1))-
     &(sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p2,Qtil)*zaQPa(p6,Qtil)*
     &zb(p5,p2)*zb(Qtil,p4)*zb(Qtil,p7)**3)/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1))-
     &(sa**2*za(p2,p6)*zaPb(Qtilsl,p2)*zaQb(p3,Qtilsl)**2*zb(p5,p2)*
     &zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)**2) - (sa**2*za(p5,p6)*zaPb(Qtilsl,p5)*
     &zaQb(p3,Qtilsl)**2*zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)**2) + (sa**2*za(p3,Qtilsl)*za(p1,p7)*
     &zaQb(p2,Qtilsl)*zaQb(p6,Qtilsl)*zb(p5,p2)
     &*zb(p7,p1)*zb(Qtilsl,p4)*
     &zb(Qtilsl,p7)**2)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)**2) - (sa**2*za(p2,p6)*za(p1,Qtilsl)*
     &zaPb(Qtilsl,p2)*zaQb(p3,Qtilsl)**2*zb(p5,p2)*zb(p3,p4)*
     &zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     & zb(Qtilsl,p1)) - (sa**2*za(p5,p6)
     &*za(p1,Qtilsl)*zaPb(Qtilsl,p5)*
     & zaQb(p3,Qtilsl)**2*zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     & zb(Qtilsl,p1)) - (sa**2*za(p2,p6)*zaPb(Qtilsl,p2)*
     & zaQb(p3,Qtilsl)**2*zaQb(p1,p7)
     &*zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)**2*
     & zb(p7,p1)*zb(Qtilsl,p1)) -
     &(sa**2*za(p5,p6)*zaPb(Qtilsl,p5)*zaQb(p3,Qtilsl)**2*zaQb(p1,p7)*
     & zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)**2*
     & zb(p7,p1)*zb(Qtilsl,p1)) +
     &(sa**2*za(p2,p6)*zaPb(p3,Qtilsl)**2
     &*zaPb(Qtilsl,p2)**2*zb(p5,p2)*
     & zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)
     &*zaPb(Qtilsl,Qtilsl)**3*zb(p7,p1)*
     &zb(Qtilsl,p2)*zb(Qtilsl,p1)) +
     &(sa**2*za(p5,p6)*zaPb(p3,Qtilsl)**2*zaPb(Qtilsl,p5)**2*zb(p5,p2)*
     &zb(p3,p4)*zb(Qtilsl,p7)**2)/
     &(2._dp*sa*sb*t(p2,p5,p6)*za(p1,p7)
     &*zaPb(Qtilsl,Qtilsl)**3*zb(p7,p1)*
     &zb(Qtilsl,p5)*zb(Qtilsl,p1))
      res=res +
     &(sa**2*za(p3,Qtilsl)**2*zaPb(p6,Qtilsl)*zb(p5,p2)*zb(Qtilsl,p4)*
     &zb(Qtilsl,p7)**2)/
     &(2._dp*sb*zaPb(p3,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*zb(Qtilsl,p1)) +
     &(sa**2*za(p3,Qtilsl)*za(p1,p7)*za(p1,Qtilsl)*zaQb(p2,Qtilsl)*
     &zaQb(p6,Qtilsl)*zb(p5,p2)
     &*zb(p7,p1)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**2)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) - (sa**2*za(p2,p1)*za(p3,Qtilsl)*za(p1,p7)*
     &zaQb(p2,Qtilsl)*zaQb(p6,Qtilsl)*zb(p5,p2)
     &*zb(p7,p1)*zb(Qtilsl,p4)*
     &zb(Qtilsl,p7)**2)/
     &(sb*zaPb(p2,Qtilsl)**2*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) - (sa**3*za(p3,Qtilsl)**2*zaPb(p6,Qtilsl)*
     &zb(Qtilsl,p5)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**2)/
     &(2._dp*sb*zaPb(p2,Qtilsl)*zaPb(p3,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*
     &zb(Qtilsl,p1)) + (sa**2*za(p3,Qtilsl)*za(p1,p7)*zaQb(p2,Qtilsl)*
     &zaQb(p6,Qtilsl)*zaQb(p1,p4)*zb(p5,p2)*zb(Qtilsl,p7)**3)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)**2*
     &zb(Qtilsl,p1)) + (2*sa**2*za(p2,p6)*zaPb(Qtilsl,p2)*
     &zaQb(p3,Qtilsl)*zaQPa(p3,Qtilsl)
     &*zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**3)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*zb(p7,p1)*
     &zb(Qtilsl,p1)) + (2*sa**2*za(p5,p6)*zaPb(Qtilsl,p5)*
     &zaQb(p3,Qtilsl)*zaQPa(p3,Qtilsl)
     &*zb(p5,p2)*zb(p3,p4)*zb(Qtilsl,p7)**3)/
     &(sb*t(p2,p5,p6)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*zb(p7,p1)*
     &zb(Qtilsl,p1)) - (sa**2*za(p3,Qtilsl)*za(p1,p7)*
     &zaPQPb(Qtilsl,p5)*zaQb(p6,Qtilsl)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**3)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) - (sa**2*za(p3,Qtilsl)*za(p1,p7)*zaQb(p6,Qtilsl)*
     &zaQPa(p2,Qtilsl)*zb(p5,p2)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**3)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) - (sa**2*za(p3,Qtilsl)*za(p1,p7)*zaQb(p2,Qtilsl)*
     &zaQPa(p6,Qtilsl)*zb(p5,p2)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**3)/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) + (sa*za(p3,p1)*za(p1,p7)*zaQb(p6,p1)*zb(p1,p4)*
     & zbQPb(p1,p5))/(sb*za(p2,p7)*zaQb(p1,p1)*zaQb(p7,p1)**2) -
     &(sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p6,Qtil)*zb(p7,p1)*zb(Qtil,p4)*
     & zb(Qtil,p7)**2*zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**3
     &*zaQb(p1,Qtil)*zb(Qtil,p1)**2)
     &-(sa**2*za(p3,Qtil)*za(p1,p7)*za(p1,Qtil)*zaQb(p6,Qtil)*zb(p7,p1)*
     &zb(Qtil,p4)*zb(Qtil,p7)**2*zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1)) +
     &(sa**2*za(p2,p1)*za(p3,Qtil)*za(p1,p7)*zaQb(p6,Qtil)*zb(p7,p1)*
     &zb(Qtil,p4)*zb(Qtil,p7)**2*zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)**2*zaPb(Qtil,Qtil)**3
     &*zaQb(p1,Qtil)*zb(Qtil,p1))
     &- (sa**2*za(p3,Qtil)*za(p1,p7)*zaQb(p6,Qtil)*zaQb(p1,p4)*
     &  zb(Qtil,p7)**3*zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**3
     &*zaQb(p1,Qtil)**2*zb(Qtil,p1))
     &+ (sa*za(p6,p3)*za(p3,Qtil)*za(p1,p7)*zb(p3,p4)*zb(Qtil,p7)**3*
     &zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**2*zaQb(p1,Qtil)*zb(Qtil,p1))+
     &(sa**2*za(p3,Qtil)*za(p1,p7)*zaQPa(p6,Qtil)*zb(Qtil,p4)*
     &zb(Qtil,p7)**3*zbQPb(Qtil,p5))/
     &(sb*zaPb(p2,Qtil)*zaPb(Qtil,Qtil)**4*zaQb(p1,Qtil)*zb(Qtil,p1))-
     &(sa**2*za(p3,Qtilsl)*za(p1,p7)*zaQb(p6,Qtilsl)*zb(p7,p1)*
     &zb(Qtilsl,p4)*zb(Qtilsl,p7)**2*zbQPb(Qtilsl,p5))/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)**2) - (sa**2*za(p3,Qtilsl)*za(p1,p7)*
     &za(p1,Qtilsl)*zaQb(p6,Qtilsl)*zb(p7,p1)*zb(Qtilsl,p4)*
     &zb(Qtilsl,p7)**2*zbQPb(Qtilsl,p5))/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) + (sa**2*za(p2,p1)*za(p3,Qtilsl)*za(p1,p7)*
     &zaQb(p6,Qtilsl)*zb(p7,p1)*zb(Qtilsl,p4)*zb(Qtilsl,p7)**2*
     &zbQPb(Qtilsl,p5))/
     &(sb*zaPb(p2,Qtilsl)**2*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1)) - (sa**2*za(p3,Qtilsl)
     &*za(p1,p7)*zaQb(p6,Qtilsl)*
     &zaQb(p1,p4)*zb(Qtilsl,p7)**3*zbQPb(Qtilsl,p5))/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*zaQb(p1,Qtilsl)**2*
     &zb(Qtilsl,p1)) + (sa**2*za(p3,Qtilsl)
     &*za(p1,p7)*zaQPa(p6,Qtilsl)*
     &zb(Qtilsl,p4)*zb(Qtilsl,p7)**3*zbQPb(Qtilsl,p5))/
     &(sb*zaPb(p2,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*zaQb(p1,Qtilsl)*
     &zb(Qtilsl,p1))

c--- multiply by overall factor (including W propagators)
      res=res/(s(p5,p6)*s(p4,p3))

c--- overall sign by hand
      bub17=-res

      return
      end

