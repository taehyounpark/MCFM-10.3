!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub234sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p2,p4,p3,p5,p6,p1,p7,v1,v2,v3,
     & Qta1,Qta2,Qta3,Qtb1,Qtb2
      integer,parameter::Qtil=8,Qtilsl=9
      real(dp):: p(mxpart,4),pin(mxpart,4),sa,sb,twopbDa,t,
     & a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,zaQb,zaPQa,zbPQb,
     & zbQPb,res,bub234sl
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)+zab(v1,Qta3,v2)
c--- <v1|Q|v2]
      zaQb(v1,v2)=(zaPb(v1,v2)+sa/sb*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2)))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2)
     &+za(v1,Qta3)*zaQb(v2,Qta3))

c--- [v1|PQ|v2]
      zbPQb(v1,v2)=(
     & zb(v1,Qta1)*zaQb(Qta1,v2)+zb(v1,Qta2)*zaQb(Qta2,v2)
     & +zb(v1,Qta3)*zaQb(Qta3,v2))

      zbQPb(v1,v2)=-zbPQb(v2,v1)


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
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamp
     &            +(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamm
     &            +(p(Qtb1,:)+p(Qtb2,:))*sa/sb

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      res=
     &(za(p6,p1)*zb(p3,p4)*(-(zb(p6,p5)*
     &(-((za(p2,p3)**2*za(p4,p6))/(za(p2,p4)*za(p2,p7)**2)) +
     &  (za(p4,p3)*za(p4,p6)*za(p3,p7))/(za(p2,p7)*za(p4,p7)**2) +
     &  (za(p2,p3)*za(p4,p6)*za(p3,p7))/(za(p2,p7)**2*za(p4,p7)) +
     &  (za(p2,p3)**2*za(p2,p6)*zaPb(p4,p4))/
     &   (za(p2,p4)*za(p2,p7)**2*zaPb(p2,p4)) +
     &  (za(p3,p7)**2*zaPb(p4,p4)*zaPb(p6,p4))/
     &   (za(p2,p7)*za(p4,p7)*zaPb(p7,p4)**2) +
     &  (za(p4,p3)*za(p3,p7)*za(p6,p7)*zaPb(p4,p4))/
     &   (za(p2,p7)*za(p4,p7)**2*zaPb(p7,p4)) +
     &  (za(p2,p3)*za(p3,p7)*za(p6,p7)*zaPb(p4,p4))/
     &   (za(p2,p7)**2*za(p4,p7)*zaPb(p7,p4)) -
     &  (sa*za(p2,p3)**2*za(p2,p6)*zb(p4,p2))/
     &   (za(p2,p7)**2*zaPb(p2,p2)*zaPb(p2,p4)) -
     &  (sa*za(p3,p7)*za(p6,p7)*zaPb(p3,p7)*zb(p7,p4))/
     &   (za(p2,p7)*zaPb(p7,p4)*zaPb(p7,p7)**2) -
     &  (sa*za(p3,p7)**2*zaPb(p6,p4)*zb(p7,p4))/
     &   (za(p2,p7)*zaPb(p7,p4)**2*zaPb(p7,p7)) -
     &  (sa*za(p2,p3)*za(p3,p7)*za(p6,p7)*zb(p7,p4))/
     &   (za(p2,p7)**2*zaPb(p7,p4)*zaPb(p7,p7)))) -
     &zb(p1,p5)*(-((za(p2,p3)**2*za(p4,p1))/(za(p2,p4)*za(p2,p7)**2)) +
     &   (za(p4,p3)*za(p4,p1)*za(p3,p7))/(za(p2,p7)*za(p4,p7)**2) +
     &   (za(p2,p3)*za(p4,p1)*za(p3,p7))/(za(p2,p7)**2*za(p4,p7)) +
     &   (za(p2,p3)**2*za(p2,p1)*zaPb(p4,p4))/
     &    (za(p2,p4)*za(p2,p7)**2*zaPb(p2,p4)) +
     &   (za(p3,p7)**2*zaPb(p4,p4)*zaPb(p1,p4))/
     &    (za(p2,p7)*za(p4,p7)*zaPb(p7,p4)**2) +
     &   (za(p4,p3)*za(p3,p7)*za(p1,p7)*zaPb(p4,p4))/
     &    (za(p2,p7)*za(p4,p7)**2*zaPb(p7,p4)) +
     &   (za(p2,p3)*za(p3,p7)*za(p1,p7)*zaPb(p4,p4))/
     &    (za(p2,p7)**2*za(p4,p7)*zaPb(p7,p4)) -
     &   (sa*za(p2,p3)**2*za(p2,p1)*zb(p4,p2))/
     &    (za(p2,p7)**2*zaPb(p2,p2)*zaPb(p2,p4)) -
     &   (sa*za(p3,p7)*za(p1,p7)*zaPb(p3,p7)*zb(p7,p4))/
     &    (za(p2,p7)*zaPb(p7,p4)*zaPb(p7,p7)**2) -
     &   (sa*za(p3,p7)**2*zaPb(p1,p4)*zb(p7,p4))/
     &    (za(p2,p7)*zaPb(p7,p4)**2*zaPb(p7,p7)) -
     &   (sa*za(p2,p3)*za(p3,p7)*za(p1,p7)*zb(p7,p4))/
     &    (za(p2,p7)**2*zaPb(p7,p4)*zaPb(p7,p7)))))/t(p5,p6,p1) +
     &  (-((za(p2,p3)*za(p2,p6)*zaPb(p3,p2))/
     &     (za(p2,p7)*zaPb(p2,p2)*zaPb(p2,p1))) +
     &  (za(p3,p1)*zaPb(p3,p1)*zaPb(p6,p1))/
     &   (zaPb(p2,p1)*zaPb(p1,p1)*zaPb(p7,p1)) +
     &  (za(p3,p7)*za(p6,p7)*zaPb(p3,p7))
     & /(za(p2,p7)*zaPb(p7,p1)*zaPb(p7,p7)))*
     &zb(p3,p4)*zb(p7,p5) - zb(p3,p4)*
     &((za(p2,p3)**2*za(p2,p6)*zaPb(p1,p5))/
     &   (za(p2,p1)*za(p2,p7)**2*zaPb(p2,p1)) +
     &  (za(p2,p3)**2*za(p6,p1)*zaPb(p1,p5))/
     &   (za(p2,p1)*za(p2,p7)**2*zaPb(p1,p1)) +
     &  (za(p3,p1)*za(p3,p7)*za(p6,p1)*zaPb(p1,p5))/
     &   (za(p2,p7)*za(p1,p7)**2*zaPb(p1,p1)) -
     &  (za(p2,p3)*za(p3,p7)*za(p6,p1)*zaPb(p1,p5))/
     &   (za(p2,p7)**2*za(p1,p7)*zaPb(p1,p1)) +
     &  (za(p3,p7)**2*zaPb(p6,p1)*zaPb(p1,p5))/
     &   (za(p2,p7)*za(p1,p7)*zaPb(p7,p1)**2) -
     &  (za(p3,p1)*za(p3,p7)*za(p6,p7)*zaPb(p1,p5))/
     &   (za(p2,p7)*za(p1,p7)**2*zaPb(p7,p1)) +
     &  (za(p2,p3)*za(p3,p7)*za(p6,p7)*zaPb(p1,p5))/
     &   (za(p2,p7)**2*za(p1,p7)*zaPb(p7,p1)) -
     &  (sa*za(p2,p3)**2*za(p2,p6)*zb(p5,p2))/
     &   (za(p2,p7)**2*zaPb(p2,p2)*zaPb(p2,p1)) -
     &  (sa*za(p3,p7)*za(p6,p7)*zaPb(p3,p7)*zb(p7,p5))/
     &   (za(p2,p7)*zaPb(p7,p1)*zaPb(p7,p7)**2) -
     &  (sa*za(p3,p7)**2*zaPb(p6,p1)*zb(p7,p5))/
     &   (za(p2,p7)*zaPb(p7,p1)**2*zaPb(p7,p7)) -
     &  (sa*za(p2,p3)*za(p3,p7)*za(p6,p7)*zb(p7,p5))/
     &   (za(p2,p7)**2*zaPb(p7,p1)*zaPb(p7,p7))) -
     &  zb(p3,p4)*zb(p7,p1)*((sa*za(p2,p3)*za(p2,p6)*za(p2,p1)*
     &   zaPb(p2,p5)*zaPb(p3,p2))/
     & (sb*za(p2,p7)*zaPb(p2,p2)*zaPb(p2,p1)*zaPQa(p2,p2)) +
     &(sa*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p3,p7)*zaPb(p7,p5))/
     & (sb*za(p2,p7)*zaPb(p7,p1)*zaPb(p7,p7)*zaPQa(p7,p7)) +
     &(za(p6,Qtil)*za(p1,Qtil)*zaPb(p3,Qtil)
     &*zaPb(Qtil,p5)*zaPQa(p3,Qtil))/
     & (sb*za(p2,Qtil)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*zaQb(Qtil,p1)) +
     &(za(p6,Qtilsl)*za(p1,Qtilsl)*zaPb(p3,Qtilsl)*zaPb(Qtilsl,p5)*
     &   zaPQa(p3,Qtilsl))/
     & (sb*za(p2,Qtilsl)*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2*
     &   zaQb(Qtilsl,p1)) + (sa*za(p3,p1)*zaPb(p3,p1)*zaPb(p6,p1)*
     &   zb(p1,p5))/(sb*zaPb(p2,p1)*zaPb(p7,p1)*zbQPb(p1,p1)))

c--- multiply by overall factor (including W propagators)
      res=res/(s(p4,p3)*s(p5,p6))

c--- overall sign by hand
      bub234sl=-res

      return
      end

