!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub127sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p2,p5,p6,p4,p3,p1,p7,p8,p9,v1,v2,v3
      real(dp):: p(mxpart,4),pin(mxpart,4),s34,s127,twopbDa,
     & a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,zaQb,zaPQa,zbPQb,zaPQPb,
     & zaQPa,bub127sl
      integer Qta1,Qta2,Qta3,Qtb1,Qtb2
      parameter(p8=8,p9=9)
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)+zab(v1,Qta3,v2)
c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)+s127/s34*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2)
     &+za(v1,Qta3)*zaQb(v2,Qta3))

      zaQPa(v1,v2)=-zaPQa(v2,v1)


c--- [v1|PQ|v2]
      zbPQb(v1,v2)=(
     & zb(v1,Qta1)*zaQb(Qta1,v2)+zb(v1,Qta2)*zaQb(Qta2,v2)
     & +zb(v1,Qta3)*zaQb(Qta3,v2))

c--- <v1|PQP|v2]
      zaPQPb(v1,v2)=(
     & -za(v1,Qta1)*zbPQb(v2,Qta1)-za(v1,Qta2)*zbPQb(v2,Qta2)
     & -za(v1,Qta3)*zbPQb(v2,Qta3)
     & )


c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P and Q
c---  P=p(Qta1)+p(Qta2)+p(Qta3), Q=P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
      Qta1=p2
      Qta2=p1
      Qta3=p7

      Qtb1=p4
      Qtb2=p3

      s34=s(p4,p3)
      s127=s(p2,p1)+s(p2,p7)+s(p1,p7)
      twopbDa=s(p4,p2)+s(p4,p1)+s(p4,p7)+s(p3,p2)+s(p3,p1)+s(p3,p7)
      a=s34
      b=twopbDa
      c=s127
      call solvequadratic(a,b,c,lamp,lamm)
      p(p8,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamp
     &          +(p(Qtb1,:)+p(Qtb2,:))*s127/s34
      p(p9,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamm
     &            +(p(Qtb1,:)+p(Qtb2,:))*s127/s34

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      bub127sl=-(
     &-((za(p5,p6)*((s127*za(p2,p3)*za(p2,p1)*zaPb(p2,p5)*zaPb(p1,p2))/
     &      (s34*za(p2,p7)*zaPb(p2,p2)*zaPQa(p2,p2)) -
     &     (s127*za(p3,p8)*zaPb(p1,p8)*zaPb(p8,p5)*
     &        zaPQa(p1,p8))/
     &      (s34*za(p7,p8)*zaPb(p8,p8)**2*zaPQa(p2,p8)) -
     &     (s127*za(p3,p9)*zaPb(p1,p9)*zaPb(p9,p5)*
     &        zaPQa(p1,p9))/
     &      (s34*za(p7,p9)*zaPb(p9,p9)**2*
     &        zaPQa(p2,p9)) -
     &     (s127*za(p3,p7)*za(p1,p7)*zaPb(p1,p7)*zaPb(p7,p5))/
     &      (s34*za(p2,p7)*zaPb(p7,p7)*zaPQa(p7,p7)))*zb(p4,p5))/
     & za(p1,p7)) + (za(p5,p6)*(-((s127*za(p3,p7)*zaPb(p1,p7)*
     &        zaPb(p7,p5))/(s34*zaPb(p7,p7)*zaPQa(p7,p7))) -
     &   (s127*za(p3,p8)*zaPb(p1,p8)*zaPb(p8,p5)*
     &      zaPQa(p1,p8))/
     &    (s34*za(p1,p8)*zaPb(p8,p8)**2*zaPQa(p7,p8)) -
     &   (s127*za(p3,p9)*zaPb(p1,p9)*zaPb(p9,p5)*
     &      zaPQa(p1,p9))/
     &    (s34*za(p1,p9)*zaPb(p9,p9)**2*zaPQa(p7,p9)))*
     & zb(p4,p5))/za(p2,p7) - (-(s127*za(p2,p6)*za(p2,p3)*
     &     zaPb(p2,p5)*zaPb(p2,p4)*zaPb(p1,p2)**2)/
     &  (2.*s34*za(p2,p7)*za(p2,p8)*zaPb(p2,p2)**2*zaPb(p2,p8)) -
     & (s127*za(p6,p8)*za(p3,p8)*zaPb(p1,p8)**2*zaPb(p8,p5)*
     &    zaPb(p8,p4))/
     &  (2.*s34*za(p2,p8)*za(p7,p8)*zaPb(p8,p8)**3) -
     & (s127*za(p6,p9)*za(p3,p9)*zaPb(p1,p9)**2*
     &    zaPb(p9,p5)*zaPb(p9,p4))/
     &  (2.*s34*za(p2,p9)*za(p7,p9)*zaPb(p9,p9)**3) +
     & (s127*za(p2,p6)*za(p2,p3)*za(p2,p1)*zaPb(p2,p4)*zaPb(p1,p2)*
     &    zaPb(p7,p5))/(s34*za(p2,p7)**2*zaPb(p2,p2)*zaPQa(p2,p2)) -
     & (s127*za(p2,p6)*za(p2,p1)*zaPb(p2,p5)*zaPb(p2,p4)*
     &    zaPb(p1,p2)*zaPQa(p2,p3))/
     &  (s34*za(p2,p7)*zaPb(p2,p2)*zaPQa(p2,p2)**2) +
     & (s127*za(p2,p6)*za(p2,p1)*za(p3,p7)*zaPb(p1,p7)*zaPb(p7,p5)*
     &    zaPb(p7,p4))/(s34*za(p2,p7)**2*zaPb(p7,p7)*zaPQa(p2,p7)) -
     & (s127*zaPb(p1,p8)*zaPb(p8,p5)*zaPb(p8,p4)*zaPQa(p2,p3)*
     &    zaPQa(p6,p8)*zaPQa(p1,p8))/
     &  (s34*za(p7,p8)*zaPb(p8,p8)**3*zaPQa(p2,p8)**2) +
     & (s127*za(p3,p8)*zaPb(p1,p8)*zaPb(p7,p5)*zaPb(p8,p4)*
     &    zaPQa(p6,p8)*zaPQa(p1,p8))/
     &  (s34*za(p7,p8)**2*zaPb(p8,p8)**3*zaPQa(p2,p8)) -
     & (s127*zaPb(p1,p9)*zaPb(p9,p5)*zaPb(p9,p4)*
     &    zaPQa(p2,p3)*zaPQa(p6,p9)*zaPQa(p1,p9))/
     &  (s34*za(p7,p9)*zaPb(p9,p9)**3*zaPQa(p2,p9)**2)
     &+ (s127*za(p3,p9)*zaPb(p1,p9)*zaPb(p7,p5)*zaPb(p9,p4)*
     &     zaPQa(p6,p9)*zaPQa(p1,p9))/
     &   (s34*za(p7,p9)**2*zaPb(p9,p9)**3*zaPQa(p2,p9))
     &- (s127*za(p3,p7)*zaPb(p1,p7)*zaPb(p7,p5)*zaPb(p7,p4)*
     &     zaPQa(p6,p7)*zaPQa(p1,p7))/
     &   (s34*zaPb(p7,p7)*zaPQa(p2,p7)*zaPQa(p7,p7)**2) -
     &  (s127*za(p6,p7)*za(p3,p7)*zaPb(p1,p7)**2*zaPb(p7,p5)*
     &     zaPb(p7,p4))/(2.*s34*za(p2,p7)*zaPb(p7,p7)**2*zaPQa(p7,p7)) -
     &  (s127*za(p3,p8)*zaPb(p1,p8)*zaPb(p8,p5)*zaPb(p8,p4)*
     &     zaPQa(p1,p8)*zaPQPb(p6,p8))/
     &   (s34*za(p7,p8)*zaPb(p8,p8)**4*zaPQa(p2,p8)) -
     &  (s127*za(p3,p9)*zaPb(p1,p9)*zaPb(p9,p5)*
     &     zaPb(p9,p4)*zaPQa(p1,p9)*zaPQPb(p6,p9))/
     &   (s34*za(p7,p9)*zaPb(p9,p9)**4*zaPQa(p2,p9)) -
     &  (s127*za(p3,p8)*zaPb(p1,p8)*zaPb(p8,p5)*zaPb(p8,p4)*
     &     zaPQa(p6,p8)*zaPQPb(p1,p8))/
     &   (s34*za(p7,p8)*zaPb(p8,p8)**4*zaPQa(p2,p8)) -
     &  (s127*za(p3,p9)*zaPb(p1,p9)*zaPb(p9,p5)*
     &     zaPb(p9,p4)*zaPQa(p6,p9)*zaPQPb(p1,p9))/
     &   (s34*za(p7,p9)*zaPb(p9,p9)**4*zaPQa(p2,p9)) +
     &  (s127**2*za(p2,p6)*za(p2,p3)*za(p2,p1)*zaPb(p2,p5)*
     &     zaPb(p1,p2)*zb(p4,p2))/
     &   (s34*za(p2,p7)*zaPb(p2,p2)**2*zaPQa(p2,p2)) +
     &  (s127**2*za(p3,p8)*zaPb(p1,p8)*zaPb(p8,p5)*
     &     zaPQa(p6,p8)*zaPQa(p1,p8)*zb(p8,p4))/
     &   (s34*za(p7,p8)*zaPb(p8,p8)**4*zaPQa(p2,p8)) +
     &  (s127**2*za(p3,p9)*zaPb(p1,p9)*zaPb(p9,p5)*
     &     zaPQa(p6,p9)*zaPQa(p1,p9)*zb(p9,p4))/
     &   (s34*za(p7,p9)*zaPb(p9,p9)**4*zaPQa(p2,p9)))/
     &za(p1,p7) + (-(s127*za(p6,p3)*za(p6,p1)*zaPb(p6,p5)*zaPb(p6,p4))/
     &   (2.*s34*za(p6,p7)*za(p6,p8)*zaPb(p6,p8)) +
     &  (s127*za(p3,p7)*za(p1,p7)*zaPb(p6,p7)**2*zaPb(p7,p5)*
     &     zaPb(p7,p4))/
     &   (2.*s34*za(p6,p7)*za(p7,p8)*zaPb(p7,p7)**2*zaPb(p7,p8)) -
     &  (s127*za(p3,p8)*za(p1,p8)*zaPb(p6,p8)**2*zaPb(p8,p5)*
     &     zaPb(p8,p4))/
     &   (2.*s34*za(p6,p8)*za(p7,p8)*zaPb(p8,p8)**3) -
     &  (s127*za(p3,p9)*za(p1,p9)*zaPb(p6,p9)**2*
     &     zaPb(p9,p5)*zaPb(p9,p4))/
     &   (2.*s34*za(p6,p9)*za(p7,p9)*zaPb(p9,p9)**3) -
     &  (s127*za(p3,p7)*zaPb(p6,p7)*zaPb(p1,p5)*zaPb(p7,p4))/
     &   (s34*zaPb(p7,p7)*zaPQa(p7,p7)) +
     &  (s127*za(p3,p8)*zaPb(p6,p8)*zaPb(p1,p5)*zaPb(p8,p4)*
     &     zaPQa(p1,p8)**2)/
     &   (s34*za(p1,p8)**2*zaPb(p8,p8)**3*zaPQa(p7,p8)) +
     &  (s127*za(p3,p9)*zaPb(p6,p9)*zaPb(p1,p5)*
     &     zaPb(p9,p4)*zaPQa(p1,p9)**2)/
     &(s34*za(p1,p9)**2*zaPb(p9,p9)**3*zaPQa(p7,p9)) -
     &  (2*s127*za(p3,p8)*zaPb(p6,p8)*zaPb(p8,p5)*zaPb(p8,p4)*
     &     zaPQa(p1,p8)*zaPQPb(p1,p8))/
     &   (s34*za(p1,p8)*zaPb(p8,p8)**4*zaPQa(p7,p8)) -
     &  (2*s127*za(p3,p9)*zaPb(p6,p9)*zaPb(p9,p5)*
     &     zaPb(p9,p4)*zaPQa(p1,p9)*zaPQPb(p1,p9))/
     &   (s34*za(p1,p9)*zaPb(p9,p9)**4*zaPQa(p7,p9)) -
     &  (s127*za(p1,p7)*zaPb(p6,p7)*zaPb(p7,p5)*zaPb(p7,p4)*
     &     zaQPa(p3,p7))/(s34*zaPb(p7,p7)*zaPQa(p7,p7)**2) +
     &  (s127*zaPb(p6,p8)*zaPb(p8,p5)*zaPb(p8,p4)*
     &     zaPQa(p1,p8)**2*zaQPa(p3,p7))/
     &   (s34*za(p1,p8)*zaPb(p8,p8)**3*zaPQa(p7,p8)**2) +
     &  (s127*zaPb(p6,p9)*zaPb(p9,p5)*zaPb(p9,p4)*
     &     zaPQa(p1,p9)**2*zaQPa(p3,p7))/
     &   (s34*za(p1,p9)*zaPb(p9,p9)**3*zaPQa(p7,p9)**2)-
     &  (s127**2*za(p3,p7)*za(p1,p7)*zaPb(p6,p7)*zaPb(p7,p5)*
     &     zb(p7,p4))/(s34*zaPb(p7,p7)**2*zaPQa(p7,p7)) +
     &  (s127**2*za(p3,p8)*zaPb(p6,p8)*zaPb(p8,p5)*
     &     zaPQa(p1,p8)**2*zb(p8,p4))/
     &   (s34*za(p1,p8)*zaPb(p8,p8)**4*zaPQa(p7,p8)) +
     &  (s127**2*za(p3,p9)*zaPb(p6,p9)*zaPb(p9,p5)*
     &     zaPQa(p1,p9)**2*zb(p9,p4))/
     &   (s34*za(p1,p9)*zaPb(p9,p9)**4*zaPQa(p7,p9)))/
     &za(p2,p7))/(s(p5,p6)*s(p4,p3))


      return
      end

