!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub12sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p2,p5,p6,p4,p3,p1,p7,Qtil,Qtilsl,Q2til,Q2tilsl,v1,v2,v3
      real(dp):: p(mxpart,4),pin(mxpart,4)
      complex(dp):: bub12sl,zab,zaPb,zaQb,zaPQa,res,zaQ2b,zaPQ2a
      integer Qta1,Qta2,Qtb1,Qtb2,Q2tb1,Q2tb2,Q2tb3
      real(dp):: sa,sb,sb2,twopbDa,twopb2Da
      real(dp):: a,b,c,lamp,lamm

      parameter(Qtil=8,Qtilsl=9,Q2til=10,Q2tilsl=11)
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)

c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)+sa/sb*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2))

c--- <v1|Q2|v2]
      zaQ2b(v1,v2)=zaPb(v1,v2)
     & +sa/sb2*(zab(v1,Q2tb1,v2)+zab(v1,Q2tb2,v2)+zab(v1,Q2tb3,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2))

c--- <v1|PQ2|v2>
      zaPQ2a(v1,v2)=(
     & za(v1,Qta1)*zaQ2b(v2,Qta1)+za(v1,Qta2)*zaQ2b(v2,Qta2))

c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P, Q and Q2
c---  P =p(Qta1)+p(Qta2)
c---  Q =P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
c---  Q2=P^2 (p(Q2tb1)+p(Q2tb2)+p(Q2tb3)) + (p(Q2tb1)+p(Q2tb2)+p(Q2tb3))^2 P
      Qta1=p2
      Qta2=p1

      Qtb1=p4
      Qtb2=p3

      Q2tb1=p4
      Q2tb2=p3
      Q2tb3=p7

c     soln1
      sb=s(Qtb1,Qtb2)
      sa=s(Qta1,Qta2)
      twopbDa=s(Qta1,Qtb1)+s(Qta2,Qtb1)+s(Qta1,Qtb2)+s(Qta2,Qtb2)
      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:)  =(p(Qta1,:)+p(Qta2,:))*lamp+(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamm+(p(Qtb1,:)+p(Qtb2,:))*sa/sb

c     soln2
      sb2=s(Q2tb1,Q2tb2)+s(Q2tb1,Q2tb3)+s(Q2tb2,Q2tb3)
      twopb2Da=s(Qta1,Q2tb1)+s(Qta1,Q2tb2)+s(Qta1,Q2tb3)
     &        +s(Qta2,Q2tb1)+s(Qta2,Q2tb2)+s(Qta2,Q2tb3)
      a=sb2
      b=twopb2Da
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Q2til,:) = (p(Qta1,:)+p(Qta2,:))*lamm
     &             +(p(Q2tb1,:)+p(Q2tb2,:)+p(Q2tb3,:))*sa/sb2
      p(Q2tilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamp
     &             +(p(Q2tb1,:)+p(Q2tb2,:)+p(Q2tb3,:))*sa/sb2

      if (swapz .eqv. .false.) then
        call spinoru(11,p,za,zb)
      else
        call spinoru(11,p,zb,za)
      endif

      res=
     &(sa**2*za(p2,p6)*za(p2,p3)*za(p2,p1)*zb(p4,p2)*zb(p1,p5))/
     & (sb*za(p2,p7)**2*zaPQa(p2,p2)*zb(p1,p2)) +
     &(sa**2*za(p2,p1)*za(p6,p7)*za(p3,p7)**2*zaPb(p7,p5)*zb(p3,p4)*
     & zb(p7,p2))/(sb*sb2*za(p2,p7)*zaPQ2a(p7,p7)*zaPQa(p7,p7))
     &+ (sa**2*za(p6,p7)*za(p3,p7)*zaPb(p7,p4)*zb(p7,p2)*zb(p7,p5))/
     & (sb*za(p2,p7)*zaPb(p7,p7)*zaPQa(p7,p7)*zb(p1,p2)) +
     &(sa**2*za(p3,p7)*za(p1,p7)*zaPb(p6,p7)*zaPb(p7,p5)*zb(p7,p4))/
     & (sb*za(p2,p7)*zaPb(p7,p7)**2*zaPQa(p7,p7)) +
     &(sa**2*za(p2,p3)*za(p6,p7)*za(p1,p7)*zaPb(p7,p5)*zb(p7,p4))/
     & (sb*za(p2,p7)**2*zaPb(p7,p7)*zaPQa(p7,p7)) +
     &(sa**2*za(p6,p7)*za(p3,p7)*za(p1,Qtil)*zaPb(p7,p5)*zb(p7,p4))/
     & (sb*za(p2,p7)*za(p7,Qtil)*zaPb(p7,p7)*zaPQa(p7,p7)) -
     &(sa**3*za(p6,p7)*za(p3,p7)*zb(p5,p2)*zb(p7,p4))/
     & (sb*za(p2,p7)*zaPb(p7,p7)*zaPQa(p7,p7)*zb(p1,p2)) +
     &(sa**2*za(p2,p1)*za(p6,Q2til)*za(p3,Q2til)**2*zaPb(Q2til,p5)*
     &   zaPb(Q2til,p7)*zb(p3,p4)*zb(Q2til,p2))/
     & (sb*sb2*za(p2,Q2til)*za(p7,Q2til)*zaPb(Q2til,Q2til)**2*
     &   zaPQa(Q2til,Q2til)) + (sa**2*za(p2,p1)*za(p6,Q2tilsl)*
     &   za(p3,Q2tilsl)**2*zaPb(Q2tilsl,p5)*zaPb(Q2tilsl,p7)*zb(p3,p4)*
     &   zb(Q2tilsl,p2))/
     & (sb*sb2*za(p2,Q2tilsl)*za(p7,Q2tilsl)*
     &   zaPb(Q2tilsl,Q2tilsl)**2*zaPQa(Q2tilsl,Q2tilsl)) +
     &(sa**2*za(p2,p1)*za(p6,Qtil)*za(p3,Qtil)**2*zaPb(Qtil,p5)*
     &   zaPb(Qtil,p7)*zb(p3,p4)*zb(Qtil,p2))/
     & (sb*sb2*za(p2,Qtil)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     &   zaPQ2a(Qtil,Qtil)) + (sa**2*za(p6,Qtil)*za(p3,Qtil)*
     &   zaPb(Qtil,p4)*zb(p7,p5)*zb(Qtil,p2))/
     & (sb*za(p2,Qtil)*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*zb(p1,p2)) -
     &(sa**2*za(p2,p1)*za(p6,p7)*za(p3,p7)*za(p7,Qtil)*zaPb(p7,p5)*
     & zb(p7,p4)*zb(Qtil,p2))/
     & (sb*za(p2,p7)*zaPb(p7,p7)*zaPQa(p7,p7)**2) +
     &(sa**2*za(p6,Qtil)*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p5)*
     &   zb(Qtil,p4))/
     & (sb*za(p2,Qtil)*za(p7,Qtil)**2*zaPb(Qtil,Qtil)**2) +
     &(sa**2*za(p2,p1)*za(p6,Qtilsl)*za(p3,Qtilsl)**2*zaPb(Qtilsl,p5)*
     &   zaPb(Qtilsl,p7)*zb(p3,p4)*zb(Qtilsl,p2))/
     & (sb*sb2*za(p2,Qtilsl)*za(p7,Qtilsl)*
     &   zaPb(Qtilsl,Qtilsl)**2*zaPQ2a(Qtilsl,Qtilsl)) +
     &(sa**2*za(p6,Qtilsl)*za(p3,Qtilsl)*zaPb(Qtilsl,p4)*zb(p7,p5)*
     &   zb(Qtilsl,p2))/
     & (sb*za(p2,Qtilsl)*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2*
     &   zb(p1,p2)) + (sa**2*za(p6,Qtilsl)*za(p3,Qtilsl)*za(p1,Qtilsl)*
     &   zaPb(Qtilsl,p5)*zb(Qtilsl,p4))/
     & (sb*za(p2,Qtilsl)*za(p7,Qtilsl)**2*zaPb(Qtilsl,Qtilsl)**2)

c--- multiply by overall factor (including W propagators)
      res=res/(s(p5,p6)*sb)

c--- overall sign by hand
      bub12sl=-res

      return
      end

