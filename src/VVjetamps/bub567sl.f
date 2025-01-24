!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub567sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: p2,p5,p6,p4,p3,p1,p7,v1,v2,v3,
     & Qta1,Qta2,Qta3,Qtb1,Qtb2
      integer, parameter:: Qtil=8,Qtilsl=9
      real(dp):: p(mxpart,4),pin(mxpart,4),sa,sb,twopbDa,t,
     & a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,res,bub567sl
      logical:: swapz

c--- statement functions
      zab(v1,v2,v3)=za(v1,v2)*zb(v2,v3)

c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)+zab(v1,Qta3,v2)

      t(v1,v2,v3)=s(v1,v2)+s(v1,v3)+s(v2,v3)

c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P and Q
c---  P=p(Qta1)+p(Qta2)+p(Qta3), Q=P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
      Qta1=p5
      Qta2=p6
      Qta3=p7

      Qtb1=p2
      Qtb2=p1
      sb=s(Qtb1,Qtb2)
      sa=s(Qta1,Qta2)+s(Qta1,Qta3)+s(Qta2,Qta3)
      twopbDa=s(Qtb1,Qta1)+s(Qtb1,Qta2)+s(Qtb1,Qta3)
     &       +s(Qtb2,Qta1)+s(Qtb2,Qta2)+s(Qtb2,Qta3)

      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamm*sb
     &          +(p(Qtb1,:)+p(Qtb2,:))*sa
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamp*sb
     &            +(p(Qtb1,:)+p(Qtb2,:))*sa

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      res=
     &-((za(p5,p6)**2*za(p3,p1)*zb(p6,p5)*zb(p4,p2)*zb(p7,p5)**2)/
     &     (t(p4,p3,p1)*zaPb(p7,p7)**2)) -
     &  (za(p3,p1)*zb(p6,p5)*((sa*za(p2,p6)**2*zb(p4,p2))/
     &  (za(p2,p7)**2*zaPb(p2,p2)) -
     & (sa*za(p2,p6)*za(p6,p7)*zb(p7,p4))/
     &  (za(p2,p7)**2*zaPb(p7,p7)) -
     & (sa*za(p5,p6)*za(p6,p7)*zb(p7,p5)*zb(p7,p4))/
     &  (za(p2,p7)*zaPb(p7,p7)**2)))/t(p4,p3,p1) -
     &  zb(p6,p5)*(-((sa**2*za(p2,p6)**2*za(p2,p3)*za(p2,p1)*zb(p4,p2))/
     &    (za(p2,p7)**2*za(p2,Qtil)*zaPb(p2,p2)*zaPb(p2,Qtil))) +
     & (sa**2*za(p6,p7)**2*za(p3,p7)*zaPb(p1,Qtil)*zb(p7,p4))/
     &  (za(p2,p7)*za(p7,Qtil)*zaPb(p7,p7)*zaPb(p7,Qtil)**2) +
     & (sa**2*za(p6,p7)**2*za(p1,p7)*zaPb(p3,p7)*zb(p7,p4))/
     &  (za(p2,p7)*za(p7,Qtil)*zaPb(p7,p7)**2*zaPb(p7,Qtil)) +
     & (sa**2*za(p6,p7)*za(p6,Qtil)*za(p3,p7)*za(p1,p7)*zb(p7,p4))/
     &  (za(p2,p7)*za(p7,Qtil)**2*zaPb(p7,p7)*zaPb(p7,Qtil)) +
     & (sa**2*za(p2,p6)*za(p6,p7)*za(p3,p7)*za(p1,p7)*zb(p7,p4))/
     &  (za(p2,p7)**2*za(p7,Qtil)*zaPb(p7,p7)*zaPb(p7,Qtil)) -
     & (sa**2*za(p6,p7)*za(p6,Qtil)*za(p3,Qtil)*za(p1,Qtil)*
     &    zb(Qtil,p4))/(za(p2,p7)*za(p7,Qtil)**2*zaPb(Qtil,Qtil)**2) -
     & (sa**2*za(p2,p6)*za(p6,Qtil)*za(p3,Qtil)*za(p1,Qtil)*
     &    zb(Qtil,p4))/(za(p2,p7)*za(p2,Qtil)
     & *za(p7,Qtil)*zaPb(Qtil,Qtil)**2)
     &- (sa**2*za(p6,p7)*za(p6,Qtilsl)*za(p3,Qtilsl)*za(p1,Qtilsl)*
     &   zb(Qtilsl,p4))/(za(p2,p7)*za(p7,Qtilsl)**2
     & *zaPb(Qtilsl,Qtilsl)**2) -
     &(sa**2*za(p2,p6)*za(p6,Qtilsl)*za(p3,Qtilsl)*za(p1,Qtilsl)*
     &   zb(Qtilsl,p4))/
     & (za(p2,p7)*za(p2,Qtilsl)*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2))

c--- multiply by overall factor (including W propagators)
      res=res/(s(p5,p6)*s(p4,p3))

c--- overall sign by hand
      bub567sl=-res

      return
      end

