!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub347sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: p2,p5,p6,p4,p3,p1,p7,v1,v2,v3,
     & Qta1,Qta2,Qta3,Qtb1,Qtb2
      integer,parameter::Qtil=8,Qtilsl=9
      real(dp):: p(mxpart,4),pin(mxpart,4),sa,sb,twopbDa,t,
     & a,b,c,lamp,lamm
      complex(dp):: zab,zaPb,res,bub347sl
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
      Qta1=p4
      Qta2=p3
      Qta3=p7

      Qtb1=p5
      Qtb2=p6

      sb=s(Qtb1,Qtb2)
      sa=+s(Qta1,Qta2)+s(Qta1,Qta3)+s(Qta2,Qta3)
      twopbDa=s(Qta1,Qtb1)+s(Qta2,Qtb1)+s(Qta3,Qtb1)
     &       +s(Qta1,Qtb2)+s(Qta2,Qtb2)+s(Qta3,Qtb2)
      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamp
     &           +(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:)+p(Qta3,:))*lamm
     &           +(p(Qtb1,:)+p(Qtb2,:))*sa/sb

      if (swapz .eqv. .false.) then
        call spinoru(9,p,za,zb)
      else
        call spinoru(9,p,zb,za)
      endif

      res=
     &(za(p6,p1)*zaPb(p3,p7)**2*zb(p5,p2)*zb(p3,p4))/
     &   (t(p2,p5,p6)*zaPb(p7,p7)**2) -
     &  (zb(p5,p2)*zb(p3,p4)*(-((sa*za(p6,p1)*za(p3,p1)**2)/
     &   (za(p1,p7)**2*zaPb(p1,p1))) -
     &(sa*za(p6,p7)*za(p3,p7)*zaPb(p3,p1)*zaPb(p1,p3))/
     & (za(p1,p7)*zaPb(p7,p3)*zaPb(p7,p1)**2) -
     &(sa*za(p6,p7)*za(p3,p7)*zaPb(p3,p3)*zaPb(p1,p3))/
     & (za(p1,p7)*zaPb(p7,p3)**2*zaPb(p7,p1)) +
     &(sa*za(p6,p1)*za(p3,p7)**2*zaPb(p1,p3))/
     & (za(p1,p7)**2*zaPb(p7,p3)*zaPb(p7,p1)) +
     &(za(p6,p1)*za(p3,p1)*zaPb(p3,p3))/(za(p1,p7)**2*zb(p1,p3)) -
     &(za(p6,p7)*zaPb(p3,p3)**3)/(za(p3,p7)*zaPb(p7,p3)**2*zb(p1,p3)) +
     &(za(p6,p7)*zaPb(p3,p3)**2*zaPb(p1,p3))/
     & (za(p1,p7)*zaPb(p7,p3)**2*zb(p1,p3)) +
     &(za(p6,p3)*zaPb(p3,p3)**2)/(za(p3,p7)*zaPb(p7,p3)*zb(p1,p3)) -
     &(za(p6,p1)*za(p3,p7)*zaPb(p3,p3)*zaPb(p1,p3))/
     & (za(p1,p7)**2*zaPb(p7,p3)*zb(p1,p3)) +
     &(sa**2*za(p3,p7)**2*zaPb(p6,p7)*zb(p7,p3))/
     & (zaPb(p7,p3)*zaPb(p7,p1)*zaPb(p7,p7)**2) +
     &(sa**2*za(p6,p7)*za(p3,p7)*zaPb(p3,p1)*zb(p7,p3))/
     & (zaPb(p7,p3)*zaPb(p7,p1)**2*zaPb(p7,p7)) +
     &(sa**2*za(p6,p7)*za(p3,p7)*zaPb(p3,p3)*zb(p7,p3))/
     & (zaPb(p7,p3)**2*zaPb(p7,p1)*zaPb(p7,p7))))/t(p2,p5,p6) +
     &  zb(p3,p4)*(-((sa**2*za(p6,p7)*za(p3,p7)*zaPb(p3,Qtil)*
     &     zaPb(p7,p5)*zb(p7,p2))/
     &   (sb*za(p7,Qtil)*zaPb(p7,p1)*zaPb(p7,p7)*zaPb(p7,Qtil)**2)) -
     &(sa**2*za(p3,p7)**2*zaPb(p6,p1)*zaPb(p7,p5)*zb(p7,p2))/
     & (sb*za(p7,Qtil)*zaPb(p7,p1)**2*zaPb(p7,p7)*zaPb(p7,Qtil)) -
     &(sa**2*za(p6,p7)*za(p3,p7)*za(p3,Qtil)*zaPb(p7,p5)*zb(p7,p2))/
     & (sb*za(p7,Qtil)**2*zaPb(p7,p1)*zaPb(p7,p7)*zaPb(p7,Qtil)) +
     &(sa**3*za(p6,p7)*za(p3,p7)**2*zb(p7,p2)*zb(p7,p5))/
     & (sb*za(p7,Qtil)*zaPb(p7,p1)*zaPb(p7,p7)**2*zaPb(p7,Qtil)) +
     &(sa**2*za(p6,p7)*za(p3,Qtil)**2*zaPb(Qtil,p5)*zb(Qtil,p2))/
     & (sb*za(p7,Qtil)**2*zaPb(p7,p1)*zaPb(Qtil,Qtil)**2) +
     &(sa**2*za(p3,Qtil)**2*zaPb(p6,p1)*zaPb(Qtil,p5)*zb(Qtil,p2))/
     & (sb*za(p7,Qtil)*zaPb(p7,p1)*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2) +
     &(sa*zaPb(p6,p1)*zaPb(p3,p1)**2*zaPb(p1,p2)*zb(p1,p5))/
     & (sb*zaPb(p1,p1)*zaPb(p7,p1)**2*zaPb(Qtil,p1)*zb(Qtil,p1)) +
     &(sa**2*za(p6,p7)*za(p3,Qtilsl)**2*zaPb(Qtilsl,p5)*
     &   zb(Qtilsl,p2))/
     & (sb*za(p7,Qtilsl)**2*zaPb(p7,p1)*zaPb(Qtilsl,Qtilsl)**2) +
     &(sa**2*za(p3,Qtilsl)**2*zaPb(p6,p1)*zaPb(Qtilsl,p5)*
     &   zb(Qtilsl,p2))/
     & (sb*za(p7,Qtilsl)*zaPb(p7,p1)*zaPb(Qtilsl,p1)*
     &   zaPb(Qtilsl,Qtilsl)**2))

c--- multiply by overall factor (including W propagators)
      res=res/(sb*s(p4,p3))

c--- no overall sign here
      bub347sl=res

      return
      end

