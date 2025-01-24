!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function bub56sl(pin,p1,p2,p3,p4,p5,p6,p7,za,zb,swapz)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: p2,p4,p3,p5,p6,p1,p7,v1,v2,v3
      real(dp):: p(mxpart,4),pin(mxpart,4),a,b,c,lamp,lamm
      complex(dp):: bub56sl,zaPb,zaQb,zaPQa,zbPQb,zaPQPb,
     & zaQPa,zbQPb,res,zaQ2b,zaPQ2a,zbPQ2b,zaPQ2Pb,zaQ2Pa,zbQ2Pb,
     & zaQ3b,zaPQ3a,zaQ3Pa,zbPQ3b,zbQ3Pb,zaPQ3Pb
      integer:: Qta1,Qta2,Qtb1,Qtb2,Q2tb1,Q2tb2,Q2tb3,Q3tb1,Q3tb2,Q3tb3
      real(dp):: sa,sb,sb2,sb3,twopbDa,twopb2Da,twopb3Da,t
      integer, parameter:: Qtil=8,Qtilsl=9,Q2til=10,Q2tilsl=11
      integer, parameter:: Q3til=12,Q3tilsl=13
      logical:: swapz

c--- statement functions
      include 'zab.f'
c--- <v1|P|v2]
      zaPb(v1,v2)=zab(v1,Qta1,v2)+zab(v1,Qta2,v2)
c--- <v1|Q|v2]
      zaQb(v1,v2)=zaPb(v1,v2)+sa/sb*(zab(v1,Qtb1,v2)+zab(v1,Qtb2,v2))
c--- <v1|Q2|v2]
      zaQ2b(v1,v2)=zaPb(v1,v2)
     & +sa/sb2*(zab(v1,Q2tb1,v2)+zab(v1,Q2tb2,v2)+zab(v1,Q2tb3,v2))
c--- <v1|Q3|v2]
      zaQ3b(v1,v2)=zaPb(v1,v2)
     & +sa/sb3*(zab(v1,Q3tb1,v2)+zab(v1,Q3tb2,v2)+zab(v1,Q3tb3,v2))

c--- <v1|PQ|v2>
      zaPQa(v1,v2)=(
     & za(v1,Qta1)*zaQb(v2,Qta1)+za(v1,Qta2)*zaQb(v2,Qta2))
c--- <v1|PQ2|v2>
      zaPQ2a(v1,v2)=(
     & za(v1,Qta1)*zaQ2b(v2,Qta1)+za(v1,Qta2)*zaQ2b(v2,Qta2))
c--- <v1|PQ3|v2>
      zaPQ3a(v1,v2)=(
     & za(v1,Qta1)*zaQ3b(v2,Qta1)+za(v1,Qta2)*zaQ3b(v2,Qta2))

      zaQPa(v1,v2)=-zaPQa(v2,v1)
      zaQ2Pa(v1,v2)=-zaPQ2a(v2,v1)
      zaQ3Pa(v1,v2)=-zaPQ3a(v2,v1)

c--- [v1|PQ|v2]
      zbPQb(v1,v2)=(
     & zb(v1,Qta1)*zaQb(Qta1,v2)+zb(v1,Qta2)*zaQb(Qta2,v2))
c--- [v1|PQ2|v2]
      zbPQ2b(v1,v2)=(
     & zb(v1,Qta1)*zaQ2b(Qta1,v2)+zb(v1,Qta2)*zaQ2b(Qta2,v2))
c--- [v1|PQ3|v2]
      zbPQ3b(v1,v2)=(
     & zb(v1,Qta1)*zaQ3b(Qta1,v2)+zb(v1,Qta2)*zaQ3b(Qta2,v2))

      zbQPb(v1,v2)=-zbPQb(v2,v1)
      zbQ2Pb(v1,v2)=-zbPQ2b(v2,v1)
      zbQ3Pb(v1,v2)=-zbPQ3b(v2,v1)

c--- <v1|PQP|v2]
      zaPQPb(v1,v2)=(
     & -za(v1,Qta1)*zbPQb(v2,Qta1)-za(v1,Qta2)*zbPQb(v2,Qta2)
     & )
c--- <v1|PQ2P|v2]
      zaPQ2Pb(v1,v2)=(
     & -za(v1,Qta1)*zbPQ2b(v2,Qta1)-za(v1,Qta2)*zbPQ2b(v2,Qta2)
     & )
c--- <v1|PQ3P|v2]
      zaPQ3Pb(v1,v2)=(
     & -za(v1,Qta1)*zbPQ3b(v2,Qta1)-za(v1,Qta2)*zbPQ3b(v2,Qta2)
     & )

      t(v1,v2,v3)=s(v1,v2)+s(v1,v3)+s(v2,v3)

c--- copy momentum array
      p(:,:)=pin(:,:)

c--- these define the momenta P, Q, Q2 and Q3
c---  P =p(Qta1)+p(Qta2)
c---  Q =P^2 (p(Qtb1)+p(Qtb2)) + (p(Qtb1)+p(Qtb2))^2 P
c---  Q2=P^2 (p(Q2tb1)+p(Q2tb2)+p(Q2tb3)) + (p(Q2tb1)+p(Q2tb2)+p(Q2tb3))^2 P
c---  Q3=P^2 (p(Q3tb1)+p(Q3tb2)+p(Q3tb3)) + (p(Q3tb1)+p(Q3tb2)+p(Q3tb3))^2 P
      Qta1=p5
      Qta2=p6

      Qtb1=p4
      Qtb2=p3

      Q2tb1=p2
      Q2tb2=p4
      Q2tb3=p3

      Q3tb1=p4
      Q3tb2=p3
      Q3tb3=p7

c     soln1
      sb=s(Qtb1,Qtb2)
      sa=s(Qta1,Qta2)
      twopbDa=s(Qta1,Qtb1)+s(Qta2,Qtb1)+s(Qta1,Qtb2)+s(Qta2,Qtb2)
      a=sb
      b=twopbDa
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Qtil,:) = (p(Qta1,:)+p(Qta2,:))*lamm+(p(Qtb1,:)+p(Qtb2,:))*sa/sb
      p(Qtilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamp+(p(Qtb1,:)+p(Qtb2,:))*sa/sb

c     soln2
      sb2=s(Q2tb1,Q2tb2)+s(Q2tb1,Q2tb3)+s(Q2tb2,Q2tb3)
      twopb2Da=s(Q2tb1,Qta1)+s(Q2tb2,Qta1)+s(Q2tb3,Qta1)
     &        +s(Q2tb1,Qta2)+s(Q2tb2,Qta2)+s(Q2tb3,Qta2)
      a=sb2
      b=twopb2Da
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Q2til,:) = (p(Qta1,:)+p(Qta2,:))*lamm
     &             +(p(Q2tb1,:)+p(Q2tb2,:)+p(Q2tb3,:))*sa/sb2
      p(Q2tilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamp
     &             +(p(Q2tb1,:)+p(Q2tb2,:)+p(Q2tb3,:))*sa/sb2

c     soln3
      sb3=+s(Q3tb1,Q3tb2)+s(Q3tb1,Q3tb3)+s(Q3tb2,Q3tb3)
      twopb3Da=s(Q3tb1,Qta1)+s(Q3tb2,Qta1)+s(Q3tb3,Qta1)
     &        +s(Q3tb1,Qta2)+s(Q3tb2,Qta2)+s(Q3tb3,Qta2)
      a=sb3
      b=twopb3Da
      c=sa
      call solvequadratic(a,b,c,lamp,lamm)
      p(Q3til,:) = (p(Qta1,:)+p(Qta2,:))*lamm
     &             +(p(Q3tb1,:)+p(Q3tb2,:)+p(Q3tb3,:))*sa/sb3
      p(Q3tilsl,:)=(p(Qta1,:)+p(Qta2,:))*lamp
     &             +(p(Q3tb1,:)+p(Q3tb2,:)+p(Q3tb3,:))*sa/sb3


      if (swapz .eqv. .false.) then
        call spinoru(13,p,za,zb)
      else
        call spinoru(13,p,zb,za)
      endif

      res=
     &((zab(p3,p2,p7) + zab(p3,p4,p7))*zb(p4,p2)*
     & ((sa*za(p6,p1)*zaPb(p1,p5))/zaPb(p1,p1) +
     & (zaPb(p6,p5)*zaPb(p1,p5))/(2._dp*zb(p1,p5))-
     & (zaPb(p6,p1)*zaPb(p1,p5)**2)/(2._dp*zaPb(p1,p1)*zb(p1,p5))))/
     & (sb2*t(p5,p6,p1)*za(p1,p7))-
     & ((zab(p1,p2,p4) + zab(p1,p7,p4))*
     & (-((sa*za(p6,p1)*zaPb(p3,p1)*zaPb(p1,p5))/zaPb(p1,p1)**2)-
     & (zaPb(p3,p5)*zaPb(p6,p5))/(2._dp*zb(p1,p5)) +
     & (zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)**2)/
     &  (2._dp*zaPb(p1,p1)**2*zb(p1,p5)) +
     & (sa**2*za(p3,p1)*za(p6,p1)*zb(p1,p5))/zaPb(p1,p1)**2))/
     & (t(p5,p6,p1)*za(p2,p7)*za(p1,p7)) +
     & ((zab(p1,p2,p4) + zab(p1,p7,p4))*za(p3,p1)*
     & (zaPb(p6,p5)/2._dp- (sa*za(p6,p1)*zb(p1,p5))/zaPb(p1,p1)))/
     & (t(p5,p6,p1)*za(p2,p7)*za(p1,p7)) +
     &(zb(p4,p2)*(-(za(p2,p3)*(-((sa*za(p6,p1)*zaPb(p1,p5)*zb(p1,p2))/
     &  (zaPb(p1,p1)*zaPb(p7,p1)))-
     & (zaPb(p6,p5)*zaPb(p1,p5)*zb(p5,p2))/
     & (2._dp*zaPb(p7,p5)*zb(p1,p5)) +
     & (zaPb(p6,p1)*zaPb(p1,p5)**2*zb(p1,p2))/
     & (2._dp*zaPb(p1,p1)*zaPb(p7,p1)*zb(p1,p5))-
     & (sa*za(p6,p1)*zaPb(p7,p2)*zb(p7,p5))/
     & (zaPb(p7,p1)*zaPb(p7,p7))-
     & (sa**2*za(p6,p7)*za(p1,p7)*zaPb(p7,p2)*zb(p7,p5)**2)/
     & (2._dp*zaPb(p7,p5)*zaPb(p7,p1)*zaPb(p7,p7)**2)))-
     &za(p4,p3)*(-((sa*za(p6,p1)*zaPb(p1,p5)*zb(p1,p4))/
     & (zaPb(p1,p1)*zaPb(p7,p1)))-
     &(zaPb(p6,p5)*zaPb(p1,p5)*zb(p5,p4))/(2._dp*zaPb(p7,p5)*zb(p1,p5))+
     &(zaPb(p6,p1)*zaPb(p1,p5)**2*zb(p1,p4))/
     & (2._dp*zaPb(p1,p1)*zaPb(p7,p1)*zb(p1,p5))-
     &(sa*za(p6,p1)*zaPb(p7,p4)*zb(p7,p5))/
     & (zaPb(p7,p1)*zaPb(p7,p7))-
     &(sa**2*za(p6,p7)*za(p1,p7)*zaPb(p7,p4)*zb(p7,p5)**2)/
     & (2._dp*zaPb(p7,p5)*zaPb(p7,p1)*zaPb(p7,p7)**2))))/
     &(sb2*za(p1,p7))- (zb(p4,p2)*
     & (za(p2,p3)*zb(p7,p2) + za(p4,p3)*zb(p7,p4))*
     & ((sa**2*za(p6,p7)*za(p1,p7)*zb(p7,p5))/
     & (sb2*zaPb(p7,p7)*zaPQ2a(p7,p7)) +
     &  (sa**2*za(p1,Q2til)*zaPQ2a(p6,Q2til)*zb(Q2til,p5))/
     & (sb2*zaPb(Q2til,Q2til)**2*zaPQ2a(p7,Q2til)) +
     &  (sa**2*za(p1,Q2tilsl)*zaPQ2a(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**2*zaPQ2a(p7,Q2tilsl))))/
     &sb2- (za(p3,p1)*zb(p3,p4)*
     & (-((sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p7,p2)*zb(p7,p5))/
     & (sb*sb2*zaPb(p7,p7)*zaPQ2a(p7,p7)*zaPQa(p7,p7)))-
     &  (sa**3*za(p3,Q2til)*za(p6,Q2til)*za(p1,Q2til)*
     & zaPb(Q2til,p2)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**2*
     & zaPQa(Q2til,Q2til))-
     &  (sa**3*za(p3,Q2tilsl)*za(p6,Q2tilsl)*za(p1,Q2tilsl)*
     & zaPb(Q2tilsl,p2)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**2*
     & zaPQa(Q2tilsl,Q2tilsl))-
     &  (sa**3*za(p3,Qtil)*za(p6,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*
     & zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     & zaPQ2a(Qtil,Qtil))-
     &  (sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*za(p1,Qtilsl)*
     & zaPb(Qtilsl,p2)*zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ2a(Qtilsl,Qtilsl))))/za(p1,p7) +
     &zb(p7,p2)*((sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p7,p4)*
     & zb(p7,p5))/
     & (sb*sb2*zaPb(p7,p7)*zaPQ2a(p7,p7)*zaPQa(p7,p7)) +
     &(sa**3*za(p3,Q2til)*za(p6,Q2til)*za(p1,Q2til)*zaPb(Q2til,p4)*
     & zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**2*
     & zaPQa(Q2til,Q2til)) + (sa**3*za(p3,Q2tilsl)*za(p6,Q2tilsl)*
     & za(p1,Q2tilsl)*zaPb(Q2tilsl,p4)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**2*
     & zaPQa(Q2tilsl,Q2tilsl)) +
     &(sa**3*za(p3,Qtil)*za(p6,Qtil)*za(p1,Qtil)*zaPb(Qtil,p4)*
     & zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     & zaPQ2a(Qtil,Qtil)) + (sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*
     & za(p1,Qtilsl)*zaPb(Qtilsl,p4)*zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ2a(Qtilsl,Qtilsl)))

      res=res +
     &(-((sa**3*za(p3,p7)*za(p1,p7)**2*zaPb(p6,Qtil)*zaPb(p7,p2)*
     & zaPb(p7,p4)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)*za(p7,Qtil)*zaPb(p7,p7)*
     & zaPb(p7,Q2til)*zaPb(p7,Qtil)**2))-
     &(sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p1,Q2til)*
     & zaPb(p7,p2)*zaPb(p7,p4)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)*za(p7,Qtil)*zaPb(p7,p7)*
     & zaPb(p7,Q2til)**2*zaPb(p7,Qtil))-
     &(sa**3*za(p3,Qtil)*za(p6,p7)*za(p1,p7)**2*zaPb(p7,p2)*
     & zaPb(p7,p4)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)*za(p7,Qtil)**2*zaPb(p7,p7)*
     & zaPb(p7,Q2til)*zaPb(p7,Qtil))-
     &(sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*za(p1,Q2til)*zaPb(p7,p2)*
     & zaPb(p7,p4)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)**2*za(p7,Qtil)*zaPb(p7,p7)*
     & zaPb(p7,Q2til)*zaPb(p7,Qtil))-
     &(sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p7,p2)*zaPb(p7,p4)*
     & zaPQ2a(p1,p7)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)**2*za(p7,Qtil)*zaPb(p7,p7)*
     & zaPb(p7,Q2til)**2*zaPb(p7,Qtil))-
     &(sa**3*za(p3,p7)*za(p1,p7)**2*zaPb(p7,p2)*zaPb(p7,p4)*
     & zaPQa(p6,p7)*zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)*za(p7,Qtil)**2*zaPb(p7,p7)*
     & zaPb(p7,Q2til)*zaPb(p7,Qtil)**2) +
     &(sa**4*za(p3,p7)*za(p6,p7)*za(p1,p7)**2*zaPb(p7,p4)*zb(p7,p2)*
     & zb(p7,p5))/
     & (sb*sb2*za(p7,Q2til)*za(p7,Qtil)*zaPb(p7,p7)**2*
     & zaPb(p7,Q2til)*zaPb(p7,Qtil))-
     &(sa**4*za(p3,p7)*za(p6,p7)*za(p1,p7)**2*zaPb(p7,p2)*
     & zaPb(p7,p4)*zb(p7,p5)**2)/
     & (2._dp*sb*sb2*za(p7,Q2til)*za(p7,Qtil)*zaPb(p7,p5)*
     & zaPb(p7,p7)**2*zaPb(p7,Q2til)*zaPb(p7,Qtil))-
     &(sa**3*za(p3,Q2til)*za(p1,Q2til)*za(Q2til,Qtil)*zaPb(p1,Qtil)*
     &zaPb(Q2til,p2)*zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**3*
     & zaPQa(Q2til,Q2til)**2) +
     &(sa**3*za(p3,Q2til)*za(p6,Q2til)*za(p1,Q2til)*zaPb(Q2til,p2)*
     & zaPb(Q2til,p4)*zaPQa(p1,Q2til)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**2*
     & zaPQa(Q2til,Q2til)**2)-
     &(sa**3*za(p3,Q2til)*za(p6,Q2til)*za(p1,Q2til)**2*zaPb(p7,p4)*
     & zaPb(Q2til,p2)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)**2*zaPb(Q2til,Q2til)**2*
     & zaPQa(Q2til,Q2til))- (sa**3*za(p3,p7)*za(p1,Q2til)**2*
     & zaPb(Q2til,p2)*zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)**2*zaPb(Q2til,Q2til)**3*
     & zaPQa(Q2til,Q2til)) + (sa**3*za(p3,Q2til)*za(p1,Q2til)*
     & za(p1,Qtil)*zaPb(Q2til,p2)*zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*
     & zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*za(Q2til,Qtil)*
     & zaPb(Q2til,Q2til)**3*zaPQa(Q2til,Q2til)) +
     &(sa**3*za(p3,Q2til)*za(p1,Q2til)**2*zaPb(Q2til,p2)*
     & zaPb(Q2til,p4)*zaPQ2Pb(p6,Q2til)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**4*
     & zaPQa(Q2til,Q2til))- (sa**4*za(p3,Q2til)*za(p1,Q2til)**2*
     & zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*zb(Q2til,p2)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**4*
     & zaPQa(Q2til,Q2til))- (sa**4*za(p3,Q2til)*za(p1,Q2til)**2*
     & zaPb(Q2til,p2)*zaPQ2a(p6,Q2til)*zb(Q2til,p4)*zb(Q2til,p5))/
     & (sb*sb2*za(p7,Q2til)*zaPb(Q2til,Q2til)**4*
     & zaPQa(Q2til,Q2til))- (sa**4*za(p3,Q2til)*za(p6,Q2til)*
     & za(p1,Q2til)**2*zaPb(Q2til,p2)*zaPb(Q2til,p4)*zb(Q2til,p5)**2)/
     & (2._dp*sb*sb2*za(p7,Q2til)*zaPb(Q2til,p5)*
     & zaPb(Q2til,Q2til)**3*zaPQa(Q2til,Q2til))-
     &(sa**3*za(p3,Q2tilsl)*za(p1,Q2tilsl)*za(Q2tilsl,Qtil)*
     & zaPb(p1,Qtil)*zaPb(Q2tilsl,p2)*zaPb(Q2tilsl,p4)*
     & zaPQ2a(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**3*
     & zaPQa(Q2tilsl,Q2tilsl)**2) +
     &(sa**3*za(p3,Q2tilsl)*za(p6,Q2tilsl)*za(p1,Q2tilsl)*
     & zaPb(Q2tilsl,p2)*zaPb(Q2tilsl,p4)*zaPQa(p1,Q2tilsl)*
     & zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**2*
     & zaPQa(Q2tilsl,Q2tilsl)**2)-
     &(sa**3*za(p3,Q2tilsl)*za(p6,Q2tilsl)*za(p1,Q2tilsl)**2*
     & zaPb(p7,p4)*zaPb(Q2tilsl,p2)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)**2*zaPb(Q2tilsl,Q2tilsl)**2*
     &zaPQa(Q2tilsl,Q2tilsl))-
     &(sa**3*za(p3,p7)*za(p1,Q2tilsl)**2*zaPb(Q2tilsl,p2)*
     & zaPb(Q2tilsl,p4)*zaPQ2a(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)**2*zaPb(Q2tilsl,Q2tilsl)**3*
     & zaPQa(Q2tilsl,Q2tilsl)) +
     &(sa**3*za(p3,Q2tilsl)*za(p1,Q2tilsl)*za(p1,Qtil)*
     & zaPb(Q2tilsl,p2)*zaPb(Q2tilsl,p4)*zaPQ2a(p6,Q2tilsl)*
     & zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*za(Q2tilsl,Qtil)*
     & zaPb(Q2tilsl,Q2tilsl)**3*zaPQa(Q2tilsl,Q2tilsl)) +
     &(sa**3*za(p3,Q2tilsl)*za(p1,Q2tilsl)**2*zaPb(Q2tilsl,p2)*
     & zaPb(Q2tilsl,p4)*zaPQ2Pb(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**4*
     & zaPQa(Q2tilsl,Q2tilsl))-
     & (sa**4*za(p3,Q2tilsl)*za(p1,Q2tilsl)**2*zaPb(Q2tilsl,p4)*
     & zaPQ2a(p6,Q2tilsl)*zb(Q2tilsl,p2)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**4*
     & zaPQa(Q2tilsl,Q2tilsl))-
     & (sa**4*za(p3,Q2tilsl)*za(p1,Q2tilsl)**2*zaPb(Q2tilsl,p2)*
     & zaPQ2a(p6,Q2tilsl)*zb(Q2tilsl,p4)*zb(Q2tilsl,p5))/
     & (sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,Q2tilsl)**4*
     & zaPQa(Q2tilsl,Q2tilsl))-
     & (sa**4*za(p3,Q2tilsl)*za(p6,Q2tilsl)*za(p1,Q2tilsl)**2*
     & zaPb(Q2tilsl,p2)*zaPb(Q2tilsl,p4)*zb(Q2tilsl,p5)**2)/
     & (2._dp*sb*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,p5)*
     & zaPb(Q2tilsl,Q2tilsl)**3*zaPQa(Q2tilsl,Q2tilsl)) +
     & (sa**2*zaPb(p3,p5)*zaPb(p6,p5)*zaPb(p1,p5)**2*zb(p5,p2)*
     & zb(p5,p4))/
     & (2._dp*sb*sb2*zaPb(p7,p5)*zaPb(Q2til,p5)*zaPb(Qtil,p5)*
     & zb(Q2til,p5)*zb(Qtil,p5)) +
     & (sa**3*za(p3,Qtil)*za(p6,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*
     & zaPb(Qtil,p4)*zaPQ2a(p1,Qtil)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**2*
     & zaPQ2a(Qtil,Qtil)**2)-
     & (sa**3*za(p3,Qtil)*za(p6,Qtil)*za(p1,Qtil)**2*zaPb(p7,p4)*
     & zaPb(Qtil,p2)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)**2*zaPb(Qtil,Qtil)**2*
     & zaPQ2a(Qtil,Qtil)) + (sa**3*za(p3,Qtil)*za(p1,Qtil)*
     & za(Q2til,Qtil)*zaPb(p1,Q2til)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zaPQa(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**3*
     & zaPQ2a(Qtil,Qtil)**2)-
     &(sa**3*za(p3,p7)*za(p1,Qtil)**2*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zaPQa(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)**2*zaPb(Qtil,Qtil)**3*
     & zaPQ2a(Qtil,Qtil))- (sa**3*za(p3,Qtil)*za(p1,Q2til)*
     & za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*zaPQa(p6,Qtil)*
     & zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*za(Q2til,Qtil)*zaPb(Qtil,Qtil)**3*
     &zaPQ2a(Qtil,Qtil)) + (sa**3*za(p3,Qtil)*za(p1,Qtil)**2*
     &zaPb(Qtil,p2)*zaPb(Qtil,p4)*zaPQPb(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**4*
     & zaPQ2a(Qtil,Qtil))- (sa**4*za(p3,Qtil)*za(p1,Qtil)**2*
     & zaPb(Qtil,p4)*zaPQa(p6,Qtil)*zb(Qtil,p2)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**4*
     & zaPQ2a(Qtil,Qtil))- (sa**4*za(p3,Qtil)*za(p1,Qtil)**2*
     & zaPb(Qtil,p2)*zaPQa(p6,Qtil)*zb(Qtil,p4)*zb(Qtil,p5))/
     & (sb*sb2*za(p7,Qtil)*zaPb(Qtil,Qtil)**4*
     & zaPQ2a(Qtil,Qtil))- (sa**4*za(p3,Qtil)*za(p6,Qtil)*
     & za(p1,Qtil)**2*zaPb(Qtil,p2)*zaPb(Qtil,p4)*zb(Qtil,p5)**2)/
     & (2._dp*sb*sb2*za(p7,Qtil)*zaPb(Qtil,p5)*
     & zaPb(Qtil,Qtil)**3*zaPQ2a(Qtil,Qtil)) +
     &(sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*za(p1,Qtilsl)*
     &zaPb(Qtilsl,p2)*zaPb(Qtilsl,p4)*zaPQ2a(p1,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ2a(Qtilsl,Qtilsl)**2)-
     &(sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*za(p1,Qtilsl)**2*
     & zaPb(p7,p4)*zaPb(Qtilsl,p2)*zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)**2*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ2a(Qtilsl,Qtilsl)) +
     &(sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*za(Q2til,Qtilsl)*
     &zaPb(p1,Q2til)*zaPb(Qtilsl,p2)*zaPb(Qtilsl,p4)*zaPQa(p6,Qtilsl)*
     & zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQ2a(Qtilsl,Qtilsl)**2)-
     &(sa**3*za(p3,p7)*za(p1,Qtilsl)**2*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb2*za(p7,Qtilsl)**2*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQ2a(Qtilsl,Qtilsl))-
     & (sa**3*za(p3,Qtilsl)*za(p1,Q2til)*za(p1,Qtilsl)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p4)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p5))/
     &  (sb*sb2*za(p7,Qtilsl)*za(Q2til,Qtilsl)*
     & zaPb(Qtilsl,Qtilsl)**3*zaPQ2a(Qtilsl,Qtilsl)) +
     & (sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)**2*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPQPb(p6,Qtilsl)*zb(Qtilsl,p5))/
     &  (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQ2a(Qtilsl,Qtilsl))-
     & (sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)**2*zaPb(Qtilsl,p4)*
     & zaPQa(p6,Qtilsl)*zb(Qtilsl,p2)*zb(Qtilsl,p5))/
     &  (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQ2a(Qtilsl,Qtilsl))-
     & (sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)**2*zaPb(Qtilsl,p2)*
     & zaPQa(p6,Qtilsl)*zb(Qtilsl,p4)*zb(Qtilsl,p5))/
     &  (sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQ2a(Qtilsl,Qtilsl))-
     & (sa**4*za(p3,Qtilsl)*za(p6,Qtilsl)*za(p1,Qtilsl)**2*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p4)*zb(Qtilsl,p5)**2)/
     &  (2._dp*sb*sb2*za(p7,Qtilsl)*zaPb(Qtilsl,p5)*
     & zaPb(Qtilsl,Qtilsl)**3*zaPQ2a(Qtilsl,Qtilsl)))/za(p1,p7) +
     &(zb(p4,p2)*(za(p2,p3)*(-((sa**2*za(p6,p7)*za(p1,p7)*zaPb(p7,p2)*
     &  zaQ2Pa(p1,p7)*zb(p7,p5))/
     & (sb2*zaPb(p7,p7)*zaPQ2a(p7,p7)**2))-
     & (sa**3*za(p6,p7)*za(p1,p7)**2*zb(p7,p2)*zb(p7,p5))/
     & (sb2*zaPb(p7,p7)**2*zaPQ2a(p7,p7)) +
     & (sa**3*za(p6,p7)*za(p1,p7)**2*zaPb(p7,p2)*zb(p7,p5)**2)/
     & (2._dp*sb2*zaPb(p7,p5)*zaPb(p7,p7)**2*zaPQ2a(p7,p7))-
     & (sa**2*za(p1,Q2til)*zaPb(Q2til,p2)*zaPQ2a(p1,Q2til)*
     & zaPQ2Pb(p6,Q2til)*zb(Q2til,p5))/
     & (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til))-
     & (sa**2*za(p1,Q2til)*zaPb(Q2til,p2)*zaPQ2a(p6,Q2til)*
     & zaPQ2Pb(p1,Q2til)*zb(Q2til,p5))/
     & (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til)) +
     & (sa**2*zaPb(Q2til,p2)*zaPQ2a(p6,Q2til)*zaPQ2a(p1,Q2til)*
     & zaQ2Pa(p1,p7)*zb(Q2til,p5))/
     & (sb2*zaPb(Q2til,Q2til)**3*zaPQ2a(p7,Q2til)**2) +
     & (sa**3*za(p1,Q2til)*zaPQ2a(p6,Q2til)*zaPQ2a(p1,Q2til)*
     & zb(Q2til,p2)*zb(Q2til,p5))/
     & (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til)) +
     & (sa**3*za(p6,Q2til)*za(p1,Q2til)**2*zaPb(Q2til,p2)*
     & zb(Q2til,p5)**2)/
     & (2._dp*sb2*za(p7,Q2til)*zaPb(Q2til,p5)*
     & zaPb(Q2til,Q2til)**3)-
     & (sa**2*za(p1,Q2tilsl)*zaPb(Q2tilsl,p2)*
     & zaPQ2a(p1,Q2tilsl)*zaPQ2Pb(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl))-
     & (sa**2*za(p1,Q2tilsl)*zaPb(Q2tilsl,p2)*
     & zaPQ2a(p6,Q2tilsl)*zaPQ2Pb(p1,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl)) +
     & (sa**2*zaPb(Q2tilsl,p2)*zaPQ2a(p6,Q2tilsl)*
     & zaPQ2a(p1,Q2tilsl)*zaQ2Pa(p1,p7)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**3*zaPQ2a(p7,Q2tilsl)**2) +
     & (sa**3*za(p1,Q2tilsl)*zaPQ2a(p6,Q2tilsl)*
     & zaPQ2a(p1,Q2tilsl)*zb(Q2tilsl,p2)*zb(Q2tilsl,p5))/
     &  (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl)) +
     & (sa**3*za(p6,Q2tilsl)*za(p1,Q2tilsl)**2*
     & zaPb(Q2tilsl,p2)*zb(Q2tilsl,p5)**2)/
     &  (2._dp*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,p5)*
     & zaPb(Q2tilsl,Q2tilsl)**3)-
     & (sa*zaPb(p5,p6)*zaPb(p1,p5)**2*zb(p5,p2))/
     &  (2._dp*sb2*zaPb(p7,p5)*zbQ2Pb(p5,p5))) +
     &za(p4,p3)*(-((sa**2*za(p6,p7)*za(p1,p7)*zaPb(p7,p4)*
     &  zaQ2Pa(p1,p7)*zb(p7,p5))/
     & (sb2*zaPb(p7,p7)*zaPQ2a(p7,p7)**2))-
     & (sa**3*za(p6,p7)*za(p1,p7)**2*zb(p7,p4)*zb(p7,p5))/
     &  (sb2*zaPb(p7,p7)**2*zaPQ2a(p7,p7)) +
     & (sa**3*za(p6,p7)*za(p1,p7)**2*zaPb(p7,p4)*zb(p7,p5)**2)/
     &  (2._dp*sb2*zaPb(p7,p5)*zaPb(p7,p7)**2*zaPQ2a(p7,p7))-
     & (sa**2*za(p1,Q2til)*zaPb(Q2til,p4)*zaPQ2a(p1,Q2til)*
     & zaPQ2Pb(p6,Q2til)*zb(Q2til,p5))/
     &  (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til))-
     & (sa**2*za(p1,Q2til)*zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*
     & zaPQ2Pb(p1,Q2til)*zb(Q2til,p5))/
     &  (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til)) +
     & (sa**2*zaPb(Q2til,p4)*zaPQ2a(p6,Q2til)*zaPQ2a(p1,Q2til)*
     & zaQ2Pa(p1,p7)*zb(Q2til,p5))/
     &  (sb2*zaPb(Q2til,Q2til)**3*zaPQ2a(p7,Q2til)**2) +
     & (sa**3*za(p1,Q2til)*zaPQ2a(p6,Q2til)*zaPQ2a(p1,Q2til)*
     & zb(Q2til,p4)*zb(Q2til,p5))/
     &  (sb2*zaPb(Q2til,Q2til)**4*zaPQ2a(p7,Q2til)) +
     & (sa**3*za(p6,Q2til)*za(p1,Q2til)**2*zaPb(Q2til,p4)*
     & zb(Q2til,p5)**2)/
     & (2._dp*sb2*za(p7,Q2til)*zaPb(Q2til,p5)*
     & zaPb(Q2til,Q2til)**3)-
     &(sa**2*za(p1,Q2tilsl)*zaPb(Q2tilsl,p4)*
     & zaPQ2a(p1,Q2tilsl)*zaPQ2Pb(p6,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl))-
     &(sa**2*za(p1,Q2tilsl)*zaPb(Q2tilsl,p4)*
     & zaPQ2a(p6,Q2tilsl)*zaPQ2Pb(p1,Q2tilsl)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl)) +
     &(sa**2*zaPb(Q2tilsl,p4)*zaPQ2a(p6,Q2tilsl)*
     & zaPQ2a(p1,Q2tilsl)*zaQ2Pa(p1,p7)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**3*zaPQ2a(p7,Q2tilsl)**2) +
     &(sa**3*za(p1,Q2tilsl)*zaPQ2a(p6,Q2tilsl)*
     & zaPQ2a(p1,Q2tilsl)*zb(Q2tilsl,p4)*zb(Q2tilsl,p5))/
     & (sb2*zaPb(Q2tilsl,Q2tilsl)**4*zaPQ2a(p7,Q2tilsl)) +
     &(sa**3*za(p6,Q2tilsl)*za(p1,Q2tilsl)**2*zaPb(Q2tilsl,p4)*
     & zb(Q2tilsl,p5)**2)/
     & (2._dp*sb2*za(p7,Q2tilsl)*zaPb(Q2tilsl,p5)*
     & zaPb(Q2tilsl,Q2tilsl)**3)-
     &(sa*zaPb(p5,p6)*zaPb(p1,p5)**2*zb(p5,p4))/
     & (2._dp*sb2*zaPb(p7,p5)*zbQ2Pb(p5,p5)))))/
     & (sb2*za(p1,p7)) + (za(p2,p1)*zb(p4,p2)*
     & ((sa**2*za(p3,p7)*za(p6,p7)*zaPb(p7,p2)*zb(p7,p5))/
     & (sb3*zaPb(p7,p1)*zaPb(p7,p7)*zaPQ3a(p7,p7)) +
     & (sa**2*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     &  zaPQ3a(p6,Q3til)*zb(Q3til,p5))/
     & (sb3*za(p7,Q3til)*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**2*
     &  zaPQ3a(p1,Q3til)) + (sa**2*za(p3,Q3tilsl)*za(p1,Q3tilsl)*
     &  zaPb(Q3tilsl,p2)*zaPQ3a(p6,Q3tilsl)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p1)*
     &  zaPb(Q3tilsl,Q3tilsl)**2*zaPQ3a(p1,Q3tilsl))-
     & (sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p1,p2))/
     & (sb3*zaPb(p1,p1)*zaPb(p7,p1)*zbQ3Pb(p1,p1))))/za(p1,p7)-
     &(-((sa**2*za(p3,p7)*za(p1,p7)*zaPb(p7,p2)*zaPb(p7,p4)*
     &  zaPQ3a(p6,p7)*zb(p7,p5))/
     & (sb3*zaPb(p7,p1)*zaPb(p7,p7)*zaPQ3a(p7,p7)**2)) +
     & (sa**3*za(p3,p7)*za(p6,p7)*za(p1,p7)*zaPb(p7,p2)*zaPb(p7,p4)*
     & zb(p7,p5)**2)/
     &  (2._dp*sb3*zaPb(p7,p5)*zaPb(p7,p1)*zaPb(p7,p7)**2*
     & zaPQ3a(p7,p7))- (sa*zaPb(p3,p5)*zaPb(p6,p5)*zaPb(p1,p5)*
     & zb(p5,p2)*zb(p5,p4))/
     &  (2._dp*sb3*zaPb(p7,p5)*zaPb(Q3til,p5)*zb(p1,p5)*zb(Q3til,p5))-
     & (sa**2*za(p1,Q3til)*zaPb(p3,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zaPQ3a(p6,Q3til)*zb(Q3til,p5))/
     &  (sb3*za(p7,Q3til)*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**4)-
     & (sa**2*za(p3,Q3til)*zaPb(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zaPQ3a(p6,Q3til)*zb(Q3til,p5))/
     &  (sb3*za(p7,Q3til)*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**4) +
     & (sa**2*za(p3,Q3til)*za(p1,Q3til)*zaPb(p7,p4)*zaPb(Q3til,p2)*
     & zaPQ3a(p6,Q3til)*zb(Q3til,p5))/
     &  (sb3*za(p7,Q3til)**2*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**3)-
     & (sa**2*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zaPQ3Pb(p6,Q3til)*zb(Q3til,p5))/
     &  (sb3*za(p7,Q3til)*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**4) +
     & (sa**3*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p4)*
     & zaPQ3a(p6,Q3til)*zb(p1,p2)*zb(Q3til,p5))/
     &  (sb3*za(p7,Q3til)*zaPb(Q3til,p1)**2*zaPb(Q3til,Q3til)**3) +
     & (sa**3*za(p3,Q3til)*za(p6,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zb(Q3til,p5)**2)/
     & (2._dp*sb3*za(p7,Q3til)*zaPb(Q3til,p5)*zaPb(Q3til,p1)*
     &  zaPb(Q3til,Q3til)**3)-
     & (sa**2*za(p1,Q3tilsl)*zaPb(p3,Q3tilsl)*zaPb(Q3tilsl,p2)*
     &  zaPb(Q3tilsl,p4)*zaPQ3a(p6,Q3tilsl)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p1)*
     &  zaPb(Q3tilsl,Q3tilsl)**4)-
     & (sa**2*za(p3,Q3tilsl)*zaPb(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     &  zaPb(Q3tilsl,p4)*zaPQ3a(p6,Q3tilsl)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p1)*
     &  zaPb(Q3tilsl,Q3tilsl)**4) +
     & (sa**2*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(p7,p4)*
     &  zaPb(Q3tilsl,p2)*zaPQ3a(p6,Q3tilsl)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)**2*zaPb(Q3tilsl,p1)*
     &  zaPb(Q3tilsl,Q3tilsl)**3)-
     & (sa**2*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     &  zaPb(Q3tilsl,p4)*zaPQ3Pb(p6,Q3tilsl)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p1)*
     &  zaPb(Q3tilsl,Q3tilsl)**4) +
     & (sa**3*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p4)*
     &  zaPQ3a(p6,Q3tilsl)*zb(p1,p2)*zb(Q3tilsl,p5))/
     & (sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p1)**2*
     &  zaPb(Q3tilsl,Q3tilsl)**3) +
     & (sa**3*za(p3,Q3tilsl)*za(p6,Q3tilsl)*za(p1,Q3tilsl)*
     &  zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p4)*zb(Q3tilsl,p5)**2)/
     &  (2._dp*sb3*za(p7,Q3tilsl)*zaPb(Q3tilsl,p5)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**3) +
     & (sa*zaPb(p3,p1)*zaPb(p1,p5)*zaPQ3Pb(p6,p1)*zb(p1,p2)*
     & zb(p1,p4))/(sb3*zaPb(p7,p1)*zbQ3Pb(p1,p1)**2)-
     & (sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)**2*zb(p1,p2)*zb(p1,p4))/
     &  (2._dp*sb3*zaPb(p1,p1)*zaPb(p7,p1)*zb(p1,p5)*zbQ3Pb(p1,p1)))/
     & za(p1,p7) + (za(p3,p1)*zb(p3,p4)*
     & (za(p2,p1)*(-((sa**2*za(p3,p1)*za(p6,p1)*zaPb(p1,p2)*
     & zb(p1,p5))/(sb*zaPb(p1,p1)**2*zaPQa(p1,p1)))-
     & (sa**2*za(p3,Qtil)*zaPb(Qtil,p2)*zaPQa(p6,Qtil)*
     &  zb(Qtil,p5))/
     &  (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2*zaPQa(p1,Qtil))-
     & (sa**2*za(p3,Qtilsl)*zaPb(Qtilsl,p2)*zaPQa(p6,Qtilsl)*
     &  zb(Qtilsl,p5))/
     &  (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**2*
     &  zaPQa(p1,Qtilsl)) +
     & (sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p1,p2))/
     &  (sb*zaPb(p1,p1)**2*zbQPb(p1,p1)))-
     &za(p1,p7)*(-((sa**2*za(p3,p1)*za(p6,p1)*zaPb(p1,p7)*
     & zb(p1,p5))/(sb*zaPb(p1,p1)**2*zaPQa(p1,p1)))-
     & (sa**2*za(p3,Qtil)*zaPb(Qtil,p7)*zaPQa(p6,Qtil)*
     &  zb(Qtil,p5))/
     &  (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2*zaPQa(p1,Qtil))-
     & (sa**2*za(p3,Qtilsl)*zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*
     &  zb(Qtilsl,p5))/
     &  (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**2*
     &  zaPQa(p1,Qtilsl))-
     & (sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p7,p1))/
     &  (sb*zaPb(p1,p1)**2*zbQPb(p1,p1)))))/(za(p2,p7)*za(p1,p7))-
     &(za(p2,p1)*(-((sa**2*za(p6,p1)*zaPb(p3,p1)*zaPb(p1,p4)*
     & zaPb(p1,p5)*zb(p1,p2))/
     & (sb*zaPb(p1,p1)**2*zaPQPb(p1,p1)))-
     &(sa**2*za(p3,p1)*za(p6,p1)*zaPb(p1,p2)*zaPb(p1,p5)*
     & zb(p1,p4))/(sb*zaPb(p1,p1)**2*zaPQa(p1,p1))-
     &(sa*zaPb(p3,p4)*zaPb(p6,p4)*zb(p4,p2)*zb(p5,p4))/
     & (2._dp*sb*zaPb(Qtil,p4)*zb(p1,p4)*zb(Qtil,p4)) +
     &(sa**3*za(p3,Qtil)*zaPb(p1,p1)*zaPb(Qtil,p2)*
     & zaPQa(p6,Qtil)*zaQb(Qtil,p5)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)**2*zaPb(Qtil,Qtil)**3*zaPQa(p1,Qtil)) +
     & (sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*
     & zaPQPb(p6,Qtil)*zaQb(Qtil,p5)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil))-
     & (sa**3*za(p1,Qtil)*zaPb(Qtil,p2)*zaPQa(p6,Qtil)*
     & zaQb(Qtil,p5)*zaQPa(p3,p1)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**3*zaPQa(p1,Qtil)**2)-
     & (sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPQa(p6,Qtil)*
     & zaQb(Qtil,p5)*zb(Qtil,p2)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil))-
     & (sa**3*za(p3,Qtil)*za(p6,Qtil)*zaPb(Qtil,p2)*
     & zaPb(Qtil,p5)*zb(Qtil,p4)**2)/
     & (2._dp*sb*zaPb(Qtil,p4)*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**3) +
     & (sa**3*za(p3,Qtilsl)*zaPb(p1,p1)*zaPb(Qtilsl,p2)*
     & zaPQa(p6,Qtilsl)*zaQb(Qtilsl,p5)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)**2*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQa(p1,Qtilsl)) +
     & (sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPQPb(p6,Qtilsl)*zaQb(Qtilsl,p5)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl))-
     & (sa**3*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*zaPQa(p6,Qtilsl)*
     & zaQb(Qtilsl,p5)*zaQPa(p3,p1)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQa(p1,Qtilsl)**2)-
     & (sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPQa(p6,Qtilsl)*
     & zaQb(Qtilsl,p5)*zb(Qtilsl,p2)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl))-
     & (sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p5)*zb(Qtilsl,p4)**2)/
     & (2._dp*sb*zaPb(Qtilsl,p4)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Qtilsl)**3)-
     & (sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*
     & zaPQa(p6,Qtil)*zb(Qtil,p4)*zbPQb(Qtil,p5))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil))-
     & (sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPQa(p6,Qtilsl)*zb(Qtilsl,p4)*zbPQb(Qtilsl,p5))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl))-
     & (sa*zaPb(p3,p1)*zaPb(p1,p4)*zaPQPb(p6,p1)*zb(p1,p2)*
     & zbPQb(p1,p5))/(sb*zaPQPb(p1,p1)*zbQPb(p1,p1)**2) +
     & (sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p4)**2*zb(p1,p2)*
     & zb(p1,p5))/(2._dp*sb*zaPb(p1,p1)**2*zb(p1,p4)*zbQPb(p1,p1)))
     &-za(p1,p7)*(-((sa**2*za(p3,p1)*za(p6,p1)*zaPb(p1,p5)*zaPb(p1,p7)*
     & zb(p1,p4))/(sb*zaPb(p1,p1)**2*zaPQa(p1,p1))) +
     &(sa**2*za(p6,p1)*zaPb(p3,p1)*zaPb(p1,p4)*zaPb(p1,p5)*
     & zb(p7,p1))/(sb*zaPb(p1,p1)**2*zaPQPb(p1,p1)) +
     &(sa*zaPb(p3,p4)*zaPb(p6,p4)*zb(p5,p4)*zb(p7,p4))/
     & (2._dp*sb*zaPb(Qtil,p4)*zb(p1,p4)*zb(Qtil,p4)) +
     &(sa**3*za(p3,Qtil)*zaPb(p1,p1)*zaPb(Qtil,p7)*
     & zaPQa(p6,Qtil)*zaQb(Qtil,p5)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)**2*zaPb(Qtil,Qtil)**3*zaPQa(p1,Qtil)) +
     &(sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p7)*
     & zaPQPb(p6,Qtil)*zaQb(Qtil,p5)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil))-
     &(sa**3*za(p1,Qtil)*zaPb(Qtil,p7)*zaPQa(p6,Qtil)*
     & zaQb(Qtil,p5)*zaQPa(p3,p1)*zb(Qtil,p4))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**3*zaPQa(p1,Qtil)**2)-
     &(sa**3*za(p3,Qtil)*za(p6,Qtil)*zaPb(Qtil,p5)*zaPb(Qtil,p7)*
     & zb(Qtil,p4)**2)/
     & (2._dp*sb*zaPb(Qtil,p4)*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**3)-
     &(sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPQa(p6,Qtil)*
     & zaQb(Qtil,p5)*zb(Qtil,p4)*zb(Qtil,p7))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil)) +
     &(sa**3*za(p3,Qtilsl)*zaPb(p1,p1)*zaPb(Qtilsl,p7)*
     & zaPQa(p6,Qtilsl)*zaQb(Qtilsl,p5)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)**2*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQa(p1,Qtilsl)) +
     &(sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p7)*
     & zaPQPb(p6,Qtilsl)*zaQb(Qtilsl,p5)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl))-
     &(sa**3*za(p1,Qtilsl)*zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*
     & zaQb(Qtilsl,p5)*zaQPa(p3,p1)*zb(Qtilsl,p4))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**3*
     & zaPQa(p1,Qtilsl)**2)-
     &(sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*zaPb(Qtilsl,p5)*
     & zaPb(Qtilsl,p7)*zb(Qtilsl,p4)**2)/
     & (2._dp*sb*zaPb(Qtilsl,p4)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Qtilsl)**3)-
     &(sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPQa(p6,Qtilsl)*
     & zaQb(Qtilsl,p5)*zb(Qtilsl,p4)*zb(Qtilsl,p7))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl))-
     &(sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p7)*
     & zaPQa(p6,Qtil)*zb(Qtil,p4)*zbPQb(Qtil,p5))/
     & (sb*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**4*zaPQa(p1,Qtil))-
     &(sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p7)*
     & zaPQa(p6,Qtilsl)*zb(Qtilsl,p4)*zbPQb(Qtilsl,p5))/
     & (sb*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**4*
     & zaPQa(p1,Qtilsl)) +
     &(sa*zaPb(p3,p1)*zaPb(p1,p4)*zaPQPb(p6,p1)*zb(p7,p1)*
     & zbPQb(p1,p5))/(sb*zaPQPb(p1,p1)*zbQPb(p1,p1)**2)-
     &(sa*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p4)**2*zb(p1,p5)*
     &zb(p7,p1))/(2._dp*sb*zaPb(p1,p1)**2*zb(p1,p4)*zbQPb(p1,p1))))/
     &(za(p2,p7)*za(p1,p7)) + (za(p2,p1)*zb(p7,p2)*
     & ((sa**3*za(p3,Q3til)*za(p6,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zb(Q3til,p5))/
     & (sb*sb3*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**2*
     & zaPQa(Q3til,Q3til)) +
     &  (sa**3*za(p3,Q3tilsl)*za(p6,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p4)*zb(Q3tilsl,p5))/
     & (sb*sb3*zaPb(Q3tilsl,p1)*zaPb(Q3tilsl,Q3tilsl)**2*
     & zaPQa(Q3tilsl,Q3tilsl)) +
     &  (sa**3*za(p3,Qtil)*za(p6,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zb(Qtil,p5))/
     & (sb*sb3*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2*
     & zaPQ3a(Qtil,Qtil)) +
     &  (sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zb(Qtilsl,p5))/
     & (sb*sb3*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ3a(Qtilsl,Qtilsl))-
     &  (sa**2*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p1,p2)*
     & zb(p1,p4))/
     & (sb*sb3*zaPb(p1,p1)*zbQ3Pb(p1,p1)*zbQPb(p1,p1))))/
     &za(p1,p7)- (za(p3,p1)*zb(p3,p4)*
     & (-((sa**3*za(p3,Q3til)*za(p6,Q3til)*zaPb(Q3til,p2)*
     &   zaPb(Q3til,p7)*zb(Q3til,p5))/
     & (sb*sb3*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**2*
     &   zaPQa(Q3til,Q3til)))-
     &  (sa**3*za(p3,Q3tilsl)*za(p6,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p7)*zb(Q3tilsl,p5))/
     & (sb*sb3*zaPb(Q3tilsl,p1)*zaPb(Q3tilsl,Q3tilsl)**2*
     & zaPQa(Q3tilsl,Q3tilsl))-
     &  (sa**3*za(p3,Qtil)*za(p6,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p7)*
     & zb(Qtil,p5))/
     & (sb*sb3*zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2*
     & zaPQ3a(Qtil,Qtil))-
     &  (sa**3*za(p3,Qtilsl)*za(p6,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p7)*zb(Qtilsl,p5))/
     & (sb*sb3*zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**2*
     & zaPQ3a(Qtilsl,Qtilsl))-
     &  (sa**2*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p1,p2)*
     & zb(p7,p1))/
     & (sb*sb3*zaPb(p1,p1)*zbQ3Pb(p1,p1)*zbQPb(p1,p1))))/
     &za(p1,p7)
      res = res + ((sa**3*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p7)*zaPb(Qtil,p4)*zaPQa(p6,Qtil)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)**2*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**2*zaPb(Qtil,Qtil)**2)-
     & (sa**3*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*zaPb(Q3til,p7)*
     & zaPb(Qtil,p4)*zaPQPb(p6,Qtil)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**2*zaPb(Q3til,Qtil)*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p6,Q3til)*za(p1,Q3til)*zaPb(p3,Qtil)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p7)*zaQ3b(Q3til,p4)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**3*zaPb(Q3til,Qtil)**2)-
     &(sa**4*za(p3,Q3til)*za(p1,Q3til)*zaPb(p6,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p7)*zaQ3b(Q3til,p4)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**4*zaPb(Q3til,Qtil))-
     &(sa**4*za(p3,Q3til)*za(p6,Q3til)*zaPb(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p7)*zaQ3b(Q3til,p4)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**4*zaPb(Q3til,Qtil))-
     &(sa**4*za(p3,Q3til)*za(p6,Q3til)*za(p1,Q3til)*zaPb(Q3til,p7)*
     & zaPb(Qtil,p2)*zaQ3b(Q3til,p4)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)**2*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**3*zaPb(Q3til,Qtil))-
     &(sa**5*za(p3,Q3til)*za(p6,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaQ3b(Q3til,p4)*zb(p7,p1)*zb(Q3til,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)**2*
     & zaPb(Q3til,Q3til)**3*zaPb(Q3til,Qtil)) +
     &(sa**4*za(p3,Q3til)*za(p6,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p4)*zaPb(Q3til,p7)*zb(Q3til,p5)**2)/
     & (2._dp*sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p5)*
     & zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**3*zaPb(Q3til,Qtil)) +
     &(sa**3*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p7)*zaPb(Qtil,p4)*zaPQa(p6,Qtil)*zb(Q3tilsl,p5))/
     & (sb*sb3*za(Q3tilsl,Qtil)**2*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**2*zaPb(Qtil,Qtil)**2)-
     &(sa**3*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p7)*zaPb(Qtil,p4)*zaPQPb(p6,Qtil)*zb(Q3tilsl,p5))/
     & (sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**2*zaPb(Q3tilsl,Qtil)*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p6,Q3tilsl)*za(p1,Q3tilsl)*zaPb(p3,Qtil)*
     & zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p7)
     &*zaQ3b(Q3tilsl,p4)*zb(Q3tilsl,p5)
     &)/(sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**3*zaPb(Q3tilsl,Qtil)**2)-
     &(sa**4*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(p6,Q3tilsl)*
     & zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p7)
     &*zaQ3b(Q3tilsl,p4)*zb(Q3tilsl,p5)
     &)/(sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**4*zaPb(Q3tilsl,Qtil))-
     &(sa**4*za(p3,Q3tilsl)*za(p6,Q3tilsl)*zaPb(p1,Q3tilsl)*
     & zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p7)
     &*zaQ3b(Q3tilsl,p4)*zb(Q3tilsl,p5)
     &)/(sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**4*zaPb(Q3tilsl,Qtil))-
     &(sa**4*za(p3,Q3tilsl)*za(p6,Q3tilsl)*za(p1,Q3tilsl)*
     & zaPb(Q3tilsl,p7)*zaPb(Qtil,p2)
     &*zaQ3b(Q3tilsl,p4)*zb(Q3tilsl,p5))/
     & (sb*sb3*za(Q3tilsl,Qtil)**2*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**3*zaPb(Q3tilsl,Qtil))-
     &(sa**5*za(p3,Q3tilsl)*za(p6,Q3tilsl)
     &*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaQ3b(Q3tilsl,p4)*zb(p7,p1)*zb(Q3tilsl,p5))/
     & (sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)**2*
     & zaPb(Q3tilsl,Q3tilsl)**3*zaPb(Q3tilsl,Qtil)) +
     &(sa**4*za(p3,Q3tilsl)*za(p6,Q3tilsl)*za(p1,Q3tilsl)*
     & zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p4)*zaPb(Q3tilsl,p7)*
     & zb(Q3tilsl,p5)**2)/
     & (2._dp*sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p5)*
     & zaPb(Q3tilsl,p1)*zaPb(Q3tilsl,Q3tilsl)**3*zaPb(Q3tilsl,Qtil)) +
     &(sa**4*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*zaPb(Q3til,p7)*
     & zaPQa(p6,Qtil)*zb(Q3til,p5)*zb(Qtil,p4))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**2*zaPb(Q3til,Qtil)*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p3,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*zaPb(Q3til,p7)*
     & zaPQPb(p6,Qtil)*zb(Q3til,p5)*zb(Qtil,p4))/
     & (sb*sb3*zaPb(Q3til,p1)*zaPb(Q3til,Q3til)**2*
     & zaPb(Q3til,Qtil)**2*zaPb(Qtil,Qtil)**2) +
     &(sa**4*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p7)*zaPQa(p6,Qtil)*zb(Q3tilsl,p5)*zb(Qtil,p4))/
     & (sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**2*zaPb(Q3tilsl,Qtil)*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p3,Q3tilsl)*za(p1,Q3tilsl)*zaPb(Q3tilsl,p2)*
     & zaPb(Q3tilsl,p7)*zaPQPb(p6,Qtil)*zb(Q3tilsl,p5)*zb(Qtil,p4))/
     & (sb*sb3*zaPb(Q3tilsl,p1)*zaPb(Q3tilsl,Q3tilsl)**2*
     & zaPb(Q3tilsl,Qtil)**2*zaPb(Qtil,Qtil)**2) +
     &(sa**2*zaPb(p3,p5)*zaPb(p6,p5)*zaPb(p1,p5)*zb(p5,p2)*zb(p5,p4)*
     & zb(p7,p5))/
     & (2._dp*sb*sb3*zaPb(Q3til,p5)*zaPb(Qtil,p5)*zb(p1,p5)*
     & zb(Q3til,p5)*zb(Qtil,p5)) +
     &(sa**3*za(p3,Qtil)*zaPb(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zaPb(Qtil,p7)*zaPQa(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p1)*zaPb(Qtil,Q3til)*
     & zaPb(Qtil,Qtil)**4) + (sa**3*za(p1,Qtil)*zaPb(p3,Q3til)*
     & zaPb(Qtil,p2)*zaPb(Qtil,p4)*zaPb(Qtil,p7)*zaPQa(p6,Qtil)*
     & zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p1)*
     & zaPb(Qtil,Q3til)**2*zaPb(Qtil,Qtil)**3)-
     &(sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Q3til,p2)*zaPb(Qtil,p4)*
     & zaPb(Qtil,p7)*zaPQa(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)**2*zaPb(Qtil,p1)*
     & zaPb(Qtil,Q3til)*zaPb(Qtil,Qtil)**3) +
     &(sa**3*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zaPb(Qtil,p7)*zaPQPb(p6,Qtil)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p1)*zaPb(Qtil,Q3til)*
     & zaPb(Qtil,Qtil)**4) + (sa**4*za(p3,Qtil)*za(p6,Q3til)*
     & za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p7)*zaQ3b(Q3til,p4)*
     & zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)**2*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtil,p1)*zaPb(Qtil,Qtil)**2) +
     &(sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPb(p6,Q3til)*zaPb(Qtil,p2)*
     & zaPb(Qtil,p7)*zaQ3b(Q3til,p4)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtil,p1)*zaPb(Qtil,Q3til)*zaPb(Qtil,Qtil)**2) +
     &(sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*
     & zaPQa(p6,Qtil)*zb(p7,p1)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p1)**2*
     & zaPb(Qtil,Q3til)*zaPb(Qtil,Qtil)**3)-
     &(sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p7)*
     & zaPQa(p6,Qtil)*zb(Qtil,p4)*zb(Qtil,p5))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p1)*zaPb(Qtil,Q3til)*
     & zaPb(Qtil,Qtil)**4)- (sa**4*za(p3,Qtil)*za(p6,Qtil)*
     & za(p1,Qtil)*zaPb(Qtil,p2)*zaPb(Qtil,p4)*zaPb(Qtil,p7)*
     & zb(Qtil,p5)**2)/
     & (2._dp*sb*sb3*za(Q3til,Qtil)*zaPb(Qtil,p5)*zaPb(Qtil,p1)*
     & zaPb(Qtil,Q3til)*zaPb(Qtil,Qtil)**3) +
     &(sa**3*za(p3,Qtilsl)*zaPb(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**4) +
     &(sa**3*za(p1,Qtilsl)*zaPb(p3,Q3til)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)**2*zaPb(Qtilsl,Qtilsl)**3)-
     &(sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Q3til,p2)*
     & zaPb(Qtilsl,p4)*zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)**2*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**3) +
     &(sa**3*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPb(Qtilsl,p7)*zaPQPb(p6,Qtilsl)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**4) +
     &(sa**4*za(p3,Qtilsl)*za(p6,Q3til)*za(p1,Qtilsl)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p7)*zaQ3b(Q3til,p4)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)**2*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtilsl,p1)*zaPb(Qtilsl,Qtilsl)**2) +
     &(sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(p6,Q3til)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p7)*zaQ3b(Q3til,p4)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtilsl,p1)*zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**2) +
     &(sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p4)*zaPQa(p6,Qtilsl)*zb(p7,p1)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p1)**2*
     & zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**3)-
     &(sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(Qtilsl,p2)*
     & zaPb(Qtilsl,p7)*zaPQa(p6,Qtilsl)*zb(Qtilsl,p4)*zb(Qtilsl,p5))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**4)-
     &(sa**4*za(p3,Qtilsl)*za(p6,Qtilsl)*za(p1,Qtilsl)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p4)
     &*zaPb(Qtilsl,p7)*zb(Qtilsl,p5)**2)/
     & (2._dp*sb*sb3*za(Q3til,Qtilsl)*zaPb(Qtilsl,p5)*
     & zaPb(Qtilsl,p1)*zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**3) +
     &(sa**4*za(p3,Q3til)*za(p6,Q3til)*za(p1,Q3til)*zaPb(Q3til,p2)*
     & zaPb(Q3til,p7)*zb(Q3til,p5)*zbPQ3b(Q3til,p4))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,p1)*
     & zaPb(Q3til,Q3til)**4*zaPb(Q3til,Qtil))-
     &(sa**4*za(p3,Qtil)*za(p1,Qtil)*zaPb(p6,Q3til)*zaPb(Qtil,p2)*
     & zaPb(Qtil,p7)*zb(Qtil,p5)*zbPQ3b(Q3til,p4))/
     & (sb*sb3*zaPb(Q3til,Q3til)**2*zaPb(Qtil,p1)*
     & zaPb(Qtil,Q3til)**2*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p3,Qtil)*za(p6,Q3til)*za(p1,Qtil)*zaPb(Qtil,p2)*
     & zaPb(Qtil,p7)*zb(Qtil,p5)*zbPQ3b(Q3til,p4))/
     & (sb*sb3*za(Q3til,Qtil)*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtil,p1)*zaPb(Qtil,Q3til)*zaPb(Qtil,Qtil)**2)-
     &(sa**4*za(p3,Qtilsl)*za(p1,Qtilsl)*zaPb(p6,Q3til)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p7)*zb(Qtilsl,p5)*zbPQ3b(Q3til,p4))/
     & (sb*sb3*zaPb(Q3til,Q3til)**2*zaPb(Qtilsl,p1)*
     & zaPb(Qtilsl,Q3til)**2*zaPb(Qtilsl,Qtilsl)**2)-
     &(sa**4*za(p3,Qtilsl)*za(p6,Q3til)*za(p1,Qtilsl)*
     & zaPb(Qtilsl,p2)*zaPb(Qtilsl,p7)*zb(Qtilsl,p5)*zbPQ3b(Q3til,p4))/
     & (sb*sb3*za(Q3til,Qtilsl)*zaPb(Q3til,Q3til)**2*
     & zaPb(Qtilsl,p1)*zaPb(Qtilsl,Q3til)*zaPb(Qtilsl,Qtilsl)**2) +
     &(sa**4*za(p3,Q3tilsl)*za(p6,Q3tilsl)*za(p1,Q3tilsl)*
     & zaPb(Q3tilsl,p2)*zaPb(Q3tilsl,p7)*zb(Q3tilsl,p5)*
     & zbPQ3b(Q3tilsl,p4))/
     & (sb*sb3*za(Q3tilsl,Qtil)*zaPb(Q3tilsl,p1)*
     & zaPb(Q3tilsl,Q3tilsl)**4*zaPb(Q3tilsl,Qtil)) +
     &(sa**2*zaPb(p3,p1)*zaPb(p1,p5)*zaPQPb(p6,p1)*zb(p1,p2)*
     & zb(p1,p4)*zb(p7,p1))/
     & (sb*sb3*zbQ3Pb(p1,p1)*zbQPb(p1,p1)**2) +
     &(sa**2*zaPb(p3,p1)*zaPb(p6,p1)*zaPb(p1,p5)*zb(p1,p2)*zb(p7,p1)*
     & zbPQ3b(p1,p4))/(sb*sb3*zbQ3Pb(p1,p1)**2*zbQPb(p1,p1))
     &- (sa**2*zaPb(p3,p1)*zaPb(p6,p1)
     &*zaPb(p1,p5)**2*zb(p1,p2)*zb(p1,p4)*
     &  zb(p7,p1))/
     & (2._dp*sb*sb3*zaPb(p1,p1)*zb(p1,p5)*zbQ3Pb(p1,p1)*
     &  zbQPb(p1,p1)))/za(p1,p7)

c--- multiply by overall factor (including W propagators)
      res=res/(sb*sa)

c--- no overall sign
      bub56sl=res

      return
      end

