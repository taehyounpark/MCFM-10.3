!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function BDK1211mp(k1,k2,k3,k4,k5,k6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: BDK1211mp
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6,j2,j1
      complex(dp):: L0,L1,Lnrat,i3m,Lsm1_2mh,zab2,zba2
      complex(dp):: T1256,T6521,box1256,box6521,
     & squarebracketex,squarebracketfl3,
     & squarebracketfl3ex,squarebracket
      real(dp):: t,DELTA3,delta
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      zba2(k1,k2,k3,k4)=za(k4,k2)*zb(k2,k1)+za(k4,k3)*zb(k3,k1)
      DELTA3(k1,k2,k3,k4,k5,k6)=
     & +s(k1,k2)**2+s(k3,k4)**2+s(k5,k6)**2-2d0*s(k1,k2)*s(k3,k4)
     & -2d0*s(k3,k4)*s(k5,k6)-2d0*s(k5,k6)*s(k1,k2)
      delta(k1,k2,k3,k4,k5,k6)=s(k1,k2)-s(k3,k4)-s(k5,k6)

      j2=k2
      j1=k1

      T1256=
     & +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
     & -DELTA3(k1,k2,k3,k4,k5,k6))
     & *(t(k1,k2,k3)-t(k1,k2,k4))*zab2(j2,k3,k4,j1)*zab2(k5,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

     & +0.5d0*s(k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))*za(j2,k5)*zb(j1,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & -zb(j1,k2)*za(k5,k6)*zab2(j2,k3,k4,k6)*zab2(k2,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(j2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)
     & *(za(k5,k2)*zb(k2,j1)*delta(k1,k2,k3,k4,k5,k6)
     & -za(k5,k6)*zb(k6,j1)*delta(k5,k6,k3,k4,k1,k2))
     & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(j2,k3)*zb(k4,k6)*zab2(k5,k2,k3,j1)/zab2(k3,k1,k2,k4)**2

     & -2d0*zb(j1,k3)*za(k4,k5)*zab2(j2,k1,k3,k4)*zab2(k3,k1,k2,k6)
     & /t(k1,k2,k3)**2/zab2(k3,k1,k2,k4)

      T6521=
     & +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
     & -DELTA3(k1,k2,k3,k4,k5,k6))
     & *(t(k1,k2,k4)-t(k1,k2,k3))*zab2(k5,k3,k4,k6)*zab2(j2,k3,k4,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

     & +0.5d0*s(k3,k4)*(t(k1,k2,k4)-t(k1,k2,k3))*za(k5,j2)*zb(k6,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & -zb(k6,k5)*za(j2,k1)*zab2(k5,k3,k4,k1)*zab2(k5,k3,k4,j1)
     & /zab2(k3,k5,k6,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(k5,k3)*zb(k4,j1)*delta(k3,k4,k1,k2,k5,k6)
     & *(za(j2,k5)*zb(k5,k6)*delta(k5,k6,k3,k4,k1,k2)
     & -za(j2,k1)*zb(k1,k6)*delta(k1,k2,k3,k4,k5,k6))
     & /zab2(k3,k5,k6,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

     & +za(k5,k3)*zb(k4,j1)*zab2(j2,k5,k3,k6)/zab2(k3,k5,k6,k4)**2

     & -2d0*zb(k6,k3)*za(k4,j2)*zab2(k5,k6,k3,k4)*zab2(k3,k6,k5,j1)
     & /t(k6,k5,k3)**2/zab2(k3,k5,k6,k4)

c      write(6,*) 'T1256',T1256
c      write(6,*) 'T6521',T6521
c      write(6,*) 'T1256+T6521',T1256+T6521
c      pause

c     T1256i=
c    &   +0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
c    & -DELTA3(k1,k2,k3,k4,k5,k6))
c    & *(t(k1,k2,k3)-t(k1,k2,k4))*zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

c    & +0.5d0*s(k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))*za(k2,k5)*zb(k1,k6)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

c    & -zb(k1,k2)*zab2(k2,k3,k4,k6)*(
c    & -zab2(k2,k3,k4,k1)*za(k5,k1)-zab2(k2,k3,k4,k2)*za(k5,k2)
c    & -zab2(k2,k3,k4,k3)*za(k5,k3)-zab2(k2,k3,k4,k4)*za(k5,k4))
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

c    & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)
c    & *(za(k5,k2)*zb(k2,k1)*delta(k1,k2,k3,k4,k5,k6)
c    & +(+za(k5,k2)*zb(k2,k1)+za(k5,k3)*zb(k3,k1)+za(k5,k4)*zb(k4,k1))
c    & *delta(k5,k6,k3,k4,k1,k2))
c    & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

c    & +za(k2,k3)*zb(k4,k6)*zab2(k5,k2,k3,k1)/zab2(k3,k1,k2,k4)**2

c    & -2d0*zb(k1,k3)*za(k4,k5)*zab2(k2,k1,k3,k4)*zab2(k3,k1,k2,k6)
c    & /t(k1,k2,k3)**2/zab2(k3,k1,k2,k4)

c     T6521i=
c    & -0.5d0*(3d0*s(k3,k4)*delta(k3,k4,k1,k2,k5,k6)
c    & -DELTA3(k1,k2,k3,k4,k5,k6))
c    & *(t(k1,k2,k4)-t(k1,k2,k3))*zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2

c    & -0.5d0*s(k3,k4)*(t(k1,k2,k4)-t(k1,k2,k3))*za(k5,k2)*zb(k6,k1)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

c    & -za(k2,k1)*(zb(k6,k1)*zab2(k1,k3,k4,k1)
c    & +zb(k6,k2)*zab2(k2,k3,k4,k1)
c    & +zb(k6,k3)*zab2(k3,k3,k4,k1)+zb(k6,k4)*zab2(k4,k3,k4,k1))
c    & *zab2(k5,k3,k4,k1)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

c    & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k1,k2,k5,k6)
c    & *(-(za(k2,k1)*zb(k1,k6)+za(k2,k3)*zb(k3,k6)+za(k2,k4)*zb(k4,k6))
c    & *delta(k5,k6,k3,k4,k1,k2)
c    & -za(k2,k1)*zb(k1,k6)*delta(k1,k2,k3,k4,k5,k6))
c    & /zab2(k3,k1,k2,k4)**2/DELTA3(k1,k2,k3,k4,k5,k6)

c    & -za(k5,k3)*zb(k4,k1)*zab2(k2,k1,k4,k6)/zab2(k3,k1,k2,k4)**2

c    & -2d0*zb(k6,k3)*za(k4,k2)*(-zab2(k5,k1,k2,k4))*zab2(k3,k2,k4,k1)
c    & /t(k6,k5,k3)**2/zab2(k3,k1,k2,k4)

c      write(6,*) 'T1256',T1256
c      write(6,*) 'T6521',T6521
      squarebracket=
     & 0.5d0*zb(k6,k4)**2*za(k4,k2)**2
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)
     & *L1(-s(k5,k6),-t(k1,k2,k3))/t(k1,k2,k3)

     & +2d0*zb(k6,k4)*za(k4,k2)*zab2(k2,k1,k3,k6)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)
     & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

     & -za(k2,k3)*za(k2,k4)*zb(k6,k4)**2*t(k1,k2,k3)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)**2
     & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

     & -0.5d0*za(k2,k3)*zb(k6,k4)*zab2(k2,k1,k3,k6)
     & /za(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)**2
     & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

     & -0.75d0*zab2(k2,k1,k3,k6)**2
     & /za(k1,k2)/zb(k5,k6)/t(k1,k2,k3)/zab2(k3,k1,k2,k4)
     & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

     & +(1.5d0*delta(k5,k6,k1,k2,k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))
     & *zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2
     & -zb(k1,k2)*za(k2,k3)*zb(k4,k6)/zab2(k3,k1,k2,k4)**2
     & /DELTA3(k1,k2,k3,k4,k5,k6)
     & *(za(k2,k5)*(t(k1,k2,k3)-t(k1,k2,k4))
     & -2d0*za(k2,k1)*zb(k1,k6)*za(k6,k5))
     & +zb(k1,k2)*za(k2,k5)*(za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
     & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & +zb(k1,k6)/zb(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & *(za(k2,k3)*zb(k3,k6)*t(k1,k2,k3)
     & -za(k2,k4)*zb(k4,k6)*t(k1,k2,k4)
     & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)*t(k1,k2,k3)
     & /zab2(k3,k1,k2,k4))
     & )*Lnrat(-s(k1,k2),-s(k3,k4))

     & -0.25d0*zb(k1,k6)*(t(k1,k2,k3)-t(k1,k2,k4))
     & *(zb(k1,k6)*delta(k3,k4,k1,k2,k5,k6)
     & -2d0*zb(k1,k2)*za(k2,k5)*zb(k5,k6))
     & /zb(k1,k2)/zb(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
     & ;

c     squarebracketi=
c    & 0.5d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
c    & -za(k3,k5)*zb(k3,k4))*zb(k6,k4)*za(k4,k2)**2
c    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)
c    & *L1(-s(k5,k6),-t(k1,k2,k3))/t(k1,k2,k3)

c    & +2d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
c    & -za(k3,k5)*zb(k3,k4))*za(k4,k2)*zab2(k2,k1,k3,k6)
c    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)
c    & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

c    & -za(k2,k3)*za(k2,k4)*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
c    & -za(k3,k5)*zb(k3,k4))
c    & *zb(k6,k4)*t(k1,k2,k3)
c    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)**2
c    & *L0(-t(k1,k2,k3),-s(k5,k6))/s(k5,k6)

c    & -0.5d0*za(k2,k3)*(-zb(k1,k4)*za(k1,k5)-zb(k2,k4)*za(k2,k5)
c    * -zb(k3,k4)*za(k3,k5))*zab2(k2,k1,k3,k6)
c    & /za(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)**2
c    & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

c    & -0.75d0*zab2(k2,k1,k3,k6)
c    & *(-zab2(k2,k1,k3,k1)*za(k1,k5)-zab2(k2,k1,k3,k2)*za(k2,k5)
c    & -zab2(k2,k1,k3,k3)*za(k3,k5)-zab2(k2,k1,k3,k4)*za(k4,k5))
c    & /za(k1,k2)/s(k5,k6)/t(k1,k2,k3)/zab2(k3,k1,k2,k4)
c    & *(Lnrat(-s(k3,k4),-s(k5,k6))+Lnrat(-t(k1,k2,k3),-s(k5,k6)))

c    & +(1.5d0*delta(k5,k6,k1,k2,k3,k4)*(t(k1,k2,k3)-t(k1,k2,k4))
c    & *zab2(k2,k3,k4,k1)*zab2(k5,k3,k4,k6)
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)**2
c    & -zb(k1,k2)*za(k2,k3)*zb(k4,k6)/zab2(k3,k1,k2,k4)**2
c    & /DELTA3(k1,k2,k3,k4,k5,k6)
c    & *(za(k2,k5)*(t(k1,k2,k3)-t(k1,k2,k4))-2d0*za(k2,k1)
c    & *(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5)))
c    & +zb(k1,k2)*za(k2,k5)*(za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
c    & /zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
c    & +(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5))
c    & /s(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)
c    & *(za(k2,k3)*zb(k3,k6)*t(k1,k2,k3)
c    & -za(k2,k4)*zb(k4,k6)*t(k1,k2,k4)
c    & +za(k2,k3)*zb(k4,k6)*delta(k3,k4,k1,k2,k5,k6)*t(k1,k2,k3)
c    & /zab2(k3,k1,k2,k4))
c    & )*Lnrat(-s(k1,k2),-s(k3,k4))

c    & -0.25d0*(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)
c    & -zb(k1,k4)*za(k4,k5))
c    & *(t(k1,k2,k3)-t(k1,k2,k4))
c    & *(zb(k1,k6)*delta(k3,k4,k1,k2,k5,k6)-2d0*zb(k1,k2)
c    & *(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6)))
c    & /zb(k1,k2)/s(k5,k6)/zab2(k3,k1,k2,k4)/DELTA3(k1,k2,k3,k4,k5,k6)

      squarebracketex=
     & 0.5d0*zb(k1,k4)**2*za(k4,k5)**2
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)
     & *L1(-s(k2,k1),-t(k6,k5,k3))/t(k6,k5,k3)

     & +2d0*zb(k1,k4)*za(k4,k5)*zab2(k5,k6,k3,k1)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)
     & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

     & -za(k5,k3)*za(k5,k4)*zb(k1,k4)**2*t(k6,k5,k3)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)**2
     & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

     & -0.5d0*za(k5,k3)*zb(k1,k4)*zab2(k5,k6,k3,k1)
     & /za(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)**2
     & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

     & -0.75d0*zab2(k5,k6,k3,k1)**2
     & /za(k6,k5)/zb(k2,k1)/t(k6,k5,k3)/zab2(k3,k6,k5,k4)
     & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

     & +(1.5d0*delta(k2,k1,k6,k5,k3,k4)*(t(k6,k5,k3)-t(k6,k5,k4))
     & *zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
     & /zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)**2
     & -zb(k6,k5)*za(k5,k3)*zb(k4,k1)/zab2(k3,k6,k5,k4)**2
     & /DELTA3(k6,k5,k3,k4,k2,k1)
     & *(za(k5,k2)*(t(k6,k5,k3)-t(k6,k5,k4))
     & -2d0*za(k5,k6)*zb(k6,k1)*za(k1,k2))
     & +zb(k6,k5)*za(k5,k2)*(za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
     & /zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)
     & +zb(k6,k1)/zb(k2,k1)/zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)
     & *(za(k5,k3)*zb(k3,k1)*t(k6,k5,k3)
     & -za(k5,k4)*zb(k4,k1)*t(k6,k5,k4)
     & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k6,k5,k2,k1)*t(k6,k5,k3)
     & /zab2(k3,k6,k5,k4))
     & )*Lnrat(-s(k6,k5),-s(k3,k4))

     & -0.25d0*zb(k6,k1)*(t(k6,k5,k3)-t(k6,k5,k4))
     & *(zb(k6,k1)*delta(k3,k4,k6,k5,k2,k1)
     & -2d0*zb(k6,k5)*za(k5,k2)*zb(k2,k1))
     & /zb(k6,k5)/zb(k2,k1)/zab2(k3,k6,k5,k4)/DELTA3(k6,k5,k3,k4,k2,k1)

c     squarebracketexi=
c    & 0.5d0*zb(k1,k4)**2*za(k4,k5)*(-za(k4,k1)*zb(k1,k6)
c    & -za(k4,k2)*zb(k2,k6)-za(k4,k3)*zb(k3,k6))
c    & /s(k5,k6)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
c    & *L1(-s(k2,k1),-t(k6,k5,k3))/t(k6,k5,k3)

c    & +2d0*zb(k1,k4)*(-za(k4,k1)*zb(k1,k6)-za(k4,k2)*zb(k2,k6)
c    & -za(k4,k3)*zb(k3,k6))*(-zab2(k5,k2,k4,k1))
c    & /s(k5,k6)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
c    & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

c    & -(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)-za(k4,k3)*zb(k4,k6))
c    & *za(k5,k4)*zb(k1,k4)**2*t(k6,k5,k3)
c    & /s(k5,k6)/zb(k2,k1)/zab2(k3,k1,k2,k4)**2
c    & *L0(-t(k6,k5,k3),-s(k2,k1))/s(k2,k1)

c    & -0.5d0*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
c    & -za(k4,k3)*zb(k4,k6))*zb(k1,k4)
c    & *(-zab2(k5,k2,k4,k1))
c    & /s(k5,k6)/zb(k2,k1)/zab2(k3,k1,k2,k4)**2
c    & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

c    & -0.75d0*zab2(k5,k2,k4,k1)
c    & *(-zab2(k1,k2,k4,k1)*zb(k1,k6)-zab2(k2,k2,k4,k1)*zb(k2,k6)
c    & -zab2(k3,k2,k4,k1)*zb(k3,k6)-zab2(k4,k2,k4,k1)*zb(k4,k6))
c    & /s(k5,k6)/zb(k2,k1)/t(k6,k5,k3)*(-1d0/zab2(k3,k1,k2,k4))
c    & *(Lnrat(-s(k3,k4),-s(k2,k1))+Lnrat(-t(k6,k5,k3),-s(k2,k1)))

c    & +(1.5d0*delta(k2,k1,k6,k5,k3,k4)*(t(k6,k5,k3)-t(k6,k5,k4))
c    & *zab2(k5,k3,k4,k6)*zab2(k2,k3,k4,k1)
c    & *(-1d0/zab2(k3,k1,k2,k4))/DELTA3(k6,k5,k3,k4,k2,k1)**2
c    & -(-zb(k6,k1)*za(k1,k3)-zb(k6,k2)*za(k2,k3)-zb(k6,k4)*za(k4,k3))
c    & *zb(k4,k1)/zab2(k3,k1,k2,k4)**2/DELTA3(k6,k5,k3,k4,k2,k1)
c    & *(za(k5,k2)*(t(k6,k5,k3)-t(k6,k5,k4))
c    & -2d0*(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)
c    & -za(k5,k4)*zb(k4,k1))*za(k1,k2))
c    & +(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))
c    & *(za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
c    & *(-1d0/zab2(k3,k1,k2,k4))/DELTA3(k6,k5,k3,k4,k2,k1)
c    & +zb(k6,k1)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
c    & /DELTA3(k6,k5,k3,k4,k2,k1)
c    & *(za(k5,k3)*zb(k3,k1)*t(k6,k5,k3)
c    & -za(k5,k4)*zb(k4,k1)*t(k6,k5,k4)
c    & +za(k5,k3)*zb(k4,k1)*delta(k3,k4,k6,k5,k2,k1)*t(k6,k5,k3)
c    & *(-1d0/zab2(k3,k1,k2,k4)))
c    & )*Lnrat(-s(k6,k5),-s(k3,k4))

c    & -0.25d0*(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)
c    & -za(k5,k4)*zb(k4,k1))
c    & *(t(k6,k5,k3)-t(k6,k5,k4))
c    & *(zb(k6,k1)*delta(k3,k4,k6,k5,k2,k1)-2
c    & *(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))
c    & *zb(k2,k1))
c    & /s(k6,k5)/zb(k2,k1)*(-1d0/zab2(k3,k1,k2,k4))
c    & /DELTA3(k6,k5,k3,k4,k2,k1)

c     squarebracketfl3i=
c    & 0.5d0*za(k2,k3)**2*zb(k3,k6)*(-zb(k3,k1)*za(k1,k5)
c    & -zb(k3,k2)*za(k2,k5)-zb(k3,k4)*za(k4,k5))
c    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
c    & *L1(-s(k1,k2),-t(k5,k6,k4))/t(k5,k6,k4)

c    & +2d0*za(k2,k3)*(-zb(k3,k1)*za(k1,k5)-zb(k3,k2)*za(k2,k5)
c    & -zb(k3,k4)*za(k4,k5))*(-zba2(k6,k1,k3,k2))
c    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
c    & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

c    & -(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)-za(k3,k5)*zb(k3,k4))
c    & *zb(k6,k3)*za(k2,k3)**2*t(k5,k6,k4)
c    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))**2
c    & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

c    & -0.5d0*(-za(k1,k5)*zb(k1,k4)-za(k2,k5)*zb(k2,k4)
c    & -za(k3,k5)*zb(k3,k4))*za(k2,k3)*(-zba2(k6,k1,k3,k2))
c    & /s(k5,k6)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))**2
c    & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

c    & -0.75d0*zba2(k6,k1,k3,k2)*(-zba2(k1,k1,k3,k2)*za(k1,k5)
c    & -zba2(k2,k1,k3,k2)*za(k2,k5)
c    & -zba2(k3,k1,k3,k2)*za(k3,k5)-zba2(k4,k1,k3,k2)*za(k4,k5))
c    & /s(k5,k6)/za(k1,k2)/t(k5,k6,k4)*(-1d0/zba2(k4,k1,k2,k3))
c    & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

c    & +(1.5d0*delta(k1,k2,k5,k6,k4,k3)*(t(k5,k6,k4)-t(k5,k6,k3))
c    & *zba2(k6,k4,k3,k5)*zba2(k1,k4,k3,k2)
c    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)**2
c    & -(-za(k5,k1)*zb(k1,k4)-za(k5,k2)*zb(k2,k4)-za(k5,k3)*zb(k3,k4))
c    & *za(k3,k2)
c    & *(-1d0/zba2(k4,k1,k2,k3))**2/DELTA3(k5,k6,k4,k3,k1,k2)
c    & *(zb(k6,k1)*(t(k5,k6,k4)-t(k5,k6,k3))-2d0*(-zb(k6,k1)*za(k1,k2)
c    & -zb(k6,k3)*za(k3,k2)-zb(k6,k4)*za(k4,k2))*zb(k2,k1))
c    & +(-za(k5,k2)*zb(k2,k1)-za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
c    & *(zb(k6,k4)*za(k4,k2)-zb(k6,k3)*za(k3,k2))
c    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)
c    & +za(k5,k2)/za(k1,k2)*(-1d0/zba2(k4,k1,k2,k3))
c    & /DELTA3(k5,k6,k4,k3,k1,k2)
c    & *(zb(k6,k4)*za(k4,k2)*t(k5,k6,k4)
c    & -zb(k6,k3)*za(k3,k2)*t(k5,k6,k3)
c    & +zb(k6,k4)*za(k3,k2)*delta(k4,k3,k5,k6,k1,k2)*t(k5,k6,k4)
c    & *(-1d0/zba2(k4,k1,k2,k3)))
c    & )*Lnrat(-s(k5,k6),-s(k4,k3))

c    & -0.25d0*(-zb(k6,k1)*za(k1,k2)-zb(k6,k3)*za(k3,k2)
c    & -zb(k6,k4)*za(k4,k2))*(t(k5,k6,k4)-t(k5,k6,k3))
c    & *(za(k5,k2)*delta(k4,k3,k5,k6,k1,k2)-2d0*(-za(k5,k2)*zb(k2,k1)
c    & -za(k5,k3)*zb(k3,k1)-za(k5,k4)*zb(k4,k1))
c    & *za(k1,k2))
c    & /s(k5,k6)/za(k1,k2)
c    & *(-1d0/zba2(k4,k1,k2,k3))/DELTA3(k5,k6,k4,k3,k1,k2)

      squarebracketfl3=
     & 0.5d0*za(k2,k3)**2*zb(k3,k6)**2
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)
     & *L1(-s(k1,k2),-t(k5,k6,k4))/t(k5,k6,k4)

     & +2d0*za(k2,k3)*zb(k3,k6)*zba2(k6,k5,k4,k2)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)
     & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

     & -zb(k6,k4)*zb(k6,k3)*za(k2,k3)**2*t(k5,k6,k4)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)**2
     & *L0(-t(k5,k6,k4),-s(k1,k2))/s(k1,k2)

     & -0.5d0*zb(k6,k4)*za(k2,k3)*zba2(k6,k5,k4,k2)
     & /zb(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)**2
     & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

     & -0.75d0*zba2(k6,k5,k4,k2)**2
     & /zb(k5,k6)/za(k1,k2)/t(k5,k6,k4)/zba2(k4,k5,k6,k3)
     & *(Lnrat(-s(k4,k3),-s(k1,k2))+Lnrat(-t(k5,k6,k4),-s(k1,k2)))

     & +(1.5d0*delta(k1,k2,k5,k6,k4,k3)*(t(k5,k6,k4)-t(k5,k6,k3))
     & *zba2(k6,k4,k3,k5)*zba2(k1,k4,k3,k2)
     & /zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)**2
     & -za(k5,k6)*zb(k6,k4)*za(k3,k2)/zba2(k4,k5,k6,k3)**2
     & /DELTA3(k5,k6,k4,k3,k1,k2)
     & *(zb(k6,k1)*(t(k5,k6,k4)-t(k5,k6,k3))
     & -2d0*zb(k6,k5)*za(k5,k2)*zb(k2,k1))
     & +za(k5,k6)*zb(k6,k1)*(zb(k6,k4)*za(k4,k2)-zb(k6,k3)*za(k3,k2))
     & /zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)
     & +za(k5,k2)/za(k1,k2)/zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)
     & *(zb(k6,k4)*za(k4,k2)*t(k5,k6,k4)
     & -zb(k6,k3)*za(k3,k2)*t(k5,k6,k3)
     & +zb(k6,k4)*za(k3,k2)*delta(k4,k3,k5,k6,k1,k2)*t(k5,k6,k4)
     & /zba2(k4,k5,k6,k3))
     & )*Lnrat(-s(k5,k6),-s(k4,k3))

     & -0.25d0*za(k5,k2)*(t(k5,k6,k4)-t(k5,k6,k3))
     & *(za(k5,k2)*delta(k4,k3,k5,k6,k1,k2)
     & -2d0*za(k5,k6)*zb(k6,k1)*za(k1,k2))
     & /za(k5,k6)/za(k1,k2)/zba2(k4,k5,k6,k3)/DELTA3(k5,k6,k4,k3,k1,k2)


      squarebracketfl3ex=
     & 0.5d0*za(k5,k3)**2*zb(k3,k1)**2/zb(k2,k1)/za(k6,k5)
     & /zba2(k4,k2,k1,k3)
     & *L1(-s(k6,k5),-t(k2,k1,k4))/t(k2,k1,k4)

     & +2d0*za(k5,k3)*zb(k3,k1)*zba2(k1,k2,k4,k5)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)
     & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

     & -zb(k1,k4)*zb(k1,k3)*za(k5,k3)**2*t(k2,k1,k4)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)**2
     & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

     & -0.5d0*zb(k1,k4)*za(k5,k3)*zba2(k1,k2,k4,k5)
     & /zb(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)**2
     & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

     & -0.75d0*zba2(k1,k2,k4,k5)**2
     & /zb(k2,k1)/za(k6,k5)/t(k2,k1,k4)/zba2(k4,k2,k1,k3)
     & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

     & +(1.5d0*delta(k6,k5,k2,k1,k4,k3)*(t(k2,k1,k4)-t(k2,k1,k3))
     & *zba2(k1,k4,k3,k2)*zba2(k6,k4,k3,k5)
     & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)**2
     & -za(k2,k1)*zb(k1,k4)*za(k3,k5)/zba2(k4,k2,k1,k3)**2
     & /DELTA3(k2,k1,k4,k3,k6,k5)
     & *(zb(k1,k6)*(t(k2,k1,k4)-t(k2,k1,k3))
     & -2d0*zb(k1,k2)*za(k2,k5)*zb(k5,k6))
     & +za(k2,k1)*zb(k1,k6)*(zb(k1,k4)*za(k4,k5)-zb(k1,k3)*za(k3,k5))
     & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
     & +za(k2,k5)/za(k6,k5)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
     & *(zb(k1,k4)*za(k4,k5)*t(k2,k1,k4)
     & -zb(k1,k3)*za(k3,k5)*t(k2,k1,k3)
     & +zb(k1,k4)*za(k3,k5)*delta(k4,k3,k2,k1,k6,k5)*t(k2,k1,k4)
     & /zba2(k4,k2,k1,k3))
     & )*Lnrat(-s(k2,k1),-s(k4,k3))

     & -0.25d0*za(k2,k5)*(t(k2,k1,k4)-t(k2,k1,k3))
     & *(za(k2,k5)*delta(k4,k3,k2,k1,k6,k5)
     & -2d0*za(k2,k1)*zb(k1,k6)*za(k6,k5))
     & /za(k2,k1)/za(k6,k5)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)


c      squarebracketfl3exi=
c    & 0.5d0*za(k5,k3)*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
c    & -za(k4,k3)*zb(k4,k6))*zb(k3,k1)**2
c    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)
c    & *L1(-s(k6,k5),-t(k2,k1,k4))/t(k2,k1,k4)

c    & +2d0*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
c    & -za(k4,k3)*zb(k4,k6))*zb(k3,k1)*zba2(k1,k2,k4,k5)
c    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)
c    & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

c    & -zb(k1,k4)*zb(k1,k3)*za(k5,k3)*(-za(k1,k3)*zb(k1,k6)
c    & -za(k2,k3)*zb(k2,k6)-za(k4,k3)*zb(k4,k6))*t(k2,k1,k4)
c    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)**2
c    & *L0(-t(k2,k1,k4),-s(k6,k5))/s(k6,k5)

c    & -0.5d0*zb(k1,k4)*(-za(k1,k3)*zb(k1,k6)-za(k2,k3)*zb(k2,k6)
c    & -za(k4,k3)*zb(k4,k6))*zba2(k1,k2,k4,k5)
c    & /zb(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)**2
c    & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

c    & -0.75d0*zba2(k1,k2,k4,k5)*(-zba2(k1,k2,k4,k1)*zb(k1,k6)
c    & -zba2(k1,k2,k4,k2)*zb(k2,k6)
c    & -zba2(k1,k2,k4,k3)*zb(k3,k6)-zba2(k1,k2,k4,k4)*zb(k4,k6))
c    & /zb(k2,k1)/s(k5,k6)/t(k2,k1,k4)/zba2(k4,k2,k1,k3)
c    & *(Lnrat(-s(k4,k3),-s(k6,k5))+Lnrat(-t(k2,k1,k4),-s(k6,k5)))

c    & +(1.5d0*delta(k6,k5,k2,k1,k4,k3)*(t(k2,k1,k4)-t(k2,k1,k3))
c    & *zba2(k1,k4,k3,k2)*zba2(k6,k4,k3,k5)
c    & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)**2
c    & -za(k2,k1)*zb(k1,k4)*za(k3,k5)/zba2(k4,k2,k1,k3)**2
c    & /DELTA3(k2,k1,k4,k3,k6,k5)
c    & *(zb(k1,k6)*(t(k2,k1,k4)-t(k2,k1,k3))-2d0*zb(k1,k2)
c    & *(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6)))
c    & +za(k2,k1)*zb(k1,k6)*(zb(k1,k4)*za(k4,k5)-zb(k1,k3)*za(k3,k5))
c    & /zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
c    & +(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)-za(k2,k4)*zb(k4,k6))
c    * /s(k5,k6)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)
c    & *(zb(k1,k4)*za(k4,k5)*t(k2,k1,k4)
c    & -zb(k1,k3)*za(k3,k5)*t(k2,k1,k3)
c    & +zb(k1,k4)*za(k3,k5)*delta(k4,k3,k2,k1,k6,k5)*t(k2,k1,k4)
c    & /zba2(k4,k2,k1,k3))
c    & )*Lnrat(-s(k2,k1),-s(k4,k3))

c    & -0.25d0*(-za(k2,k1)*zb(k1,k6)-za(k2,k3)*zb(k3,k6)
c    & -za(k2,k4)*zb(k4,k6))*(t(k2,k1,k4)-t(k2,k1,k3))
c    & *(za(k2,k5)*delta(k4,k3,k2,k1,k6,k5)-2d0*za(k2,k1)
c    & *(-zb(k1,k2)*za(k2,k5)-zb(k1,k3)*za(k3,k5)-zb(k1,k4)*za(k4,k5)))
c    & /za(k2,k1)/s(k5,k6)/zba2(k4,k2,k1,k3)/DELTA3(k2,k1,k4,k3,k6,k5)


      box1256=(za(k4,k5)*zb(k1,k3))**2
     & /(zb(k1,k2)*za(k5,k6)*t(k1,k2,k3)*zab2(k4,k1,k2,k3))
     & -(zab2(k3,k1,k2,k6)*zab2(k2,k1,k3,k4))**2
     & /(za(k1,k2)*zb(k5,k6)*t(k1,k2,k3)*zab2(k3,k1,k2,k4)**3)

      box6521=(za(k4,k2)*zb(k6,k3))**2
     & /(zb(k6,k5)*za(k2,k1)*t(k6,k5,k3)*zab2(k4,k6,k5,k3))
     & -(zab2(k3,k6,k5,k1)*zab2(k5,k6,k3,k4))**2
     & /(za(k6,k5)*zb(k2,k1)*t(k6,k5,k3)*zab2(k3,k6,k5,k4)**3)

      BDK1211mp=
     & +Lsm1_2mh(s(k3,k4),t(k1,k2,k3),s(k1,k2),s(k5,k6))*box1256
     & +Lsm1_2mh(s(k3,k4),t(k6,k5,k3),s(k6,k5),s(k2,k1))*box6521
     & +i3m(s(k1,k2),s(k3,k4),s(k5,k6))*(T1256+T6521)
     & +squarebracket+squarebracketex
     & -squarebracketfl3-squarebracketfl3ex
      return
      end
