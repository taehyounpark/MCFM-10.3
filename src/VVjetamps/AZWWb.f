!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZWWb1(Alo,Av,j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7
      real(dp):: s123
      complex(dp):: l12,l23,L0,L1,Lsm1,Alo,Anlo(-2:0)
      complex(dp):: lnrat,ZWWcurr_ab,boxcoeff,bubcoeff,Av
      integer p1,p2,p3,p4,p5,p6,p7
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

      p1=j3
      p2=j1
      p3=j4
      p4=j5
      p5=j6
      p6=j7
      p7=j2

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
c    -i * A5tree rewritten in a suitable way
c      A5lom =-za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
      Alo =-(za(j3,j1)*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &      +za(j3,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb))
     &     /(za(j1,j2)*za(j2,j3)*s123)

      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))

c--leading N
c      Vcc=
c     & -(epinv**2+epinv*l12+0.5d0*l12**2)
c     & -(epinv**2+epinv*l23+0.5d0*l23**2)
c     & -2d0*(epinv+l23)-4d0
      Anlo(-2)=-2d0*Alo
      Anlo(-1)=-(l12+l23+2d0)*Alo
      Anlo( 0)=(-0.5d0*(l12**2+l23**2)-2d0*l23-4d0)*Alo

c      Fcc=za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
c     & *(Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
c     & -2d0*za(j3,j1)*zb(j1,j5)*za(j5,j4)/za(j3,j4)
c     &   *L0(-s(j2,j3),-s(j4,j5))/s(j4,j5))

      Anlo(0)=Anlo(0)
     & +(za(j3,j1)*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  +za(j3,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb))
     &  /(za(j1,j2)*za(j2,j3)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & -2d0*za(j1,j3)/(za(j1,j2)*za(j2,j3))
     &  *ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  *L0(-s(j2,j3),-s123)/s123

c      Vsc =0.5d0*(epinv+l23)+1d0
      Anlo(-1)=Anlo(-1)+0.5d0*Alo
      Anlo( 0)=Anlo( 0)+(0.5d0*l23+1d0)*Alo

c      Fsc =za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4)
c     & /(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-s(j2,j3),-s123)/s123
c     & +0.5d0*(za(j3,j1)*zb(j1,j5))**2*za(j4,j5)
c     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

      Anlo(0)=Anlo(0)
     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  *L0(-s(j2,j3),-s123)/s123
     & +0.5d0*za(j1,j3)**2/(za(j1,j2)*za(j2,j3))
     &  *(zb(j1,j2)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &   +zb(j1,j3)*ZWWcurr_ab(j3,j4,j1,j5,za,zb))
     &  *L1(-s(j2,j3),-s123)/s123**2

c      A51=(Vcc+Vsc)*Alo+Fcc+Fsc

      Av=Anlo(-2)*epinv*epinv2+Anlo(-1)*epinv+Anlo(0)

c      write(6,*) 'AZWWb1 Alo',Alo
c      write(6,'(a15,g24.15)') 'AZWWb1 abs(Alo)=',cdabs(Alo)

c      write(6,*) 'AZWWb1 Anlo(-2)',Anlo(-2)
c      write(6,*) 'AZWWb1 Anlo(-1)',Anlo(-1)
c      write(6,*) 'AZWWb1 Anlo( 0)',Anlo( 0)
c      write(6,'(a19,2g24.15)') 'AZWWb1 Anlo(-2)/Alo=',Anlo(-2)/Alo
c      write(6,'(a19,2g24.15)') 'AZWWb1 Anlo(-1)/Alo=',Anlo(-1)/Alo
c      write(6,'(a19,2g24.15)') 'AZWWb1 Anlo( 0)/Alo=',Anlo( 0)/Alo

      boxcoeff=
     & +(za(j3,j1)*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  +za(j3,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb))
     &  /(za(j1,j2)*za(j2,j3)*s123)
     &  *s(j1,j2)*s(j2,j3)/2d0

      bubcoeff=
     & -2d0*za(j1,j3)/(za(j1,j2)*za(j2,j3))
     &  *ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  /(1d0-s(j2,j3)/s123)/s123

c--- old version
      bubcoeff=bubcoeff
     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
     &  /(1d0-s(j2,j3)/s123)/s123
     & +0.5d0*za(j1,j3)**2/(za(j1,j2)*za(j2,j3))
     &  *(zb(j1,j2)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &   +zb(j1,j3)*ZWWcurr_ab(j3,j4,j1,j5,za,zb))
     &  /(1d0-s(j2,j3)/s123)**2/s123**2

c--- new version
c      bubcoeff=bubcoeff
c     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
c     &  /(1d0-s(j2,j3)/s123)/s123
c     & +0.5d0*za(j3,j1)/(za(j1,j2)*za(j2,j3))
c     &       *(za(j1,j2)*zb(j2,j1)+za(j1,j3)*zb(j3,j1))
c     &  *ZWWcurr_ab(j3,j4,j1,j5)
c     &  /(1d0-s(j2,j3)/s123)**2/s123**2
c     & -0.5d0*za(j3,j1)/(za(j1,j2)*za(j2,j3))
c     &       *za(j2,j3)*zb(j1,j2)
c     &  *ZWWcurr_ab(j1,j4,j1,j5)
c     &  /(1d0-s(j2,j3)/s123)**2/s123**2


c      write(6,*) 'boxcoeff/Alo',boxcoeff/Alo
c      write(6,*) 'bubcoeff1/Alo',bubcoeff/Alo
c      write(6,*) 'bubcoeff2/Alo',(-bubcoeff-3d0/2d0*Alo)/Alo

c      A51=(Vcc+Vsc)*Alo+Fcc+Fsc

      return
      end






      subroutine AZWWb2(Alo,Av,j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7
      real(dp):: s123
      complex(dp):: l12,l123,L0,L1,Lsm1,Alo,Anlo(-2:0)
      complex(dp):: lnrat,ZWWcurr_ab,boxcoeff1,boxcoeff2,bubcoeff,Av,
     & bubcoeffa,bubcoeffb,bubcoeffc,bubcoeff12,bubcoeff17,bubcoeff127
      integer p1,p2,p3,p4,p5,p6,p7
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

      p1=j2
      p2=j1
      p3=j4
      p4=j5
      p5=j6
      p6=j7
      p7=j3

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
c    -i * A5tree rewritten in a suitable way
c      Alom=za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))
      Alo =-(za(j2,j1)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &      +za(j2,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     &     /(za(j1,j3)*za(j2,j3)*s123)

      l12=lnrat(musq,-s(j1,j2))
      l123=lnrat(musq,-s123)

c      Vcc=-(epinv**2+epinv*l12+0.5d0*l12**2)
c     & -2d0*(epinv+l123)-4d0
      Anlo(-2)=-Alo
      Anlo(-1)=-(l12+2d0)*Alo
      Anlo( 0)=-(0.5d0*l12**2+2d0*l123+4d0)*Alo

c--subleading N
      Anlo(0)=Anlo(0)
     & -(za(j2,j1)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &  +za(j2,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j3,j1)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j1,j3), -s123)
     & +(za(j1,j2)*za(j3,j1)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &  +za(j1,j2)*za(j3,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  -za(j1,j2)*za(j2,j3)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  -za(j1,j3)*za(j2,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j1,j3)**2*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & +2d0*zb(j1,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     & *L0(-s(j2,j3),-s123)/s123

c     & +2d0*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5))
c     & *L0(-s(j2,j3),-s123)/s123

c      Vsc=0.5d0*(epinv+l123)+0.5d0
      Anlo(-1)=Anlo(-1)+0.5d0*Alo
      Anlo( 0)=Anlo( 0)+0.5d0*(l123+1d0)*Alo

      Anlo(0)=Anlo(0)
     & +(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     &  *za(j2,j3)/(za(j1,j3)**3*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & -0.5d0*zb(j1,j3)**2*za(j2,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     & *L1(-s123,-s(j2,j3))/s(j2,j3)**2

     & +za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     & *L0(-s123,-s(j2,j3))/s(j2,j3)

     & +za(j2,j1)*zb(j1,j3)*ZWWcurr_ab(j3,j4,j3,j5,za,zb)/za(j1,j3)
     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2

     & -za(j2,j1)*zb(j1,j3)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j3,j4,j3,j5,za,zb))
     & *L0(-s123,-s(j1,j2))/s(j1,j2)

     & +0.5d0*zb(j1,j3)
     & *(zb(j3,j1)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +zb(j3,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

     & +0.5d0*zb(j2,j3)
     & *(zb(j3,j1)*ZWWcurr_ab(j1,j4,j1,j5,za,zb)
     &  +zb(j3,j2)*ZWWcurr_ab(j2,j4,j1,j5,za,zb))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

c     & -0.5d0*(za(j4,j1)*zb(j1,j3))**2*za(j2,j3)/(za(j1,j3)*za(j4,j5))
c     & *L1(-s123,-s(j2,j3))/s(j2,j3)**2

c     & +za(j1,j4)**2*za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*za(j4,j5))
c     & *L0(-s123,-s(j2,j3))/s(j2,j3)

c     & -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3)
c     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2

c     & -za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(za(j1,j3)**2*za(j4,j5))
c     & *L0(-s123,-s(j1,j2))/s(j1,j2)

c     & -0.5d0*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5))
c     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5))

c      A52=(Vcc+Vsc)*A5lom+Fcc+Fsc

      Av=-(Anlo(-2)*epinv*epinv2+Anlo(-1)*epinv+Anlo(0))

c      write(6,*)
c      write(6,*) 'AZWWb2 Alo',Alo
c      write(6,'(a15,g24.15)') 'AZWWb2 abs(Alo)=',cdabs(Alo)

c      write(6,*) 'AZWWb2 Anlo(-2)',Anlo(-2)
c      write(6,*) 'AZWWb2 Anlo(-1)',Anlo(-1)
c      write(6,*) 'AZWWb2 Anlo( 0)',Anlo( 0)
c      write(6,'(a19,2g24.15)') 'AZWWb2 Anlo(-2)/Alo=',Anlo(-2)/Alo
c      write(6,'(a19,2g24.15)') 'AZWWb2 Anlo(-1)/Alo=',Anlo(-1)/Alo
c      write(6,'(a19,2g24.15)') 'AZWWb2 Anlo( 0)/Alo=',Anlo( 0)/Alo

      return












      boxcoeff1=
     & -(za(j2,j1)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &  +za(j2,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j3,j1)*s123)
     &  *s(j1,j2)*s(j1,j3)/2d0
c--- this version is direct from paper
      boxcoeff2=(
     & +(za(j1,j2)*za(j3,j1)*ZWWcurr_ab(j2,j4,j1,j5,za,zb)
     &  +za(j1,j2)*za(j3,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  -za(j1,j2)*za(j2,j3)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  -za(j1,j3)*za(j2,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j1,j3)**2*s123)
     & +(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     &  *za(j2,j3)/(za(j1,j3)**3*s123)
     &  )*s(j1,j2)*s(j2,j3)/2d0
c--- this version has been simplified by hand
c      boxcoeff2=za(j1,j2)**2
c     & *(za(j3,j1)*ZWWcurr_ab(j3,j4,j1,j5,za,zb)
c     &  +za(j3,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb))
c     & /(za(j1,j3)**3*za(j2,j3)*s123)
c     & *s(j1,j2)*s(j2,j3)/2d0

      write(6,*) 'boxcoeff1/Alo',boxcoeff1/Alo,cdabs(boxcoeff1/Alo)
      write(6,*) 'boxcoeff2/Alo',boxcoeff2/Alo,cdabs(boxcoeff2/Alo)

      bubcoeff=
     & +za(j2,j1)*zb(j1,j3)*ZWWcurr_ab(j3,j4,j3,j5,za,zb)/za(j1,j3)
     & /(1d0-s123/s(j1,j2))**2/s(j1,j2)**2

     & -za(j2,j1)*zb(j1,j3)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j3,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j3,j4,j3,j5,za,zb))
     & /(1d0-s123/s(j1,j2))/s(j1,j2)

      write(6,*) 's12=s12 bubcoeff/Alo',bubcoeff/Alo,cdabs(bubcoeff/Alo)
      bubcoeff12=bubcoeff

      bubcoeffa=
     & +2d0*zb(j1,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j2,j4,j3,j5,za,zb))
     & /(1d0-s(j2,j3)/s123)/s123

      bubcoeffb=
     & -0.5d0*zb(j1,j3)**2*za(j2,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     & /(1d0-s123/s(j2,j3))**2/s(j2,j3)**2

      bubcoeffc=
     & +za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
     & /(1d0-s123/s(j2,j3))/s(j2,j3)

      bubcoeff=-bubcoeffa+bubcoeffb+bubcoeffc
c      write(6,*) 's23 bubcoeffa/Alo',bubcoeffa/Alo,cdabs(bubcoeffa/Alo)
c      write(6,*) 's23 bubcoeffb/Alo',bubcoeffb/Alo,cdabs(bubcoeffb/Alo)
c      write(6,*) 's23 bubcoeffc/Alo',bubcoeffc/Alo,cdabs(bubcoeffc/Alo)
      write(6,*) 's23=s17 bubcoeff/Alo',bubcoeff/Alo,cdabs(bubcoeff/Alo)
      bubcoeff17=bubcoeff
      bubcoeff127=-1.5d0*Alo-bubcoeff17-bubcoeff12
      write(6,*) 'bubcoeff127/Alo',bubcoeff127/Alo

c      ratbit=
c     & -0.5d0*zb(j1,j3)**2*za(j2,j3)/(za(j1,j3)*s123)
c     & *(za(j1,j2)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
c     &  +za(j1,j3)*ZWWcurr_ab(j1,j4,j3,j5,za,zb))
c     & /(1d0-s123/s(j2,j3))/s(j2,j3)**2 ! from L1

c     & +za(j2,j1)*zb(j1,j3)*ZWWcurr_ab(j3,j4,j3,j5)/za(j1,j3)
c     & /(1d0-s123/s(j1,j2))/s(j1,j2)**2 ! from L1

c     & +0.5d0*zb(j1,j3)
c     & *(zb(j3,j1)*ZWWcurr_ab(j1,j4,j2,j5,za,zb)
c     &  +zb(j3,j2)*ZWWcurr_ab(j2,j4,j2,j5,za,zb))
c     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

c     & +0.5d0*zb(j2,j3)
c     & *(zb(j3,j1)*ZWWcurr_ab(j1,j4,j1,j5,za,zb)
c     &  +zb(j3,j2)*ZWWcurr_ab(j2,j4,j1,j5,za,zb))
c     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)
c      ratbit=ratbit-0.5d0*Alo

c      write(6,*) 'ratbit/Alo',ratbit/Alo,cdabs(ratbit/Alo)


c      write(6,*) '0 bit',(Anlo(0)-ratbit)/Alo


c      Anlo(0)=
c     & +boxcoeff1*Lsm1(-s(j1,j2),-s123,-s(j1,j3), -s123)
c     &           /(s(j1,j2)*s(j1,j3)/2d0)
c     & +boxcoeff2*Lsm1(-s(j1,j2),-s123,-s(j2,j3), -s123)
c     &           /(s(j1,j2)*s(j2,j3)/2d0)
c     & +bubcoeff12*lnrat(musq,-s(j1,j2))
c     & +bubcoeff17*lnrat(musq,-s(j2,j3))
c     & +bubcoeff127*lnrat(musq,-s123)
c     & -(0.5d0*l12**2+2d0*l123*0d0+4d0)*Alo
c     & +0.5d0*(l123*0d0+1d0)*Alo
c     & +0.5d0*Alo ! added by hand
c           write(6,*) 'New Anlo(0)/Alo',Anlo(0)/Alo

      return
      end
