!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWZsr1(Alo,Av,j1,j2,j3,j4,j5,j6,j7,he,swapen,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      logical:: swapen
      integer:: j1,j2,j3,j4,j5,j6,j7,he
      real(dp):: s123
      complex(dp):: l12,l23,L0,L1,Lsm1,Alo,Anlo(-2:0)
      complex(dp):: lnrat,WWZcurr_sr,Av
      integer:: p1,p2,p3,p4,p5,p6,p7
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
      Alo =-(za(j3,j1)*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
     &      +za(j3,j2)*WWZcurr_sr(he,swapen,j3,j4,j2,j5,za,zb))
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
     & +(za(j3,j1)*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
     &  +za(j3,j2)*WWZcurr_sr(he,swapen,j3,j4,j2,j5,za,zb))
     &  /(za(j1,j2)*za(j2,j3)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & -2d0*za(j1,j3)/(za(j1,j2)*za(j2,j3))
     &  *WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
     &  *L0(-s(j2,j3),-s123)/s123

c      Vsc =0.5d0*(epinv+l23)+1d0
      Anlo(-1)=Anlo(-1)+0.5d0*Alo
      Anlo( 0)=Anlo( 0)+(0.5d0*l23+1d0)*Alo

c      Fsc =za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4)
c     & /(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-s(j2,j3),-s123)/s123
c     & +0.5d0*(za(j3,j1)*zb(j1,j5))**2*za(j4,j5)
c     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

      Anlo(0)=Anlo(0)
     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
     &  *L0(-s(j2,j3),-s123)/s123
     & +0.5d0*za(j1,j3)**2/(za(j1,j2)*za(j2,j3))
     &  *(zb(j1,j2)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb)
     &   +zb(j1,j3)*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb))
     &  *L1(-s(j2,j3),-s123)/s123**2

c      A51=(Vcc+Vsc)*Alo+Fcc+Fsc

      Av=Anlo(-2)*epinv*epinv2+Anlo(-1)*epinv+Anlo(0)

      return
      end


c!! c      write(6,*) 'AWWZsr1 Alo',Alo
c!! c      write(6,'(a15,g24.15)') 'AWWZsr1 abs(Alo)=',cdabs(Alo)
c!!
c!! c      write(6,*) 'AWWZsr1 Anlo(-2)',Anlo(-2)
c!! c      write(6,*) 'AWWZsr1 Anlo(-1)',Anlo(-1)
c!! c      write(6,*) 'AWWZsr1 Anlo( 0)',Anlo( 0)
c!! c      write(6,'(a19,2g24.15)') 'AWWZsr1 Anlo(-2)/Alo=',Anlo(-2)/Alo
c!! c      write(6,'(a19,2g24.15)') 'AWWZsr1 Anlo(-1)/Alo=',Anlo(-1)/Alo
c!! c      write(6,'(a19,2g24.15)') 'AWWZsr1 Anlo( 0)/Alo=',Anlo( 0)/Alo
c!!
c!!       boxcoeff=
c!!      & +(za(j3,j1)*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
c!!      &  +za(j3,j2)*WWZcurr_sr(he,swapen,j3,j4,j2,j5,za,zb))
c!!      &  /(za(j1,j2)*za(j2,j3)*s123)
c!!      &  *s(j1,j2)*s(j2,j3)/2d0
c!!
c!!       bubcoeff=
c!!      & -2d0*za(j1,j3)/(za(j1,j2)*za(j2,j3))
c!!      &  *WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
c!!      &  /(1d0-s(j2,j3)/s123)/s123
c!!
c!! c--- old version
c!!       bubcoeff=bubcoeff
c!!      & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb)
c!!      &  /(1d0-s(j2,j3)/s123)/s123
c!!      & +0.5d0*za(j1,j3)**2/(za(j1,j2)*za(j2,j3))
c!!      &  *(zb(j1,j2)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb)
c!!      &   +zb(j1,j3)*WWZcurr_sr(he,swapen,j3,j4,j1,j5,za,zb))
c!!      &  /(1d0-s(j2,j3)/s123)**2/s123**2
c!!
c!! c--- new version
c!! c      bubcoeff=bubcoeff
c!! c     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*WWZcurr_sr(he,swapen,j3,j4,j1,j5)
c!! c     &  /(1d0-s(j2,j3)/s123)/s123
c!! c     & +0.5d0*za(j3,j1)/(za(j1,j2)*za(j2,j3))
c!! c     &       *(za(j1,j2)*zb(j2,j1)+za(j1,j3)*zb(j3,j1))
c!! c     &  *WWZcurr_sr(he,swapen,j3,j4,j1,j5)
c!! c     &  /(1d0-s(j2,j3)/s123)**2/s123**2
c!! c     & -0.5d0*za(j3,j1)/(za(j1,j2)*za(j2,j3))
c!! c     &       *za(j2,j3)*zb(j1,j2)
c!! c     &  *WWZcurr_sr(he,swapen,j1,j4,j1,j5)
c!! c     &  /(1d0-s(j2,j3)/s123)**2/s123**2
c!!
c!!
c!! c      write(6,*) 'boxcoeff/Alo',boxcoeff/Alo
c!! c      write(6,*) 'bubcoeff1/Alo',bubcoeff/Alo
c!! c      write(6,*) 'bubcoeff2/Alo',(-bubcoeff-3d0/2d0*Alo)/Alo
c!!
c!! c      A51=(Vcc+Vsc)*Alo+Fcc+Fsc
c!!
c!!       return
c!!       end






      subroutine AWWZsr2(Alo,Av,j1,j2,j3,j4,j5,j6,j7,he,swapen,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      logical:: swapen
      integer j1,j2,j3,j4,j5,j6,j7,he
      real(dp):: s123
      complex(dp):: l12,l123,L0,L1,Lsm1,Alo,Anlo(-2:0)
      complex(dp):: lnrat,WWZcurr_sr,Av
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
      Alo =-(za(j2,j1)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb)
     &      +za(j2,j3)*WWZcurr_sr(he,swapen,j2,j4,j3,j5,za,zb))
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
     & -(za(j2,j1)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb)
     &  +za(j2,j3)*WWZcurr_sr(he,swapen,j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j3,j1)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j1,j3), -s123)
     & +(za(j1,j2)*za(j3,j1)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb)
     &  +za(j1,j2)*za(j3,j2)*WWZcurr_sr(he,swapen,j2,j4,j2,j5,za,zb)
     &  -za(j1,j2)*za(j2,j3)*WWZcurr_sr(he,swapen,j2,j4,j2,j5,za,zb)
     &  -za(j1,j3)*za(j2,j3)*WWZcurr_sr(he,swapen,j2,j4,j3,j5,za,zb))
     &  /(za(j2,j3)*za(j1,j3)**2*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & +2d0*zb(j1,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*WWZcurr_sr(he,swapen,j2,j4,j2,j5,za,zb)
     &  +za(j1,j3)*WWZcurr_sr(he,swapen,j2,j4,j3,j5,za,zb))
     & *L0(-s(j2,j3),-s123)/s123

c     & +2d0*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5))
c     & *L0(-s(j2,j3),-s123)/s123

c      Vsc=0.5d0*(epinv+l123)+0.5d0
      Anlo(-1)=Anlo(-1)+0.5d0*Alo
      Anlo( 0)=Anlo( 0)+0.5d0*(l123+1d0)*Alo

      Anlo(0)=Anlo(0)
     & +(za(j1,j2)*WWZcurr_sr(he,swapen,j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*WWZcurr_sr(he,swapen,j1,j4,j3,j5,za,zb))
     &  *za(j2,j3)/(za(j1,j3)**3*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & -0.5d0*zb(j1,j3)**2*za(j2,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*WWZcurr_sr(he,swapen,j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*WWZcurr_sr(he,swapen,j1,j4,j3,j5,za,zb))
     & *L1(-s123,-s(j2,j3))/s(j2,j3)**2

     & +za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*WWZcurr_sr(he,swapen,j1,j4,j2,j5,za,zb)
     &  +za(j1,j3)*WWZcurr_sr(he,swapen,j1,j4,j3,j5,za,zb))
     & *L0(-s123,-s(j2,j3))/s(j2,j3)

     & +za(j2,j1)*zb(j1,j3)*WWZcurr_sr(he,swapen,j3,j4,j3,j5,za,zb)/za(j1,j3)
     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2

     & -za(j2,j1)*zb(j1,j3)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*WWZcurr_sr(he,swapen,j3,j4,j2,j5,za,zb)
     &  +za(j1,j3)*WWZcurr_sr(he,swapen,j3,j4,j3,j5,za,zb))
     & *L0(-s123,-s(j1,j2))/s(j1,j2)

     & +0.5d0*zb(j1,j3)
     & *(zb(j3,j1)*WWZcurr_sr(he,swapen,j1,j4,j2,j5,za,zb)
     &  +zb(j3,j2)*WWZcurr_sr(he,swapen,j2,j4,j2,j5,za,zb))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

     & +0.5d0*zb(j2,j3)
     & *(zb(j3,j1)*WWZcurr_sr(he,swapen,j1,j4,j1,j5,za,zb)
     &  +zb(j3,j2)*WWZcurr_sr(he,swapen,j2,j4,j1,j5,za,zb))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

      Av=-(Anlo(-2)*epinv*epinv2+Anlo(-1)*epinv+Anlo(0))

      return
      end
