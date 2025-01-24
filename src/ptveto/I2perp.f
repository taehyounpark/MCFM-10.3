!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine I2perp(x,Lperp,LQ,I2)
      implicit none
!     written by writeI2perp.frm
!     using results from https://arxiv.org/abs/2207.07037
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-1,dmax=2,delt=-1,plus=0,lpls=1,rglr=2
      real(dp):: I2(0:6,dmin:dmax)
c      I2, index one transition flavor
c      I2, index2,distribution type

      real(dp):: x,Lperp,LQ,Li2,Li3,omx,opx,lnx,lnomx,lnopx,
     & Li2x,Li2mx,Li2omx,Li3x,Li3mx,Li3omx,Li3xonopx
      integer,parameter::n1=-1,n2=1
      integer::nf
      real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2)

! Explicit expressions for hplog (in terms of Li2, Li3) are now included
! below, but the following may be uncommented for comparison
! Using Hplog program, hep-ph/0107173

      omx=one-x
      opx=one+x
      lnx=log(x); lnomx=log(omx); lnopx=log(opx)
      Li2x=Li2(x); Li2mx=Li2(-x); Li2omx=Li2(omx)
      Li3x=Li3(x); Li3mx=Li3(-x); Li3omx=Li3(omx)
      Li3xonopx=Li3(x/(opx))

!      we need Hr1(-1),Hr1( 0),Hr1(+1)
      Hr1( 0)=lnx
      Hr1(-1)=lnopx
      Hr1(+1)=-lnomx

      Hr2(0,+1)=(Li2x)
      Hr2(0,-1)=(-Li2mx)

!      we need Hr3( 0, 0,-1),Hr3( 0, 0,+1),Hr3( 0,-1,-1),Hr3( 0,+1,+1)
      Hr3(0,0,+1)=+Li3x
      Hr3(0,0,-1)=-Li3mx
      Hr3(0,+1,+1)=0.5d0*lnx*lnomx**2+lnomx*Li2omx-Li3omx+zeta3
      Hr3(0,-1,-1)=-lnopx**3/6d0-lnopx*Li2mx+Li3mx+Li3xonopx

      nf=nfl
c     initialize to zero
      I2(:,:)=0

      I2(gg,delt)= - 224.D0/27.D0*CA*TF*nf*LQ - 160.D0/9.D0*CA*TF*nf*LQ
     & *Lperp - 32.D0/3.D0*CA*TF*nf*LQ*Lperp**2 + 808.D0/27.D0*CA**2*LQ
     &  + 536.D0/9.D0*CA**2*LQ*Lperp + 88.D0/3.D0*CA**2*LQ*Lperp**2 + 
     & 32*CA**2*LQ**2*Lperp**2 - 16.D0/3.D0*CA**2*pisq*Lperp**2 - 8.D0/
     & 3.D0*CA**2*pisq*LQ*Lperp - 28*CA**2*zeta3*LQ

      I2(gg,plus)=224.D0/27.D0*CA*TF*nf + 160.D0/9.D0*CA*TF*nf*Lperp + 
     & 32.D0/3.D0*CA*TF*nf*Lperp**2 - 808.D0/27.D0*CA**2 - 536.D0/9.D0*
     & CA**2*Lperp - 88.D0/3.D0*CA**2*Lperp**2 - 64*CA**2*LQ*Lperp**2
     &  + 8.D0/3.D0*CA**2*pisq*Lperp + 28*CA**2*zeta3

      I2(gg,lpls)=64*CA**2*Lperp**2

      I2(gg,rglr)= - 8.D0/3.D0*CF*TF*nf*x**(-1) - 32.D0/3.D0*CF*TF*nf*
     & x**(-1)*Lperp + 64.D0/3.D0*CF*TF*nf*x**(-1)*Lperp**2 + 64*CF*TF*
     & nf + 112*CF*TF*nf*Lperp + 16*CF*TF*nf*Lperp**2 - 64*CF*TF*nf*x
     &  - 80*CF*TF*nf*x*Lperp - 16*CF*TF*nf*x*Lperp**2 + 8.D0/3.D0*CF*
     & TF*nf*x**2 - 64.D0/3.D0*CF*TF*nf*x**2*Lperp - 64.D0/3.D0*CF*TF*
     & nf*x**2*Lperp**2 + 368.D0/9.D0*CA*TF*nf*x**(-1)*opx**(-1)*Lperp
     &  + 484.D0/27.D0*CA*TF*nf*x**(-1) + 32.D0/3.D0*CA*TF*nf*x**(-1)*
     & Lperp**2 - 32.D0/3.D0*CA*TF*nf*opx**(-1)*Lperp - 664.D0/27.D0*CA
     & *TF*nf - 64.D0/3.D0*CA*TF*nf*Lperp**2 - 160.D0/9.D0*CA*TF*nf*x*
     & opx**(-1)*Lperp + 440.D0/27.D0*CA*TF*nf*x + 32.D0/3.D0*CA*TF*nf*
     & x*Lperp**2 - 64.D0/9.D0*CA*TF*nf*x**2*opx**(-1)*Lperp - 556.D0/
     & 27.D0*CA*TF*nf*x**2 - 32.D0/3.D0*CA*TF*nf*x**2*Lperp**2 - 368.D0/
     & 9.D0*CA*TF*nf*x**3*opx**(-1)*Lperp - 3160.D0/27.D0*CA**2*x**(-1)
     &  - 440.D0/3.D0*CA**2*x**(-1)*Lperp**2 - 64*CA**2*x**(-1)*LQ*
     & Lperp**2 + 44.D0/9.D0*CA**2*x**(-1)*pisq + 24*CA**2*x**(-1)*
     & zeta3*opx**(-1)
      I2(gg,rglr) = I2(gg,rglr) + 24*CA**2*x**(-1)*zeta3 + 3704.D0/27.D0
     & *CA**2 + 100.D0/9.D0*CA**2*Lperp + 464.D0/3.D0*CA**2*Lperp**2 + 
     & 128*CA**2*LQ*Lperp**2 - 8*CA**2*pisq*opx**(-1)*Lperp - 16.D0/3.D0
     & *CA**2*pisq + 24*CA**2*zeta3*omx**(-1) - 42*CA**2*zeta3*
     & opx**(-1) - 42*CA**2*zeta3 - 3040.D0/27.D0*CA**2*x + 436.D0/9.D0
     & *CA**2*x*Lperp - 376.D0/3.D0*CA**2*x*Lperp**2 - 64*CA**2*x*LQ*
     & Lperp**2 - 32.D0/3.D0*CA**2*x*pisq*opx**(-1)*Lperp + 16.D0/3.D0*
     & CA**2*x*pisq + 10*CA**2*x*zeta3*opx**(-1) - 32*CA**2*x*zeta3 + 
     & 3340.D0/27.D0*CA**2*x**2 + 440.D0/3.D0*CA**2*x**2*Lperp**2 + 64*
     & CA**2*x**2*LQ*Lperp**2 - 16.D0/3.D0*CA**2*x**2*pisq*opx**(-1)*
     & Lperp - 44.D0/9.D0*CA**2*x**2*pisq + 28*CA**2*x**2*zeta3*
     & opx**(-1) - 4*CA**2*x**2*zeta3 - 16.D0/3.D0*CA**2*x**3*pisq*
     & opx**(-1)*Lperp - 24*CA**2*x**3*zeta3*opx**(-1) - 28*CA**2*x**3*
     & zeta3 + 28*CA**2*x**4*zeta3*opx**(-1) + 16*Hr3(0,-1,-1)*CA**2*
     & x**(-1)*opx**(-1) + 32*Hr3(0,-1,-1)*CA**2*opx**(-1) + 48*Hr3(0,
     & -1,-1)*CA**2*x*opx**(-1)
      I2(gg,rglr) = I2(gg,rglr) + 32*Hr3(0,-1,-1)*CA**2*x**2*opx**(-1)
     &  + 16*Hr3(0,-1,-1)*CA**2*x**3*opx**(-1) + 24*Hr3(0,0,-1)*CA**2*
     & x**(-1)*opx**(-1) + 48*Hr3(0,0,-1)*CA**2*opx**(-1) + 72*Hr3(0,0,
     & -1)*CA**2*x*opx**(-1) + 48*Hr3(0,0,-1)*CA**2*x**2*opx**(-1) + 24
     & *Hr3(0,0,-1)*CA**2*x**3*opx**(-1) - 20*Hr3(0,0,1)*CA**2*x**(-1)*
     & opx**(-1) - 20*Hr3(0,0,1)*CA**2*x**(-1) - 24*Hr3(0,0,1)*CA**2*
     & omx**(-1) + 4*Hr3(0,0,1)*CA**2*opx**(-1) + 8*Hr3(0,0,1)*CA**2 - 
     & 20*Hr3(0,0,1)*CA**2*x*opx**(-1) - 12*Hr3(0,0,1)*CA**2*x - 4*Hr3(
     & 0,0,1)*CA**2*x**2*opx**(-1) - 16*Hr3(0,0,1)*CA**2*x**2 + 20*Hr3(
     & 0,0,1)*CA**2*x**3*opx**(-1) + 4*Hr3(0,0,1)*CA**2*x**3 - 4*Hr3(0,
     & 0,1)*CA**2*x**4*opx**(-1) + 32*Hr2(0,-1)*CA**2*x**(-1)*opx**(-1)
     & *Lperp + 64*Hr2(0,-1)*CA**2*opx**(-1)*Lperp + 96*Hr2(0,-1)*CA**2
     & *x*opx**(-1)*Lperp + 64*Hr2(0,-1)*CA**2*x**2*opx**(-1)*Lperp + 
     & 32*Hr2(0,-1)*CA**2*x**3*opx**(-1)*Lperp - 16*Hr2(0,-1)*Hr1(-1)*
     & CA**2*x**(-1)*opx**(-1) - 32*Hr2(0,-1)*Hr1(-1)*CA**2*opx**(-1)
     &  - 48*Hr2(0,-1)*Hr1(-1)*CA**2*x*opx**(-1)
      I2(gg,rglr) = I2(gg,rglr) - 32*Hr2(0,-1)*Hr1(-1)*CA**2*x**2*
     & opx**(-1) - 16*Hr2(0,-1)*Hr1(-1)*CA**2*x**3*opx**(-1) - 8*Hr2(0,
     & -1)*Hr1(0)*CA**2*x**(-1)*opx**(-1) - 16*Hr2(0,-1)*Hr1(0)*CA**2*
     & opx**(-1) - 24*Hr2(0,-1)*Hr1(0)*CA**2*x*opx**(-1) - 16*Hr2(0,-1)
     & *Hr1(0)*CA**2*x**2*opx**(-1) - 8*Hr2(0,-1)*Hr1(0)*CA**2*x**3*
     & opx**(-1) - 88.D0/3.D0*Hr2(0,1)*CA**2*x**(-1) + 32*Hr2(0,1)*
     & CA**2 - 32*Hr2(0,1)*CA**2*x + 88.D0/3.D0*Hr2(0,1)*CA**2*x**2 + 
     & 12*Hr2(0,1)*Hr1(0)*CA**2*x**(-1)*opx**(-1) + 12*Hr2(0,1)*Hr1(0)*
     & CA**2*x**(-1) + 16*Hr2(0,1)*Hr1(0)*CA**2*omx**(-1) - 4*Hr2(0,1)*
     & Hr1(0)*CA**2*opx**(-1) - 8*Hr2(0,1)*Hr1(0)*CA**2 + 12*Hr2(0,1)*
     & Hr1(0)*CA**2*x*opx**(-1) + 4*Hr2(0,1)*Hr1(0)*CA**2*x + 4*Hr2(0,1
     & )*Hr1(0)*CA**2*x**2*opx**(-1) + 8*Hr2(0,1)*Hr1(0)*CA**2*x**2 - 
     & 12*Hr2(0,1)*Hr1(0)*CA**2*x**3*opx**(-1) - 4*Hr2(0,1)*Hr1(0)*
     & CA**2*x**3 + 4*Hr2(0,1)*Hr1(0)*CA**2*x**4*opx**(-1) + 4.D0/3.D0*
     & Hr1(-1)*CA**2*x**(-1)*pisq*opx**(-1) + 8.D0/3.D0*Hr1(-1)*CA**2*
     & pisq*opx**(-1)
      I2(gg,rglr) = I2(gg,rglr) + 4*Hr1(-1)*CA**2*x*pisq*opx**(-1) + 8.D
     & 0/3.D0*Hr1(-1)*CA**2*x**2*pisq*opx**(-1) + 4.D0/3.D0*Hr1(-1)*
     & CA**2*x**3*pisq*opx**(-1) + 8*Hr1(-1)**2*Hr1(0)*CA**2*x**(-1)*
     & opx**(-1) + 16*Hr1(-1)**2*Hr1(0)*CA**2*opx**(-1) + 24*Hr1(-1)**2
     & *Hr1(0)*CA**2*x*opx**(-1) + 16*Hr1(-1)**2*Hr1(0)*CA**2*x**2*
     & opx**(-1) + 8*Hr1(-1)**2*Hr1(0)*CA**2*x**3*opx**(-1) - 32*Hr1(-1
     & )*Hr1(0)*CA**2*x**(-1)*opx**(-1)*Lperp - 64*Hr1(-1)*Hr1(0)*CA**2
     & *opx**(-1)*Lperp - 96*Hr1(-1)*Hr1(0)*CA**2*x*opx**(-1)*Lperp - 
     & 64*Hr1(-1)*Hr1(0)*CA**2*x**2*opx**(-1)*Lperp - 32*Hr1(-1)*Hr1(0)
     & *CA**2*x**3*opx**(-1)*Lperp - 4*Hr1(-1)*Hr1(0)**2*CA**2*x**(-1)*
     & opx**(-1) - 8*Hr1(-1)*Hr1(0)**2*CA**2*opx**(-1) - 12*Hr1(-1)*
     & Hr1(0)**2*CA**2*x*opx**(-1) - 8*Hr1(-1)*Hr1(0)**2*CA**2*x**2*
     & opx**(-1) - 4*Hr1(-1)*Hr1(0)**2*CA**2*x**3*opx**(-1) + 24*Hr1(0)
     & *CF*TF*nf + 48*Hr1(0)*CF*TF*nf*Lperp + 32*Hr1(0)*CF*TF*nf*
     & Lperp**2 + 24*Hr1(0)*CF*TF*nf*x + 48*Hr1(0)*CF*TF*nf*x*Lperp + 
     & 32*Hr1(0)*CF*TF*nf*x*Lperp**2
      I2(gg,rglr) = I2(gg,rglr) + 52.D0/9.D0*Hr1(0)*CA*TF*nf + 32.D0/3.D
     & 0*Hr1(0)*CA*TF*nf*Lperp + 40.D0/9.D0*Hr1(0)*CA*TF*nf*x + 32.D0/3.
     & D0*Hr1(0)*CA*TF*nf*x*Lperp - 32*Hr1(0)*CA**2*x**(-1)*Lperp**2 - 
     & 32*Hr1(0)*CA**2*omx**(-1)*Lperp**2 - 701.D0/9.D0*Hr1(0)*CA**2 + 
     & 200.D0/3.D0*Hr1(0)*CA**2*Lperp - 149.D0/9.D0*Hr1(0)*CA**2*x - 88.
     & D0/3.D0*Hr1(0)*CA**2*x*Lperp - 96*Hr1(0)*CA**2*x*Lperp**2 - 536.D
     & 0/9.D0*Hr1(0)*CA**2*x**2 + 352.D0/3.D0*Hr1(0)*CA**2*x**2*Lperp
     &  + 32*Hr1(0)*CA**2*x**2*Lperp**2 + 6*Hr1(0)**2*CF*TF*nf + 16*
     & Hr1(0)**2*CF*TF*nf*Lperp + 2*Hr1(0)**2*CF*TF*nf*x + 16*Hr1(0)**2
     & *CF*TF*nf*x*Lperp + 4.D0/3.D0*Hr1(0)**2*CA*TF*nf + 4.D0/3.D0*
     & Hr1(0)**2*CA*TF*nf*x - 8*Hr1(0)**2*CA**2*omx**(-1)*Lperp - 8*
     & Hr1(0)**2*CA**2*opx**(-1)*Lperp + 25.D0/3.D0*Hr1(0)**2*CA**2 - 
     & 16*Hr1(0)**2*CA**2*x*opx**(-1)*Lperp - 11.D0/3.D0*Hr1(0)**2*
     & CA**2*x - 16*Hr1(0)**2*CA**2*x*Lperp + 8*Hr1(0)**2*CA**2*x**2*
     & opx**(-1)*Lperp + 44.D0/3.D0*Hr1(0)**2*CA**2*x**2 - 8*Hr1(0)**2*
     & CA**2*x**2*Lperp
      I2(gg,rglr) = I2(gg,rglr) + 16*Hr1(0)**2*CA**2*x**3*opx**(-1)*
     & Lperp + 8*Hr1(0)**2*CA**2*x**3*Lperp - 8*Hr1(0)**2*CA**2*x**4*
     & opx**(-1)*Lperp + 4.D0/3.D0*Hr1(0)**3*CF*TF*nf + 4.D0/3.D0*Hr1(0
     & )**3*CF*TF*nf*x - 2.D0/3.D0*Hr1(0)**3*CA**2*omx**(-1) - 2.D0/3.D0
     & *Hr1(0)**3*CA**2*opx**(-1) - 4.D0/3.D0*Hr1(0)**3*CA**2*x*
     & opx**(-1) - 4.D0/3.D0*Hr1(0)**3*CA**2*x + 2.D0/3.D0*Hr1(0)**3*
     & CA**2*x**2*opx**(-1) - 2.D0/3.D0*Hr1(0)**3*CA**2*x**2 + 4.D0/3.D0
     & *Hr1(0)**3*CA**2*x**3*opx**(-1) + 2.D0/3.D0*Hr1(0)**3*CA**2*x**3
     &  - 2.D0/3.D0*Hr1(0)**3*CA**2*x**4*opx**(-1) - 4*Hr1(0)**2*Hr1(1)
     & *CA**2*x**(-1) - 4*Hr1(0)**2*Hr1(1)*CA**2*omx**(-1) + 8*Hr1(0)**
     & 2*Hr1(1)*CA**2 - 4*Hr1(0)**2*Hr1(1)*CA**2*x + 4*Hr1(0)**2*Hr1(1)
     & *CA**2*x**2 + 88.D0/3.D0*Hr1(0)*Hr1(1)*CA**2*x**(-1) - 32*Hr1(0)
     & *Hr1(1)*CA**2*x**(-1)*Lperp - 32*Hr1(0)*Hr1(1)*CA**2*omx**(-1)*
     & Lperp - 32*Hr1(0)*Hr1(1)*CA**2 + 64*Hr1(0)*Hr1(1)*CA**2*Lperp + 
     & 32*Hr1(0)*Hr1(1)*CA**2*x - 32*Hr1(0)*Hr1(1)*CA**2*x*Lperp - 88.D0
     & /3.D0*Hr1(0)*Hr1(1)*CA**2*x**2
      I2(gg,rglr) = I2(gg,rglr) + 32*Hr1(0)*Hr1(1)*CA**2*x**2*Lperp + 4
     & *Hr1(0)*Hr1(1)**2*CA**2*x**(-1) + 4*Hr1(0)*Hr1(1)**2*CA**2*
     & omx**(-1) - 8*Hr1(0)*Hr1(1)**2*CA**2 + 4*Hr1(0)*Hr1(1)**2*CA**2*
     & x - 4*Hr1(0)*Hr1(1)**2*CA**2*x**2 + 4.D0/3.D0*Hr1(1)*CA*TF*nf*x
     &  - 64*Hr1(1)*CA**2*x**(-1)*Lperp**2 + 128*Hr1(1)*CA**2*Lperp**2
     &  - 2.D0/3.D0*Hr1(1)*CA**2*x - 64*Hr1(1)*CA**2*x*Lperp**2 + 64*
     & Hr1(1)*CA**2*x**2*Lperp**2


      I2(qq,delt)= - 224.D0/27.D0*CF*TF*nf*LQ - 160.D0/9.D0*CF*TF*nf*LQ
     & *Lperp - 32.D0/3.D0*CF*TF*nf*LQ*Lperp**2 + 32*CF**2*LQ**2*
     & Lperp**2 - 16.D0/3.D0*CF**2*pisq*Lperp**2 + 808.D0/27.D0*CA*CF*
     & LQ + 536.D0/9.D0*CA*CF*LQ*Lperp + 88.D0/3.D0*CA*CF*LQ*Lperp**2
     &  - 8.D0/3.D0*CA*CF*pisq*LQ*Lperp - 28*CA*CF*zeta3*LQ

      I2(qq,plus)=224.D0/27.D0*CF*TF*nf + 160.D0/9.D0*CF*TF*nf*Lperp + 
     & 32.D0/3.D0*CF*TF*nf*Lperp**2 - 64*CF**2*LQ*Lperp**2 - 808.D0/27.D
     & 0*CA*CF - 536.D0/9.D0*CA*CF*Lperp - 88.D0/3.D0*CA*CF*Lperp**2 + 
     & 8.D0/3.D0*CA*CF*pisq*Lperp + 28*CA*CF*zeta3

      I2(qq,lpls)=64*CF**2*Lperp**2

      I2(qq,rglr)=344.D0/27.D0*CF*TF*x**(-1) - 208.D0/9.D0*CF*TF*
     & x**(-1)*Lperp + 32.D0/3.D0*CF*TF*x**(-1)*Lperp**2 - 8.D0/9.D0*CF
     & *TF*x**(-1)*pisq - 70.D0/3.D0*CF*TF + 32*CF*TF*Lperp + 8*CF*TF*
     & Lperp**2 + 4.D0/3.D0*CF*TF*pisq + 62.D0/3.D0*CF*TF*x - 48*CF*TF*
     & x*Lperp - 8*CF*TF*x*Lperp**2 - 4.D0/3.D0*CF*TF*x*pisq - 272.D0/
     & 27.D0*CF*TF*x**2 + 352.D0/9.D0*CF*TF*x**2*Lperp - 32.D0/3.D0*CF*
     & TF*x**2*Lperp**2 + 8.D0/9.D0*CF*TF*x**2*pisq - 148.D0/27.D0*CF*
     & TF*nf - 32.D0/9.D0*CF*TF*nf*Lperp - 16.D0/3.D0*CF*TF*nf*Lperp**2
     &  - 76.D0/27.D0*CF*TF*nf*x - 128.D0/9.D0*CF*TF*nf*x*Lperp - 16.D0/
     & 3.D0*CF*TF*nf*x*Lperp**2 - 22*CF**2 + 56*CF**2*Lperp - 16*CF**2*
     & Lperp**2 + 16*CF**2*LQ*Lperp + 32*CF**2*LQ*Lperp**2 + 4.D0/3.D0*
     & CF**2*pisq + 48*CF**2*zeta3*omx**(-1) - 24*CF**2*zeta3 + 22*
     & CF**2*x - 56*CF**2*x*Lperp + 16*CF**2*x*Lperp**2 - 16*CF**2*x*LQ
     & *Lperp + 32*CF**2*x*LQ*Lperp**2 - 4.D0/3.D0*CF**2*x*pisq - 24*
     & CF**2*x*zeta3 + 800.D0/27.D0*CA*CF - 80.D0/9.D0*CA*CF*Lperp + 44.
     & D0/3.D0*CA*CF*Lperp**2
      I2(qq,rglr) = I2(qq,rglr) - CA*CF*pisq - 4.D0/3.D0*CA*CF*pisq*
     & Lperp - 24*CA*CF*zeta3*omx**(-1) - 2*CA*CF*zeta3 + 8.D0/27.D0*CA
     & *CF*x + 616.D0/9.D0*CA*CF*x*Lperp + 44.D0/3.D0*CA*CF*x*Lperp**2
     &  + CA*CF*x*pisq - 4.D0/3.D0*CA*CF*x*pisq*Lperp - 2*CA*CF*x*zeta3
     &  - 40*Hr3(0,0,1)*CF**2*omx**(-1) + 20*Hr3(0,0,1)*CF**2 + 20*Hr3(
     & 0,0,1)*CF**2*x + 16*Hr3(0,0,1)*CA*CF*omx**(-1) - 8*Hr3(0,0,1)*CA
     & *CF - 8*Hr3(0,0,1)*CA*CF*x - 8*Hr3(0,1,1)*CF**2*omx**(-1) + 4*
     & Hr3(0,1,1)*CF**2 + 4*Hr3(0,1,1)*CF**2*x + 8*Hr3(0,1,1)*CA*CF*
     & omx**(-1) - 4*Hr3(0,1,1)*CA*CF - 4*Hr3(0,1,1)*CA*CF*x + 16.D0/3.D
     & 0*Hr2(0,1)*CF*TF*x**(-1) - 8*Hr2(0,1)*CF*TF + 8*Hr2(0,1)*CF*TF*x
     &  - 16.D0/3.D0*Hr2(0,1)*CF*TF*x**2 - 8*Hr2(0,1)*CF**2 + 8*Hr2(0,1
     & )*CF**2*x + 4*Hr2(0,1)*CA*CF - 4*Hr2(0,1)*CA*CF*x + 24*Hr2(0,1)*
     & Hr1(0)*CF**2*omx**(-1) - 12*Hr2(0,1)*Hr1(0)*CF**2 - 12*Hr2(0,1)*
     & Hr1(0)*CF**2*x - 8*Hr2(0,1)*Hr1(0)*CA*CF*omx**(-1) + 4*Hr2(0,1)*
     & Hr1(0)*CA*CF + 4*Hr2(0,1)*Hr1(0)*CA*CF*x + 28.D0/3.D0*Hr1(0)*CF*
     & TF
      I2(qq,rglr) = I2(qq,rglr) - 8*Hr1(0)*CF*TF*Lperp + 16*Hr1(0)*CF*
     & TF*Lperp**2 - 40.D0/3.D0*Hr1(0)*CF*TF*x - 24*Hr1(0)*CF*TF*x*
     & Lperp + 16*Hr1(0)*CF*TF*x*Lperp**2 + 128.D0/9.D0*Hr1(0)*CF*TF*
     & x**2 - 64.D0/3.D0*Hr1(0)*CF*TF*x**2*Lperp + 40.D0/9.D0*Hr1(0)*CF
     & *TF*nf*omx**(-1) + 32.D0/3.D0*Hr1(0)*CF*TF*nf*omx**(-1)*Lperp - 
     & 20.D0/9.D0*Hr1(0)*CF*TF*nf - 16.D0/3.D0*Hr1(0)*CF*TF*nf*Lperp - 
     & 20.D0/9.D0*Hr1(0)*CF*TF*nf*x - 16.D0/3.D0*Hr1(0)*CF*TF*nf*x*
     & Lperp + 16*Hr1(0)*CF**2*omx**(-1) + 24*Hr1(0)*CF**2*omx**(-1)*
     & Lperp - 32*Hr1(0)*CF**2*omx**(-1)*Lperp**2 - 6*Hr1(0)*CF**2 + 8*
     & Hr1(0)*CF**2*Lperp + 24*Hr1(0)*CF**2*Lperp**2 - 32*Hr1(0)*CF**2*
     & x + 8*Hr1(0)*CF**2*x*Lperp + 24*Hr1(0)*CF**2*x*Lperp**2 - 152.D0/
     & 9.D0*Hr1(0)*CA*CF*omx**(-1) - 88.D0/3.D0*Hr1(0)*CA*CF*omx**(-1)*
     & Lperp + 94.D0/9.D0*Hr1(0)*CA*CF + 20.D0/3.D0*Hr1(0)*CA*CF*Lperp
     &  + 166.D0/9.D0*Hr1(0)*CA*CF*x + 20.D0/3.D0*Hr1(0)*CA*CF*x*Lperp
     &  - Hr1(0)**2*CF*TF + 8*Hr1(0)**2*CF*TF*Lperp - Hr1(0)**2*CF*TF*x
     &  + 8*Hr1(0)**2*CF*TF*x*Lperp
      I2(qq,rglr) = I2(qq,rglr) - 8.D0/3.D0*Hr1(0)**2*CF*TF*x**2 + 4.D0/
     & 3.D0*Hr1(0)**2*CF*TF*nf*omx**(-1) - 2.D0/3.D0*Hr1(0)**2*CF*TF*nf
     &  - 2.D0/3.D0*Hr1(0)**2*CF*TF*nf*x + 3*Hr1(0)**2*CF**2*omx**(-1)
     &  + 4*Hr1(0)**2*CF**2*Lperp + 2*Hr1(0)**2*CF**2*x + 4*Hr1(0)**2*
     & CF**2*x*Lperp - 11.D0/3.D0*Hr1(0)**2*CA*CF*omx**(-1) - 8*Hr1(0)
     & **2*CA*CF*omx**(-1)*Lperp + 11.D0/6.D0*Hr1(0)**2*CA*CF + 4*Hr1(0
     & )**2*CA*CF*Lperp - 1.D0/6.D0*Hr1(0)**2*CA*CF*x + 4*Hr1(0)**2*CA*
     & CF*x*Lperp + 2.D0/3.D0*Hr1(0)**3*CF*TF + 2.D0/3.D0*Hr1(0)**3*CF*
     & TF*x + 1.D0/3.D0*Hr1(0)**3*CF**2 + 1.D0/3.D0*Hr1(0)**3*CF**2*x
     &  - 2.D0/3.D0*Hr1(0)**3*CA*CF*omx**(-1) + 1.D0/3.D0*Hr1(0)**3*CA*
     & CF + 1.D0/3.D0*Hr1(0)**3*CA*CF*x - 4*Hr1(0)**2*Hr1(1)*CF**2*
     & omx**(-1) + 2*Hr1(0)**2*Hr1(1)*CF**2 + 2*Hr1(0)**2*Hr1(1)*CF**2*
     & x - 16.D0/3.D0*Hr1(0)*Hr1(1)*CF*TF*x**(-1) + 8*Hr1(0)*Hr1(1)*CF*
     & TF - 8*Hr1(0)*Hr1(1)*CF*TF*x + 16.D0/3.D0*Hr1(0)*Hr1(1)*CF*TF*
     & x**2 - 32*Hr1(0)*Hr1(1)*CF**2*omx**(-1)*Lperp + 12*Hr1(0)*Hr1(1)
     & *CF**2
      I2(qq,rglr) = I2(qq,rglr) + 16*Hr1(0)*Hr1(1)*CF**2*Lperp - 12*
     & Hr1(0)*Hr1(1)*CF**2*x + 16*Hr1(0)*Hr1(1)*CF**2*x*Lperp - 4*Hr1(0
     & )*Hr1(1)*CA*CF + 4*Hr1(0)*Hr1(1)*CA*CF*x + 8*Hr1(0)*Hr1(1)**2*
     & CF**2*omx**(-1) - 4*Hr1(0)*Hr1(1)**2*CF**2 - 4*Hr1(0)*Hr1(1)**2*
     & CF**2*x - 4*Hr1(0)*Hr1(1)**2*CA*CF*omx**(-1) + 2*Hr1(0)*Hr1(1)**
     & 2*CA*CF + 2*Hr1(0)*Hr1(1)**2*CA*CF*x + 16*Hr1(1)*CF**2*Lperp + 
     & 32*Hr1(1)*CF**2*Lperp**2 + 2*Hr1(1)*CF**2*x - 16*Hr1(1)*CF**2*x*
     & Lperp + 32*Hr1(1)*CF**2*x*Lperp**2 - 2*Hr1(1)*CA*CF*x


      I2(qbq,delt)= 0

      I2(qbq,plus)= 0

      I2(qbq,lpls)= 0

      I2(qbq,rglr)=344.D0/27.D0*CF*TF*x**(-1) - 208.D0/9.D0*CF*TF*
     & x**(-1)*Lperp + 32.D0/3.D0*CF*TF*x**(-1)*Lperp**2 - 8.D0/9.D0*CF
     & *TF*x**(-1)*pisq - 70.D0/3.D0*CF*TF + 32*CF*TF*Lperp + 8*CF*TF*
     & Lperp**2 + 4.D0/3.D0*CF*TF*pisq + 62.D0/3.D0*CF*TF*x - 48*CF*TF*
     & x*Lperp - 8*CF*TF*x*Lperp**2 - 4.D0/3.D0*CF*TF*x*pisq - 272.D0/
     & 27.D0*CF*TF*x**2 + 352.D0/9.D0*CF*TF*x**2*Lperp - 32.D0/3.D0*CF*
     & TF*x**2*Lperp**2 + 8.D0/9.D0*CF*TF*x**2*pisq + 30*CF**2*
     & opx**(-1) - 32*CF**2*opx**(-1)*Lperp + 8.D0/3.D0*CF**2*pisq*
     & opx**(-1)*Lperp - 2*CF**2*pisq + 4*CF**2*zeta3*opx**(-1) + 2.D0/
     & 3.D0*CF**2*x*pisq - 30*CF**2*x**2*opx**(-1) + 32*CF**2*x**2*
     & opx**(-1)*Lperp + 8.D0/3.D0*CF**2*x**2*pisq*opx**(-1)*Lperp + 4*
     & CF**2*x**2*zeta3*opx**(-1) - 15*CA*CF*opx**(-1) + 16*CA*CF*
     & opx**(-1)*Lperp - 4.D0/3.D0*CA*CF*pisq*opx**(-1)*Lperp + CA*CF*
     & pisq - 2*CA*CF*zeta3*opx**(-1) - 1.D0/3.D0*CA*CF*x*pisq + 15*CA*
     & CF*x**2*opx**(-1) - 16*CA*CF*x**2*opx**(-1)*Lperp - 4.D0/3.D0*CA
     & *CF*x**2*pisq*opx**(-1)*Lperp
      I2(qbq,rglr) = I2(qbq,rglr) - 2*CA*CF*x**2*zeta3*opx**(-1) - 16*
     & Hr3(0,-1,-1)*CF**2*opx**(-1) - 16*Hr3(0,-1,-1)*CF**2*x**2*
     & opx**(-1) + 8*Hr3(0,-1,-1)*CA*CF*opx**(-1) + 8*Hr3(0,-1,-1)*CA*
     & CF*x**2*opx**(-1) - 24*Hr3(0,0,-1)*CF**2*opx**(-1) - 24*Hr3(0,0,
     & -1)*CF**2*x**2*opx**(-1) + 12*Hr3(0,0,-1)*CA*CF*opx**(-1) + 12*
     & Hr3(0,0,-1)*CA*CF*x**2*opx**(-1) + 16*Hr3(0,0,1)*CF**2*opx**(-1)
     &  + 16*Hr3(0,0,1)*CF**2*x**2*opx**(-1) - 8*Hr3(0,0,1)*CA*CF*
     & opx**(-1) - 8*Hr3(0,0,1)*CA*CF*x**2*opx**(-1) - 32*Hr2(0,-1)*
     & CF**2*opx**(-1)*Lperp + 8*Hr2(0,-1)*CF**2 + 8*Hr2(0,-1)*CF**2*x
     &  - 32*Hr2(0,-1)*CF**2*x**2*opx**(-1)*Lperp + 16*Hr2(0,-1)*CA*CF*
     & opx**(-1)*Lperp - 4*Hr2(0,-1)*CA*CF - 4*Hr2(0,-1)*CA*CF*x + 16*
     & Hr2(0,-1)*CA*CF*x**2*opx**(-1)*Lperp + 16*Hr2(0,-1)*Hr1(-1)*
     & CF**2*opx**(-1) + 16*Hr2(0,-1)*Hr1(-1)*CF**2*x**2*opx**(-1) - 8*
     & Hr2(0,-1)*Hr1(-1)*CA*CF*opx**(-1) - 8*Hr2(0,-1)*Hr1(-1)*CA*CF*
     & x**2*opx**(-1) + 8*Hr2(0,-1)*Hr1(0)*CF**2*opx**(-1) + 8*Hr2(0,-1
     & )*Hr1(0)*CF**2*x**2*opx**(-1)
      I2(qbq,rglr) = I2(qbq,rglr) - 4*Hr2(0,-1)*Hr1(0)*CA*CF*opx**(-1)
     &  - 4*Hr2(0,-1)*Hr1(0)*CA*CF*x**2*opx**(-1) + 16.D0/3.D0*Hr2(0,1)
     & *CF*TF*x**(-1) - 8*Hr2(0,1)*CF*TF + 8*Hr2(0,1)*CF*TF*x - 16.D0/3.
     & D0*Hr2(0,1)*CF*TF*x**2 + 8*Hr2(0,1)*CF**2 - 8*Hr2(0,1)*CF**2*x
     &  - 4*Hr2(0,1)*CA*CF + 4*Hr2(0,1)*CA*CF*x - 8*Hr2(0,1)*Hr1(0)*
     & CF**2*opx**(-1) - 8*Hr2(0,1)*Hr1(0)*CF**2*x**2*opx**(-1) + 4*
     & Hr2(0,1)*Hr1(0)*CA*CF*opx**(-1) + 4*Hr2(0,1)*Hr1(0)*CA*CF*x**2*
     & opx**(-1) - 4.D0/3.D0*Hr1(-1)*CF**2*pisq*opx**(-1) - 4.D0/3.D0*
     & Hr1(-1)*CF**2*x**2*pisq*opx**(-1) + 2.D0/3.D0*Hr1(-1)*CA*CF*pisq
     & *opx**(-1) + 2.D0/3.D0*Hr1(-1)*CA*CF*x**2*pisq*opx**(-1) - 8*
     & Hr1(-1)**2*Hr1(0)*CF**2*opx**(-1) - 8*Hr1(-1)**2*Hr1(0)*CF**2*
     & x**2*opx**(-1) + 4*Hr1(-1)**2*Hr1(0)*CA*CF*opx**(-1) + 4*Hr1(-1)
     & **2*Hr1(0)*CA*CF*x**2*opx**(-1) + 32*Hr1(-1)*Hr1(0)*CF**2*
     & opx**(-1)*Lperp - 8*Hr1(-1)*Hr1(0)*CF**2 - 8*Hr1(-1)*Hr1(0)*
     & CF**2*x + 32*Hr1(-1)*Hr1(0)*CF**2*x**2*opx**(-1)*Lperp - 16*Hr1(
     & -1)*Hr1(0)*CA*CF*opx**(-1)*Lperp
      I2(qbq,rglr) = I2(qbq,rglr) + 4*Hr1(-1)*Hr1(0)*CA*CF + 4*Hr1(-1)*
     & Hr1(0)*CA*CF*x - 16*Hr1(-1)*Hr1(0)*CA*CF*x**2*opx**(-1)*Lperp + 
     & 4*Hr1(-1)*Hr1(0)**2*CF**2*opx**(-1) + 4*Hr1(-1)*Hr1(0)**2*CF**2*
     & x**2*opx**(-1) - 2*Hr1(-1)*Hr1(0)**2*CA*CF*opx**(-1) - 2*Hr1(-1)
     & *Hr1(0)**2*CA*CF*x**2*opx**(-1) + 28.D0/3.D0*Hr1(0)*CF*TF - 8*
     & Hr1(0)*CF*TF*Lperp + 16*Hr1(0)*CF*TF*Lperp**2 - 40.D0/3.D0*Hr1(0
     & )*CF*TF*x - 24*Hr1(0)*CF*TF*x*Lperp + 16*Hr1(0)*CF*TF*x*Lperp**2
     &  + 128.D0/9.D0*Hr1(0)*CF*TF*x**2 - 64.D0/3.D0*Hr1(0)*CF*TF*x**2*
     & Lperp + 6*Hr1(0)*CF**2 - 16*Hr1(0)*CF**2*Lperp + 22*Hr1(0)*CF**2
     & *x - 16*Hr1(0)*CF**2*x*Lperp - 3*Hr1(0)*CA*CF + 8*Hr1(0)*CA*CF*
     & Lperp - 11*Hr1(0)*CA*CF*x + 8*Hr1(0)*CA*CF*x*Lperp - Hr1(0)**2*
     & CF*TF + 8*Hr1(0)**2*CF*TF*Lperp - Hr1(0)**2*CF*TF*x + 8*Hr1(0)**
     & 2*CF*TF*x*Lperp - 8.D0/3.D0*Hr1(0)**2*CF*TF*x**2 - 8*Hr1(0)**2*
     & CF**2*opx**(-1)*Lperp - 8*Hr1(0)**2*CF**2*x**2*opx**(-1)*Lperp
     &  + 4*Hr1(0)**2*CA*CF*opx**(-1)*Lperp + 4*Hr1(0)**2*CA*CF*x**2*
     & opx**(-1)*Lperp
      I2(qbq,rglr) = I2(qbq,rglr) + 2.D0/3.D0*Hr1(0)**3*CF*TF + 2.D0/3.D
     & 0*Hr1(0)**3*CF*TF*x - 2.D0/3.D0*Hr1(0)**3*CF**2*opx**(-1) - 2.D0/
     & 3.D0*Hr1(0)**3*CF**2*x**2*opx**(-1) + 1.D0/3.D0*Hr1(0)**3*CA*CF*
     & opx**(-1) + 1.D0/3.D0*Hr1(0)**3*CA*CF*x**2*opx**(-1) - 16.D0/3.D0
     & *Hr1(0)*Hr1(1)*CF*TF*x**(-1) + 8*Hr1(0)*Hr1(1)*CF*TF - 8*Hr1(0)*
     & Hr1(1)*CF*TF*x + 16.D0/3.D0*Hr1(0)*Hr1(1)*CF*TF*x**2 - 8*Hr1(0)*
     & Hr1(1)*CF**2 + 8*Hr1(0)*Hr1(1)*CF**2*x + 4*Hr1(0)*Hr1(1)*CA*CF
     &  - 4*Hr1(0)*Hr1(1)*CA*CF*x


      I2(qpq,delt)= 0

      I2(qpq,plus)= 0

      I2(qpq,lpls)= 0

      I2(qpq,rglr)=344.D0/27.D0*CF*TF*x**(-1) - 208.D0/9.D0*CF*TF*
     & x**(-1)*Lperp + 32.D0/3.D0*CF*TF*x**(-1)*Lperp**2 - 8.D0/9.D0*CF
     & *TF*x**(-1)*pisq - 70.D0/3.D0*CF*TF + 32*CF*TF*Lperp + 8*CF*TF*
     & Lperp**2 + 4.D0/3.D0*CF*TF*pisq + 62.D0/3.D0*CF*TF*x - 48*CF*TF*
     & x*Lperp - 8*CF*TF*x*Lperp**2 - 4.D0/3.D0*CF*TF*x*pisq - 272.D0/
     & 27.D0*CF*TF*x**2 + 352.D0/9.D0*CF*TF*x**2*Lperp - 32.D0/3.D0*CF*
     & TF*x**2*Lperp**2 + 8.D0/9.D0*CF*TF*x**2*pisq + 16.D0/3.D0*Hr2(0,
     & 1)*CF*TF*x**(-1) - 8*Hr2(0,1)*CF*TF + 8*Hr2(0,1)*CF*TF*x - 16.D0/
     & 3.D0*Hr2(0,1)*CF*TF*x**2 + 28.D0/3.D0*Hr1(0)*CF*TF - 8*Hr1(0)*CF
     & *TF*Lperp + 16*Hr1(0)*CF*TF*Lperp**2 - 40.D0/3.D0*Hr1(0)*CF*TF*x
     &  - 24*Hr1(0)*CF*TF*x*Lperp + 16*Hr1(0)*CF*TF*x*Lperp**2 + 128.D0/
     & 9.D0*Hr1(0)*CF*TF*x**2 - 64.D0/3.D0*Hr1(0)*CF*TF*x**2*Lperp - 
     & Hr1(0)**2*CF*TF + 8*Hr1(0)**2*CF*TF*Lperp - Hr1(0)**2*CF*TF*x + 
     & 8*Hr1(0)**2*CF*TF*x*Lperp - 8.D0/3.D0*Hr1(0)**2*CF*TF*x**2 + 2.D0
     & /3.D0*Hr1(0)**3*CF*TF + 2.D0/3.D0*Hr1(0)**3*CF*TF*x - 16.D0/3.D0
     & *Hr1(0)*Hr1(1)*CF*TF*x**(-1)
      I2(qpq,rglr) = I2(qpq,rglr) + 8*Hr1(0)*Hr1(1)*CF*TF - 8*Hr1(0)*
     & Hr1(1)*CF*TF*x + 16.D0/3.D0*Hr1(0)*Hr1(1)*CF*TF*x**2


      I2(qg,delt)= 0

      I2(qg,plus)= 0

      I2(qg,lpls)= 0

      I2(qg,rglr)= - 13*CF*TF - 40*CF*TF*Lperp - 28*CF*TF*Lperp**2 - 32
     & *CF*TF*LQ*Lperp**2 + 8.D0/3.D0*CF*TF*pisq*Lperp + 28*CF*TF*zeta3
     &  + 75*CF*TF*x + 132*CF*TF*x*Lperp + 64*CF*TF*x*Lperp**2 + 32*CF*
     & TF*x*LQ*Lperp + 64*CF*TF*x*LQ*Lperp**2 - 4.D0/3.D0*CF*TF*x*pisq
     &  - 16.D0/3.D0*CF*TF*x*pisq*Lperp - 56*CF*TF*x*zeta3 - 72*CF*TF*
     & x**2 - 112*CF*TF*x**2*Lperp - 48*CF*TF*x**2*Lperp**2 - 32*CF*TF*
     & x**2*LQ*Lperp - 64*CF*TF*x**2*LQ*Lperp**2 + 4.D0/3.D0*CF*TF*x**2
     & *pisq + 16.D0/3.D0*CF*TF*x**2*pisq*Lperp + 56*CF*TF*x**2*zeta3
     &  + 344.D0/27.D0*CA*TF*x**(-1) - 208.D0/9.D0*CA*TF*x**(-1)*Lperp
     &  + 32.D0/3.D0*CA*TF*x**(-1)*Lperp**2 - 8.D0/9.D0*CA*TF*x**(-1)*
     & pisq - 70.D0/3.D0*CA*TF + 32*CA*TF*Lperp + 8*CA*TF*Lperp**2 + 4.D
     & 0/3.D0*CA*TF*pisq + 86.D0/3.D0*CA*TF*x - 120*CA*TF*x*Lperp + 64*
     & CA*TF*x*Lperp**2 - 4*CA*TF*x*pisq + 16.D0/3.D0*CA*TF*x*pisq*
     & Lperp + 8*CA*TF*x*zeta3 - 596.D0/27.D0*CA*TF*x**2 + 928.D0/9.D0*
     & CA*TF*x**2*Lperp - 248.D0/3.D0*CA*TF*x**2*Lperp**2 + 44.D0/9.D0*
     & CA*TF*x**2*pisq
      I2(qg,rglr) = I2(qg,rglr) - 8*Hr3(0,-1,-1)*CA*TF - 16*Hr3(0,-1,-1
     & )*CA*TF*x - 16*Hr3(0,-1,-1)*CA*TF*x**2 - 12*Hr3(0,0,-1)*CA*TF - 
     & 24*Hr3(0,0,-1)*CA*TF*x - 24*Hr3(0,0,-1)*CA*TF*x**2 - 4*Hr3(0,0,1
     & )*CF*TF + 8*Hr3(0,0,1)*CF*TF*x - 8*Hr3(0,0,1)*CF*TF*x**2 + 32*
     & Hr3(0,0,1)*CA*TF*x + 4*Hr3(0,1,1)*CF*TF - 8*Hr3(0,1,1)*CF*TF*x
     &  + 8*Hr3(0,1,1)*CF*TF*x**2 - 4*Hr3(0,1,1)*CA*TF + 8*Hr3(0,1,1)*
     & CA*TF*x - 8*Hr3(0,1,1)*CA*TF*x**2 - 16*Hr2(0,-1)*CA*TF*Lperp - 8
     & *Hr2(0,-1)*CA*TF*x - 32*Hr2(0,-1)*CA*TF*x*Lperp - 8*Hr2(0,-1)*CA
     & *TF*x**2 - 32*Hr2(0,-1)*CA*TF*x**2*Lperp + 8*Hr2(0,-1)*Hr1(-1)*
     & CA*TF + 16*Hr2(0,-1)*Hr1(-1)*CA*TF*x + 16*Hr2(0,-1)*Hr1(-1)*CA*
     & TF*x**2 + 4*Hr2(0,-1)*Hr1(0)*CA*TF + 8*Hr2(0,-1)*Hr1(0)*CA*TF*x
     &  + 8*Hr2(0,-1)*Hr1(0)*CA*TF*x**2 + 16.D0/3.D0*Hr2(0,1)*CA*TF*
     & x**(-1) - 8*Hr2(0,1)*CA*TF + 32*Hr2(0,1)*CA*TF*x - 88.D0/3.D0*
     & Hr2(0,1)*CA*TF*x**2 + 4*Hr2(0,1)*Hr1(0)*CF*TF - 8*Hr2(0,1)*Hr1(0
     & )*CF*TF*x + 8*Hr2(0,1)*Hr1(0)*CF*TF*x**2 - 16*Hr2(0,1)*Hr1(0)*CA
     & *TF*x
      I2(qg,rglr) = I2(qg,rglr) - 2.D0/3.D0*Hr1(-1)*CA*TF*pisq - 4.D0/3.
     & D0*Hr1(-1)*CA*TF*x*pisq - 4.D0/3.D0*Hr1(-1)*CA*TF*x**2*pisq - 4*
     & Hr1(-1)**2*Hr1(0)*CA*TF - 8*Hr1(-1)**2*Hr1(0)*CA*TF*x - 8*Hr1(-1
     & )**2*Hr1(0)*CA*TF*x**2 + 16*Hr1(-1)*Hr1(0)*CA*TF*Lperp + 8*Hr1(
     & -1)*Hr1(0)*CA*TF*x + 32*Hr1(-1)*Hr1(0)*CA*TF*x*Lperp + 8*Hr1(-1)
     & *Hr1(0)*CA*TF*x**2 + 32*Hr1(-1)*Hr1(0)*CA*TF*x**2*Lperp + 2*Hr1(
     & -1)*Hr1(0)**2*CA*TF + 4*Hr1(-1)*Hr1(0)**2*CA*TF*x + 4*Hr1(-1)*
     & Hr1(0)**2*CA*TF*x**2 + 8*Hr1(0)*CF*TF - 4*Hr1(0)*CF*TF*Lperp - 8
     & *Hr1(0)*CF*TF*Lperp**2 + 15*Hr1(0)*CF*TF*x + 32*Hr1(0)*CF*TF*x*
     & Lperp + 16*Hr1(0)*CF*TF*x*Lperp**2 - 8*Hr1(0)*CF*TF*x**2 - 32*
     & Hr1(0)*CF*TF*x**2*Lperp - 32*Hr1(0)*CF*TF*x**2*Lperp**2 + 28.D0/
     & 3.D0*Hr1(0)*CA*TF - 8*Hr1(0)*CA*TF*Lperp + 16*Hr1(0)*CA*TF*
     & Lperp**2 - 40.D0/3.D0*Hr1(0)*CA*TF*x + 64*Hr1(0)*CA*TF*x*
     & Lperp**2 + 272.D0/9.D0*Hr1(0)*CA*TF*x**2 - 352.D0/3.D0*Hr1(0)*CA
     & *TF*x**2*Lperp + 1.D0/2.D0*Hr1(0)**2*CF*TF - 4*Hr1(0)**2*CF*TF*
     & Lperp
      I2(qg,rglr) = I2(qg,rglr) + 6*Hr1(0)**2*CF*TF*x + 8*Hr1(0)**2*CF*
     & TF*x*Lperp - 4*Hr1(0)**2*CF*TF*x**2 - 16*Hr1(0)**2*CF*TF*x**2*
     & Lperp - Hr1(0)**2*CA*TF + 8*Hr1(0)**2*CA*TF*Lperp + 4*Hr1(0)**2*
     & CA*TF*x + 16*Hr1(0)**2*CA*TF*x*Lperp - 44.D0/3.D0*Hr1(0)**2*CA*
     & TF*x**2 - 1.D0/3.D0*Hr1(0)**3*CF*TF + 2.D0/3.D0*Hr1(0)**3*CF*TF*
     & x - 4.D0/3.D0*Hr1(0)**3*CF*TF*x**2 + 2.D0/3.D0*Hr1(0)**3*CA*TF
     &  + 4.D0/3.D0*Hr1(0)**3*CA*TF*x - 2*Hr1(0)**2*Hr1(1)*CF*TF + 4*
     & Hr1(0)**2*Hr1(1)*CF*TF*x - 4*Hr1(0)**2*Hr1(1)*CF*TF*x**2 - 16*
     & Hr1(0)*Hr1(1)*CF*TF*Lperp + 8*Hr1(0)*Hr1(1)*CF*TF*x + 32*Hr1(0)*
     & Hr1(1)*CF*TF*x*Lperp - 8*Hr1(0)*Hr1(1)*CF*TF*x**2 - 32*Hr1(0)*
     & Hr1(1)*CF*TF*x**2*Lperp - 16.D0/3.D0*Hr1(0)*Hr1(1)*CA*TF*x**(-1)
     &  + 8*Hr1(0)*Hr1(1)*CA*TF - 32*Hr1(0)*Hr1(1)*CA*TF*x + 88.D0/3.D0
     & *Hr1(0)*Hr1(1)*CA*TF*x**2 + 2*Hr1(0)*Hr1(1)**2*CA*TF - 4*Hr1(0)*
     & Hr1(1)**2*CA*TF*x + 4*Hr1(0)*Hr1(1)**2*CA*TF*x**2 - 16*Hr1(1)*CF
     & *TF*Lperp**2 + 6*Hr1(1)*CF*TF*x + 32*Hr1(1)*CF*TF*x*Lperp + 32*
     & Hr1(1)*CF*TF*x*Lperp**2
      I2(qg,rglr) = I2(qg,rglr) - 8*Hr1(1)*CF*TF*x**2 - 32*Hr1(1)*CF*TF
     & *x**2*Lperp - 32*Hr1(1)*CF*TF*x**2*Lperp**2 - 16*Hr1(1)*CA*TF*
     & Lperp**2 - 6*Hr1(1)*CA*TF*x + 32*Hr1(1)*CA*TF*x*Lperp**2 + 8*
     & Hr1(1)*CA*TF*x**2 - 32*Hr1(1)*CA*TF*x**2*Lperp**2 - 8*Hr1(1)**2*
     & CF*TF*Lperp + 4*Hr1(1)**2*CF*TF*x + 16*Hr1(1)**2*CF*TF*x*Lperp
     &  - 4*Hr1(1)**2*CF*TF*x**2 - 16*Hr1(1)**2*CF*TF*x**2*Lperp + 8*
     & Hr1(1)**2*CA*TF*Lperp - 4*Hr1(1)**2*CA*TF*x - 16*Hr1(1)**2*CA*TF
     & *x*Lperp + 4*Hr1(1)**2*CA*TF*x**2 + 16*Hr1(1)**2*CA*TF*x**2*
     & Lperp + 2.D0/3.D0*Hr1(1)**3*CF*TF - 4.D0/3.D0*Hr1(1)**3*CF*TF*x
     &  + 4.D0/3.D0*Hr1(1)**3*CF*TF*x**2 - 2.D0/3.D0*Hr1(1)**3*CA*TF + 
     & 4.D0/3.D0*Hr1(1)**3*CA*TF*x - 4.D0/3.D0*Hr1(1)**3*CA*TF*x**2


      I2(gq,delt)= 0

      I2(gq,plus)= 0

      I2(gq,lpls)= 0

      I2(gq,rglr)=448.D0/27.D0*CF*TF*nf*x**(-1) + 320.D0/9.D0*CF*TF*nf*
     & x**(-1)*Lperp + 64.D0/3.D0*CF*TF*nf*x**(-1)*Lperp**2 - 448.D0/27.
     & D0*CF*TF*nf - 320.D0/9.D0*CF*TF*nf*Lperp - 64.D0/3.D0*CF*TF*nf*
     & Lperp**2 + 104.D0/27.D0*CF*TF*nf*x + 160.D0/9.D0*CF*TF*nf*x*
     & Lperp + 32.D0/3.D0*CF*TF*nf*x*Lperp**2 + 10*CF**2 + 12*CF**2*
     & Lperp + 16*CF**2*Lperp**2 - CF**2*x + 24*CF**2*x*Lperp - 4*CF**2
     & *x*Lperp**2 - 3160.D0/27.D0*CA*CF*x**(-1) - 8*CA*CF*x**(-1)*
     & Lperp - 424.D0/3.D0*CA*CF*x**(-1)*Lperp**2 - 64*CA*CF*x**(-1)*LQ
     & *Lperp**2 + 44.D0/9.D0*CA*CF*x**(-1)*pisq + 48*CA*CF*x**(-1)*
     & zeta3 + 3164.D0/27.D0*CA*CF - 152.D0/9.D0*CA*CF*Lperp + 368.D0/3.
     & D0*CA*CF*Lperp**2 + 64*CA*CF*LQ*Lperp**2 - 16.D0/3.D0*CA*CF*pisq
     &  - 16.D0/3.D0*CA*CF*pisq*Lperp - 56*CA*CF*zeta3 - 1072.D0/27.D0*
     & CA*CF*x - 32.D0/9.D0*CA*CF*x*Lperp - 64.D0/3.D0*CA*CF*x*Lperp**2
     &  + 16*CA*CF*x*LQ*Lperp - 32*CA*CF*x*LQ*Lperp**2 + 4.D0/3.D0*CA*
     & CF*x*pisq + 24*CA*CF*x*zeta3 + 608.D0/27.D0*CA*CF*x**2 - 352.D0/
     & 9.D0*CA*CF*x**2*Lperp
      I2(gq,rglr) = I2(gq,rglr) + 32.D0/3.D0*CA*CF*x**2*Lperp**2 - 8.D0/
     & 9.D0*CA*CF*x**2*pisq + 16*Hr3(0,-1,-1)*CA*CF*x**(-1) + 16*Hr3(0,
     & -1,-1)*CA*CF + 8*Hr3(0,-1,-1)*CA*CF*x + 24*Hr3(0,0,-1)*CA*CF*
     & x**(-1) + 24*Hr3(0,0,-1)*CA*CF + 12*Hr3(0,0,-1)*CA*CF*x - 40*
     & Hr3(0,0,1)*CA*CF*x**(-1) + 8*Hr3(0,0,1)*CA*CF - 20*Hr3(0,0,1)*CA
     & *CF*x + 32*Hr2(0,-1)*CA*CF*x**(-1)*Lperp + 32*Hr2(0,-1)*CA*CF*
     & Lperp - 4*Hr2(0,-1)*CA*CF*x + 16*Hr2(0,-1)*CA*CF*x*Lperp - 16*
     & Hr2(0,-1)*Hr1(-1)*CA*CF*x**(-1) - 16*Hr2(0,-1)*Hr1(-1)*CA*CF - 8
     & *Hr2(0,-1)*Hr1(-1)*CA*CF*x - 8*Hr2(0,-1)*Hr1(0)*CA*CF*x**(-1) - 
     & 8*Hr2(0,-1)*Hr1(0)*CA*CF - 4*Hr2(0,-1)*Hr1(0)*CA*CF*x - 88.D0/3.D
     & 0*Hr2(0,1)*CA*CF*x**(-1) + 32*Hr2(0,1)*CA*CF - 8*Hr2(0,1)*CA*CF*
     & x + 16.D0/3.D0*Hr2(0,1)*CA*CF*x**2 + 24*Hr2(0,1)*Hr1(0)*CA*CF*
     & x**(-1) - 8*Hr2(0,1)*Hr1(0)*CA*CF + 12*Hr2(0,1)*Hr1(0)*CA*CF*x
     &  + 4.D0/3.D0*Hr1(-1)*CA*CF*x**(-1)*pisq + 4.D0/3.D0*Hr1(-1)*CA*
     & CF*pisq + 2.D0/3.D0*Hr1(-1)*CA*CF*x*pisq + 8*Hr1(-1)**2*Hr1(0)*
     & CA*CF*x**(-1)
      I2(gq,rglr) = I2(gq,rglr) + 8*Hr1(-1)**2*Hr1(0)*CA*CF + 4*Hr1(-1)
     & **2*Hr1(0)*CA*CF*x - 32*Hr1(-1)*Hr1(0)*CA*CF*x**(-1)*Lperp - 32*
     & Hr1(-1)*Hr1(0)*CA*CF*Lperp + 4*Hr1(-1)*Hr1(0)*CA*CF*x - 16*Hr1(
     & -1)*Hr1(0)*CA*CF*x*Lperp - 4*Hr1(-1)*Hr1(0)**2*CA*CF*x**(-1) - 4
     & *Hr1(-1)*Hr1(0)**2*CA*CF - 2*Hr1(-1)*Hr1(0)**2*CA*CF*x - 15*Hr1(
     & 0)*CF**2 - 16*Hr1(0)*CF**2*Lperp + 16*Hr1(0)*CF**2*Lperp**2 + 5*
     & Hr1(0)*CF**2*x - 20*Hr1(0)*CF**2*x*Lperp - 8*Hr1(0)*CF**2*x*
     & Lperp**2 - 32*Hr1(0)*CA*CF*x**(-1)*Lperp**2 - 166.D0/3.D0*Hr1(0)
     & *CA*CF + 96*Hr1(0)*CA*CF*Lperp - 32*Hr1(0)*CA*CF*Lperp**2 + 4.D0/
     & 3.D0*Hr1(0)*CA*CF*x + 40*Hr1(0)*CA*CF*x*Lperp - 32*Hr1(0)*CA*CF*
     & x*Lperp**2 - 176.D0/9.D0*Hr1(0)*CA*CF*x**2 + 64.D0/3.D0*Hr1(0)*
     & CA*CF*x**2*Lperp - 2*Hr1(0)**2*CF**2 + 8*Hr1(0)**2*CF**2*Lperp
     &  - 3.D0/2.D0*Hr1(0)**2*CF**2*x - 4*Hr1(0)**2*CF**2*x*Lperp + 12*
     & Hr1(0)**2*CA*CF - 16*Hr1(0)**2*CA*CF*Lperp + 3*Hr1(0)**2*CA*CF*x
     &  - 8*Hr1(0)**2*CA*CF*x*Lperp + 8.D0/3.D0*Hr1(0)**2*CA*CF*x**2 + 
     & 2.D0/3.D0*Hr1(0)**3*CF**2
      I2(gq,rglr) = I2(gq,rglr) - 1.D0/3.D0*Hr1(0)**3*CF**2*x - 4.D0/3.D
     & 0*Hr1(0)**3*CA*CF - 2.D0/3.D0*Hr1(0)**3*CA*CF*x - 4*Hr1(0)**2*
     & Hr1(1)*CA*CF*x**(-1) + 4*Hr1(0)**2*Hr1(1)*CA*CF - 2*Hr1(0)**2*
     & Hr1(1)*CA*CF*x + 88.D0/3.D0*Hr1(0)*Hr1(1)*CA*CF*x**(-1) - 32*
     & Hr1(0)*Hr1(1)*CA*CF*x**(-1)*Lperp - 32*Hr1(0)*Hr1(1)*CA*CF + 32*
     & Hr1(0)*Hr1(1)*CA*CF*Lperp + 12*Hr1(0)*Hr1(1)*CA*CF*x - 16*Hr1(0)
     & *Hr1(1)*CA*CF*x*Lperp - 16.D0/3.D0*Hr1(0)*Hr1(1)*CA*CF*x**2 + 4*
     & Hr1(0)*Hr1(1)**2*CA*CF*x**(-1) - 4*Hr1(0)*Hr1(1)**2*CA*CF + 2*
     & Hr1(0)*Hr1(1)**2*CA*CF*x - 80.D0/9.D0*Hr1(1)*CF*TF*nf*x**(-1) - 
     & 64.D0/3.D0*Hr1(1)*CF*TF*nf*x**(-1)*Lperp + 80.D0/9.D0*Hr1(1)*CF*
     & TF*nf + 64.D0/3.D0*Hr1(1)*CF*TF*nf*Lperp - 16.D0/9.D0*Hr1(1)*CF*
     & TF*nf*x - 32.D0/3.D0*Hr1(1)*CF*TF*nf*x*Lperp - 32*Hr1(1)*CF**2*
     & x**(-1) - 48*Hr1(1)*CF**2*x**(-1)*Lperp - 32*Hr1(1)*CF**2*
     & x**(-1)*Lperp**2 + 32*Hr1(1)*CF**2 + 48*Hr1(1)*CF**2*Lperp + 32*
     & Hr1(1)*CF**2*Lperp**2 - 10*Hr1(1)*CF**2*x - 24*Hr1(1)*CF**2*x*
     & Lperp
      I2(gq,rglr) = I2(gq,rglr) - 16*Hr1(1)*CF**2*x*Lperp**2 + 304.D0/9.
     & D0*Hr1(1)*CA*CF*x**(-1) + 176.D0/3.D0*Hr1(1)*CA*CF*x**(-1)*Lperp
     &  - 32*Hr1(1)*CA*CF*x**(-1)*Lperp**2 - 304.D0/9.D0*Hr1(1)*CA*CF
     &  - 176.D0/3.D0*Hr1(1)*CA*CF*Lperp + 32*Hr1(1)*CA*CF*Lperp**2 + 
     & 86.D0/9.D0*Hr1(1)*CA*CF*x + 136.D0/3.D0*Hr1(1)*CA*CF*x*Lperp - 
     & 16*Hr1(1)*CA*CF*x*Lperp**2 + 8.D0/3.D0*Hr1(1)**2*CF*TF*nf*
     & x**(-1) - 8.D0/3.D0*Hr1(1)**2*CF*TF*nf + 4.D0/3.D0*Hr1(1)**2*CF*
     & TF*nf*x + 6*Hr1(1)**2*CF**2*x**(-1) + 16*Hr1(1)**2*CF**2*x**(-1)
     & *Lperp - 6*Hr1(1)**2*CF**2 - 16*Hr1(1)**2*CF**2*Lperp + Hr1(1)**
     & 2*CF**2*x + 8*Hr1(1)**2*CF**2*x*Lperp - 22.D0/3.D0*Hr1(1)**2*CA*
     & CF*x**(-1) - 16*Hr1(1)**2*CA*CF*x**(-1)*Lperp + 22.D0/3.D0*Hr1(1
     & )**2*CA*CF + 16*Hr1(1)**2*CA*CF*Lperp - 5.D0/3.D0*Hr1(1)**2*CA*
     & CF*x - 8*Hr1(1)**2*CA*CF*x*Lperp - 4.D0/3.D0*Hr1(1)**3*CF**2*
     & x**(-1) + 4.D0/3.D0*Hr1(1)**3*CF**2 - 2.D0/3.D0*Hr1(1)**3*CF**2*
     & x + 4.D0/3.D0*Hr1(1)**3*CA*CF*x**(-1) - 4.D0/3.D0*Hr1(1)**3*CA*
     & CF
      I2(gq,rglr) = I2(gq,rglr) + 2.D0/3.D0*Hr1(1)**3*CA*CF*x


      I2(qbpq,delt)= 0

      I2(qbpq,plus)= 0

      I2(qbpq,lpls)= 0

      I2(qbpq,rglr)=344.D0/27.D0*CF*TF*x**(-1) - 208.D0/9.D0*CF*TF*
     & x**(-1)*Lperp + 32.D0/3.D0*CF*TF*x**(-1)*Lperp**2 - 8.D0/9.D0*CF
     & *TF*x**(-1)*pisq - 70.D0/3.D0*CF*TF + 32*CF*TF*Lperp + 8*CF*TF*
     & Lperp**2 + 4.D0/3.D0*CF*TF*pisq + 62.D0/3.D0*CF*TF*x - 48*CF*TF*
     & x*Lperp - 8*CF*TF*x*Lperp**2 - 4.D0/3.D0*CF*TF*x*pisq - 272.D0/
     & 27.D0*CF*TF*x**2 + 352.D0/9.D0*CF*TF*x**2*Lperp - 32.D0/3.D0*CF*
     & TF*x**2*Lperp**2 + 8.D0/9.D0*CF*TF*x**2*pisq + 16.D0/3.D0*Hr2(0,
     & 1)*CF*TF*x**(-1) - 8*Hr2(0,1)*CF*TF + 8*Hr2(0,1)*CF*TF*x - 16.D0/
     & 3.D0*Hr2(0,1)*CF*TF*x**2 + 28.D0/3.D0*Hr1(0)*CF*TF - 8*Hr1(0)*CF
     & *TF*Lperp + 16*Hr1(0)*CF*TF*Lperp**2 - 40.D0/3.D0*Hr1(0)*CF*TF*x
     &  - 24*Hr1(0)*CF*TF*x*Lperp + 16*Hr1(0)*CF*TF*x*Lperp**2 + 128.D0/
     & 9.D0*Hr1(0)*CF*TF*x**2 - 64.D0/3.D0*Hr1(0)*CF*TF*x**2*Lperp - 
     & Hr1(0)**2*CF*TF + 8*Hr1(0)**2*CF*TF*Lperp - Hr1(0)**2*CF*TF*x + 
     & 8*Hr1(0)**2*CF*TF*x*Lperp - 8.D0/3.D0*Hr1(0)**2*CF*TF*x**2 + 2.D0
     & /3.D0*Hr1(0)**3*CF*TF + 2.D0/3.D0*Hr1(0)**3*CF*TF*x - 16.D0/3.D0
     & *Hr1(0)*Hr1(1)*CF*TF*x**(-1)
      I2(qbpq,rglr) = I2(qbpq,rglr) + 8*Hr1(0)*Hr1(1)*CF*TF - 8*Hr1(0)*
     & Hr1(1)*CF*TF*x + 16.D0/3.D0*Hr1(0)*Hr1(1)*CF*TF*x**2


      return
      end
