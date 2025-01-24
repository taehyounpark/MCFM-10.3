!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine calI2(z,cI2)
      implicit none
c     1909.00811v2, Eq.3.26
      include 'types.f'
      include 'Lw.f'
      include 'constants.f'
      include 'qtconstants.f'
      include 'scale.f'
      include 'facscale.f'
c     First argument gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     Second argument type of distribution
c---  -2 distribution regular piece
c---  -1 the coefficient of delta(1-z)
c---  =0 the coefficient of L0(1-z)
c---  =1 the coefficient of L1(1-z)
c     Third argument power of Lb
      real(dp):: cI2(0:6,dmin:dmax,0:2),
     & tI1(0:6,dmin:dmax),tI2(0:6,dmin:dmax),
     & P0(0:6,dmin:dmax),P0P0(0:6,dmin:dmax),P1(0:6,dmin:dmax),
     & tI1xP0(0:6,dmin:dmax)
      integer i,ii
      integer,parameter::x(0:6)=[0,1,1,1,1,0,1]
      real(dp)::z
c     Our generic decomposition of the quark kernels adequate for two loops
c     K(qi,qj)=K(qbi,qbj)=delta_(ij)*KqqV +KqqS
c     K(qi,qbj)=K(qbi,qj)=delta_(ij)*KqqbV+KqqS (K(qi,qbj) is zero at order alphas)

      cI2(:,:,:)=0._dp

      call tildeI1(z,tI1)
      call tildeI2(z,tI2)
      call AP0(z,P0)
      call AP1(z,P1)
      call APP0xP0(z,P0P0)
      call APtI1xP0(z,tI1xP0)

c     1909.00811v2, Eq 3.26
      do i=0,6
      ii=x(i)
c     delta function piece
      if (i < 2) then   ! only for i=j (gg, qqV) contributions
      cI2(i,delt,2)=Lw**2*Gamma0(ii)**2/two
     & +Lw*Gamma0(ii)/two*(beta0+tgammaB0(ii))
     & +(beta0+tgammaB0(ii)/two)*tgammaB0(ii)/4._dp
      cI2(i,delt,1)=-Lw**2*Gamma0(ii)*tgamman0(ii)/two
     & +Lw*(-(beta0+tgammaB0(ii)/two)*tgamman0(ii)/two+Gamma1(ii))
     & +tgammaB1(ii)/two
      cI2(i,delt,0)=Lw**2*tgamman0(ii)**2/8._dp-Lw*tgamman1(ii)/two
      endif

c     all other z-pieces
c Note extra term for factorization scale dependence below (check?)
      cI2(i,:,2)=cI2(i,:,2)
     & -(Lw*Gamma0(ii)+beta0/two+tgammaB0(ii)/two)*P0(i,:)
     & +P0P0(i,:)/two
      cI2(i,:,1)=cI2(i,:,1)+Lw*tgamman0(ii)/two*P0(i,:)-P1(i,:)
     & +(Lw*Gamma0(ii)+beta0+tgammaB0(ii)/two)*tI1(i,:)-tI1xP0(i,:)
     & +2*((Lw*Gamma0(ii)+0*beta0+tgammaB0(ii)/two)*P0(i,:)
     &     -P0P0(i,:))*log(scale/facscale)
      cI2(i,:,0)=cI2(i,:,0)-Lw*tgamman0(ii)/two*tI1(i,:)+tI2(i,:)
     & +2*(-Lw*tgamman0(ii)/two*P0(i,:)
     &     +P1(i,:)+tI1xP0(i,:))*log(scale/facscale)
     & +2*(beta0*P0(i,:)+P0P0(i,:))*log(scale/facscale)**2
      enddo
      return
      end

