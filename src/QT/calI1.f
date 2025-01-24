!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine calI1(z,cI1)
      implicit none
c     1909.00811v2, Eq.3.26
      include 'types.f'
      include 'Lw.f'
      include 'qtconstants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'singletlabels.f'
c     First argument gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6

      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1

      real(dp):: cI1(0:6,dmin:dmax,0:1)
c     Second argument type of distribution
c---  -2 distribution regular piece
c---  -1 the coefficient of delta(1-z)
c---  =0 the coefficient of L0(1-z)
c---  =1 the coefficient of L1(1-z)
c     Third argument power of Lb
      real(dp)::tI1(0:6,dmin:dmax),P0(0:6,dmin:dmax)
      integer i
      real(dp)::z

      cI1(:,:,:)=0._dp

      call tildeI1(z,tI1)
      call AP0(z,P0)
c     1909.00811v2, Eq. 3.26
      do i=gg,qqV
      cI1(i,delt,0)=-0.5_dp*Lw*tgamman0(i)         !delta function piece Lb^0
      cI1(i,delt,1)=Lw*Gamma0(i)+tgammaB0(i)/2._dp !delta function piece Lb^1
      enddo

c Note extra term for factorization scale dependence below
      cI1(:,:,1)=cI1(:,:,1)-P0(:,:)
      cI1(:,:,0)=cI1(:,:,0)+tI1(:,:)+2*P0(:,:)*log(scale/facscale)
      return
      end

