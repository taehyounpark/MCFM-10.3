!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tildeI2(z,tI2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'singletlabels.f'
      include 'qtconstants.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     delta=-1,0=L0(z),1=L1(z),-2=finite
      real(dp):: z,omz,tI2(0:6,dmin:dmax),I2slsh(0:6,dmin:dmax)
      tI2(:,:)=0
      omz=1._dp-z
c     1602.01829v3, Eq.6.18
      tI2(gg,delt)=0.5_dp*(CA**2*pi**4/36._dp-S2(gg))
      tI2(gq,rglr)=CF*CA*pisq/3._dp*z
c     1602.01829v3, Eq.6.19
      tI2(qqV,delt)=0.5_dp*(CF**2*pi**4/36._dp-S2(qqV))
      tI2(qqV,rglr)=CF**2*pisq/3._dp*omz
      tI2(qg,rglr)=2*CF*TF*z*omz*pisq/3._dp

      call I2slash(z,I2slsh)
      tI2(:,:)=tI2(:,:)+I2slsh(:,:)
      return
      end
