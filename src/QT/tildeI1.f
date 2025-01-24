!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tildeI1(z,tI1)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6

      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     -2=regular,delta=-1,0=L0(z),1=L1(z)
      real(dp):: z,omz,tI1(0:6,dmin:dmax)
      tI1(:,:)=0
      omz=1._dp-z
c     1602.01829v3
c     Eq.6.14
      tI1(gq,rglr)=2*CF*z
c     Eq.6.15
      tI1(qqV,rglr)=2*CF*omz
      tI1(qg,rglr)=4*TF*z*omz
      return
      end
