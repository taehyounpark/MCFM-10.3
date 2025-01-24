!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine I1bar(x,lnmuonpt,I1)
      implicit none
!     using results from 1412.8408 Eqs(A.1-A.3)
      include 'types.f'
      include 'constants.f'
      include 'transitionlabels.f'
      include 'distributions.f'
      include 'nfl.f'
c     dmin=-1,dmax=2,delt=-1,plus=0,lpls=1,rglr=2
      real(dp)::I1(0:6,dmin:dmax)
c      I1, index one transition flavor
c      I2, index2,distribution type
      real(dp):: x,lnmuonpt,omx,beta0

      beta0=(11*CA-4*TF*nfl)/3d0

      I1(:,:)=0
      omx=1._dp-x
      I1(gg,plus)=-4*CA*2d0*lnmuonpt
      I1(gg,rglr)=-4*CA*(2._dp/x-4._dp+2._dp*x*omx)*lnmuonpt
      I1(gg,delt)=-CA*pisq/6._dp-2._dp*beta0*lnmuonpt

! Remember [(1+z^2)/(1-z)]_+ = 2/(1-z)_+ - (1+z) + 3/2*delta(1-z)
      I1(qq,plus)=-4*CF*2._dp*lnmuonpt
      I1(qq,rglr)=-4*CF*(-1._dp-x)*lnmuonpt+2*CF*omx
      I1(qq,delt)=-CF*pisq/6._dp-6._dp*CF*lnmuonpt

      I1(gq,rglr)=-4*CF*(1._dp+omx**2)/x*lnmuonpt+2._dp*CF*x
      I1(qg,rglr)=-4*TF*(x**2+omx**2)*lnmuonpt+4._dp*TF*x*omx

      return
      end
