!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
       subroutine AP0(z,p0)
       implicit none
c      See, for example, 1401.5478v2, Eq.(A.11,A.12)
       include 'types.f'
       include 'constants.f'
       include 'qtconstants.f'
       include 'singletlabels.f'
c      gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6

      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     general definition of functions
c     -1 delta function
c      0 coefficient of L0(z)
c      1 coefficient of L1(z)
c      -2 distribution-free piece
       real(dp), intent(in) :: z
       real(dp),intent(out) :: p0(0:6,dmin:dmax)


       p0(:,:)=0

       p0(gg,delt)=0.5_dp*beta0  !delta function
       p0(gg, plus)= 2*CA  !coeff of 1/(1-z)
       p0(gg,rglr)= 2*CA*(-1._dp+(1._dp-z)/z+z*(1._dp-z)) ! distribution-free piece

       p0(gq,rglr)=CF*(1._dp+(1._dp-z)**2)/z
       p0(qg,rglr)=((1._dp-z)**2+z**2)*TR

       p0(qqV,delt) = 1.5_dp*CF !delta function
       p0(qqV,plus) = 2*CF  !coeff of 1/(1-z)_+
       p0(qqV,rglr) = -CF*(1._dp+z) ! distribution-free piece

c  correct to as/4/pi
       p0(:,:)=2*p0(:,:)
       return
       end
