!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine APP0xP0(z,p0p0)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'qtconstants.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     regular=-2,delta=-1,0=L0(z),1=L1(z)
      real(dp), intent(in) :: z
      real(dp) :: opz,omz,p0xp0(0:6,0:6,dmin:dmax),P0P0(0:6,dmin:dmax),
     & p0gg,p0gq,p0qg,p0qq

      p0xp0(:,:,:)=0
      p0p0(:,:)=0

      omz=1._dp-z
      opz=1._dp+z
      p0qq=(1._dp+z**2)/omz
      p0gq=(1.0_dp+omz**2)/z
      p0qg=(z**2+omz**2)
      p0gg=2._dp*(omz+z**2)**2/z

c  qqV qqV
c      See, for example, 1401.5478v2, Eq.(A.17)
      p0xp0(qqV,qqV,delt) = -(9/4.0_dp+2*pisq/3.0_dp) + 3*3/2.0_dp
      p0xp0(qqV,qqV,plus) = 6._dp
      p0xp0(qqV,qqV,lpls) = 8._dp
      if (z == 1._dp) then
      p0xp0(qqV,qqV,rglr)=4._dp
      else
      p0xp0(qqV,qqV,rglr)=-2*omz+(1._dp+z-2*p0qq)*log(z)
      p0xp0(qqV,qqV,rglr)=p0xp0(qqV,qqV,rglr)-3*opz-4*opz*log(omz)
      endif
      p0xp0(qqV,qqV,:)=p0xp0(qqV,qqV,:)*CF**2

c  qg gq
      p0xp0(qg,gq,rglr)=2*opz*log(z)+4/3.0_dp/z+omz-4/3.0_dp*z**2
      p0xp0(qg,gq,rglr)=p0xp0(qg,gq,rglr)*CF*TR

c  qg gg
c This is an implementation of the full convolution of
c P(0)qg x P(0)gg (z), including all color factors and
c endpoint contributions.
c It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
c factor CA*TR and adding the endpoint explicit in Eq. (A.5)
      p0xp0(qg,gg,rglr) = 2*p0qg*log(omz)+2*(1._dp+4*z)*log(z)
     &        +4/3.0_dp/z+1._dp+8*z-31/3.0_dp*z**2
     &        +beta0/two*p0qg/CA

      p0xp0(qg,gg,rglr)=p0xp0(qg,gg,rglr)*CA*TR

c  qq qg
c      See, for example, 1405.1044v2, Eq.(A.17)
      p0xp0(qqV,qg,rglr)=2*p0qg*log(omz/z)+(1._dp-2*z)*log(z)-0.5_dp+2*z
      p0xp0(qqV,qg,rglr)=p0xp0(qqV,qg,rglr)*CF*TR

c  gg gg
c This is an implementation of the full convolution of
c P(0)gg x P(0)gg (z), including all color factors and
c endpoint contributions.
c It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
c factor CA^2 and adding the endpoint explicit in Eq. (A.5)
      p0xp0(gg,gg,delt) = -2*pisq/3.0_dp + beta0**2/4.0_dp/CA**2
      p0xp0(gg,gg,lpls)=8._dp
      p0xp0(gg,gg,plus)=beta0*2.0_dp/CA
      if (z == 1._dp) then
      p0xp0(gg,gg,rglr)=0
      else
      p0xp0(gg,gg,rglr) = - 2*(p0gg/omz+4*opz)*log(z)
     &            - 44/3.0_dp/z + 12*omz+44/3.0_dp*z**2
     &            +8*(1._dp-2*z+z**2-z**3)/z*log(omz)
     &            +(-2._dp+1._dp/z+z-z**2)*beta0*2/CA
c      p0xp0(gg,gg,rglr)=p0xp0(gg,gg,rglr)+8._dp*omz/z*log(omz)
c      p0xp0(gg,gg,rglr)=p0xp0(gg,gg,rglr)+beta0*2.0_dp/CA
      endif
      p0xp0(gg,gg,:)=p0xp0(gg,gg,:)*CA**2

c  gq qg
      p0xp0(gq,qg,rglr) = 2*opz*log(z)+4/3.0_dp/z+omz-4/3.0_dp*z**2
      p0xp0(gq,qg,rglr)=p0xp0(gq,qg,rglr)*CF*TR

c  gq qqV

      p0xp0(gq,qqV,rglr)=2*p0gq*log(omz)+(2._dp-z)*log(z)+2._dp-z/2.0_dp
      p0xp0(gq,qqV,rglr)=p0xp0(gq,qqV,rglr)*CF**2


c This is an implementation of the full convolution of
c P(0)gg x P(0)gq (z), including all color factors and
c endpoint contributions.
c It is obtained by dressing Eq. (A.11) of 1405.1044 with a color
c factor CA*CF and adding the endpoint explicit in Eq. (A.5)
 !  gg gq
      p0xp0(gg,gq,rglr) = 2*p0gq*log(omz/z)-2*(4._dp+z)*log(z)
     &        -31/3.0_dp/z+8._dp+z+4/3.0_dp*z**2
     &         +beta0/two*p0gq/CA

      p0xp0(gg,gq,rglr)=p0xp0(gg,gq,rglr)*CA*CF

c   Implementation of 1909.00811v2, Eq. D13
c   for the special case of low order kernels P(qqS)=0,P(qbV)=0
      p0p0(gg,:)=p0xp0(gg,gg,:)+2*nfl*p0xp0(gq,qg,:)
      p0p0(qg,:)=p0xp0(qqV,qg,:)+p0xp0(qg,gg,:)
      p0p0(gq,:)=p0xp0(gq,qqV,:)+p0xp0(gg,gq,:)
      p0p0(qqV,:)=p0xp0(qqV,qqV,:)
      p0p0(qqS,:)=p0xp0(qg,gq,:)

c     correct overall normalization for (as/4/pi)
      p0p0(:,:)=4*p0p0(:,:)

      return
      end
