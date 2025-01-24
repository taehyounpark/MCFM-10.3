!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine APtI1xP0(z,tI1xP0)
      implicit none
c     Iqq=2*CF*(1-y)
c     Iqg=4*TR*y*(1-y)
c     Igq=2*CF*y
c     Igg=0

c     Pqq=2*CF*(2/(1-x)_+ -1-x +3/2*delta(1-x)
c     Pqg=2*TR*(x^2+(1-x)^2)
c     Pgq=2*CF*(1+(1-x)^2)/x
c     Pgg=4*CA*(1/(1-x)_+ -1+(1-x)/x+x*(1-x))+beta0*delta(1-x)

cId,agg=CF*TR*2*nfl*(4*omz*(2*z+1)+8*lnz*z);
cId,agq=CF**2*(4+2*z-4*lnz*z+8*lnomz*z);
cId,aqg=TR*4*z*omz*beta0
c      +CA*TR*(-8*omz*(17*z^2+2*z-1)/3/z-32*lnz*z+16*lnomz*z*omz)
c      +CF*TR*(-4*omz*(z+2)-4*lnz-8*lnz*z);
cId,aqqV=CF**2*omz*(-2+4*(2*lnomz-lnz));
cId,aqqS=CF*TR*(-8*omz*(2*z^2+2*z-1)/(3*z)-8*lnz*z);

      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'qtconstants.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6

      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     -2=regular,delta=-1,0=L0(z),1=L1(z)
      real(dp), intent(in) :: z
      real(dp),intent(out) :: tI1xP0(0:6,dmin:dmax)
      real(dp):: lnz,omz,lnomz

c     initialize to zero
      tI1xP0(:,:)=0

      lnz=log(z)
      omz=1._dp-z
      lnomz=log(omz)
      tI1xP0(gg,rglr)=CF*TR*2*nfl*(4*omz*(2*z+1._dp)+8*lnz*z)
      tI1xP0(gq,rglr)=CF**2*(4._dp+2*z-4*lnz*z+8*lnomz*z)
      tI1xP0(qg,rglr)=TR*4*z*omz*beta0+CA*TR*(-8*omz*(17*z**2+2*z-1._dp)
     & /(3._dp*z)-32*lnz*z+16*lnomz*z*omz)
     & +CF*TR*(-4*omz*(z+2._dp)-4*lnz-8*lnz*z)
      tI1xP0(qqV,rglr)=CF**2*omz*(-2._dp+4*(2*lnomz-lnz))
      tI1xP0(qqS,rglr)=CF*TR*(-8*omz*(2*z**2+2*z-1._dp)
     & /(3._dp*z)-8*lnz*z)
      return
      end

