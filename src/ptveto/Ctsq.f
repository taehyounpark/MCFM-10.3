!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine Ctsq(order,musq,coeff)
      use LHAPDF, only: getalphas
      implicit none
! Ct(mt^2,mut^2) given in Eq. (12) of 0809.4283
! and evolution to Ct(mt^2,mu^2) given in Eq. (20)
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nfl.f'
      include 'facscale.f'
      integer,intent(in)::order
      real(dp),intent(in)::musq
      real(dp),intent(out)::coeff(2)
      real(dp):: Ct1,Ct2,alphasmu,ason4pi,logmt2omu2,alphasptveto,alphasmut,beta0,beta1,beta2

! Compute alphas at the resummation scale (facscale, close to ptveto)
      alphasmu=getalphas(facscale)
      ason4pi=alphasmu/fourpi

! For Ct we fix the scale mut = mt so that the log is zero
      logmt2omu2=zip ! 2*log(mt/mut)
!      logmt2omu2=log(mt**2/abs(musq))   ! debug for 1-step/2-step comparison
      Ct1=5*CA-3*CF
      Ct2=27._dp/2._dp*CF**2+(11*logmt2omu2-100._dp/3._dp)*CF*CA
     &   -(7*logmt2omu2-1063._dp/36._dp)*CA**2
     &   -4._dp/3._dp*CF*TF-5._dp/6._dp*CA*TF
     &   -(8*logmt2omu2+5)*CF*TF*nfl-47._dp/9._dp*CA*TF*nfl
! Expansion of Ct^2, c.f. Eq. (2) of 1205.3806
      beta0=11._dp-2._dp*nfl/3._dp
      beta1=102._dp-38*nfl/3._dp
      beta2=2857/2._dp-5033/18._dp*nfl+325/54._dp*nfl**2
      alphasmut=getalphas(mt)
!      alphasmut=alphasmu ! debug for numerical study
      coeff(:)=0._dp
      if (order > 0) then
        coeff(1)=alphasmut/fourpi*2*Ct1
     &          +(alphasmu-alphasmut)/fourpi*2*beta1/beta0
      endif
      if (order > 1) then
        coeff(2)=(alphasmut/fourpi)**2*(Ct1**2+2*Ct2)
     &           +alphasmut*(alphasmu-alphasmut)/fourpi**2
     &              *(4*beta2+4*beta1*Ct1-2*beta1**2/beta0)/beta0
     &           +(alphasmu-alphasmut)**2/fourpi**2
     &              *(2*beta2+beta1**2/beta0)/beta0
      endif

      return
      end

