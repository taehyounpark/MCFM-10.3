!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_dm_monojet(p,msq)
      implicit none
      include 'types.f'

c===== gluon produced DM, shamelessly adapted from Higgs routnes
c===== with bb=>dm
c---- matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
c----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(-->  b(p3)+b~(p4))+g(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'dm_params.f'
      integer:: j,k,iglue
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ss,tt,uu,mhsq,hdecay,s(mxpart,mxpart)
      real(dp):: Asq,fac
      parameter(iglue=5)

      msq(:,:)=0._dp

      call dotem(iglue,p,s)
      call dmsdecay(p,3,4,hdecay)


      Asq=hdecay/dm_lam**6*as**2*16d0
c--   calculate propagators
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      mhsq=ss+tt+uu

      if(effective_th) then
         fac=gsq*Asq
      else
         write(6,*) 'Full theory not implemented for GG operator'
         write(6,*) 'Setting EFT to true : We are back in EFT'
         fac=gsq*Asq
         effective_th=.true.
      endif

      msq(0,0)=
     & avegg*fac*V*xn*(mhsq**4+ss**4+tt**4+uu**4)/(ss*tt*uu)
      msq(1,-1)=+aveqq*fac*V/2d0*(tt**2+uu**2)/ss
      msq(0,+1)=-aveqg*fac*V/2d0*(ss**2+tt**2)/uu
      msq(+1,0)=-aveqg*fac*V/2d0*(ss**2+uu**2)/tt

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      if ((k == -j) .and. (j  /=  0)) then
      msq(j,k)=msq(1,-1)
      elseif ((j == 0) .and. (k  /=  0)) then
      msq(j,k)=msq(0,1)
      elseif ((j  /=  0) .and. (k == 0)) then
      msq(j,k)=msq(1,0)
      endif
      enddo
      enddo

      return
      end


      subroutine dmsdecay(p,ib,ibb,msq)
      implicit none
      include 'types.f'
c**** scalar decay to dm taken from H=>bb routine

      include 'mxpart.f'
      include 'dm_params.f'
      integer:: ib,ibb
      real(dp):: p(mxpart,4),s56,msq

      s56=2d0*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))


      msq=4d0*(s56/2d0-xmass**2)

      return
      end





