!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_gamgam_vbis(p,msq,order)
c------ routine for qqb_gamgam using hard coeff
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple_cms.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      integer:: order
      real(dp) :: hard(2),qqb,qbq
      real(dp)::statfac,fac
      integer:: j,k
      parameter(statfac=0.5_dp)

      fac=8._dp*xn*abs(zesq)**2*statfac

c======= compute Matrix elements
      call gamgamampsq(order,p,1,2,3,4,qqb,hard)

c---- apply overall factor
      qqb=fac*aveqq*qqb
      qbq=qqb

c----- higher order corrections
      if(order==1) then
         qqb=qqb*ason2pi*hard(1)
         qbq=qbq*ason2pi*hard(1)
      endif

      msq(:,:)=0._dp
      do j=-nf,nf
         k=-j
         if (j*k >= 0) cycle    ! skip gluons, qq, aa

         if ((j > 0) .and. (k < 0)) then
            msq(j,k)=qqb*Q(j)**4
         elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=qbq*Q(k)**4
         endif


      enddo
      return
      end
