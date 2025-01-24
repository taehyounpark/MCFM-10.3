!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_gamgam(p,msq)
      implicit none
      include 'types.f'
c=====
c-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+gamma(p4)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zcouple_cms.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qa,aq,gg,facgg,Qsum
      real(dp),parameter::statfac=0.5_dp

c--set msq=0 to initalize
      msq(:,:)=0._dp

      call dotem(3,p,s)

      fac=8._dp*xn*abs(zesq)**2*statfac

      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      facgg=4._dp*abs(zesq)*gsq/(16._dp*pisq)*Qsum

      gg=zip
c JC removing gg contribution
c      if (omitgg) then
c        gg=0._dp
c      else
c        gg=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
c      endif

      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa

      do j=-nf,nf
      k=-j
c--qa
      if (j > 0) then
        msq(j,k)=Q(j)**4*qa
c--aq
      elseif (k > 0) then
        msq(j,k)=Q(k)**4*aq
c--gg
      elseif (j == 0) then
        msq(j,k)=gg
      endif

      enddo

      return
      end

      subroutine gg_gamgam(p,msq)
      implicit none
      include 'types.f'
c=====
c-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+gamma(p4)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zcouple_cms.f'
      real(dp), intent(out) :: msq
      real(dp), intent(in) :: p(mxpart,4)

      real(dp) :: facgg,msqgggaga,Qsum
      real(dp), parameter :: statfac=0.5_dp

      call dotem(3,p,s)

      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      facgg=4*abs(zesq)*gsq/(16*pisq)*Qsum

      msq = avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac

      end
