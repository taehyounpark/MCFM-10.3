!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hard_gamgam(p,order,msq0,msq1,msq2)
      implicit none
      include 'types.f'
c---- LO * (hard function) for diphoton processes, in units of as/4/pi
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      real(dp),intent(out):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      integer order
      real(dp):: p(mxpart,4),qqb,tlrem,msq0(-nf:nf,-nf:nf),hard(2),msqextra(-nf:nf,-nf:nf)

      call gamgamampsq_new(order,p,1,2,3,4,qqb,hard,tlrem)
c get back to coefficients of (as/4/pi) from (as/2/pi) returned by this routine
      hard(1)=hard(1)*2._dp
      hard(2)=hard(2)*4._dp

      call qqb_gamgam(p,msq0)

      msq1(:,:)=msq0(:,:)*hard(1)
      msq2(:,:)=msq0(:,:)*hard(2)

c--- add contributions that do not factorize over LO
      call qqb_gamgamas2(p,tlrem,msqextra)
      msq2(:,:)=msq2(:,:)+msqextra(:,:)/ason4pi**2

      return
      end


