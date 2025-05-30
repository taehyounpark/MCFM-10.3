!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zh_ww_v(p,msqv)
      implicit none
      include 'types.f'

      integer:: j,k
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      real(dp):: p(mxpart,4),
     & msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     & dot,virt,xl12

      xl12=log(two*dot(p,1,2)/musq)
      scheme='dred'
c---  calculate lowest order matrix element
      call qqb_zh_ww(p,msq)
c---calculate the multiple of the lowest order
      virt=ason2pi*cf*(-2._dp*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &                 -3._dp*(epinv-xl12)
     &                 +pisq-7._dp)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=virt*msq(j,k)
      enddo
      enddo
      end

