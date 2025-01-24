!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_w_v(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      integer:: j,k
      real(dp):: p(mxpart,4),
     & msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     & dot,virt,xl12

      xl12=log(two*dot(p,1,2)/musq)

c---  calculate lowest order matrix element
      call qqb_w(p,msq)
c---calculate the multiple of the lowest order
      scheme='dred'
      virt=ason2pi*cf*(-two*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &                 -three*(epinv-xl12)
     &                 +pisq-seven)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=virt*msq(j,k)
      enddo
      enddo
      end

