!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine storedip(msq_dip,msq_dipv,dsub,dsubv,
     &                    sub_dip,sub_dipv,n)
      implicit none
      include 'types.f'
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration

      include 'nf.f'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      integer:: i,j,k,n
      real(dp):: msq_dip(6,0:2,-nf:nf,-nf:nf),dsub(4),sub_dip(6,4),
     &           msq_dipv(6,0:2,-nf:nf,-nf:nf),dsubv,sub_dipv(6)

      msq_dip(n,:,:,:) = msq_cs(:,:,:)
      msq_dipv(n,:,:,:) = msqv_cs(:,:,:)

      sub_dip(n,:) = dsub(:)
      sub_dipv(n)=dsubv

      return
      end
