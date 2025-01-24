!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hard_Vgam(p,order,msq0,msq1,msq2)
      implicit none
      include 'types.f'
c---- LO * (hard function) for Wgamma and Zgamma processes, in units of as/4/pi
      include 'kprocess.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      real(dp),intent(out):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      integer:: order
      real(dp)::p(mxpart,4),msq0(-nf:nf,-nf:nf),msq(-nf:nf, -nf:nf, 0:2)

      if     (kcase==kWgamma) then
        call wgam_mat(p,msq)
      elseif (kcase==kZgamma) then
        call zgam_mat(p,msq)
      else
        write(6,*) 'Unexpected process in hard_Vgam: ',kcase
        stop
      endif

      msq0(:,:)=msq(:,:,0)
      msq1(:,:)=msq(:,:,1)/ason4pi
      msq2(:,:)=msq(:,:,2)/ason4pi**2

      return
      end


