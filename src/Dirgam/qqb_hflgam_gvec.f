!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_hflgam_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     January, 2013.                                                   *
c     Matrix element for gamma production                              *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p4)              *
c     q(-p1)+qbar(-p2)--> gamma(p3)+ b/c(p4)                           *
c***********************************************************************
      include 'nf.f'
      include 'mxpart.f'
      include 'heavyflav.f'
      integer:: in
c--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),msqa(-nf:nf,-nf:nf),
     & p(mxpart,4),n(4),nDn

      msq(:,:)=0._dp

      call qqb_hflgam(p,msqa)
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2

      call checkndotp(p,n,in)

      if (in == 1) then
        msq(0,flav)=-0.5_dp*nDn*msqa(0,flav)
      elseif (in == 2) then
        msq(flav,0)=-0.5_dp*nDn*msqa(flav,0)
      elseif (in == 4) then
        write(6,*) 'Check code in qqb_hflgam_gvec.f'
        stop
      endif

      return
      end
