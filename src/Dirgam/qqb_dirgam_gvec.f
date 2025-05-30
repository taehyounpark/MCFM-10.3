!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dirgam_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     October, 2002.                                                   *
c     Matrix element for gamma production                              *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p4)              *
c     q(-p1)+qbar(-p2)--> gamma(p3)+ g(p4)                             *
c***********************************************************************
      include 'nf.f'
      include 'mxpart.f'
      integer:: j,in
c--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),msqa(-nf:nf,-nf:nf),
     & p(mxpart,4),n(4),nDn

      msq(:,:)=0._dp

      call qqb_dirgam(p,msqa)
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2

      call checkndotp(p,n,in)

      do j=-nf,nf

      if (in == 1) then
        msq(0,j)=-0.5_dp*nDn*msqa(0,j)
      elseif (in == 2) then
        msq(j,0)=-0.5_dp*nDn*msqa(j,0)
      elseif (in == 4) then
        msq(j,-j)=-0.5_dp*nDn*msqa(j,-j)
      endif

      enddo

      return
      end
