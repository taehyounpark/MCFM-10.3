!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine VV_Hqq_gs(p,msqc)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion : sums up WW and VV contributions
c***********************************************************************
c     Author: J. M. Campbell                                           *
c     July, 2002.                                                      *
c                                                                      *
c     Weak Boson Fusion : sums up WW and ZZ contributions              *
c     This routine calculates the dipole subtraction terms             *
c     for the process:                                                 *
c     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
c                                                                      *
c***********************************************************************
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf),
     & msqc_ww(maxd,-nf:nf,-nf:nf),msqc_zz(maxd,-nf:nf,-nf:nf)
      integer:: j,k,nd

      call WW_Hqq_gs(p,msqc_ww)
      call ZZ_Hqq_gs(p,msqc_zz)

      do nd=1,ndmax
      do j=-nf,nf
      do k=-nf,nf
        msqc(nd,j,k)=msqc_ww(nd,j,k)+msqc_zz(nd,j,k)
      enddo
      enddo
      enddo

      return
      end

