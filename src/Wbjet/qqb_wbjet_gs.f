!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wbjet_gs(p,msqc)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J. M. Campbell                                           *
c     January, 2004.                                                   *
c                                                                      *
c     This is merely a wrapper routine to qqb_w(m/p)bjet_gs            *
c***********************************************************************

      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf)

      if     (nwz == +1) then
        call qqb_wpbjet_gs(p,msqc)
      elseif (nwz == -1) then
        call qqb_wmbjet_gs(p,msqc)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_wbjet_gs.f'
      endif

      return
      end

