!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_w2jet_g(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J. M. Campbell                                           *
c     August, 2002.                                                    *
c                                                                      *
c     This is merely a wrapper routine to qqb_w(m/p)2jet_g             *
c***********************************************************************
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)

      if     (nwz == +1) then
        call qqb_wp2jet_g(p,msq)
c      write(*,*) msq(-1,2),p(7,:)
      elseif (nwz == -1) then
        call qqb_wm2jet_g(p,msq)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_w2jet_g.f'
      endif
      return
      end

