!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine spinoruorz(N,p,za,zb)
      implicit none
      include 'types.f'
c--- This routine just provides an easy way of switching between
c--- spinoru (normal MCFM running) and spinorz (checks of virtual)
      include 'mxpart.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4)
      integer:: N

      call spinoru(N,p,za,zb)
c      call spinorz(N,p,za,zb)

      return
      end
