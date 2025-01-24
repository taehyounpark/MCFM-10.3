!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine iperm(IHEL,PERM,IHELX,NGLUONS)
      implicit none
      integer::J,NGLUONS,PERM(NGLUONS),IHEL(NGLUONS),IHELX(NGLUONS)
c permutes helicities to match momenta
      do j=1,ngluons
      ihelx(j)=IHEL(PERM(j))
      enddo
      return
      end
