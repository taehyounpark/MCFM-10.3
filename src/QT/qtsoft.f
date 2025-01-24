!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c#################################################################
c##### a factor of alphas/4/pi has been extracted out
c#################################################################
      subroutine qtsoft(i,order,tSb1,tSb2)
      implicit none
      include 'types.f'
      integer :: i,order
c     i=0 for gluon, i=1 for quark
c     tSb1 -- order as/4/pi contribution to soft function
c     tSb2 -- order (as/4/pi)^2 contribution to soft function
c     index (0:2),(0:4) encodes power of Lb

      real(dp) :: tSb1(0:2),tSb2(0:4)

      if (order < 1) return
      call tildeSb1(i,tSb1)

      if (order < 2) return
      call tildeSb2(i,tSb2)

      return
      end


