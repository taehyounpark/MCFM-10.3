!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tildeSb1(i,tSb1)
c     first argument is gg or qq
c     second argument is power of Lb
c     1909.00811v2, Eq 3.16
      implicit none
      include 'types.f'
      include 'qtconstants.f'
      include 'Lnu.f'
      integer i
      real(dp)::tSb1(0:2)
      tSb1(2)=-Gamma0(i)/2._dp
      tSb1(1)=Lnu*2*Gamma0(i)+0.5_dp*(tgammaS0(i)+tgamman0(i))
      tSb1(0)=-Lnu*tgamman0(i)+tildes1(i)
      return
      end
