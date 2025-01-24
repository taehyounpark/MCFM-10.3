!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A6texact(s23,mt2)
      use loopI2_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scalarselect.f'
      complex(dp)::A6texact
      real(dp)::s23,mt2,musq

c--- musq is irrelevant, set it to some value
      musq=mt2 !abs(s23)

c--- this is the form of the top loop vacuum polarization obtained
c--- after renormalization by subtraction at zero momentum transfer
c--- note: typo in Eq. (B12) of arXiv:1610.02189, sign of final term is wrong
      A6texact=-two/three*(
     &  (one+two*mt2/s23)*(loopI2(s23,mt2,mt2,musq,0)-loopI2(0._dp,mt2,mt2,musq,0))
     &   -one/three)

      return
      end
