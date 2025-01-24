!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Bdiff(s34,s12,msq)
      use loopI2_generic
      implicit none
      include 'types.f'
      complex(dp):: Bdiff
      include 'scale.f'
      include 'scalarselect.f'
      real(dp):: s34,s12,msq

      Bdiff=loopI2(s34,msq,msq,musq,0)-loopI2(s12,msq,msq,musq,0)
      return
      end
