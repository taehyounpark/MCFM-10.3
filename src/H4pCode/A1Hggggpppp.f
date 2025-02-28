!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function A1Hggggpppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1Hggggpppp
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: A1ggggmmmmCC,A1ggggmmmmNCC
c--- Use 0704.3914v3 Eqs. (2.4) and (2.6) to relate pppp amplitude
c--- to mmmm one; c.c. is equivalent to interchanging za and zb

c--- Expresssion of Eq. (9) of Badger and Glover, hep-ph/0607139v2
      A1Hggggpppp=A1ggggmmmmCC(j1,j2,j3,j4,zb,za)
     &           +A1ggggmmmmNCC(j1,j2,j3,j4,zb,za)

      return
      end
