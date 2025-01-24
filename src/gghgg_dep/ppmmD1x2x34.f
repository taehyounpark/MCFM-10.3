!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmD1x2x34_generic
      implicit none
      public ppmmD1x2x34,ppmmD1x2x34_qp

      interface ppmmD1x2x34
      module procedure ppmmD1x2x34,ppmmD1x2x34_qp
      end interface

      contains

      function ppmmD1x2x34(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppmmD1x2x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmD1x2x34_inc.f'
      end function ppmmD1x2x34

      function ppmmD1x2x34_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppmmD1x2x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmD1x2x34_inc.f'
      end function ppmmD1x2x34_qp

      end module ppmmD1x2x34_generic

