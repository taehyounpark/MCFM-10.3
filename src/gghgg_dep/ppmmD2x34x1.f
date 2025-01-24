!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppmmD2x34x1_generic
      implicit none
      public ppmmD2x34x1,ppmmD2x34x1_qp

      interface ppmmD2x34x1
      module procedure ppmmD2x34x1,ppmmD2x34x1_qp
      end interface

      contains

      function ppmmD2x34x1(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppmmD2x34x1_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppmmD2x34x1_inc.f'
      end function ppmmD2x34x1

      function ppmmD2x34x1_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppmmD2x34x1_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppmmD2x34x1_inc.f'
      end function ppmmD2x34x1_qp

      end module ppmmD2x34x1_generic

