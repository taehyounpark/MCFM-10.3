!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmD1x2x3_generic
      implicit none
      public pppmD1x2x3,pppmD1x2x3_qp

      interface pppmD1x2x3
      module procedure pppmD1x2x3,pppmD1x2x3_qp
      end interface

      contains

      function pppmD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD1x2x3_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmD1x2x3_inc.f'
      end function pppmD1x2x3

      function pppmD1x2x3_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD1x2x3_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmD1x2x3_inc.f'
      end function pppmD1x2x3_qp

      end module pppmD1x2x3_generic

