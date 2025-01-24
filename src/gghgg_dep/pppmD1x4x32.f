!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmD1x4x32_generic
      implicit none
      public pppmD1x4x32,pppmD1x4x32_qp

      interface pppmD1x4x32
      module procedure pppmD1x4x32,pppmD1x4x32_qp
      end interface

      contains

      function pppmD1x4x32(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD1x4x32_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmD1x4x32_inc.f'
      end function pppmD1x4x32

      function pppmD1x4x32_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD1x4x32_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmD1x4x32_inc.f'
      end function pppmD1x4x32_qp

      end module pppmD1x4x32_generic

