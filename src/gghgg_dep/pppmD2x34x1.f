!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmD2x34x1_generic
      implicit none
      public pppmD2x34x1,pppmD2x34x1_qp

      interface pppmD2x34x1
      module procedure pppmD2x34x1,pppmD2x34x1_qp
      end interface

      contains

      function pppmD2x34x1(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD2x34x1_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmD2x34x1_inc.f'
      end function pppmD2x34x1

      function pppmD2x34x1_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD2x34x1_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmD2x34x1_inc.f'
      end function pppmD2x34x1_qp

      end module pppmD2x34x1_generic

