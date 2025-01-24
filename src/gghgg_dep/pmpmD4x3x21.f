!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmD4x3x21_generic
      implicit none
      public pmpmD4x3x21,pmpmD4x3x21_qp

      interface pmpmD4x3x21
      module procedure pmpmD4x3x21,pmpmD4x3x21_qp
      end interface

      contains

      function pmpmD4x3x21(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pmpmD4x3x21_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmD4x3x21_inc.f'
      end function pmpmD4x3x21

      function pmpmD4x3x21_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pmpmD4x3x21_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmD4x3x21_inc.f'
      end function pmpmD4x3x21_qp

      end module pmpmD4x3x21_generic

