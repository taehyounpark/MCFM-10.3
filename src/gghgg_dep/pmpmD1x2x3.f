!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pmpmD1x2x3_generic
      implicit none
      public pmpmD1x2x3,pmpmD1x2x3_qp

      interface pmpmD1x2x3
      module procedure pmpmD1x2x3,pmpmD1x2x3_qp
      end interface

      contains

      function pmpmD1x2x3(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pmpmD1x2x3_res)
      use double_precision
      use sprod_dp
      include 'Inc/pmpmD1x2x3_inc.f'
      end function pmpmD1x2x3

      function pmpmD1x2x3_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pmpmD1x2x3_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pmpmD1x2x3_inc.f'
      end function pmpmD1x2x3_qp

      end module pmpmD1x2x3_generic

