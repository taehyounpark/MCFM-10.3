!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pppmD2x1x43_generic
      implicit none
      public pppmD2x1x43,pppmD2x1x43_qp

      interface pppmD2x1x43
      module procedure pppmD2x1x43,pppmD2x1x43_qp
      end interface

      contains

      function pppmD2x1x43(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD2x1x43_res)
      use double_precision
      use sprod_dp
      include 'Inc/pppmD2x1x43_inc.f'
      end function pppmD2x1x43

      function pppmD2x1x43_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(pppmD2x1x43_res)
      use quad_precision
      use sprod_qp
      include 'Inc/pppmD2x1x43_inc.f'
      end function pppmD2x1x43_qp

      end module pppmD2x1x43_generic

