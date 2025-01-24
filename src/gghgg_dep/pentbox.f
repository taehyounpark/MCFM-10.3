!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module pentbox_generic
      implicit none
      public pentbox,pentbox_qp

      interface pentbox
      module procedure pentbox,pentbox_qp
      end interface

      contains

      subroutine pentbox(p1,p2,p3,p4,mtsq,C)
      use double_precision
      use sprod_dp
      include 'Inc/pentbox_inc.f'
      end subroutine pentbox

      subroutine pentbox_qp(p1,p2,p3,p4,mtsq,C)
      use quad_precision
      use sprod_qp
      include 'Inc/pentbox_inc.f'
      end subroutine pentbox_qp

      end module pentbox_generic

