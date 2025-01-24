!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_integralfill_generic
      implicit none
      public hgggg_integralfill,hgggg_integralfill_qp

      interface hgggg_integralfill
      module procedure hgggg_integralfill,hgggg_integralfill_qp
      end interface

      contains

      subroutine hgggg_integralfill(p1,p2,p3,p4,mtsq,Dint,Cint,Bint)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_integralfill_inc.f'
      end subroutine hgggg_integralfill

      subroutine hgggg_integralfill_qp(p1,p2,p3,p4,mtsq,Dint,Cint,Bint)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_integralfillqp_inc.f'
      end subroutine hgggg_integralfill_qp

      end module hgggg_integralfill_generic

