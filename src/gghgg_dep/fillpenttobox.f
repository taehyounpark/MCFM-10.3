!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module fillpenttobox_generic
      implicit none
      public fillpenttobox,fillpenttobox_qp

      interface fillpenttobox
      module procedure fillpenttobox,fillpenttobox_qp
      end interface

      contains

      subroutine fillpenttobox(mtsq,Cred,I5to4)
      use double_precision
      use sprod_dp
      include 'Inc/fillpenttobox_inc.f'
      end subroutine fillpenttobox

      subroutine fillpenttobox_qp(mtsq,Cred,I5to4)
      use quad_precision
      use sprod_qp
      include 'Inc/fillpenttobox_inc.f'
      end subroutine fillpenttobox_qp

      end module fillpenttobox_generic

