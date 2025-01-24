!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module types
          implicit none

          public

          integer, parameter :: sp = selected_real_kind(6)
          integer, parameter :: dp = selected_real_kind(15)
          integer, parameter :: ex = selected_real_kind(18)
          integer, parameter :: qp = selected_real_kind(33)
      end module
