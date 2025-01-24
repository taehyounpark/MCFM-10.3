!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module gghgg_dep_params
      implicit none
      real(selected_real_kind(precision(1.0d0))), save::
     &  mt_dp,mb_dp,mt_yuk_dp,mb_yuk_dp,as_dp,vevsq_dp
!$omp threadprivate(mt_dp,mb_dp,mt_yuk_dp,mb_yuk_dp,as_dp,vevsq_dp)
      endmodule
