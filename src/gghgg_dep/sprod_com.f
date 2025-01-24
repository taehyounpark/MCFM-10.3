!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module sprod_dp
      use double_precision
      implicit none
      integer,parameter:: mxpart=14
      real(dp),save:: s(mxpart,mxpart)
!$omp threadprivate(s)
      endmodule

      module sprod_qp
      use quad_precision
      implicit none
      integer,parameter:: mxpart=14
      real(dp),save:: s(mxpart,mxpart)
!$omp threadprivate(s)
      endmodule
