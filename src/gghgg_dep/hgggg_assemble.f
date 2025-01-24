!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module hgggg_assemble_generic
      implicit none
      public hgggg_assemble,hgggg_assemble_qp

      interface hgggg_assemble
      module procedure hgggg_assemble,hgggg_assemble_qp
      end interface

      contains

      function hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
     & result(hgggg_assemble_res)
      use double_precision
      use sprod_dp
      include 'Inc/hgggg_assemble_inc.f'
      end function hgggg_assemble

      function hgggg_assemble_qp(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
     & result(hgggg_assemble_res)
      use quad_precision
      use sprod_qp
      include 'Inc/hgggg_assemble_inc.f'
      end function hgggg_assemble_qp

      end module hgggg_assemble_generic

