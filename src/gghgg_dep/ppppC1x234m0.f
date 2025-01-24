!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppppC1x234m0_generic
      implicit none
      public ppppC1x234m0,ppppC1x234m0_qp

      interface ppppC1x234m0
      module procedure ppppC1x234m0,ppppC1x234m0_qp
      end interface

      contains

      function ppppC1x234m0(p1,p2,p3,p4,za,zb) result(ppppC1x234m0_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppppC1x234m0_inc.f'
      end function ppppC1x234m0

      function ppppC1x234m0_qp(p1,p2,p3,p4,za,zb) result(ppppC1x234m0_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppppC1x234m0_inc.f'
      end function ppppC1x234m0_qp

      end module ppppC1x234m0_generic

