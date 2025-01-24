!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppppC1x234m2_generic
      implicit none
      public ppppC1x234m2,ppppC1x234m2_qp

      interface ppppC1x234m2
      module procedure ppppC1x234m2,ppppC1x234m2_qp
      end interface

      contains

      function ppppC1x234m2(p1,p2,p3,p4,za,zb) result(ppppC1x234m2_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppppC1x234m2_inc.f'
      end function ppppC1x234m2

      function ppppC1x234m2_qp(p1,p2,p3,p4,za,zb) result(ppppC1x234m2_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppppC1x234m2_inc.f'
      end function ppppC1x234m2_qp

      end module ppppC1x234m2_generic

