!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module ppppD1x2x34_generic
      implicit none
      public ppppD1x2x34,ppppD1x2x34_qp

      interface ppppD1x2x34
      module procedure ppppD1x2x34,ppppD1x2x34_qp
      end interface

      contains

      function ppppD1x2x34(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppppD1x2x34_res)
      use double_precision
      use sprod_dp
      include 'Inc/ppppD1x2x34_inc.f'
      end function ppppD1x2x34

      function ppppD1x2x34_qp(p1,p2,p3,p4,mtsq,za,zb,Cred,I5to4) result(ppppD1x2x34_res)
      use quad_precision
      use sprod_qp
      include 'Inc/ppppD1x2x34_inc.f'
      end function ppppD1x2x34_qp

      end module ppppD1x2x34_generic

