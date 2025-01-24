!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module spinor
      implicit none
      public spinoru,spinoru_qp

      interface spinoru
      module procedure spinoru,spinoru_qp
      end interface

      contains

      subroutine spinoru(N,p,za,zb)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      use double_precision
      use sprod_dp
      implicit none
      include 'Inc/spinor_inc.f'
      end subroutine spinoru

      subroutine spinoru_qp(N,p,za,zb)
      use quad_precision
      use sprod_qp
      implicit none
      include 'Inc/spinor_inc.f'
      end subroutine spinoru_qp

      end module spinor

