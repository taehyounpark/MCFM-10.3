!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine Qformfill
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nflav.f'
      include 'zeta.f'
      include 'Qform.f'

      Qform(0)=cone
      Qform(1)=CF*(-4._dp+pisq/2._dp-1.5_dp*impi)
      Qform(2)=CF**2*(255._dp/32._dp-pisq*23._dp/12._dp
     & +pi**4*23._dp/360._dp-zeta3*15._dp/2._dp
     & +impi*(45._dp/8._dp-pisq/4._dp-6._dp*zeta3))
     & +NC*CF*(-1535._dp/96._dp+pisq*215._dp/108._dp-(pi**4)/120._dp
     & +zeta3*59._dp/9._dp
     & +impi*(-2545._dp/216._dp-pisq*11._dp/72._dp+zeta3*13._dp/2._dp))
     & +nflav*CF*(127._dp/48._dp-pisq*19._dp/54._dp+zeta3*4._dp/9._dp
     & +impi*(209._dp/108._dp+pisq/36._dp))
      end
