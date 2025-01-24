!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqbAjF(order,AjF)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'Qform.f'
      complex(dp):: AjF(1:10)
      integer:: order

      AjF(1:6)=czip
      AjF(7)=cmplx(2._dp,kind=dp)
      AjF(8)=cmplx(2._dp,kind=dp)
      AjF(9)=-cmplx(1._dp,kind=dp)
      AjF(10)=-cmplx(1._dp,kind=dp)

      AjF(7:10)=Qform(order)*AjF(7:10)
      end
