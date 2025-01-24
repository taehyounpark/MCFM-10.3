!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Li22(x1,x2)
      use handyG
      implicit none
      include 'types.f'
      real(dp):: Li22
      complex(dp):: zLi22
      real(dp):: x1,x2

      zLi22=G([ 0._dp, (1._dp/x1), 0._dp, (1._dp/x1/x2), 1.0_dp])
      Li22=real(zLi22,kind=dp)
      return
      end
