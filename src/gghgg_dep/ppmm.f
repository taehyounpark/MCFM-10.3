!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ppmm(i,c1,c2)
      implicit none
      integer ppmm,i,c1,c2

      if ((i==c1) .or. (i==c2)) then
        ppmm=1
      else
        ppmm=2
      endif

      return
      end
