!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function pppm(i,c)
      implicit none
      integer pppm,i,c

      if (i==c) then
        pppm=1
      else
        pppm=2
      endif

      return
      end

