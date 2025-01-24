!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module jettagging
        use types
        implicit none
        private

        include 'mxpart.f'

        integer, save, public :: jetcontent(mxpart,0:40)
!$omp threadprivate(jetcontent)

        integer, parameter, public :: bitflag(14) = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192]

      end module
