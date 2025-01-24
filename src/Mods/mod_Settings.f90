!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module MCFMSettings
      implicit none

      logical, public, save :: newStyleHistograms = .false.

      ! set this to true to calculate linear power corrections with the resexp part
      ! then resexprange should be two equal numbers, equal to the qt cutoff
      logical, public, save :: resexp_linPC = .false.
end module
