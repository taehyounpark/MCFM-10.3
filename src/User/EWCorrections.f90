!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module EWCorrections
      use types
      implicit none

      public

      contains

      function deltaZEW_central(qt)
          implicit none
          real(dp) :: deltaZEW_central
          real(dp), intent(in) :: qt
      end function

end module
