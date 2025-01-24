!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module ResummationTransition
      use types
      implicit none

      public :: transition

      private


      contains

      function slog(x,l,r,red)
          implicit none
          real(dp) :: slog

          real(dp), intent(in) :: x,l,r,red

          real(dp) :: m, w

          m = (r+l)/2._dp
          w = (r-l)/2._dp

          slog = 1._dp/(1._dp + exp(log((1-red)/red)*(x-m)/w))
      end function

      ! below ptmin exactly 1
      ! above ptmax exactly 0
      function transition(pt,ptmin,ptmax,red)
          implicit none
          real(dp) :: transition

          real(dp), intent(in) :: pt,ptmin,ptmax
          real(dp), intent(in), optional :: red
          real(dp) :: red_use

          if (present(red)) then
              red_use = red
          else
              red_use = 0.001d0
          endif

          if (pt < ptmin) then
              transition = 1d0
              return
          else
              transition = slog(pt,ptmin,ptmax,red)/slog(ptmin,ptmin,ptmax,red)
              if (transition < 1d-3) then
                  transition = 0d0
              endif
          endif

      end function

end module
