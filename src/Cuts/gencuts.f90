!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module m_gencuts
    implicit none

    public :: gencuts
    public :: gencuts_user
    public :: reweight_user

    logical, public :: enable_reweight_user = .false.

    interface
        module logical function gencuts_user(pjet, njets)
            use types
            implicit none
            include 'mxpart.f'
            real(dp), intent(in) :: pjet(mxpart,4)
            integer, intent(in) :: njets
        end function

        module function reweight_user(pjet)
            use types
            implicit none
            include 'mxpart.f'
            real(dp), intent(in) :: pjet(mxpart,4)
            real(dp) :: reweight_user
        end function
    end interface

    private

    contains

    logical function gencuts(pjet, njets)
        use types
        implicit none
        include 'first.f'
        include 'lhcb.f'
        include 'mxpart.f'
    
        real(dp), intent(in) :: pjet(mxpart,4)
        integer, intent(in) :: njets
    
        logical :: lhcb_cuts
        logical :: gencuts_input
    
        gencuts = .false.
    
        if (first) then
            first = .false.
        endif
    
        ! LHCb cuts
        if (cut_mode > 0) then
          lhcb_pass(1) = .false.
          lhcb_pass(2) = .false.
          if (dir_mode == 0) then
             gencuts = lhcb_cuts(pjet, njets, 1)
             gencuts = lhcb_cuts(pjet, njets, 2) .and. gencuts
          else
             gencuts = lhcb_cuts(pjet, njets, dir_mode)
          endif
          return
        endif

        if (gencuts_user(pjet,njets) .eqv. .true.) then
            gencuts = .true.
            return
        endif
    
        gencuts = gencuts_input(pjet,njets)
    
    end function
end module
