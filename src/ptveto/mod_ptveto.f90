!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

module ptveto
    use types
    implicit none

    logical, public, save :: usept

    ! use BNR version of d3veto
    logical, public, save :: useBNR

    ! value of ptveto
    real(dp), public, save :: jetptveto

    ! use Banfi et al. version of d3veto
    logical, public, save :: useBanfid3veto

    ! kappa (BNR)
    real(dp), public, save :: d3vetokappa

    ! R0 (Banfi et al.)
    real(dp), public, save :: d3vetoR0
    
    logical, public, save :: gghsinglestep
    
    logical, public, save :: timelikemusq

! Value of x above which regular term is replaced by series expansion
    real(dp), public, parameter :: xlim = 0.9_dp
    
    private

end module
