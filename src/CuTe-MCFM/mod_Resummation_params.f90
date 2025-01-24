!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module qtResummation_params
    use types
    implicit none

    real(dp), save, public :: qtcutoff = 1.0d0

    real(dp), save, public :: qtmaxRes = 80d0
    real(dp), save, public :: qtminRes = 0d0

    real(dp), save, public :: qtmaxResexp = 80d0
    real(dp), save, public :: qtminResexp = 1.0d0

    real(dp), save, public :: transitionSwitch = 0.4d0

    logical, save, public :: enable_fixed_y = .false.
    real(dp), save, public :: fixed_y = 0.0d0
    
    logical, save, public :: enable_dsigma_dQ = .false.

    integer, save, public :: generations = 5

    logical, save, public ::  fix_alphas_nf5 = .false.

    logical, save, public :: scalevar_rapidity = .false.

    integer, save, public :: scalevar_rapidity_i = 0
!$omp threadprivate(scalevar_rapidity_i)

    real(dp), public, parameter :: &
        scalevar_rapidity_mult(1:2) = [6d0, 0.167d0]

    real(dp), public, parameter :: &
        scalevar_rapidity_mult_higgs(1:2) = [2d0, 0.5d0]
end module
