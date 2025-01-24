!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module Scalevar
    use Integration
    use types
    implicit none

    integer, public, save :: maxscalevar = -1
    logical, public, save :: doScalevar = .false.

    ! for variation about a dynamic scale and also jetptveto
    logical, public, save :: vetoScalevar = .false.

    ! for additional scale variations (resummation)
    integer, public, save :: extrascalevar = 0

    logical, public, save :: foundpow(maxParts)
!$omp threadprivate(foundpow)
    integer, public, save :: alphaspow(maxParts)
!$omp threadprivate(alphaspow)
    real(dp), public, save :: scalereweight(17)
!$omp threadprivate(scalereweight)

    real(dp), public, parameter :: &
        scalevarmult(1:8) = [2d0, 0.5d0, 2d0, 0.5d0, 1d0, 1d0, 2d0, 0.5d0], &
        facscalevarmult(1:8) = [2d0, 0.5d0, 1d0, 1d0, 2d0, 0.5d0, 0.5d0, 2d0]

    include 'initialscales.f'

    public :: usescales, updateAlphas

    public :: scaleset_pt3pt4d2

    private

    contains

    ! forces an update of all variables related to alphas
    subroutine updateAlphas(renscale_in)
        use constants
        use ieee_arithmetic
        implicit none
        include 'qcdcouple.f'
        include 'couple.f'
        include 'nlooprun.f'
        real(dp), intent(in) :: renscale_in
        real(dp) :: alphas
        real(dp) :: renscale_use

        if (renscale_in < 2d0) then
            !write (*,*) "WARNING: renscale_in = ", renscale_in
            renscale_use = 2d0
        else
            renscale_use = renscale_in
        endif

        if (.not. ieee_is_finite(renscale_in)) then
            write (*,*) "WARNING: renscale_in = ", renscale_in
            renscale_use = 2d0
        endif

        if (ieee_is_nan(renscale_in)) then
            write (*,*) "WARNING: renscale_in = ", renscale_in
            renscale_use = 2d0
        endif

        as = alphas(renscale_use,amz,nlooprun)
        ason2pi=as/(2._dp*pi)
        ason4pi=as/(4._dp*pi)
        gsq=4._dp*pi*as

    end subroutine

    subroutine usescales(renscale_in,facscale_in)
        use constants
        use ptveto, only: timelikemusq
        implicit none
        real(dp), intent(in) :: renscale_in, facscale_in
        include 'nf.f'
        include 'qcdcouple.f'
        include 'stopscales.f'
        include 'facscale.f'
        include 'scale.f'
        include 'couple.f'
        include 'nlooprun.f'

        logical :: mustUpdateAlphas

        if (scale == renscale_in .and. facscale == facscale_in) then
            return
        endif

        mustUpdateAlphas = .false.

        if (scale /= renscale_in) then
            mustUpdateAlphas = .true.
        endif
      
        scale=renscale_in
        facscale=facscale_in

        if  (scale > 60000._dp) scale=60000._dp
        if  (facscale > 60000._dp) facscale=60000._dp
        if  (scale < 1._dp) scale=1._dp
        if  (facscale < 1._dp) facscale=1._dp

        ! these are additional scales used in the t-channel single top + b routines
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale

        musq=scale**2
        if (timelikemusq) musq=-musq

        if (mustUpdateAlphas) then
            call updateAlphas(scale)
        endif

    end subroutine

    ! (pt3+pt4)/2
    function scaleset_pt3pt4d2(p)
        use types
        implicit none
        include 'mxpart.f'
        include 'constants.f'
        real(dp), intent(in) :: p(mxpart,4)
        real(dp) :: scaleset_pt3pt4d2

        real(dp) :: pt

        scaleset_pt3pt4d2 = max(2._dp, (pt(3,p) + pt(4,p))/2._dp)

    end function scaleset_pt3pt4d2
      

end module
