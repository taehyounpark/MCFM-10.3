!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module SCET
    use types
    implicit none

    public :: getdynamictau
    public :: softqq

    private

    ! whether to actually do multitau
    logical, public, save :: doMultitaucut = .false.

    logical, public, save :: useQT = .false.
    real(dp), public, save, allocatable :: tcutarray(:)
    real(dp), public, save :: smallestTaucut

    real(dp), public, save, allocatable :: scetreweight(:)
!$omp threadprivate(scetreweight)

    ! this flag specifies for the above cut pieces the result of
    ! maketaucut for taucut; for each dipole contribution
    logical, public, save :: includeTaucutgrid(0:40)
!$omp threadprivate(includeTaucutgrid)

    ! used to transfor current dipole contribution information
    ! from nplotter to bookplot without changing the bookplot
    ! interface
    integer, public, save :: currentNd
!$omp threadprivate(currentNd)

    contains

    function getdynamictau(p,taucut_in)
        implicit none
        include 'types.f'
        include 'kprocess.f'
        include 'mxpart.f'
        include 'taucut.f'
        include 'nproc.f'
        include 'limits.f'

        real(dp):: getdynamictau
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in), optional :: taucut_in

        real(dp) :: puremass, pt, pttwo

        real(dp) :: taucut_use

        if (present(taucut_in)) then
            taucut_use = taucut_in
        else
            taucut_use = taucut
        endif

! Return normal value of tau if not(dynamictau)
        if (dynamictau .eqv. .false.) then
          getdynamictau=taucut_use
          return
        endif

        if (any(nprocbelow == [11,16,41]) .and. (pt34min > 0._dp)) then
            ! This branch for matching corrections in N4LL qt calculations, where pt34min > 0
            ! The definition for W+jet and Z+jet divides by 1050 to match previous non-dynamic results.
            getdynamictau = sqrt(pttwo(3,4,p)**2+puremass(p(3,:)+p(4,:))**2)/1050._dp
        elseif (any(nprocbelow == [11,16,41,204,210]))  then
            getdynamictau = sqrt(pt(5,p)**2+puremass(p(3,:)+p(4,:))**2)
        elseif ( (kcase == kW_only) .or. (kcase == kW_1jet) &
             .or.(kcase == kZ_only) .or. (kcase == kZ_1jet) &
             .or.(kcase == kggfus0) .or. (kcase == kggfus1) &
             .or.(kcase == kggfus2) &
             .or.(kcase == kgamgam) .or. (kcase == kgmgmjt) ) then
          getdynamictau=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2 &
                            -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
        elseif ( (kcase == kWgamma) .or. (kcase == kWgajet) .or. (kcase==kWgajew) &
             .or. (kcase == kZgamma) .or. (kcase == kZgajet) ) then
          getdynamictau=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2 &
                            -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
        elseif ( (kcase == kWHbbar) .or. (kcase == kWH1jet)  &
             .or.(kcase == kZHbbar) .or. (kcase == kZH1jet)  &
             .or.(kcase == kWHgaga) .or. (kcase == kZHgaga)  &
             .or.(kcase == kWWqqbr) .or. (kcase == kWW_jet)  &
             .or.(kcase == kWZbbar) .or. (kcase == kWZ_jet)  &
             .or.(kcase == kZZlept) .or. (kcase == kZZ_jet) ) then
          getdynamictau=sqrt((p(3,4)+p(4,4)+p(5,4)+p(6,4))**2-(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2 &
                            -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2-(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2)
        elseif (kcase == kbq_tpq .or. kcase==kbq_tpq_jet) then
            getdynamictau = puremass(+p(1,:)-p(6,:))
        else
          write(6,*) 'Invalid process for dynamic taucut!'
          stop
        endif
          
        getdynamictau=getdynamictau*taucut_use
    end function
      
    subroutine softqq(order,soft1,soft2,logy_in)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'scale.f'
        include 'nf.f'
        include 'scet_const.f'
        integer, intent(in) :: order
        real(dp), intent(out) :: soft1(-1:1),soft2(-1:3)
        real(dp), intent(in), optional :: logy_in
        real(dp) :: logy
        include 'scet_beta.f'

        if (present(logy_in)) then
            logy = logy_in
        else
            logy = 0._dp
        endif


        soft1(-1) = CF*(pisqo6 - logy**2)
        soft1(0) = 4._dp*CF*logy
        soft1(1) = -8*CF

        soft2(-1) =  CF*be0 * ( &
          - 5._dp/27._dp                          &
          - 14._dp/9._dp*logy                     &
          - 5._dp/6._dp*logy**2                   &
          - 1._dp/6._dp*logy**3                   &
          + 29._dp/6._dp*zeta3                    &
          - 37._dp/12._dp*zeta2                   &
          + zeta2*logy )                          &
                  + CF*CA * (  &
          - 160._dp/27._dp                        &
          - 16._dp/9._dp*logy                     &
          - 2._dp/3._dp*logy**2                   &
          + 7._dp*zeta3*logy                      &
          + 2._dp*zeta2                           &
          + zeta2*logy**2                         &
          + 22._dp/5._dp*zeta2**2 )               &
                  + CF**2 * (  &
          + 1._dp/2._dp*logy**4                   &
          - 32._dp*zeta3*logy                     &
          - 9._dp*zeta2*logy**2                   &
          - 27._dp/10._dp*zeta2**2 )
        soft2(0) = &
         + 28._dp/9._dp*CF*be0 &
         + 10._dp/3._dp*CF*be0*logy &
         + CF*be0*logy**2 &
         - 2._dp*CF*zeta2*be0 &
         + 32._dp/9._dp*CF*CA &
         + 8._dp/3._dp*CF*CA*logy &
         - 14._dp*CF*CA*zeta3 &
         - 4._dp*CF*CA*zeta2*logy &
         - 4._dp*CF**2*logy**3 &
         + 64._dp*CF**2*zeta3 &
         + 36._dp*CF**2*zeta2*logy
        soft2(1) = &
          - 20._dp/3._dp*CF*be0 &
          - 4._dp*CF*be0*logy &
          - 16._dp/3._dp*CF*CA &
          + 8._dp*CF*CA*zeta2 &
          + 24._dp*CF**2*logy**2 &
          - 72._dp*CF**2*zeta2
        soft2(2) = &
          + 4._dp*CF*be0 &
          - 48._dp*CF**2*logy
        soft2(3) = 32*CF**2

    end subroutine

end module
