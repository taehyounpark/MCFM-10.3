!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop2_scet_heavy_prod
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use singletop2_scet_light, only: singletop2_scet_tree, &
          singletop2_scet_tree_ub, singletop2_scet_tree_bu
      use singletop2_nnlo_vars
    use LHAPDF
    implicit none
    private

    public :: lumxmsq_singletop_prod
    public :: passed_taucut_heavyprod

    public :: qqb_tbb_g_heavy_all
    public :: singletop2_scet_heavy_prod_gs_all

    public :: singletop2_scet_virt_heavy_prod_all
    public :: singletop2_scet_heavy_prod_z

    public :: virtqqb_heavy

    contains

    subroutine singletop2_scet_heavy_prod_gs_all(p,ndmx,msqall)
        use singletop2_scale_m
        use singletop2_nnlo_vars
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'nwz.f'
        include 'ptilde.f'
        include 'qqgg.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ndmx
        real(dp), intent(out) :: msqall(ndmx, -nf:nf, -nf:nf, max_bcontrib, max_corr_on_beam)

        real(dp) :: dummyv(-nf:nf,-nf:nf), dsubv
        real(dp) :: msqborn(-nf:nf,-nf:nf), sub(4)
        integer, parameter :: noglue(4) = [-3,-1,2,4]

        external donothing_gvec

        ndmax = 6

        if (nwz == -1) then
            error stop "nwz = -1 will need adjustments in singletop2_scet_heavy_prod_gs"
        endif

        msqall = 0._dp
        msqborn = 0

        qqproc = .true.

        corr_islight = .false.
        corr_beam1 = .false.
        corr_on_beam = 2

        !if (corr_on_beam == 2) then

        call dips_mass(1,p,2,5,3,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_ub,donothing_gvec)
        msqall(1,noglue,5, 1,2) = 2._dp*cf*sub(qq)*msqborn(noglue,5)

        call dips_mass(2,p,3,5,2,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_ub,donothing_gvec)
        msqall(2,noglue,5, 1,2) = 2._dp*cf*sub(qq)*msqborn(noglue,5)

        ! binned as final state b~ contrib
        call dips_mass(5,p,2,5,1,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_ub,donothing_gvec)
        msqall(5,noglue,0, 3,2) = 2._dp*tr*sub(qg)*msqborn(noglue,5)

        !else

        corr_islight = .false.
        corr_beam1 = .true.
        corr_on_beam = 1

        call dips_mass(3,p,1,5,3,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_bu,donothing_gvec)
        msqall(3,5,noglue, 1,1) = 2._dp*cf*sub(qq)*msqborn(5,noglue)

        call dips_mass(4,p,3,5,1,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_bu,donothing_gvec)
        msqall(4,5,noglue, 1,1) = 2._dp*cf*sub(qq)*msqborn(5,noglue)

        ! binned as final state b~ contrib
        call dips_mass(6,p,1,5,2,sub,dsubv,msqborn,dummyv,singletop2_scet_tree_bu,donothing_gvec)
        msqall(6,0,noglue, 3,1) = 2._dp*tr*sub(qg)*msqborn(5,noglue)


        !endif

    end subroutine

    subroutine lumxmsq_singletop_prod(p,xx,z1,z2,QB,order,xmsq,central)
        use types
        use SCET
        use SCET_Beamfunctions
        use SCET_Jet
        use singletop2_scet_heavy_decay, only : softfun, hardfun
        use singletop2_nnlo_vars, only: maxbeams, beams_enabled
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'taucut.f'
        include 'masses.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'scet_const.f'
        include 'beamtype.f'
    
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in) :: xx(2), z1, z2, QB(2)
        integer, intent(in) :: order
        real(dp), intent(out) :: xmsq
        logical, intent(in) :: central
    
        real(dp) :: beama0_light(-5:5), beama0_heavy(-5:5)
        real(dp) :: beamb0_light(-5:5), beamb0_heavy(-5:5)
        real(dp) :: beama1(-5:5, -1:1)
        real(dp) :: beama2(-5:5, -1:3)
        real(dp) :: beamb1(-5:5, -1:1)
        real(dp) :: beamb2(-5:5, -1:3)

        real(dp) :: soft(2,0:4), jet(2,0:4)
        real(dp) :: hard_ub(0:2), hard_bu(0:2), hard_ubarb(0:2), hard_bubar(0:2)

        ! functions
        real(dp) :: massvec
        real(dp) :: x, tauc

        real(dp), allocatable :: taucuts(:)
        real(dp), allocatable :: xmsqs(:)

        integer :: k,m

        xmsq = 0._dp


        tauc = taucut
        if (central .and. doMultitaucut) then
            allocate(taucuts(1+size(tcutarray)))
            allocate(xmsqs(1+size(tcutarray)))
            do m=1,size(tcutarray)
                !taucuts(m+1) = getdynamictau(p, tcutarray(m))
                taucuts(m+1) = tcutarray(m)
            enddo
        else
            allocate(taucuts(1))
            allocate(xmsqs(1))
        endif

        taucuts(1) = tauc

        xmsqs(:) = 0._dp

        soft = softfun(order)

        if (any(beams_enabled(1:maxbeams) == 2)) then
            call fdist(ih1,xx(1),facscale_beam1_islight_onheavy,beama0_light,1)
            call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy,beamb0_heavy,2)

            call xbeam1bis_new(ih2,z2,xx(2),QB(2),beamb1, 2, renscale_beam2_isheavy_onheavy, &
                    facscale_beam2_isheavy_onheavy,disttau=.false.)
            if (order >= 2) then
                call xbeam2bis_new(ih2,z2,xx(2),QB(2),beamb2, 2, &
                        renscale_beam2_isheavy_onheavy, &
                        facscale_beam2_isheavy_onheavy,disttau=.false.)
            endif
        endif

        if (any(beams_enabled(1:maxbeams) == 1)) then
            call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy,beama0_heavy,1)
            call fdist(ih2,xx(2),facscale_beam2_islight_onheavy,beamb0_light,2)

            call xbeam1bis_new(ih1,z1,xx(1),QB(1),beama1, 1, renscale_beam1_isheavy_onheavy, &
                facscale_beam1_isheavy_onheavy,disttau=.false.)
            if (order >= 2) then
                call xbeam2bis_new(ih1,z1,xx(1),QB(1),beama2, 1, &
                    renscale_beam1_isheavy_onheavy, &
                    facscale_beam1_isheavy_onheavy,disttau=.false.)
            endif
        endif


        ! 1) Carefully compared at one-loop against existing routines
        ! remember to set epinv=epinv2=0 and add a pi^2/6 term when comparing.
        ! 2) Numerically confirmed that at 1-loop there is no scale-dependence left
        ! when using the assemble_production_pieces routine.

        if (any(beams_enabled(1:maxbeams) == 2)) then

! ub channel
        x = massvec(-p(1,:)-p(6,:)) / mt**2
        hard_ub    = real(hardfun(order, renscale_beam2_isheavy_onheavy, x, p, 1,2,3,4,5,6, production=.true.),dp)
        hard_ubarb = real(hardfun(order, renscale_beam2_isheavy_onheavy, x, p, 6,2,3,4,5,1, production=.true.),dp)

        !hard_ub(1) = hard_ub(1) - pi**2/6._dp * cf * hard_ub(0)

        k = 5
        ! convert to Laplace-space log-notation
        jet(1,0) = beamb1(k,-1)
        jet(1,1) = beamb1(k,0)
        jet(1,2) = beamb1(k,1)/2._dp

        jet(2,0) = beamb2(k,-1)
        jet(2,1) = beamb2(k,0)
        jet(2,2) = beamb2(k,1)/2._dp
        jet(2,3) = beamb2(k,2)/3._dp
        jet(2,4) = beamb2(k,3)/4._dp
        xmsqs = xmsqs + (beama0_light(2) + beama0_light(4))*hard_ub(0)* &
            assemble_production_pieces(order, x, &
            renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
            hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

        xmsqs = xmsqs + (beama0_light(-1) + beama0_light(-3))*hard_ubarb(0)* &
            assemble_production_pieces(order, x, &
            renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
            hard_ubarb/hard_ubarb(0), jet, soft, beamb0_heavy(k))

        endif

        ! debug for ug channel only; for mu=mt; log(1-x) is fixed by taucut-independence of ub channel
        !xmsqs = beama0_light(2)*hard_ub(0)*(as/4._dp/pi)*(jet(1,0) + log(1-x)*jet(1,1) + log(taucuts)*jet(1,1))

! bu channel
        if (any(beams_enabled(1:maxbeams) == 1)) then

        x = massvec(-p(2,:)-p(6,:)) / mt**2
        hard_bu    = real(hardfun(order, renscale_beam1_isheavy_onheavy, x, p, 2,1,3,4,5,6, production=.true.),dp)
        hard_bubar = real(hardfun(order, renscale_beam1_isheavy_onheavy, x, p, 6,1,3,4,5,2, production=.true.),dp)

        k = 5
        ! convert to Laplace-space log-notation
        jet(1,0) = beama1(k,-1)
        jet(1,1) = beama1(k,0)
        jet(1,2) = beama1(k,1)/2._dp

        jet(2,0) = beama2(k,-1)
        jet(2,1) = beama2(k,0)
        jet(2,2) = beama2(k,1)/2._dp
        jet(2,3) = beama2(k,2)/3._dp
        jet(2,4) = beama2(k,3)/4._dp
        xmsqs = xmsqs + (beamb0_light(2) + beamb0_light(4))*hard_bu(0)* &
            assemble_production_pieces(order, x, &
            renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucuts, &
            hard_bu/hard_bu(0), jet, soft, beama0_heavy(5))

        xmsqs = xmsqs + (beamb0_light(-1) + beamb0_light(-3))*hard_bubar(0)* &
            assemble_production_pieces(order, x, &
            renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucuts, &
            hard_bubar/hard_bubar(0), jet, soft, beama0_heavy(5))

        endif

        ! TEST SCALE DEPENDENCE

!       ! only vary renormalizations cale
!       newscale = (1d0 + 1d-8)*renscale_beam2_isheavy_onheavy
!       call xbeam1bis_new(ih2,z2,xx(2),QB(2),beamb1, 2, newscale, &
!               facscale_beam2_isheavy_onheavy,disttau=.false.)
!       call xbeam2bis_new(ih2,z2,xx(2),QB(2),beamb2, 2, newscale, &
!               facscale_beam2_isheavy_onheavy,disttau=.false.)

!       call fdist(ih1,xx(1),facscale_beam2_isheavy_onheavy,beama0_light,1)

!       x = massvec(-p(1,:)-p(6,:)) / mt**2
!       hard_ub = real(hardfun(order, newscale, x, p, 1,2,3,4,5,6, production=.true.),dp)
!       k = 5
!       jet(1,0) = beamb1(k,-1)
!       jet(1,1) = beamb1(k,0)
!       jet(1,2) = beamb1(k,1)/2._dp ! convert to log-notation

!       jet(2,0) = beamb2(k,-1)
!       jet(2,1) = beamb2(k,0)
!       jet(2,2) = beamb2(k,1)/2._dp
!       jet(2,3) = beamb2(k,2)/3._dp
!       jet(2,4) = beamb2(k,3)/4._dp
!       xmsqs = beama0_light(2)*hard_ub(0)*assemble_production_pieces(order, x, &
!           newscale, as_heavy_beam2, taucuts, &
!           hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

!       test2 = xmsqs(1)

!       write (*,*) test2/test1, test2-test1, test2, test1
!       ! this returns 1.0, 0.0, ..

!       ! compute just order 1 piece for NNLO check
!       xmsqs = beama0_light(2)*hard_ub(0)*assemble_production_pieces(1, x, &
!           newscale, as_heavy_beam2, taucuts, &
!           hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

!       write (*,*) "NNLO CHECK", 1d8*(test2-test1)/(xmsqs(1)*(as_heavy_beam2/4._dp/pi)*be0*2)
!       pause
!       ! END TEST SCALE DEPENDENCE


!       ! this checks the one-loop hard function
!       epinv = 0._dp
!       epinv2 = 0._dp
!       call singletop2_scet_virt_heavy_prod(p,msqv)
!       msqv(2,5) = msqv(2,5) + pi**2/6._dp * cf * as_heavy_beam2/4._dp/pi * hard_ub(0)
!       write (*,*) "DEBUG"
!       write (*,*) "hard_ub(1), msqv(2,5) ", hard_ub(1), msqv(2,5)
!       write (*,*) "ratio, diff", as_heavy_beam2/4/pi*hard_ub(1)/msqv(2,5), &
!           (as_heavy_beam2/4._dp/pi*hard_ub(1) - msqv(2,5))/ hard_ub(0)
!       ! this return 1.0, 0.0
!       pause

        xmsq = xmsqs(1)

        if (central .and. doMultitaucut) then
            scetreweight(:) = 0._dp
            if (xmsq /= 0._dp) then
                scetreweight(:) = xmsqs(2:) / xmsq
            endif
        endif

    end subroutine

    ! this is the same as assemble_decay_pieces but with
    ! leading term of jet function to be j0 instead of 1
    ! so we can feed it the leading beam function term, i.e. PDF
    function assemble_production_pieces(order, x, scale, as, taucut, h, j, s, j0)
        use types
        use constants
        implicit none
        include 'masses.f'
        include 'kpart.f'

        integer, intent(in) :: order
        real(dp), intent(in) :: x, scale, as
        real(dp), intent(in) :: taucut(:)
        real(dp), intent(in) :: h(0:2)
        real(dp), intent(in) :: j(2,0:4)
        real(dp), intent(in) :: s(2,0:4)
        real(dp), intent(in) :: j0
        real(dp) :: assemble_production_pieces(size(taucut(:)))

        real(dp) :: LogMtMu
        real(dp) :: full(0:2,0:4)
        real(dp) :: LSX, LJX

        LogMtMu = log(mt/scale)

        ! expressing everything in terms of LJX/LSX + log(tau) allows
        ! us to repurpose this assembly function for the production process
        LJX = 2*LogMtMu + log(1-x)
        LSX = LogMtMu

        full(:,:) = 0._dp

        if (h(0) /= 1._dp) then
            write (*,*) "WARNING: bad hard function normalization: ", h(0)
        endif

        if (coeffonly) then
            full(0,0) = 0._dp
        else
            full(0,0) = j0
        endif

        if ((order == 1) .or. (order >= 1 .and. coeffonly .eqv. .false.)) then
            full(1,2) = j(1,2) + j0*s(1,2)
            full(1,1) = j(1,1) + 2*LJX*j(1,2) + j0*(s(1,1) + 2*LSX*s(1,2))
            full(1,0) = j(1,0) + LJX*(j(1,1) + LJX*j(1,2)) + &
                j0*(h(1) + s(1,0) + LSX*(s(1,1) + LSX*s(1,2)))
        endif

        ! The additional LogMtMu terms in this routine are from
        ! the flavor scheme conversion.
        ! In the MM notebook we convert the hard function to 4-flavor scheme
        ! and renormalize everything consistently in it. Structurally
        ! this is H(Nf4) * S(Nf4) * J(Nf4). We re-use the Nf5 beam function,
        ! so only the renormalized hard and soft functions obtain a
        ! scheme conversion, which leads to the LogMtMu below.

        ! Fixed: this term must is zero!

        LogMtMu = 0._dp
        if (order >= 2) then
            full(2,4) = j(2,4) + j(1,2)*s(1,2) + j0*s(2,4)
            full(2,3) = j(2,3) + 4*LJX*j(2,4) + j(1,2)*s(1,1) + j(1,1)*s(1,2) + &
                2*LJX*j(1,2)*s(1,2) + 2*LSX*j(1,2)*s(1,2) + j0*s(2,3) + 4*j0*LSX*s(2,4)
            full(2,2) = j(2,2) + 3*LJX*j(2,3) + 6*LJX**2*j(2,4) + j(1,2)*s(1,0) + &
                j(1,1)*s(1,1) + 2*LJX*j(1,2)*s(1,1) + LSX*j(1,2)*s(1,1) + &
                (4*j0*LogMtMu*s(1,2))/3._dp + j(1,0)*s(1,2) + LJX*j(1,1)*s(1,2) + &
                2*LSX*j(1,1)*s(1,2) + LJX**2*j(1,2)*s(1,2) + 4*LJX*LSX*j(1,2)*s(1,2) + &
                LSX**2*j(1,2)*s(1,2) - (2*Pi**2*j(1,2)*s(1,2))/3._dp + &
                h(1)*(j(1,2) + j0*s(1,2)) + j0*s(2,2) + 3*j0*LSX*s(2,3) + &
                6*j0*LSX**2*s(2,4)
            full(2,1) = j(2,1) + 2*LJX*j(2,2) + 3*LJX**2*j(2,3) + 4*LJX**3*j(2,4) + &
                j(1,1)*s(1,0) + 2*LJX*j(1,2)*s(1,0) + (4*j0*LogMtMu*s(1,1))/3._dp + &
                j(1,0)*s(1,1) + LJX*j(1,1)*s(1,1) + LSX*j(1,1)*s(1,1) + &
                LJX**2*j(1,2)*s(1,1) + 2*LJX*LSX*j(1,2)*s(1,1) - &
                (Pi**2*j(1,2)*s(1,1))/3._dp + (8*j0*LogMtMu*LSX*s(1,2))/3._dp + &
                2*LSX*j(1,0)*s(1,2) + 2*LJX*LSX*j(1,1)*s(1,2) + LSX**2*j(1,1)*s(1,2) - &
                (Pi**2*j(1,1)*s(1,2))/3._dp + 2*LJX**2*LSX*j(1,2)*s(1,2) + &
                2*LJX*LSX**2*j(1,2)*s(1,2) - (2*LJX*Pi**2*j(1,2)*s(1,2))/3._dp - &
                (2*LSX*Pi**2*j(1,2)*s(1,2))/3._dp + 8*zeta3*j(1,2)*s(1,2) + &
                h(1)*(j(1,1) + 2*LJX*j(1,2) + j0*s(1,1) + 2*j0*LSX*s(1,2)) + &
                j0*s(2,1) + 2*j0*LSX*s(2,2) + 3*j0*LSX**2*s(2,3) + 4*j0*LSX**3*s(2,4)
            full(2,0) = h(1)*(j(1,0) + LJX*(j(1,1) + LJX*j(1,2))) + j(2,0) + LJX*j(2,1) + &
                LJX**2*j(2,2) + LJX**3*j(2,3) + LJX**4*j(2,4) + j(1,0)*s(1,0) + &
                LJX*j(1,1)*s(1,0) + LJX**2*j(1,2)*s(1,0) + LSX*j(1,0)*s(1,1) + &
                LJX*LSX*j(1,1)*s(1,1) - (Pi**2*j(1,1)*s(1,1))/6._dp + &
                LJX**2*LSX*j(1,2)*s(1,1) - (LJX*Pi**2*j(1,2)*s(1,1))/3._dp + &
                2*zeta3*j(1,2)*s(1,1) + LSX**2*j(1,0)*s(1,2) + &
                LJX*LSX**2*j(1,1)*s(1,2) - (LSX*Pi**2*j(1,1)*s(1,2))/3._dp + &
                2*zeta3*j(1,1)*s(1,2) + LJX**2*LSX**2*j(1,2)*s(1,2) - &
                (2*LJX*LSX*Pi**2*j(1,2)*s(1,2))/3._dp - (Pi**4*j(1,2)*s(1,2))/90._dp + &
                4*LJX*zeta3*j(1,2)*s(1,2) + 4*LSX*zeta3*j(1,2)*s(1,2) + &
                (4*j0*LogMtMu*(h(1) + s(1,0) + LSX*(s(1,1) + LSX*s(1,2))))/3._dp + &
                j0*(h(2) + h(1)*(s(1,0) + LSX*(s(1,1) + LSX*s(1,2))) + s(2,0) + &
                   LSX*s(2,1) + LSX**2*s(2,2) + LSX**3*s(2,3) + LSX**4*s(2,4))

        endif

        assemble_production_pieces = 0._dp

        if (order == 1 .or. (order > 1 .and. coeffonly .eqv. .false.)) then
            assemble_production_pieces = full(0,0) + &
                as/4._dp/pi*( &
                    full(1,0) + full(1,1)*log(taucut) + full(1,2)*log(taucut)**2 )
        endif

        if (order >= 2) then
            assemble_production_pieces = assemble_production_pieces + (as/4._dp/pi)**2 * ( &
                full(2,0) + full(2,1)*log(taucut) + full(2,2)*log(taucut)**2 &
                + full(2,3)*log(taucut)**3 + full(2,4)*log(taucut)**4 )
        endif

    end function

    subroutine singletop2_scet_heavy_prod_z(p,z)
        use singletop2_nnlo_vars
      implicit none

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      include 'nwz.f'
      include 'masses.f'

      real(dp), intent(in) :: p(mxpart,4), z
      real(dp) :: dot,ii_qg,if_mqq, fi_mqq

      integer :: is
      real(dp) :: mbar15, mbar25, xl25, xl15

      xl25=log((-two*dot(p,1,6)+mt**2)/renscale_beam2_isheavy_onheavy**2)
      xl15=log((-two*dot(p,2,6)+mt**2)/renscale_beam1_isheavy_onheavy**2)
      mbar15=mt/sqrt(-two*dot(p,2,6)+mt**2)
      mbar25=mt/sqrt(-two*dot(p,1,6)+mt**2)

      Q1 = zip
      Q2 = zip
      B1 = zip
      B2 = zip

      if (nwz == -1) then
          error stop "nwz == -1, might need adjustments in singletop2_scet_heavy_prod_z"
      endif

      do is=1,3
         !if (corr_on_beam == 2) then
             B2(b,b,q,is) = as_heavy_beam2/2/pi * cf*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
             B2(q,g,q,is) = as_heavy_beam2/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam2_isheavy_onheavy**2),is)
         !else
             B1(b,b,q,is) = as_heavy_beam1/2/pi * cf*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
             B1(q,g,q,is) = as_heavy_beam1/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam1_isheavy_onheavy**2),is)
         !endif
      enddo

    end subroutine

    function passed_taucut_heavyprod(pparton,scetreweight_local,taucut_in)
        use SCET
        use singletop2_nnlo_vars
        implicit none
        include 'constants.f'
        include 'mxpart.f'
        include 'npart.f'
        include 'nqcdjets.f'
        include 'taucut.f'
        include 'plabel.f'
        include 'kprocess.f'
        include 'masses.f'

        logical :: passed_taucut_heavyprod
        real(dp), intent(in) :: pparton(mxpart,4)
        real(dp), intent(inout), optional :: scetreweight_local(:)
        real(dp), intent(in), optional :: taucut_in

        integer j
        real(dp) :: tau,tauc
        real(dp) :: qsq

        real(dp) :: dotvec, massvec, nn(4)

        logical :: bin
        common/bin/bin

        passed_taucut_heavyprod = .false.

        if (present(taucut_in)) then
            tauc = taucut_in
        else
            tauc = taucut
        endif

        if (corr_on_beam == 1) then
            qsq = massvec(-pparton(2,:)-pparton(6,:))
            nn(4) = 1._dp
            nn(1:3) = -pparton(1,1:3)/pparton(1,4)
            nn(:) = -dotvec(pparton(8,:)+pparton(7,:)+pparton(1,:),nn(:))*pparton(1,:)/pparton(1,4) / two
        else
            qsq = massvec(-pparton(1,:)-pparton(6,:))
            nn(4) = 1._dp
            nn(1:3) = -pparton(2,1:3)/pparton(2,4)
            nn(:) = -dotvec(pparton(8,:)+pparton(7,:)+pparton(2,:),nn(:))*pparton(2,:)/pparton(2,4) / two
        endif

        tau = 2._dp * dotvec(pparton(7,:) + pparton(8,:), nn) / (mt**2 - qsq)

!---    heck to make sure no NaN
        if (ieee_is_nan(tau)) then
          !call writeout(pparton)
          !write(6,*) __FILE__//' tau=',tau
          !stop
        endif

        if (bin .and. doMultitaucut .and. present(scetreweight_local)) then
            scetreweight_local(:) = 0._dp

            ! no sampled value for taus will pass cuts, if smaller than the smallest taucut
            if (tau < smallestTaucut*(tauc/taucut)) then
                scetreweight_local(:) = 0._dp
                return
            endif

            ! otherwise compute "weights" for other taucuts
            do j=1,size(tcutarray)
                if (tau < tcutarray(j)*(tauc/taucut)) then
                    scetreweight_local(j) = 0._dp
                else
                    scetreweight_local(j) = 1._dp
                endif
            enddo

            ! and postpone this cut for later in the *int routines
            ! and in the plotting
            if (tau < tauc) then
                !includeTaucutgrid(nd) = .false.
                return
            endif
        else
            if (tau < tauc) return
        endif

        passed_taucut_heavyprod = .true.

    end function passed_taucut_heavyprod

    subroutine singletop2_scet_virt_heavy_prod_all(p,msqv)
        use singletop2_nnlo_vars
    implicit none
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'ewcouple.f'
    include 'ckm.f'
    include 'nwz.f'
    include 'zprods_com.f'
    include 'scheme.f'

    real(dp), intent(in) :: p(mxpart,4)
    real(dp), intent(out) :: msqv(-nf:nf,-nf:nf, max_bcontrib,max_corr_on_beam)
    real(dp):: fac

    scheme = 'tH-V'

    call spinoru(6,p,za,zb)

    if (nwz == +1) then
        fac=aveqq*gw**8*xn**2*(as_heavy_beam2/2._dp/pi*cf)
        msqv([2,4],5, 1,2) = fac*virtqqb_heavy(1,2,3,4,5,6, renscale_beam2_isheavy_onheavy**2)
        msqv([-1,-3],5, 1,2) = fac*virtqqb_heavy(6,2,3,4,5,1, renscale_beam2_isheavy_onheavy**2)

        fac=aveqq*gw**8*xn**2*(as_heavy_beam1/2._dp/pi*cf)
        msqv(5,[2,4], 1,1) = fac*virtqqb_heavy(2,1,3,4,5,6, renscale_beam1_isheavy_onheavy**2)
        msqv(5,[-1,-3], 1,1) = fac*virtqqb_heavy(6,1,3,4,5,2, renscale_beam1_isheavy_onheavy**2)
    elseif (nwz == -1) then
      !ub=fac*virtqqb_heavy(6,2,4,3,5,1)
      !bu=fac*virtqqb_heavy(6,1,4,3,5,2)
      !bubar=fac*virtqqb_heavy(2,1,4,3,5,6)
      !ubarb=fac*virtqqb_heavy(1,2,4,3,5,6)
    endif

    end



    function virtqqb_heavy(ju,jb,jn,je,jc,jd,musq)
    implicit none
    real(dp):: virtqqb_heavy

    integer:: ju,jd,jn,je,jc,jb
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'masses.f'
    include 'sprods_com.f'
    include 'zprods_com.f'
    real(dp):: snec,prop,cv,cv0,mtsq
    complex(dp):: c1,amp,ampho

    real(dp), intent(in) :: musq


    mtsq=mt**2
    snec=+s(jn,je)+s(je,jc)+s(jc,jn)

    call coefs_new(s(ju,jd),mtsq,cv0,cv,c1,musq)

    if (s(ju,jd) < 0._dp) then
    prop=(s(ju,jd)-wmass**2)**2
    else
    prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
    endif
    prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
    prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

    ! only heavy line corrections for now
    cv0 = 0._dp

    amp=za(jc,jn)*zb(ju,jb) &
     *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
    ampho=za(jc,jn)*zb(ju,jb) &
     *(cplx1(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)) &
     +c1*chalf*zb(je,jb)*za(jb,jd))

    virtqqb_heavy=real(amp*conjg(ampho))/prop

    end


      subroutine coefs_new(s12,mtsq,cv0,cv,c1,musq)
      implicit none
      include 'types.f'
!-----In this routine:-
!-----cv0 is the results for all massless vertex function
!-----cv and c1 are  is the results for one-mass vertex function
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      real(dp):: cv,cv0,Li2la
      real(dp):: s12,mtsq,taucs,ddilog,eta,la,oml
      complex(dp):: lnrat,logoml,logla,xl12,logsca,Kfun,c1

      real(dp), intent(in) :: musq

      if (scheme =='dred') then
!------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
!------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

!**********************************************************************
!   Massless case
!   Taken from
!   %\cite{Altarelli:1979ub}
!   \bibitem{Altarelli:1979ub}
!   G.~Altarelli, R.~K.~Ellis and G.~Martinelli,
!   %``Large Perturbative Corrections To The Drell-Yan Process In QCD,''
!   Nucl.\ Phys.\ B {\bf 157}, 461 (1979).
!   %%CITATION = NUPHA,B157,461;%%
!   Using Eqn(58) with normalization changed to 
!   as/2/pi*cf*(4*pi)^ep/Gamma(1-ep) 
!   Taking account that Gamma(1-ep)^2/Gamma(1-2*ep)=1-ep^2*pi^2/6
!**********************************************************************
      xl12=lnrat(-s12,musq) 

!-----2/22/2012
!-----This appears to be the correction to the vertex in units of as/4/pi*cf
!-----despite the above comment
      cv0=-2._dp*(epinv2*epinv-epinv*real(xl12))-real(xl12**2) &
                 -3._dp*(epinv-real(xl12))-7._dp-eta

!---- this routine has been constructed following closely 
!---- the notation of
!---- %\cite{Gottschalk:1980rv}
!---- \bibitem{Gottschalk:1980rv}
!---- T.~Gottschalk,
!---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
!---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
!---- %%CITATION = PHRVA,D23,56;%%
!----- Adapted from Eqs.(A8,A9)

! NB  s12=-Q^2, -taucs=mtsq+Q^2
      taucs=s12-mtsq
      la=-s12/(mtsq-s12)
      oml=1._dp-la
!-----oml=mtsq/(mtsq-s12)
      logoml=-lnrat(-taucs,mtsq)
      logsca=lnrat(-taucs,musq)
      Kfun=cplx1(oml/la)*logoml

!--- Minus sign relative to Gottschalk since incoming b has momentum
!--- vector reversed for the t-channel process
!--- s-channel process follows by crossing
      c1=-ctwo*Kfun

      if (la < 1._dp) then
      Li2la=ddilog(la)
      else
      logla=lnrat(-s12,-taucs)
      Li2la=pisqo6-ddilog(oml)-real(logla*logoml)
      endif
!-----Again from A8 and A9 these are in units of alpha_s/4/pi*CF
      cv=-epinv2*epinv &
       -epinv*(2.5_dp+real(logoml-logsca))  &
       -0.5_dp*(11._dp+eta)-pisqo6+2._dp*Li2la-real(Kfun) &
        -0.5_dp*real(logoml*(cone-logoml)) &
        +2.5_dp*real(logsca)+real(logsca*logoml)-0.5_dp*real(logsca**2)

! answer from gotts.mac
!   ans:
!   -1/ep^2;
!   -2.5/ep -(log(oml)/ep-+log(sca)/ep)
!  -6-pisqo6+2*Li2-Ka
!  -0.5*log(oml)*(1-log(oml))
!  +2.5*log(sca)+log(oml)*log(sca)-log(sca)^2/2
      return

      end subroutine coefs_new

      ! chopped down from qqb_tbb_g
      subroutine qqb_tbb_g_heavy_all(p,msqall)
        use singletop2_nnlo_vars, only: maxbeams, beams_enabled
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'stopbmass.f'
      include 'taucut.f'! for usescet

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp):: ubtdg_h,fac

      call spinoru(7,p,za,zb)

      msqall = 0._dp
      
      if (nwz == +1) then
            if (any(beams_enabled(1:maxbeams) == 2)) then
                fac=2._dp*(4*pi*as_heavy_beam2)*cf*gw**8*xn**2
                msqall([2,4],5, 1,2) = aveqq*fac*(ubtdg_h(1,2,3,4,5,6,7,p))
                msqall([-1,-3],5, 1,2) = aveqq*fac*(ubtdg_h(6,2,3,4,5,1,7,p))
                msqall([2,4],0, 3,2) = aveqg*fac*(ubtdg_h(1,7,3,4,5,6,2,p))
                msqall([-1,-3],0, 3,2) = aveqg*fac*(ubtdg_h(6,7,3,4,5,1,2,p))
            endif

            if (any(beams_enabled(1:maxbeams) == 1)) then
                fac=2._dp*(4*pi*as_heavy_beam1)*cf*gw**8*xn**2
                msqall(5,[2,4], 1,1) = aveqq*fac*(ubtdg_h(2,1,3,4,5,6,7,p))
                msqall(5,[-1,-3], 1,1) = aveqq*fac*(ubtdg_h(6,1,3,4,5,2,7,p))
                msqall(0,[2,4], 3,1) = aveqg*fac*(ubtdg_h(2,7,3,4,5,6,1,p))
                msqall(0,[-1,-3], 3,1) = aveqg*fac*(ubtdg_h(6,7,3,4,5,2,1,p))
            endif
      elseif (nwz == -1) then
          error stop "nwz = -1 not implemented in qqb_tbb_g_heavy"
      endif
      
      end subroutine

end module
