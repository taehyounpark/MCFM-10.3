!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop2_scet_heavy_decay
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use singletop2_scet_light, only: singletop2_scet_tree, &
          singletop2_scet_tree_ub, &
          singletop2_scet_tree_bu
      use singletop2_nnlo_vars
    use LHAPDF
    implicit none
    private

    public :: passed_taucut_decay

    public :: lumxmsq_singletop_decay_jetmass
    public :: singletop2_heavy_decay_g_all
    public :: singletop2_heavy_decay_gs_all
    public :: singletop2_heavy_decay_gs_all_new
    public :: singletop2_scet_virt_heavy_decay_all

    public :: dips_fi_mt, dips_fi_mt_gg
    public :: dips_ff_ident
    public :: intdip_fi_mt, intdip_fi_mt_gg

    public :: softfun, jetfun, hardfun

    public :: qqbtbbargd

    contains

      subroutine singletop2_scet_virt_heavy_decay_all(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'scheme.f'
      real(dp):: msqv(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam),p(mxpart,4), &
       fac,bu,ub,bubar,ubarb

      scheme='tH-V'

      call spinoru(6,p,za,zb)

      msqv = 0._dp

      if (nwz == +1) then
            fac=(as_heavy_beam2/2._dp/pi*cf)*aveqq*gw**8*xn**2
            ub=fac*virtqqbdk(p, 1,2,3,4,5,6,renscale_beam2_isheavy_onheavy**2)
            ubarb=fac*virtqqbdk(p, 6,2,3,4,5,1,renscale_beam2_isheavy_onheavy**2)
            msqv([2,4],5, 1,2) = ub
            msqv([-1,-3],5, 1,2) = ubarb

            fac=(as_heavy_beam1/2._dp/pi*cf)*aveqq*gw**8*xn**2
            bubar=fac*virtqqbdk(p, 6,1,3,4,5,2,renscale_beam1_isheavy_onheavy**2)
            bu=fac*virtqqbdk(p, 2,1,3,4,5,6,renscale_beam1_isheavy_onheavy**2)
            msqv(5,[2,4], 1,1) = bu
            msqv(5,[-1,-3], 1,1) = bubar
      elseif (nwz == -1) then
        !ub=fac*virtqqbdk(6,2,4,3,5,1)
        !bu=fac*virtqqbdk(6,1,4,3,5,2)
        !bubar=fac*virtqqbdk(2,1,4,3,5,6)
        !ubarb=fac*virtqqbdk(1,2,4,3,5,6)
      endif

      end


      function virtqqbdk(p, ju,jb,jn,je,jc,jd,musq)
      implicit none
      include 'types.f'
      real(dp):: virtqqbdk
      

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      real(dp):: snec,prop,mtsq,cv,ct,c1
      complex(dp):: amp,ampho

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: ju,jd,jn,je,jc,jb
      real(dp), intent(in) :: musq

      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      call coefsdk(s(jn,je),mtsq,ct,cv,c1,musq)
      if (s(ju,jd) < 0._dp) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      amp=za(jc,jn)*zb(ju,jb) &
       *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

      ! disable integrated counterterm ct for scet calculation
      !ct = 0._dp
      !ct = intdip_fi_mt(p)

      ! coefficients are in terms of alphas/4/pi
      ampho=za(jc,jn)*zb(ju,jb) &
       *(cplx1(ct+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)) &
       +cplx1(0.5_dp*c1)*zb(je,jc)*za(jc,jd))

      ! but factor of two is left out here
      virtqqbdk=real(amp*conjg(ampho))/prop
      return
      end

 
      subroutine coefsdk(s12,mtsq,ct,cv,c1,musq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'alfacut.f'
      !include 'includect.f'
      real(dp):: cv,ct,c1,s12,mtsq,ddilog,rsq,omrsq,eta,Kfun, &
       wlog,rlog,mulog

      real(dp), intent(in) :: musq

      if (scheme =='dred') then
         eta=0._dp
      elseif (scheme == 'tH-V') then
         eta=1._dp
      endif

      rsq=s12/mtsq
      omrsq=1._dp-rsq
      wlog=log(omrsq)
      rlog=log(rsq)

!--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
!----  ct is the integrated counter-term, including alpha-dependence,
!----  see Eq. (10) of arXiv:1102.1967
      mulog=log(musq/mtsq)
      !if (includect) then
      ct=(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2) &
       +(epinv+mulog)*(2.5_dp-2._dp*wlog) & 
       +25._dp/4._dp+0.5_dp*(1._dp/omrsq**2-8._dp/omrsq+7._dp)*rlog &
       +0.5_dp/omrsq+2._dp*ddilog(omrsq)-5._dp*pisqo6 &
       -5._dp*wlog+2._dp*wlog**2+eta/2._dp &
       -2._dp*log(aif)**2-(3.5_dp-4._dp*aif+aif**2/2._dp)*log(aif) &
       +2._dp*(1._dp-aif)*rsq/omrsq*log(rsq/(1._dp-aif+rsq*aif))
      !else
      !ct=0._dp
      !endif

!---- this routine has been constructed from 
!---- %\cite{Gottschalk:1980rv}
!---- \bibitem{Gottschalk:1980rv}
!---- T.~Gottschalk,
!---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
!---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
!---- %%CITATION = PHRVA,D23,56;%%
!----- Adapted from Eqs.(A8,A9)

      Kfun=1._dp/rsq*wlog
      c1=2._dp*Kfun      
      cv=-(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2) &
       -(epinv+mulog)*(2.5_dp-2._dp*wlog) &
       -0.5_dp*(11._dp+eta)-pisqo6-2._dp*ddilog(rsq) &
        +3._dp*wlog-2._dp*wlog**2-Kfun
      
      return
      end

    subroutine lumxmsq_singletop_decay_jetmass(p,xx,z1,z2,order,xmsq,central)
        use types
        use SCET
        use SCET_Jet
        implicit none
        include 'constants.f'
        include 'masses.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'taucut.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'beamtype.f'
    
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in) :: xx(2), z1, z2
        integer, intent(in) :: order
        real(dp), intent(out) :: xmsq
        logical, intent(in) :: central

        real(dp) :: beama0(-5:5), beamb0(-5:5)

        real(dp) :: jet(2,0:4), soft(2,0:4)
        real(dp) :: hard_ub(0:2), hard_bu(0:2), hard_ubarb(0:2), hard_bubar(0:2)

        integer :: m
        real(dp) :: x

        real(dp), allocatable :: taucuts(:)
        real(dp), allocatable :: xmsqs(:)

        real(dp) :: puremass

        real(dp) :: tauc

        ! provides access to renormalization and factorization scales
        ! on heavy and light lines separately
        ! now called previously.
        !call singletop2_scale_setup(p)

        xmsq = 0._dp

        x = puremass(p(3,:)+p(4,:))**2 / mt**2

        !tauc = getdynamictau(p, taucut)
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

        soft = softfun(order)
        jet = jetfun(order)

        xmsqs(:) = 0._dp

        ! Numerically confirmed that at 1-loop there is no scale dependence left

        ! There is a lot of duplicate calculation here,
        ! but since this is just the below cut, I don't care
        hard_ub = real(hardfun(order, renscale_beam2_isheavy_onheavy, x, p, 1,2,3,4,5,6),dp)
        hard_ubarb = real(hardfun(order, renscale_beam2_isheavy_onheavy, x, p, 6,2,3,4,5,1),dp)

        call fdist(ih1,xx(1),facscale_beam1_islight_onheavy,beama0,1)
        call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy,beamb0,2)

        xmsqs = xmsqs + hard_ub(0)*assemble_decay_pieces(order, x, &
            renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
            hard_ub/hard_ub(0), jet, soft)*(beama0(2) + beama0(4))*beamb0(5)

        xmsqs = xmsqs + hard_ubarb(0)*assemble_decay_pieces(order, x, &
            renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
            hard_ubarb/hard_ubarb(0), jet, soft)*(beama0(-1) + beama0(-3))*beamb0(5)


        hard_bu = real(hardfun(order, renscale_beam1_isheavy_onheavy, x, p, 2,1,3,4,5,6),dp)
        hard_bubar = real(hardfun(order, renscale_beam1_isheavy_onheavy, x, p, 6,1,3,4,5,2),dp)

        call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy,beama0,1)
        call fdist(ih2,xx(2),facscale_beam2_islight_onheavy,beamb0,2)

!       !!! DEBUG
!       debug: block
!       use NNTopDec
!       real(dp) :: new0loop, new1loop, new2loop

!       write (*,*) "=== LUMXMSQ CALLED ==="
!       call initloop
!       write (*,*) "scale = ", renscale_beam2_isheavy_onheavy
!       write (*,*) "g11, g21+g31", 2*g11, (2*g21 + 2*g31)/mt
!       write (*,*) "g12, g22+g32", 2*g12, (2*g22 + 2*g32)/mt
!       write (*,*) "g11sq", g11**2/ g11sq
!       write (*,*) "g21sq", g21**2/ g21sq
!       write (*,*) "g31sq", g31**2/ g31sq
!       write (*,*) "g11g21", g11g21/ (g11*g21)
!       write (*,*) "g21g31", g21g31/ (g21*g31)

!       call top22treesqamp(p,new0loop)
!       write (*,*) "Born", new0loop/hard_ub(0)
!       !! to run with snloV, set order = 2 in nntopdec.f
!       !call top22loopsqamp(p,as_heavy_beam2,1,new1loop)
!       !write (*,*) "1-Loop", new1loop/(hard_ub(0)*assemble_decay_pieces(order, x, &
!       !    renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
!       !    hard_ub/hard_ub(0), jet, soft))

!       !! to run with nnloVV, set order = 3 in nntopdec.f
!       call top22loopsqamp(p,as_heavy_beam2,2,new1loop)
!       write (*,*) "2-Loop", new1loop/(hard_ub(0)*assemble_decay_pieces(order, x, &
!           renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucuts, &
!           hard_ub/hard_ub(0), jet, soft))
!       pause
!       end block debug
!       !!! END DEBUG

        xmsqs = xmsqs + hard_bu(0)*assemble_decay_pieces(order, x, &
            renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucuts, &
            hard_bu/hard_bu(0), jet, soft)*beama0(5)*(beamb0(2) + beamb0(4))

        xmsqs = xmsqs + hard_bubar(0)*assemble_decay_pieces(order, x, &
            renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucuts, &
            hard_bubar/hard_bubar(0), jet, soft)*beama0(5)*(beamb0(-1) + beamb0(-3))

        xmsq = xmsqs(1)

        if (central .and. doMultitaucut) then
            scetreweight(:) = 0._dp
            if (xmsq /= 0._dp) then
                scetreweight(:) = xmsqs(2:) / xmsq
            endif
        endif

    end subroutine

    function jetfun(order)
        use types
        use constants
        implicit none
        include 'masses.f'
        include 'nf.f'

        integer, intent(in) :: order
        real(dp) :: jetfun(2,0:4)

        real(dp), parameter :: beta0 = 11/3._dp*CA - 4/3._dp*tf*nf

        jetfun(:,:) = 0._dp

        ! Becher, Neubert, hep-ph/0603140

        if (order >= 1) then
            jetfun(1,2) = 2*CF
            jetfun(1,1) = -3*CF
            jetfun(1,0) = 7*CF - CF*Pi**2
        endif

        if (order >= 2) then
            jetfun(2,4) = 2*CF**2
            jetfun(2,3) = (-2*beta0*CF)/3._dp - 6*CF**2
            jetfun(2,2) = (3*beta0*CF)/2._dp + (134*CA*CF)/9._dp + (37*CF**2)/2._dp - (2*CA*CF*Pi**2)/3._dp - &
                (10*CF**2*Pi**2)/3._dp - (40*CF*NF*TF)/9._dp
            jetfun(2,1) = -7*beta0*CF - (1769*CA*CF)/54._dp - (45*CF**2)/2._dp + beta0*CF*Pi**2 - &
                (11*CA*CF*Pi**2)/9._dp + 7*CF**2*Pi**2 + (242*CF*NF*TF)/27._dp + &
                (4*CF*NF*Pi**2*TF)/9._dp + 40*CA*CF*zeta3 - 8*CF**2*zeta3
            jetfun(2,0) = (53129*CA*CF)/648._dp + (205*CF**2)/8._dp - (208*CA*CF*Pi**2)/27._dp - &
                (67*CF**2*Pi**2)/6._dp - (17*CA*CF*Pi**4)/180._dp + (14*CF**2*Pi**4)/15._dp - &
                (4057*CF*NF*TF)/162._dp + (68*CF*NF*Pi**2*TF)/27._dp - (206*CA*CF*zeta3)/9._dp - &
                18*CF**2*zeta3 + (16*CF*NF*TF*zeta3)/9._dp
        endif

    end function

    ! checked again 2020-02-18
    ! Numerically agrees with HardFunctionNew.nb
    ! for complex tree and struc and arguments 0 < x < 1 for decay
    ! and x < 0 for production.
    ! tree and struc expressions tested through agreement at NLO.
    function hardfun(order, scale, x, p, ju,jb,jn,je,jc,jd, production)
        use types
        use constants
        implicit none
        include 'masses.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ewcouple.f'

        integer, intent(in) :: order
        real(dp), intent(in) :: scale
        real(dp), intent(inout) :: x
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ju,jd,jn,je,jc,jb
        logical, intent(in), optional :: production
        complex(dp) :: hardfun(0:2)

        real(dp) :: Li2, Li3
        real(dp) :: LSX, LJX

        real(dp) :: snec, prop, fac
        complex(dp) :: tree, struc, treeCC, strucCC
        real(dp) :: LogMtMu

        complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp) :: s(mxpart,mxpart)

        complex(dp) :: w(127)
        real(dp) :: omx

        real(dp), parameter :: c0mb = -60.2493267_dp

        ! n1,n2 range of indices (0,1) or (-1,0) or (-1,1)
        ! nw, maximum weight
        integer, parameter :: n1=-1, n2=1, nw = 4

        complex(dp):: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
                  Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
        real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), &
                  Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
        real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), &
                  Hi4(n1:n2,n1:n2,n1:n2,n1:n2)

        call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

        call spinoru_dp_s(6,p,za,zb,s)

        fac = (1._dp/(4._dp*9._dp))*gw**8*ca**2

        snec = +s(jn,je)+s(je,jc)+s(jc,jn)
        if (s(ju,jd) < 0._dp) then
            prop=(s(ju,jd)-wmass**2)**2
        else
            prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
        endif
            prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
            prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

        tree = za(jc,jn)*zb(ju,jb) *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

        treeCC = conjg(tree)

        LogMtMu = log(mt/scale)
        omx = 1._dp - x

        if (present(production)) then
            if (production .eqv. .true.) then
                struc = -za(jc,jn)*zb(ju,jb)*zb(je,jb)*za(jb,jd) * 0.5_dp * mt
                LJX = 2*LogMtMu + log(omx)
                LSX = LogMtMu
            else
                ! decay
                struc = za(jc,jn)*zb(ju,jb)*zb(je,jc)*za(jc,jd) * 0.5_dp * mt
                LJX = 2*LogMtMu
                LSX = LogMtMu - log(omx)
            endif
        else
            ! decay
            struc = za(jc,jn)*zb(ju,jb)*zb(je,jc)*za(jc,jd) * 0.5_dp * mt
            LJX = 2*LogMtMu
            LSX = LogMtMu - log(omx)
        endif

        strucCC = conjg(struc)


        ! This assembly assumes for simplification that the form factors are
        ! always real which holds for mW^/mt^2 = x < 1, so as long as the W in decay
        ! doesn't go off-shell too far, this is ok.

        w = 0._dp

!       !!! DEBUG
!       tree = 1._dp + 2*im
!       treeCC = 1._dp - 2*im

!       struc = 3._dp + 4*im
!       strucCC = 3._dp - 4*im

!       mt = 172.5_dp
!       LogMtMu = log(172.5_dp/400._dp)
!       x = -200._dp

!       omx = 1._dp - x
!       LJX = 2*LogMtMu
!       LSX = LogMtMu - log(omx)

!       call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

!       !!! END DEBUG

        if (x > 1._dp) then
            write (*,*) "WARNING: x evaluated > 1, x = ", x
        endif

        !!! DEBUG for nntopdec comparison
        !struc = 1._dp
        !strucCC = 1._dp
        !tree = 1._dp
        !treeCC = 1._dp
        !!! END DEBUG for nntopdec comparison

        if (order < 2) then
#include "hardfun_nlo.f90.inc"
        endif

        if (order >= 2) then
#include "hardfun_nnlo.f90.inc"
        endif

!       write (*,*) "nf = ", nf
!       write (*,*) "hardfun(0)", hardfun(0)
!       write (*,*) "hardfun(1)", hardfun(1)
!       write (*,*) "hardfun(2)", hardfun(2)
        !write (*,*) CF*LogMtMu*hardfun(1) / hardfun(2)

        !hardfun(2) = hardfun(2) + CF*LogMtMu*hardfun(1)


!       pause


!       !!! DEBUG
!       debug: block
!       use NNTopDec
!       real(dp) :: new0loop, new1loop, new2loop

!       write (*,*) "=== HARDFUN CALLED ==="
!       call initloop
!       ! for comparison set struc, strucCC, tree, treeCC to 1
!       ! above before the evaluation of hardfun
!       write (*,*) "g11, g21+g31", 2*g11, (2*g21 + 2*g31)/mt
!       write (*,*) "g12, g22+g32", 2*g12, (2*g22 + 2*g32)/mt
!       write (*,*) "g11sq", g11**2/ g11sq
!       write (*,*) "g21sq", g21**2/ g21sq
!       write (*,*) "g31sq", g31**2/ g31sq
!       write (*,*) "g11g21", g11g21/ (g11*g21)
!       write (*,*) "g21g31", g21g31/ (g21*g31)

!       ! these are only valid for Born amplitudes set to one
!       write (*,*) "h(1): ", ((2*g21 + 2*g31)/mt + 2*g11) / hardfun(1)
!       write (*,*) "h(2): ", ((g21sq + 2*g21g31 + g31sq)/mt**2 + &
!           (2*g11g21 + 2*g22 + 2*g11*g31 + 2*g32)/mt + &
!           (g11sq + 2*g12)) / hardfun(2)
!       pause
!       end block debug
        !!! END DEBUG

        ! factors of as/4/pi are attached in assembly routine
        hardfun(:) = fac * hardfun(:) / prop

    end function

    ! checked again 2020-02-18
    function softfun(order)
        use types
        use constants
        implicit none
        include 'masses.f'
        include 'nf.f'

        integer, intent(in) :: order
        real(dp) :: softfun(2,0:4)

        real(dp), parameter :: beta0 = 11/3._dp*CA - 4/3._dp*tf*nf

        softfun(:,:) = 0._dp

        ! Becher, Neubert, hep-ph/0512208

        if (order >= 1) then
            softfun(1,2) = -4*CF
            softfun(1,1) = -4*CF
            softfun(1,0) = -(CF*Pi**2)/6._dp
        endif

        if (order >= 2) then
            softfun(2,4) = 8*CF**2
            softfun(2,3) = (8*beta0*CF)/3._dp + 16*CF**2
            softfun(2,2) = 4*beta0*CF - (268*CA*CF)/9._dp + 8*CF**2 + (4*CA*CF*Pi**2)/3._dp - &
                (14*CF**2*Pi**2)/3._dp + (80*CF*NF*TF)/9._dp
            softfun(2,1) = (220*CA*CF)/27._dp + (beta0*CF*Pi**2)/3._dp + (CA*CF*Pi**2)/9._dp - &
                (14*CF**2*Pi**2)/3._dp + (16*CF*NF*TF)/27._dp + (4*CF*NF*Pi**2*TF)/9._dp - &
                36*CA*CF*zeta3 + 64*CF**2*zeta3
            softfun(2,0) = (-326*CA*CF)/81._dp - (427*CA*CF*Pi**2)/108._dp - (4*CF**2*Pi**2)/3._dp + &
                (67*CA*CF*Pi**4)/180._dp - (3*CF**2*Pi**4)/40._dp - (8*CF*NF*TF)/81._dp + &
                (5*CF*NF*Pi**2*TF)/27._dp - (107*CA*CF*zeta3)/9._dp + 32*CF**2*zeta3 - &
                (20*CF*NF*TF*zeta3)/9._dp
        endif

    end function

    function assemble_decay_pieces(order, x, scale, as, taucut, h, j, s)
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
        real(dp) :: assemble_decay_pieces(size(taucut(:)))

        real(dp) :: LogMtMu
        !complex(dp) :: treeCC, strucCC
        real(dp) :: full(0:2,0:4)
        real(dp) :: LSX, LJX

        !treeCC = conjg(tree)
        !strucCC = conjg(struc)
        LogMtMu = log(mt/scale)

        ! expressing everything in terms of LJX/LSX + log(tau) allows
        ! us to repurpose this assembly function for the production process
        LJX = 2*LogMtMu
        LSX = LogMtMu - log(1-x)

        full(:,:) = 0._dp

        if (h(0) /= 1._dp) then
            write (*,*) "WARNING: bad hard function normalization!"
        endif

        if (coeffonly) then
            full(0,0) = 0._dp
        else
            full(0,0) = 1._dp
        endif

        if (order == 1 .or. (order > 1 .and. coeffonly .eqv. .false.)) then
            full(1,2) = j(1,2) + s(1,2)
            full(1,1) = j(1,1) + 2*LJX*j(1,2) + s(1,1) + 2*LSX*s(1,2)
            full(1,0) = h(1) + j(1,0) + LJX*j(1,1) + LJX**2*j(1,2) + s(1,0) + LSX*s(1,1) + &
               LSX**2*s(1,2)
        endif

        if (order >= 2) then
            full(2,4) = j(2,4) + j(1,2)*s(1,2) + s(2,4)
            full(2,3) = j(2,3) + 4*LJX*j(2,4) + j(1,2)*s(1,1) + j(1,1)*s(1,2) + &
                2*LJX*j(1,2)*s(1,2) + 2*LSX*j(1,2)*s(1,2) + s(2,3) + 4*LSX*s(2,4)
            full(2,2) = j(2,2) + 3*LJX*j(2,3) + 6*LJX**2*j(2,4) + j(1,2)*s(1,0) + &
                j(1,1)*s(1,1) + 2*LJX*j(1,2)*s(1,1) + LSX*j(1,2)*s(1,1) + &
                j(1,0)*s(1,2) + LJX*j(1,1)*s(1,2) + 2*LSX*j(1,1)*s(1,2) + &
                LJX**2*j(1,2)*s(1,2) + 4*LJX*LSX*j(1,2)*s(1,2) + LSX**2*j(1,2)*s(1,2) - &
                (2*Pi**2*j(1,2)*s(1,2))/3._dp + h(1)*(j(1,2) + s(1,2)) + s(2,2) + &
                3*LSX*s(2,3) + 6*LSX**2*s(2,4)
            full(2,1) = j(2,1) + 2*LJX*j(2,2) + 3*LJX**2*j(2,3) + 4*LJX**3*j(2,4) + &
                j(1,1)*s(1,0) + 2*LJX*j(1,2)*s(1,0) + j(1,0)*s(1,1) + &
                LJX*j(1,1)*s(1,1) + LSX*j(1,1)*s(1,1) + LJX**2*j(1,2)*s(1,1) + &
                2*LJX*LSX*j(1,2)*s(1,1) - (Pi**2*j(1,2)*s(1,1))/3._dp + &
                2*LSX*j(1,0)*s(1,2) + 2*LJX*LSX*j(1,1)*s(1,2) + LSX**2*j(1,1)*s(1,2) - &
                (Pi**2*j(1,1)*s(1,2))/3._dp + 2*LJX**2*LSX*j(1,2)*s(1,2) + &
                2*LJX*LSX**2*j(1,2)*s(1,2) - (2*LJX*Pi**2*j(1,2)*s(1,2))/3._dp - &
                (2*LSX*Pi**2*j(1,2)*s(1,2))/3._dp + 8*zeta3*j(1,2)*s(1,2) + &
                h(1)*(j(1,1) + 2*LJX*j(1,2) + s(1,1) + 2*LSX*s(1,2)) + s(2,1) + &
                2*LSX*s(2,2) + 3*LSX**2*s(2,3) + 4*LSX**3*s(2,4)
            full(2,0) = h(2) + j(2,0) + LJX*j(2,1) + LJX**2*j(2,2) + LJX**3*j(2,3) + &
                LJX**4*j(2,4) + j(1,0)*s(1,0) + LJX*j(1,1)*s(1,0) + &
                LJX**2*j(1,2)*s(1,0) + LSX*j(1,0)*s(1,1) + LJX*LSX*j(1,1)*s(1,1) - &
                (Pi**2*j(1,1)*s(1,1))/6._dp + LJX**2*LSX*j(1,2)*s(1,1) - &
                (LJX*Pi**2*j(1,2)*s(1,1))/3._dp + 2*zeta3*j(1,2)*s(1,1) + &
                LSX**2*j(1,0)*s(1,2) + LJX*LSX**2*j(1,1)*s(1,2) - &
                (LSX*Pi**2*j(1,1)*s(1,2))/3._dp + 2*zeta3*j(1,1)*s(1,2) + &
                LJX**2*LSX**2*j(1,2)*s(1,2) - (2*LJX*LSX*Pi**2*j(1,2)*s(1,2))/3._dp - &
                (Pi**4*j(1,2)*s(1,2))/90._dp + 4*LJX*zeta3*j(1,2)*s(1,2) + &
                4*LSX*zeta3*j(1,2)*s(1,2) + &
                h(1)*(j(1,0) + LJX*j(1,1) + LJX**2*j(1,2) + s(1,0) + LSX*s(1,1) + &
                   LSX**2*s(1,2)) + s(2,0) + LSX*s(2,1) + LSX**2*s(2,2) + &
                LSX**3*s(2,3) + LSX**4*s(2,4)

        endif

        assemble_decay_pieces = 0._dp

        if (order == 1 .or. (order > 1 .and. coeffonly .eqv. .false.)) then
            assemble_decay_pieces = full(0,0) + &
                as/4._dp/pi*( &
                    full(1,0) + full(1,1)*log(taucut) + full(1,2)*log(taucut)**2 )
        endif

        if (order >= 2) then
            assemble_decay_pieces = assemble_decay_pieces + (as/4._dp/pi)**2 * ( &
                full(2,0) + full(2,1)*log(taucut) + full(2,2)*log(taucut)**2 &
                + full(2,3)*log(taucut)**3 + full(2,4)*log(taucut)**4 )
        endif

    end function


    ! this doesn't decompose into jet, and soft functions
    ! only for debugging of one-loop expression
    function assemble_beam2_decay_full(order, taucut, beama0, beamb0, hard)
        use types 
        use constants
        implicit none
        include 'kpart.f'
        include 'masses.f'

        real(dp) :: assemble_beam2_decay_full
        integer, intent(in) :: order
        real(dp), intent(in) :: taucut
        real(dp), intent(in) :: beama0, beamb0
        real(dp), intent(in) :: hard(2)

        real(dp) :: lcut, full(-1:3)
        real(dp) :: scale, x

        lcut = log(taucut)
        scale = renscale_beam2_isheavy_onheavy
        x = wmass**2/mt**2

        error stop "assemble_beam2_decay_full is deprecated."

        if (coeffonly) then
            assemble_beam2_decay_full = 0._dp
        else
            assemble_beam2_decay_full = 1._dp
        endif

        full(:) = 0._dp

        full(-1) = (4*log(scale/mt)**2 -8*log(1-x)*log(scale/mt) + &
            10*log(scale/mt) - 4*log(1-x)**2 + 4*log(1-x) + 7 &
            - 7*pi**2/6._dp)*cf
        full(0) = +(8*log(1-x) - 7)*cf
        full(1) = -4*cf

        ! the additional pi^2 term is from the 1/Gamma(1-eps) expansion

        assemble_beam2_decay_full = assemble_beam2_decay_full + as_heavy_beam2/4._dp/pi*( &
            full(-1) + full(0)*lcut + full(1)*lcut**2/2._dp) + hard(1) &
                + pi**2/6._dp*as_heavy_beam2/4._dp/pi*cf
        assemble_beam2_decay_full = assemble_beam2_decay_full * beama0 * beamb0

    end function

      function passed_taucut_decay(pparton,scetreweight_local,fvec)
          use SCET
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'kprocess.f'
      include 'masses.f'

      logical :: passed_taucut_decay
      real(dp), intent(in) :: pparton(mxpart,4)
      real(dp), intent(inout) :: scetreweight_local(:)
      real(dp), intent(in), optional :: fvec(4)

      integer j
      real(dp) :: tau,tauc,massvec

      logical :: bin
      common/bin/bin

      scetreweight_local(:) = 0._dp

      passed_taucut_decay = .false.

      !tauc = getdynamictau(pparton,taucut)
      tauc = taucut

      if (present(fvec)) then
          tau = massvec(fvec) / mt**2
      else
          tau = massvec(pparton(5,:) + pparton(7,:) + pparton(8,:)) / mt**2
      endif

      if ((tau < 0._dp)) then
          if (tau > -1d-8) then
              tau = 0._dp
          else
              write (*,*) "WARNING, abnormal tau = ", tau
              tau = 0._dp
          endif
      endif

!--- check to make sure no NaN
!     if (ieee_is_nan(tau)) then
!       write(*,*) __FILE__//"maketaucut.f:  tau=",tau
!       write (*,*) "p5 = ", pparton(5,:)
!       write (*,*) "p7 = ", pparton(7,:)
!       write (*,*) "p8 = ", pparton(8,:)
!       write (*,*) "(p5+p7+p8)^2 = ", &
!           massvec(pparton(5,:)+pparton(7,:)+pparton(8,:))
!       write (*,*) "----"
!     endif

      if (bin .and. doMultitaucut) then
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

      passed_taucut_decay = .true.

      return
      end function passed_taucut_decay

    subroutine singletop2_heavy_decay_g_all(p,msq)
        use singletop2_nnlo_vars
    implicit none
    include 'types.f'
          
!     Matrix element for t-bbar production
!      b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
!     averaged(summed) over initial(final) colours and spins
!--NB average over spins only -- colour factors cancel
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'ewcouple.f'
    include 'ckm.f'
    include 'nwz.f'
    include 'zprods_com.f'
    real(dp):: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam),p(mxpart,4)
    real(dp):: fac,ub,bu,bubar,ubarb
    integer:: ib

    call spinoru(7,p,za,zb)
    ib=5*nwz

    msq(:,:,:,:) = 0._dp

    if (nwz == +1) then
            fac=2._dp*(4*pi*as_heavy_beam2)*cf*aveqq*gw**8*xn**2
            ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
            ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
            msq([2,4],5, 1,2) = ub
            msq([-1,-3],5, 1,2) = ubarb

            fac=2._dp*(4*pi*as_heavy_beam1)*cf*aveqq*gw**8*xn**2
            bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
            bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
            msq(5,[2,4], 1,1) = bu
            msq(5,[-1,-3], 1,1) = bubar
    elseif (nwz == -1) then
        error stop "to do"
        ub=fac*qqbtbbargd(6,2,4,3,5,1,7,p)
        bu=fac*qqbtbbargd(6,1,4,3,5,2,7,p)
        ubarb=fac*qqbtbbargd(1,2,4,3,5,6,7,p)
        bubar=fac*qqbtbbargd(2,1,4,3,5,6,7,p)
    endif

    end subroutine singletop2_heavy_decay_g_all


    function qqbtbbargd(ju,jb,jn,je,jc,jd,jg,p)
    implicit none
    include 'types.f'
    real(dp):: qqbtbbargd
!     Matrix element squared for single top production with gluon
!     radiation in decay (radiation from final line)
!      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc)+g(jg))+d(jd)
!     masses of b quarks c.c=b.b=0
           
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'masses.f'
    include 'zprods_com.f'
    include 'sprods_com.f'
    integer:: ju,jb,jn,je,jc,jd,jg,nu
    real(dp):: p(mxpart,4),pt(4),ptDpt
    real(dp):: sne,sdu,prop,twoptg
    complex(dp):: ampf(2),amp0

    do nu=1,4
        pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)+p(jg,nu)
    enddo
    ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
    twoptg= &
    &  2._dp*(pt(4)*p(jg,4)-pt(1)*p(jg,1)-pt(2)*p(jg,2)-pt(3)*p(jg,3))

    sne=s(jn,je)
    sdu=s(jd,ju)
    if (sdu < 0._dp) then
        prop=(sdu-wmass**2)**2
    else
        prop=(sdu-wmass**2)**2+(wmass*wwidth)**2
    endif
    prop=((sne-wmass**2)**2+(wmass*wwidth)**2) &
    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop
          
!  -Lefthanded gluon tb-line
!      ampf(1)=(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))
!     & *(+zb(jc,je)*za(je,jd)+zb(jc,jn)*za(jn,jd)+zb(jc,jg)*za(jg,jd))
!     & +mt**2*zb(je,jc)*za(jg,jd)
!      ampf(1)=-za(jc,jn)*zb(ju,jb)/(twoptg*zb(jg,jc))*ampf(1)
!      ampf(1)=ampf(1)-za(jg,jn)*zb(ju,jb)/zb(jg,jc)
!     & *(+zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))
!      write(6,*) 'ampf(1)',ampf(1)


    amp0=(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd)) &
    *za(jc,jn)*zb(ju,jb)
!---eikonal form
    ampf(1)=-amp0 &
    *(za(jg,je)*zb(je,jc)+za(jg,jn)*zb(jn,jc))/zb(jg,jc)/twoptg
    ampf(1)=ampf(1)-za(jg,jn)*zb(ju,jb)/zb(jg,jc) &
    *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))
     
!  -Righthanded gluon tb-line
!      ampf(2)=-(zb(je,jn)*za(jn,jc)
!     & *(+zb(jg,jc)*za(jc,jd)+zb(jg,je)*za(je,jd)+zb(jg,jn)*za(jn,jd))
!     & +mt**2*zb(je,jg)*za(jc,jd))
!      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jc,jg)/twoptg*ampf(2)
!      write(6,*) 'ampf(2)',ampf(2)

!---eikonal form
    ampf(2)=-amp0 &
    *(zb(jg,je)*za(je,jc)+zb(jg,jn)*za(jn,jc))/za(jc,jg)/twoptg &
    -za(jc,jn)*zb(ju,jb)/twoptg*zb(je,jg) &
    *(zb(jg,jc)*za(jc,jd)+zb(jg,je)*za(je,jd)+zb(jg,jn)*za(jn,jd))


    qqbtbbargd=(abs(ampf(1))**2+abs(ampf(2))**2)/prop

    return
    end function qqbtbbargd

    subroutine singletop2_heavy_decay_gs_all_new(p,ndmx,msqc)
        use singletop2_nnlo_vars, only: corr_on_beam
        use singletop2_scale_m
        use types
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'qqgg.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ndmx
        real(dp), intent(out) :: msqc(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

        real(dp) :: msqborn(-nf:nf,-nf:nf)

        integer, parameter :: noglue(4) = [-3,-1,2,4]

        external donothing_gvec

        msqc = 0._dp

        qqproc = .true.

        ! shouldn't be necessary for decay, but let's have them anyway:
        corr_islight = .false.
        corr_beam1 = .false.
        corr_on_beam = 2

        call dips_fi_mt(1, p, singletop2_scet_tree_ub, msqborn, 5,7)
        msqc(1,noglue,5, 1,2) = cf*msqborn(noglue,5)

        corr_islight = .false.
        corr_beam1 = .true.
        corr_on_beam = 1

        call dips_fi_mt(2, p, singletop2_scet_tree_bu, msqborn, 5,7)
        msqc(2,5,noglue, 1,1) = cf*msqborn(5,noglue)

    end subroutine singletop2_heavy_decay_gs_all_new

    ! Campbell, Ellis, Tramontano, PR D70 094012, p.7
    subroutine ltrans(pt,pw,vec,pout)
        use types
        implicit none

        real(dp), intent(in) :: pt(4), pw(4), vec(4)
        real(dp), intent(out) :: pout(4)

        real(dp) :: dotvec,massvec
        
        real(dp) :: sinhx, coshx

        sinhx = (-(massvec(pt) - massvec(pw))*dotvec(pt,pw) + (massvec(pt) + massvec(pw)) * &
                    sqrt(dotvec(pt,pw)**2 - massvec(pw)*massvec(pt)) ) / 2._dp / massvec(pt) / massvec(pw)

        coshx = ((massvec(pt)+massvec(pw))*dotvec(pt,pw) - (massvec(pt)-massvec(pw))* &
                    sqrt(dotvec(pt,pw)**2 - massvec(pw)*massvec(pt))) / 2._dp / massvec(pt) / massvec(pw)

        pout = vec + sinhx*(pt*dotvec(pw,vec) - pw*dotvec(pt,vec)) / &
                sqrt(dotvec(pt,pw)**2 - massvec(pw)*massvec(pt)) &
                + (coshx - 1._dp)*( dotvec(pt,pw)*(pt*dotvec(pw,vec) + pw*dotvec(pt,vec)) &
                    - massvec(pw)*pt*dotvec(pt,vec) - massvec(pt)*pw*dotvec(pw,vec) ) &
                    / ( dotvec(pt,pw)**2 - massvec(pw)*massvec(pt) )

    end subroutine

    function intdip_fi_mt_gg(p, recoiladd)
        use types
        use singletop2_nnlo_vars
        use singletop2_scale_m
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'masses.f'
        include 'constants.f'
        include 'alfacut.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'scheme.f'

        real(dp) :: intdip_fi_mt_gg
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in), optional :: recoiladd

        real(dp) :: puremass, ddilog

        real(dp) :: r, r1, rsq, logmu
        real(dp) :: pw(4)

        pw = p(3,:) + p(4,:)
        if (present(recoiladd)) then
            pw = pw + p(recoiladd,:)
        endif

        if (scheme == 'dred') then
            error stop __FILE__//': dred not supported in intdip_fi_mt!'
        endif

        r = puremass(pw)/mt
        rsq = r**2
        r1 = 1._dp - r**2

        if (corr_on_beam == 1) then
            logmu = log(renscale_beam2_isheavy_onheavy**2/mt**2)
        else
            logmu = log(renscale_beam1_isheavy_onheavy**2/mt**2)
        endif

!       ! Melnikov, Scharf, Schulze 2011 (1111.4991), eq. 22
        intdip_fi_mt_gg = epinv2*epinv/2._dp + 17._dp/(12._dp)*epinv + &
        (-8*afi**9*(-1 + r**2)**5 + 6*afi**8*(-1 + r**2)**4*(-5 + 7*r**2) - &
           afi**7*(-1 + r**2)**3*(67 - 162*r**2 + 115*r**4) + &
           afi**6*(-1 + r**2)**2*(-141 + 448*r**2 - 483*r**4 + 216*r**6) + &
           afi**5*(241 - 1219*r**2 + 2516*r**4 - 2524*r**6 + 1291*r**8 - &
              305*r**10) + afi**3*&
            (205 - 917*r**2 + 2078*r**4 - 1062*r**6 + 393*r**8 - 97*r**10) - &
           afi**2*(67 - 358*r**2 + 272*r**4 + 602*r**6 - 263*r**8 + &
              40*r**10) + afi**4*&
            (-283 + 1360*r**2 - 3050*r**4 + 2650*r**6 - 1195*r**8 + 278*r**10) &
            - 10*(-1 + r**2)**2*(-91 + 186*r**2 - 91*r**4 + &
              10*Pi**2*(-1 + r**2)**2) - &
           10*afi*(91 - 451*r**2 + 952*r**4 - 952*r**6 + 463*r**8 - &
              91*r**10 + 10*Pi**2*(-1 + r**2)**5) + &
           240*(-1 + r**2)**4*(1 + afi*(-1 + r**2))*ddilog(r1))/ &
         (240._dp*(1 + afi*(-1 + r**2))*r1**4) + &
        ((-23 + 24*afi - 3*afi**2 + 2*afi**3)*Log(afi))/12._dp - &
        Log(afi)**2 - (r**2*(-6*afi**3*(r**2 + r**4) + &
             3*afi**2*(r**2 + 5*r**4) + &
             12*afi*(1 - 3*r**2 + 5*r**4 - 4*r**6 + r**8) + &
             (12 - 21*r**2 + 11*r**4)*r1**2)*Log(r))/(6._dp*r1**5) - &
        (2*Log(1 + r))/3._dp + (-2.5_dp - epinv)*Log(r1) + Log(r1)**2 + &
        (r**2*(-4 + 13*r**2 - 23*r**4 + 16*r**6 - 4*r**8 - &
             2*afi**3*(r**2 + r**4) + afi**2*(r**2 + 5*r**4) + &
             4*afi*(1 - 3*r**2 + 5*r**4 - 4*r**6 + r**8))*Log(1 - afi*r1))/ &
         (4._dp*r1**5)
        intdip_fi_mt_gg = intdip_fi_mt_gg + logmu**2/4._dp + &
            logmu*(17._dp - 12._dp*Log(r1))/12._dp  + logmu/2._dp*epinv

        ! this fixes a typo in my MM notebook
        intdip_fi_mt_gg = intdip_fi_mt_gg - &
            (- 8*Log(1 + r) + 4*Log(r1))/12._dp


    end function

    ! Campbell, Ellis, Tramontano 2004, eq. 38, without CF
    ! Melnikov, Scharf, Schulze 2011 (1111.4991), eq. 8
    function intdip_fi_mt(p, recoiladd)
        use types
        use singletop2_nnlo_vars
        use singletop2_scale_m
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'masses.f'
        include 'constants.f'
        include 'alfacut.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'scheme.f'

        real(dp) :: intdip_fi_mt
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in), optional :: recoiladd

        real(dp) :: puremass, ddilog

        real(dp) :: r, omrsq, rsq
        real(dp) :: pw(4)
        real(dp) :: mulog, wlog, rlog

        pw = p(3,:) + p(4,:)
        if (present(recoiladd)) then
            pw = pw + p(recoiladd,:)
        endif

        r = puremass(pw)/mt
        rsq = r**2
        omrsq = 1._dp - r**2
        wlog = log(omrsq)
        rlog = log(rsq)

        if (corr_on_beam == 1) then
            mulog=log(renscale_beam2_isheavy_onheavy**2/mt**2)
        else
            mulog=log(renscale_beam1_isheavy_onheavy**2/mt**2)
        endif

        intdip_fi_mt = 0._dp

        ! fixed in tH-V scheme

        if (scheme == 'dred') then
            error stop __FILE__//': dred not supported in intdip_fi_mt!'
        endif

        intdip_fi_mt = (epinv2*epinv+epinv*mulog+0.5_dp*mulog**2) &
         +(epinv+mulog)*(2.5_dp-2._dp*wlog) & 
         +25._dp/4._dp+0.5_dp*(1._dp/omrsq**2-8._dp/omrsq+7._dp)*rlog &
         +0.5_dp/omrsq+2._dp*ddilog(omrsq)-5._dp*pisqo6 &
         -5._dp*wlog+2._dp*wlog**2+1._dp/2._dp &
         -2._dp*log(aif)**2-(3.5_dp-4._dp*aif+aif**2/2._dp)*log(aif) &
         +2._dp*(1._dp-aif)*rsq/omrsq*log(rsq/(1._dp-aif+rsq*aif))

    end function

    ! only for use in g -> b b~ splitting
    subroutine dips_ff_ident(nd,p,ip,jp,kp,sub,subv,msq,msqv,subr_born,subr_corr)
        use types
        use singletop2_scale_m
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'kprocess.f'
      include 'incldip.f'

      real(dp), intent(in) :: p(mxpart,4)

      real(dp):: ptrans(mxpart,4),sub(4),subv
      real(dp) :: ptrans_save(4)
      real(dp):: z,omz,y,omy,sij,sik,sjk,dot,vec(4)
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)
      real(dp) :: gsq
      integer:: nd,ip,jp,kp,nu
      external subr_born,subr_corr
            
      sub = 0._dp
      subv=0._dp
      msq = 0._dp
      msqv = 0._dp
      incldip(nd)=.true.
      
      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

      y=sij/(sij+sjk+sik)

        if (y > aff) then
          incldip(nd)=.false.
          return
        endif

        z=sik/(sjk+sik)
        omz=one-z
        omy=one-y

        call transform(p,ptrans,y,ip,jp,kp)

        ! hack for g -> b b~ splitting
        if (ip == 5 .and. jp == 8 .and. kp == 7) then
            ptrans_save(:) = ptrans(5,:)
            ptrans(5,:) = ptrans(7,:)
            ptrans(7,:) = ptrans_save(:)
        endif

        call storeptilde(nd,ptrans)

        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo

        call singletop2_set_dipscale(nd, ptrans, gsq)
        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,7,msqv)
         
        sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
        sub(gq)=gsq/sij
        sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
        subv   =+4._dp*gsq/sij/sij
      
    end
      


    subroutine dips_fi_mt_gg(ndip, p, sub_tree, sub_gvec, msqc, emitter, emitted, recoiladd)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'alfacut.f'
        include 'maxd.f'
        include 'incldip.f'
        include 'masses.f'
        include 'npart.f'

        integer, intent(in) :: ndip
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msqc(-nf:nf,-nf:nf)

        integer, intent(in) :: emitter, emitted
        integer, intent(in), optional :: recoiladd

        interface
            subroutine sub_tree(p,msq)
                use types
                implicit none
                include 'nf.f'
                include 'mxpart.f'
                real(dp), intent(in) :: p(mxpart,4)
                real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
            end subroutine

            subroutine sub_gvec(p,n,ip,msq)
                use types
                implicit none
                include 'nf.f'
                include 'mxpart.f'
                real(dp), intent(in) :: p(mxpart,4)
                real(dp), intent(in) :: n(4)
                integer, intent(in) :: ip
                real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
            end subroutine
        end interface

        real(dp) :: dotvec, massvec

        real(dp) :: pt(4), pw(4)
        real(dp) :: p3tilde(4), p4tilde(4), prectilde(4), p12tilde(4), pwtilde(4)
        real(dp) :: q(mxpart,4), msq(-nf:nf,-nf:nf)

        real(dp) :: pimu(4), amu(4)

        real(dp) :: omz, z, pwsq, xr, ymax, y, gsq, fyz
        real(dp) :: ptDpg1, ptDpg2, pgDpg

        msqc = 0._dp

        pt = -p(1,:) - p(2,:) - p(6,:)
        pw = p(3,:) + p(4,:)
        if (present(recoiladd)) then
            pw = pw + p(recoiladd,:)
        endif

        call ltrans(pt,pw, p(3,:), p3tilde)
        call ltrans(pt,pw, p(4,:), p4tilde)
        pwtilde = p3tilde + p4tilde

        if (present(recoiladd)) then
            call ltrans(pt,pw, p(recoiladd,:), prectilde)
            pwtilde = pwtilde + prectilde
        endif

        p12tilde = pt - pwtilde

        q(:,:) = 0._dp
        q(1,:) = p(1,:)
        q(2,:) = p(2,:)
        q(6,:) = p(6,:)
        q(3,:) = p3tilde
        q(4,:) = p4tilde
        q(5,:) = prectilde
        q(emitted,:) = 0._dp
        if (present(recoiladd)) then
            q(7,:) = p12tilde
        endif

        pgDpg = dotvec(p(emitter,:), p(emitted,:))
        ptDpg1 = dotvec(pt,p(emitter,:))
        ptDpg2 = dotvec(pt,p(emitted,:))

        omz = 2._dp*ptDpg1/(mt**2 - massvec(pw))

        z=1._dp-omz
        pwsq = massvec(pwtilde)

        xr=sqrt(pwsq/mt**2)
        ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
        y=2._dp*pgDpg/mt**2/(1._dp-xr)**2
        if ((z < 1._dp-afi) .AND. (y > afi*ymax)) then
            incldip(ndip) = .false. 
            return
        endif

        incldip(ndip) = .true.
        call storeptilde(ndip, q)

        call singletop2_set_dipscale(ndip,q,gsq)

        ! see Melnikov, Scharf, Schulze '11; 1111.4991
        amu = 2._dp/(mt**2*(1-xr**2))*(ptDpg2*p(emitter,:) - ptDpg1*p(emitted,:))
        !amu = z*p(emitter,:) - (1-z)*p(emitted,:)
        pimu = 1._dp/dotvec(pt,p12tilde)*(dotvec(pt,p12tilde)*amu - dotvec(p12tilde,amu)*pt)


!       ! checks
!       write (*,*) "CHECKS"
!       write (*,*) "XR z", xr, z
!       write (*,*) ptDpg1, mt**2/2._dp*(1-xr**2)*(1-z)
!       write (*,*) pgDpg, mt**2/2._dp*(1-xr)**2*y
!       write (*,*) dotvec(pt,p12tilde), mt**2/2._dp*(1-xr**2)
!       write (*,*) "RATIO", ptDpg1/ptDpg2, (1-z)/z

        ! g_mu_nu piece
        call sub_tree(q,msq)
        msqc(:,:) = gsq/(2._dp*pgDpg)*(-(z/omz - mt**2/4._dp*2*pgDpg/ptDpg1**2)*msq)

        ! pi_mu pi_nu piece
        call sub_gvec(q,pimu,7,msq)
        fyz = 4._dp/mt**4*(dotvec(pt,pw)**2 - xr**2*mt**4)/(1-xr**2)**2
        msqc(:,:) = msqc(:,:) - gsq/(2._dp*pgDpg)*( fyz/2._dp/pgDpg*msq )

    end subroutine

    subroutine dips_fi_mt(ndip, p, sub_tree, msqc, emitter, emitted, recoiladd)
        use types
        implicit none
        include 'nf.f'
        include 'mxpart.f'
        include 'alfacut.f'
        include 'maxd.f'
        include 'incldip.f'
        include 'masses.f'
        include 'npart.f'

        integer, intent(in) :: ndip
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msqc(-nf:nf,-nf:nf)

        integer, intent(in) :: emitter, emitted
        integer, intent(in), optional :: recoiladd

        interface
            subroutine sub_tree(p,msq)
                use types
                implicit none
                include 'nf.f'
                include 'mxpart.f'
                real(dp), intent(in) :: p(mxpart,4)
                real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
            end subroutine
        end interface

        real(dp) :: q(mxpart,4)
        real(dp) :: pbDpg, ptDpg, ptDpb
        real(dp) :: omz, pwsq, xr, ymax, z, y, fac, gsq

        real(dp) :: pt(4), pw(4), p3tilde(4), p4tilde(4), prectilde(4), pwtilde(4)
        real(dp) :: massvec,dotvec


        msqc = 0._dp

        pt = -p(1,:) - p(2,:) - p(6,:)
        pw = p(3,:) + p(4,:)
        if (present(recoiladd)) then
            pw = pw + p(recoiladd,:)
        endif

        pbDpg = dotvec(p(emitter,:),p(emitted,:))
        ptDpg = dotvec(pt,p(emitted,:))
        ptDpb = dotvec(pt,p(emitter,:))

        call ltrans(pt,pw, p(3,:), p3tilde)
        call ltrans(pt,pw, p(4,:), p4tilde)
        pwtilde = p3tilde + p4tilde

        if (present(recoiladd)) then
            call ltrans(pt,pw, p(recoiladd,:), prectilde)
            pwtilde = pwtilde + prectilde
        endif

        q(1:npart+2,:) = 0._dp
        q(1,:) = p(1,:)
        q(2,:) = p(2,:)
        q(6,:) = p(6,:)
        q(3,:) = p3tilde
        q(4,:) = p4tilde
        q(5,:) = pt - pwtilde
        q(emitted,:) = 0._dp
        if (present(recoiladd)) then
            q(7,:) = prectilde
        endif

        omz = 2._dp*ptDpg/(mt**2 - massvec(pw))
        !omz=ptDpg/(ptDpb+ptDpg-pbDpg)

        z=1._dp-omz
        pwsq = massvec(pwtilde)

        xr=sqrt(pwsq/mt**2)
        ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
        y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
        if ((z < 1._dp-aif) .AND. (y > aif*ymax)) then
            incldip(ndip) = .false. 
            return
        endif

        incldip(ndip) = .true.
        call storeptilde(ndip, q)

        call singletop2_set_dipscale(ndip,q,gsq)
        call sub_tree(q,msqc)

        fac=gsq*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)

        !msqc(:,:) = msqc(:,:) * fac

        msqc(5,[-3,-1,2,4]) = msqc(5,[-3,-1,2,4]) * fac
        msqc([-3,-1,2,4],5) = msqc([-3,-1,2,4],5) * fac

    end subroutine

    ! don't use this routine, ptilde is stored via wtransform
    ! only for dipole 1
    subroutine singletop2_heavy_decay_gs_all(p,ndmx,msqc)
    implicit none
    include 'types.f'
          

!     Matrix element for t-bbar production with radiation in decay
!      b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
!     averaged(summed) over initial(final) colours and spins
!--- g(p7) represents a gluon


    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'masses.f'
    include 'ptilde.f'
    include 'alfacut.f'
    include 'incldip.f'

    real(dp), intent(in) :: p(mxpart,4)
    integer, intent(in) :: ndmx
    real(dp), intent(out) :: msqc(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

    real(dp) :: msq(-nf:nf,-nf:nf)

    real(dp) :: q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb,pwsq,xr
    real(dp) :: y,ymax, gsq

    msqc = 0._dp

    ndmax=2
    incldip([1,2])= .TRUE. 

    call wtransform(p,q,pbDpg,ptDpg,ptDpb)
    omz=ptDpg/(ptDpb+ptDpg-pbDpg)
    z=1._dp-omz
    pwsq=2._dp*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))
    xr=sqrt(pwsq/mt**2)
    ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
    y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
    if ((z < 1._dp-afi) .AND. (y > afi*ymax)) then
        incldip([1,2])= .FALSE. 
        return
    endif

    corr_islight = .false.

    corr_beam1 = .false.
    call singletop2_set_dipscale(1,q,gsq)
    call singletop2_scet_tree_ub(q,msq)
    ! Campbell, Ellis, Tramontano, PR D70 094012, p.7
    fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)
    msqc(1, [2,4],5, 1,2) = fac*msq([2,4],5)
    msqc(1, [-1,-3],5, 1,2) = fac*msq([-1,-3],5)

    corr_beam1 = .true.
    call singletop2_set_dipscale(2,q,gsq)
    call singletop2_scet_tree_bu(q,msq)
    fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)
    msqc(2, 5,[2,4], 1,1) = fac*msq(5,[2,4])
    msqc(2, 5,[-1,-3], 1,1) = fac*msq(5,[-1,-3])

    end subroutine singletop2_heavy_decay_gs_all

end module
