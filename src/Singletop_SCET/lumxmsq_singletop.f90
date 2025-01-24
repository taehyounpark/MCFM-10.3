!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop2_scet_light
    use ieee_arithmetic
    use types
    use singletop2_scale_m
    use LHAPDF

    public :: singletop2_scet_tree
    public :: singletop2_scet_tree_ub
    public :: singletop2_scet_tree_bu

    public :: passed_taucut_light
    public :: lumxmsq_singletop
    public :: lumxmsq_singletop_new

    public :: singletop2_scet_virt_light
    public :: singletop2_scet_virt_light_all
    public :: singletop2_scet_z

    !public :: singletop2_scet_gs
    public :: singletop2_scet_gsall

    ! internal tree piece
    private :: qqbtbbar

    ! internal virt pieces
    private :: virtqqb_light

    private

#define DYNAMICTAU 1
    

    contains

    subroutine singletop2_scet_virt_light(p,msqv)
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
    real(dp), intent(out) :: msqv(-nf:nf,-nf:nf)
    real(dp):: fac

    scheme = 'tH-V'

    call spinoru(6,p,za,zb)

    msqv = 0._dp

    if (nwz == +1) then
        if (corr_on_beam == 1) then
            fac = (as_light_beam1/2._dp/pi*cf) * aveqq*gw**8*xn**2
            msqv([2,4],5) = fac*virtqqb_light(1,2,3,4,5,6, renscale_beam1_islight_onlight**2)
            msqv([-1,-3],5) = fac*virtqqb_light(6,2,3,4,5,1, renscale_beam1_islight_onlight**2)
        endif

        if (corr_on_beam == 2) then
            fac = (as_light_beam2/2._dp/pi*cf) * aveqq*gw**8*xn**2
            msqv(5,[2,4]) = fac*virtqqb_light(2,1,3,4,5,6, renscale_beam2_islight_onlight**2)
            msqv(5,[-1,-3]) = fac*virtqqb_light(6,1,3,4,5,2, renscale_beam2_islight_onlight**2)
        endif
    endif

    ! t~
    ! ub=fac*virtqqb_light(6,2,4,3,5,1)
    ! bu=fac*virtqqb_light(6,1,4,3,5,2)
    ! bubar=fac*virtqqb_light(2,1,4,3,5,6)
    ! ubarb=fac*virtqqb_light(1,2,4,3,5,6)

    end

    subroutine singletop2_scet_virt_light_all(p,msqvall)
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
    real(dp), intent(out) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
    real(dp):: fac

    scheme = 'tH-V'

    call spinoru(6,p,za,zb)

    msqvall = 0._dp

    if (nwz == +1) then
        fac = (as_light_beam1/2._dp/pi*cf) * aveqq*gw**8*xn**2
        msqvall([2,4],5, 1, 1) = fac*virtqqb_light(1,2,3,4,5,6, renscale_beam1_islight_onlight**2)
        msqvall([-1,-3],5, 1, 1) = fac*virtqqb_light(6,2,3,4,5,1, renscale_beam1_islight_onlight**2)

        fac = (as_light_beam2/2._dp/pi*cf) * aveqq*gw**8*xn**2
        msqvall(5,[2,4], 1, 2) = fac*virtqqb_light(2,1,3,4,5,6, renscale_beam2_islight_onlight**2)
        msqvall(5,[-1,-3], 1, 2) = fac*virtqqb_light(6,1,3,4,5,2, renscale_beam2_islight_onlight**2)
    endif

    ! t~
    ! ub=fac*virtqqb_light(6,2,4,3,5,1)
    ! bu=fac*virtqqb_light(6,1,4,3,5,2)
    ! bubar=fac*virtqqb_light(2,1,4,3,5,6)
    ! ubarb=fac*virtqqb_light(1,2,4,3,5,6)

    end

    function virtqqb_light(ju,jb,jn,je,jc,jd,musq)
    implicit none
    real(dp):: virtqqb_light

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

    ! only light line corrections for now
    cv = 0._dp
    c1 = 0._dp

    amp=za(jc,jn)*zb(ju,jb) &
     *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
    ampho=za(jc,jn)*zb(ju,jb) &
     *(cplx1(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)) &
     +c1*chalf*zb(je,jb)*za(jb,jd))

    virtqqb_light=real(amp*conjg(ampho))/prop

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
      end



      function qqbtbbar(ju,jb,jn,je,jc,jd)
          implicit none
          include 'types.f'
          real(dp):: qqbtbbar
          
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'cplx.h'
          include 'sprods_com.f'
          include 'masses.f'
          integer:: ju,jb,jn,je,jc,jd
          real(dp):: st,prop

          st=s(jn,je)+s(je,jc)+s(jn,jc)
          if (s(ju,jd) < 0._dp) then
          prop=(s(ju,jd)-wmass**2)**2
          else
          prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
          endif
          prop=prop &
              *((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2) &
              *((st-mt**2)**2+(mt*twidth)**2)
          qqbtbbar=s(jc,jn)*s(ju,jb) &
           *(-(s(ju,jd)+s(jb,jd))*(s(jn,je)+s(jc,je))-st*s(je,jd))/prop
      end

      subroutine singletop2_scet_tree_ub(p,msq)
          use singletop2_nnlo_vars
          !use NNTopDec
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'cplx.h'
          include 'ewcouple.f'
          include 'ckm.f'
          include 'nwz.f'
          include 'zprods_com.f'
          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
          real(dp):: fac,bu,ub,bubar,ubarb
          integer:: ib

          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2

          if     (nwz == +1) then
            ub=fac*qqbtbbar(1,2,3,4,5,6)
            ubarb=fac*qqbtbbar(6,2,3,4,5,1)
          elseif (nwz == -1) then
            ub=fac*qqbtbbar(6,2,4,3,5,1)
            bu=fac*qqbtbbar(6,1,4,3,5,2)
            ubarb=fac*qqbtbbar(1,2,4,3,5,6)
            bubar=fac*qqbtbbar(2,1,4,3,5,6)
          endif

          ib=5*nwz

          msq(:,:) = 0._dp

          if (nwz == +1) then
              msq([2,4],5) = ub
              msq([-1,-3],5) = ubarb
          elseif (nwz == -1) then
              write(6,*) 'Abort in singletop2_scet_tree_ub'
      stop
          endif

          !call top22treesqamp(p, newub)
          !write (*,*) "comparison = ", newub / msq(2,5)
          !pause
          
      end subroutine singletop2_scet_tree_ub

      subroutine singletop2_scet_tree_bu(p,msq)
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
          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
          real(dp):: fac,bu,ub,bubar,ubarb
          integer:: ib

          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2

          if     (nwz == +1) then
            bu=fac*qqbtbbar(2,1,3,4,5,6)
            bubar=fac*qqbtbbar(6,1,3,4,5,2)
          elseif (nwz == -1) then
            ub=fac*qqbtbbar(6,2,4,3,5,1)
            bu=fac*qqbtbbar(6,1,4,3,5,2)
            ubarb=fac*qqbtbbar(1,2,4,3,5,6)
            bubar=fac*qqbtbbar(2,1,4,3,5,6)
          endif

          ib=5*nwz

          msq(:,:) = 0._dp

          if (nwz == +1) then
              msq(5,[2,4]) = bu
              msq(5,[-1,-3]) = bubar
          elseif (nwz == -1) then
              stop 'Abort in singletop2_scet_tree_bu'
          endif
          
      end subroutine singletop2_scet_tree_bu

      subroutine singletop2_scet_tree(p,msq)
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
          real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
          real(dp):: fac,bu,ub,bubar,ubarb
          integer:: ib


          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2

          ub = 0._dp
          bu = 0._dp
          ubarb = 0._dp
          bubar = 0._dp

          if     (nwz == +1) then
            if (any(beams_enabled(1:maxbeams) == 1)) then
                ub=fac*qqbtbbar(1,2,3,4,5,6)
                ubarb=fac*qqbtbbar(6,2,3,4,5,1)
            endif

            if (any(beams_enabled(1:maxbeams) == 2)) then
                bu=fac*qqbtbbar(2,1,3,4,5,6)
                bubar=fac*qqbtbbar(6,1,3,4,5,2)
            endif
          elseif (nwz == -1) then
            ub=fac*qqbtbbar(6,2,4,3,5,1)
            bu=fac*qqbtbbar(6,1,4,3,5,2)
            ubarb=fac*qqbtbbar(1,2,4,3,5,6)
            bubar=fac*qqbtbbar(2,1,4,3,5,6)
          endif

          ib=5*nwz

          msq(:,:) = 0._dp

          if (nwz == +1) then
              msq([2,4],5) = ub
              msq([-1,-3],5) = ubarb
              msq(5,[2,4]) = bu
              msq(5,[-1,-3]) = bubar
          elseif (nwz == -1) then
              write(6,*) 'Abort in singletop2_scet_tree'
              stop
          endif
          
      end subroutine singletop2_scet_tree


    subroutine singletop2_scet_z(p,z)
        use singletop2_nnlo_vars
      implicit none

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      include 'nwz.f'

      real(dp), intent(in) :: p(mxpart,4), z
      real(dp) :: dot, if_qq, fi_qq, ii_qg

      integer :: is
      real(dp) :: xl16,xl26

      xl16 = log(-2*dot(p,1,6)/renscale_beam1_islight_onlight**2)
      xl26 = log(-2*dot(p,2,6)/renscale_beam2_islight_onlight**2)

      Q1 = zip
      Q2 = zip
      B1 = zip
      B2 = zip

      do is=1,3
          ! corr_on_beam = 1
          B1(q,q,b,is) = as_light_beam1/2/pi * cf*(if_qq(z,xl16,is) + fi_qq(z,xl16,is))
          Q1(q,g,q,is) = as_light_beam1/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam1_islight_onlight**2),is)

          ! corr_on_beam = 2
          B2(q,q,b,is) = as_light_beam2/2/pi * cf*(if_qq(z,xl26,is) + fi_qq(z,xl26,is))
          Q2(q,g,q,is) = as_light_beam2/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam2_islight_onlight**2),is)
      enddo

    end subroutine

!   subroutine singletop2_scet_gs(p,msq)
!       use singletop2_nnlo_vars
!       implicit none
!       include 'nf.f'
!       include 'mxpart.f'
!       include 'ptilde.f'
!       real(dp), intent(in) :: p(mxpart,4)
!       real(dp), intent(out) :: msq(maxd,-nf:nf,-nf:nf)

!       real(dp) :: msqall(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

!       call singletop2_scet_gsall(p,4,msqall)

!       !if (corr_on_beam == 1) then
!           msq = sum(msqall(:,:,:,:,1),4)
!       !else
!           !msq = sum(msqall(:,:,:,:,2),4)
!       !endif

!   end subroutine

!   subroutine singletop2_scet_gs(p,msq)
!       use singletop2_scale_m
!       use singletop2_nnlo_vars
!       implicit none
!       include 'constants.f'
!       include 'nf.f'
!       include 'mxpart.f'
!       include 'nwz.f'
!       include 'ptilde.f'
!       include 'qqgg.f'

!       real(dp), intent(in) :: p(mxpart,4)
!       real(dp), intent(out) :: msq(maxd,-nf:nf,-nf:nf)

!       real(dp) :: dummyv(-nf:nf,-nf:nf), dummy(-nf:nf,-nf:nf), dsubv
!       real(dp) :: msq17_6(-nf:nf,-nf:nf), msq67_1(-nf:nf,-nf:nf)
!       real(dp) :: msq27_6(-nf:nf,-nf:nf), msq67_2(-nf:nf,-nf:nf)
!       real(dp) :: msq17_2_light(-nf:nf,-nf:nf), msq16_2(-nf:nf,-nf:nf)
!       real(dp) :: msq27_1_light(-nf:nf,-nf:nf), msq26_1(-nf:nf,-nf:nf)
!       real(dp) :: sub17_6(4), sub67_1(4)
!       real(dp) :: sub27_6(4), sub67_2(4)
!       real(dp) :: sub17_2_light(4), sub16_2(4)
!       real(dp) :: sub27_1_light(4), sub26_1(4)

!       integer :: ib

!       integer :: noglue(4) = [-1,-3,2,4]
!       integer :: i,j

!       external donothing_gvec

!       msq17_6 = 0._dp
!       msq67_1 = 0._dp
!       sub17_6 = 0._dp
!       sub67_1 = 0._dp

!       msq17_2_light = 0._dp
!       sub17_2_light = 0._dp
!       msq16_2 = 0._dp
!       sub16_2 = 0._dp

!       msq27_6 = 0._dp
!       sub27_6 = 0._dp
!       msq67_2 = 0._dp
!       sub67_2 = 0._dp

!       msq27_1_light = 0._dp
!       sub27_1_light = 0._dp
!       msq26_1 = 0._dp
!       sub26_1 = 0._dp

!       ndmax = 4

!       msq = 0._dp
!       ib = 5*nwz
!       
!       corr_islight = .true.
!       
!       msq = 0._dp

!       if (corr_on_beam == 1) then
!           corr_beam1 = .true.
!           call dips(1,p, 1,7,6,sub17_6,dsubv,msq17_6,dummyv,singletop2_scet_tree_ub,donothing_gvec)
!           call dips(2,p, 6,7,1,sub67_1,dsubv,msq67_1,dummyv,singletop2_scet_tree_ub,donothing_gvec)

!           call dips(3,p, 1,7,2,sub17_2_light,dsubv,msq17_2_light,dummyv,singletop2_scet_tree_ub,donothing_gvec)
!           call dips(4,p, 1,6,2,sub16_2,dsubv,msq16_2,dummyv,singletop2_scet_tree_ub,donothing_gvec)

!           msq(1,noglue,5) = 2._dp*cf*sub17_6(qq)*msq17_6(noglue,5)
!           msq(2,noglue,5) = 2._dp*cf*sub67_1(qq)*msq67_1(noglue,5)

!           msq(3,0,5) = 2._dp*tr*sub17_2_light(qg) * sum(msq17_2_light(1:5,5))
!           msq(4,0,5) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1, 5))
!       else
!           corr_beam1 = .false.
!           call dips(1,p, 2,7,6,sub27_6,dsubv,msq27_6,dummyv,singletop2_scet_tree_bu,donothing_gvec)
!           call dips(2,p, 6,7,2,sub67_2,dsubv,msq67_2,dummyv,singletop2_scet_tree_bu,donothing_gvec)

!           call dips(3,p, 2,7,1,sub27_1_light,dsubv,msq27_1_light,dummyv,singletop2_scet_tree_bu,donothing_gvec)
!           call dips(4,p, 2,6,1,sub26_1,dsubv,msq26_1,dummyv,singletop2_scet_tree_bu,donothing_gvec)

!           msq(1,5,noglue) = 2._dp*cf*sub27_6(qq)*msq27_6(5,noglue)
!           msq(2,5,noglue) = 2._dp*cf*sub67_2(qq)*msq67_2(5,noglue)

!           msq(3,5,0) = 2._dp*tr*sub27_1_light(qg) * sum(msq27_1_light(5,1:5))
!           msq(4,5,0) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(5,-5:-1))
!       endif


!       !       msq = 0._dp
!       
!       !       ib = 5*nwz
!       !       do i=-nf,nf; do j=-nf,nf
!       !           if (i /= 0 .and. j == ib) then
!       !             msq(1,i,j) = 2._dp*cf*sub17_6(qq)*msq17_6(i,j) ! light
!       !             msq(2,i,j) = 2._dp*cf*sub67_1(qq)*msq67_1(i,j) ! light
!       !           elseif ((i == ib) .and. (j /= 0)) then
!       !             msq(3,i,j) = 2._dp*cf*sub27_6(qq)*msq27_6(i,j) ! light
!       !             msq(4,i,j) = 2._dp*cf*sub67_2(qq)*msq67_2(i,j) ! light
!       !           elseif ((i == 0) .and. (j == ib)) then
!       !             msq(5,i,j) = 2._dp*tr*sub17_2_light(qg) * sum(msq17_2_light(1:5,j)) ! light
!       !             msq(6,i,j) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1,j)) ! light
!       !           elseif ((i == ib) .and. (j==0)) then
!       !             msq(7,i,j) = 2._dp*tr*sub27_1_light(qg) * sum(msq27_1_light(i,1:5)) ! light
!       !             msq(8,i,j) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(i,-5:-1)) ! light            
!       !           endif
!       !       enddo; enddo

!   end subroutine

    subroutine singletop2_scet_gsall(p,ndmx,msq)
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
        real(dp), intent(inout) :: msq(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

        real(dp) :: dummyv(-nf:nf,-nf:nf),dsubv
        real(dp) :: msq17_6(-nf:nf,-nf:nf), msq67_1(-nf:nf,-nf:nf)
        real(dp) :: msq27_6(-nf:nf,-nf:nf), msq67_2(-nf:nf,-nf:nf)
        real(dp) :: msq17_2_light(-nf:nf,-nf:nf), msq16_2(-nf:nf,-nf:nf)
        real(dp) :: msq27_1_light(-nf:nf,-nf:nf), msq26_1(-nf:nf,-nf:nf)
        real(dp) :: sub17_6(4), sub67_1(4)
        real(dp) :: sub27_6(4), sub67_2(4)
        real(dp) :: sub17_2_light(4), sub16_2(4)
        real(dp) :: sub27_1_light(4), sub26_1(4)

        integer :: ib

        integer :: noglue(4) = [-1,-3,2,4]

        external donothing_gvec

        msq17_6 = 0._dp
        msq67_1 = 0._dp
        sub17_6 = 0._dp
        sub67_1 = 0._dp

        ndmax = 4

        msq = 0._dp
        ib = 5*nwz
        
        corr_islight = .true.

        !if (corr_on_beam == 1) then
            corr_on_beam = 1
            corr_beam1 = .true.
            call dips(1,p, 1,7,6,sub17_6,dsubv,msq17_6,dummyv,singletop2_scet_tree_ub,donothing_gvec)
            call dips(2,p, 6,7,1,sub67_1,dsubv,msq67_1,dummyv,singletop2_scet_tree_ub,donothing_gvec)

            call dips(3,p, 1,7,2,sub17_2_light,dsubv,msq17_2_light,dummyv,singletop2_scet_tree_ub,donothing_gvec)
            call dips(4,p, 1,6,2,sub16_2,dsubv,msq16_2,dummyv,singletop2_scet_tree_ub,donothing_gvec)

            msq(1,noglue,5, 1,1) = 2._dp*cf*sub17_6(qq)*msq17_6(noglue,5)
            msq(2,noglue,5, 1,1) = 2._dp*cf*sub67_1(qq)*msq67_1(noglue,5)

            msq(3,0,5, 1,1) = 2._dp*tr*sub17_2_light(qg) * sum(msq17_2_light(1:5,5))
            msq(4,0,5, 1,1) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1, 5))
        !else
            corr_on_beam = 2
            corr_beam1 = .false.
            call dips(5,p, 2,7,6,sub27_6,dsubv,msq27_6,dummyv,singletop2_scet_tree_bu,donothing_gvec)
            call dips(6,p, 6,7,2,sub67_2,dsubv,msq67_2,dummyv,singletop2_scet_tree_bu,donothing_gvec)

            call dips(7,p, 2,7,1,sub27_1_light,dsubv,msq27_1_light,dummyv,singletop2_scet_tree_bu,donothing_gvec)
            call dips(8,p, 2,6,1,sub26_1,dsubv,msq26_1,dummyv,singletop2_scet_tree_bu,donothing_gvec)

            msq(5,5,noglue, 1,2) = 2._dp*cf*sub27_6(qq)*msq27_6(5,noglue)
            msq(6,5,noglue, 1,2) = 2._dp*cf*sub67_2(qq)*msq67_2(5,noglue)

            msq(7,5,0, 1,2) = 2._dp*tr*sub27_1_light(qg) * sum(msq27_1_light(5,1:5))
            msq(8,5,0, 1,2) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(5,-5:-1))
        !endif


!       msq = 0._dp

!       ib = 5*nwz
!       do i=-nf,nf; do j=-nf,nf
!           if (i /= 0 .and. j == ib) then
!             msq(1,i,j) = 2._dp*cf*sub17_6(qq)*msq17_6(i,j) ! light
!             msq(2,i,j) = 2._dp*cf*sub67_1(qq)*msq67_1(i,j) ! light
!           elseif ((i == ib) .and. (j /= 0)) then
!             msq(3,i,j) = 2._dp*cf*sub27_6(qq)*msq27_6(i,j) ! light
!             msq(4,i,j) = 2._dp*cf*sub67_2(qq)*msq67_2(i,j) ! light
!           elseif ((i == 0) .and. (j == ib)) then
!             msq(5,i,j) = 2._dp*tr*sub17_2_light(qg) * sum(msq17_2_light(1:5,j)) ! light
!             msq(6,i,j) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1,j)) ! light
!           elseif ((i == ib) .and. (j==0)) then
!             msq(7,i,j) = 2._dp*tr*sub27_1_light(qg) * sum(msq27_1_light(i,1:5)) ! light
!             msq(8,i,j) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(i,-5:-1)) ! light            
!           endif
!       enddo; enddo


    end subroutine


    ! like the routine without "_new", but does use assembly routine
    ! that directly takes coefficients in log(tau)
    subroutine lumxmsq_singletop_new(p,xx,z1,z2,QB,order,xmsq,central)
        use types
        use SCET
        use SCET_Jet
        use SCET_Beamfunctions
        use singletop2_nnlo_vars, only: maxbeams, beams_enabled
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'taucut.f'
        include 'beamtype.f'
    
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in) :: xx(2), z1, z2, QB(2)
        integer, intent(in) :: order
        real(dp), intent(out) :: xmsq
        logical, intent(in) :: central

        real(dp) :: msqlo(-nf:nf, -nf:nf)

        real(dp) :: beama0(-5:5), beamb0(-5:5)
        real(dp) :: beama1(-5:5, -1:1), beamb1(-5:5, -1:1)
        real(dp) :: beama2(-5:5, -1:3), beamb2(-5:5, -1:3)

        real(dp) :: soft1(-1:1), soft2(-1:3)
        real(dp) :: jet1q(-1:1), jet2q(-1:3)
        real(dp) :: hard(0:2)

        real(dp) :: soft(2,0:4), jet(2,0:4), beam(2,0:4)

        integer :: i,j,k

        real(dp) :: dot, dotvec
        real(dp) :: logy16, Qj, Qsq

        real(dp), allocatable :: taucutsx(:), coeffsx(:)

        !integer, allocatable :: xjx(:)
        integer, parameter :: xjx(4) = [-3,-1,2,4]

        ! on heavy and light lines separately

        Qj = 2._dp*p(6,4)

        xmsq = 0._dp

        call singletop2_scet_tree(p,msqlo)

        if (central .and. doMultitaucut) then
            allocate(taucutsx(size(tcutarray)+1))
            allocate(coeffsx(size(tcutarray)+1))
            scetreweight(:) = 0._dp
        else
            allocate(taucutsx(1))
            allocate(coeffsx(1))
        endif
        taucutsx(1) = taucut

        if (size(taucutsx) > 1) then
            taucutsx(2:) = tcutarray(:)
        endif

        !allocate(xjx, source=[2,4,-1,-3])

        if (maxbeams > 1) then
            error stop "light line below piece only with one beam at a time"
        endif

#if DYNAMICTAU == 1
        if (beams_enabled(1) == 1) then
            Qsq = -dotvec(p(2,:)+p(3,:)+p(4,:)+p(5,:), &
                          p(2,:)+p(3,:)+p(4,:)+p(5,:))
        else
            Qsq = -dotvec(p(1,:)+p(3,:)+p(4,:)+p(5,:), &
                          p(1,:)+p(3,:)+p(4,:)+p(5,:))
        endif

        taucutsx(:) = taucutsx(:) * sqrt(Qsq) / 80._dp
#endif

        if (any(beams_enabled(1:maxbeams) == 1)) then

        ! COMPUTE ALL SCALE DEPENDENT INGREDIENTS
        ! THIS WILL BE REPEATED BELOW FOR CROSSING
        logy16 = log(-2._dp*dot(p,1,6)/((-2._dp*p(1,4))*2._dp*p(6,4)))
        call softqq(order,soft1,soft2,logy16)

        ! convert to Laplace-space log notation
        ! multiply by 2 because soft1 is in 2pi normalization
        soft = 0._dp
        soft(1,0) = soft1(-1)*2._dp
        soft(1,1) = soft1(0)*2._dp
        soft(1,2) = soft1(1)

        ! divide by 4 because soft1 is in 2pi normalization
        soft(2,0) = soft2(-1)*4._dp
        soft(2,1) = soft2(0)*4._dp
        soft(2,2) = soft2(1)*2._dp
        soft(2,3) = soft2(2)/3._dp*4._dp
        soft(2,4) = soft2(3)

        call jetq(order,Qj,jet1q,jet2q, renscale_beam1_islight_onlight, disttau=.false.)

        jet = 0._dp
        jet(1,0) = jet1q(-1)
        jet(1,1) = jet1q(0)
        jet(1,2) = jet1q(1)/2._dp

        jet(2,0) = jet2q(-1)
        jet(2,1) = jet2q(0)
        jet(2,2) = jet2q(1)/2._dp
        jet(2,3) = jet2q(2)/3._dp
        jet(2,4) = jet2q(3)/4._dp
    
        if (order >= 0) then
            call fdist(ih1,xx(1),facscale_beam1_islight_onlight,beama0,1)
            call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0,2)
        endif
    
        if (order >= 1) then
            call xbeam1bis_new(ih1,z1,xx(1),QB(1),beama1, 1, renscale_beam1_islight_onlight, &
                facscale_beam1_islight_onlight, disttau=.false.)
        endif

        if (order >= 2) then
            call xbeam2bis_new(ih1,z1,xx(1),QB(1),beama2, 1, renscale_beam1_islight_onlight, &
                facscale_beam1_islight_onlight, disttau=.false.)
        endif

        call hardqq(2._dp*dot(p,1,6),renscale_beam1_islight_onlight**2,hard(1:2))
        hard(0) = 1._dp
        hard(1) = hard(1)*2._dp
        hard(2) = hard(2)*4._dp
        ! END COMPUTE ALL SCALE DEPENDENT INGREDIENTS

        ! b on second beam
        do i=1,size(xjx)
            j = xjx(i)
            k = 5

            beam = 0._dp
            beam(1,0) = beama1(j, -1)
            beam(1,1) = beama1(j, 0)
            beam(1,2) = beama1(j, 1)/2._dp

            beam(2,0) = beama2(j, -1)
            beam(2,1) = beama2(j, 0)
            beam(2,2) = beama2(j, 1)/2._dp
            beam(2,3) = beama2(j, 2)/3._dp
            beam(2,4) = beama2(j, 3)/4._dp

            coeffsx(:) = assemble_light_pieces(order, Qj, QB(1), renscale_beam1_islight_onlight, as_light_beam1, &
                taucutsx, hard, jet, soft, beam, beama0(j)) * beamb0(k) * msqlo(j,k)

            !write (*,*) "beama0", beama0(j)
            !write (*,*) "Coeffsx ",j,k, coeffsx(1)

            xmsq = xmsq + coeffsx(1)
            if (size(taucutsx) > 1) then
                scetreweight(:) = scetreweight(:) + coeffsx(2:)
            endif
        enddo

        endif

        if (any(beams_enabled(1:maxbeams) == 2)) then

        ! b on first beam, recompute ingredients
        ! notice the beama beamb exchange for the do_assemble_beam1 routine

        ! COMPUTE ALL SCALE DEPENDENT INGREDIENTS, SEE ABOVE
        logy16 = log(-2._dp*dot(p,2,6)/((-2._dp*p(2,4))*2._dp*p(6,4)))
        call softqq(order,soft1,soft2,logy16)
        ! convert to Laplace-space log notation
        ! multiply by 2 because soft1 is in 2pi normalization
        soft = 0._dp
        soft(1,0) = soft1(-1)*2._dp
        soft(1,1) = soft1(0)*2._dp
        soft(1,2) = soft1(1)

        ! divide by 4 because soft1 is in 2pi normalization
        soft(2,0) = soft2(-1)*4._dp
        soft(2,1) = soft2(0)*4._dp
        soft(2,2) = soft2(1)*2._dp
        soft(2,3) = soft2(2)/3._dp*4._dp
        soft(2,4) = soft2(3)

        call jetq(order,2._dp*p(6,4),jet1q,jet2q, renscale_beam2_islight_onlight, disttau=.false.)
        jet = 0._dp
        jet(1,0) = jet1q(-1)
        jet(1,1) = jet1q(0)
        jet(1,2) = jet1q(1)/2._dp

        jet(2,0) = jet2q(-1)
        jet(2,1) = jet2q(0)
        jet(2,2) = jet2q(1)/2._dp
        jet(2,3) = jet2q(2)/3._dp
        jet(2,4) = jet2q(3)/4._dp
    
        if (order >= 0) then
            call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0,1)
            call fdist(ih2,xx(2),facscale_beam2_islight_onlight,beamb0,2)
        endif

        if (order >= 1) then
            call xbeam1bis_new(ih2,z2,xx(2),QB(2),beamb1, 2, renscale_beam2_islight_onlight, &
                facscale_beam2_islight_onlight, disttau=.false.)
        endif

        if (order >= 2) then
            call xbeam2bis_new(ih2,z2,xx(2),QB(2),beamb2, 2, renscale_beam2_islight_onlight, &
                facscale_beam2_islight_onlight, disttau=.false.)
        endif


        call hardqq(2._dp*dot(p,2,6),renscale_beam2_islight_onlight**2,hard(1:2))
        hard(0) = 1._dp
        hard(1) = hard(1)*2._dp
        hard(2) = hard(2)*4._dp
        ! END COMPUTE ALL SCALE DEPENDENT INGREDIENTS

        do i=1,size(xjx)
            j = 5
            k = xjx(i)

            beam = 0._dp
            beam(1,0) = beamb1(k, -1)
            beam(1,1) = beamb1(k, 0)
            beam(1,2) = beamb1(k, 1)/2._dp

            beam(2,0) = beamb2(k, -1)
            beam(2,1) = beamb2(k, 0)
            beam(2,2) = beamb2(k, 1)/2._dp
            beam(2,3) = beamb2(k, 2)/3._dp
            beam(2,4) = beamb2(k, 3)/4._dp

            coeffsx(:) = assemble_light_pieces(order, Qj, QB(2), renscale_beam2_islight_onlight, as_light_beam2, &
                taucutsx, hard, jet, soft, beam, beamb0(k)) * beama0(j) * msqlo(j,k)

            xmsq = xmsq + coeffsx(1)
            if (size(taucutsx) > 1) then
                scetreweight(:) = scetreweight(:) + coeffsx(2:)
            endif
        enddo

        endif

        if (size(taucutsx) > 1 .and. xmsq /= 0._dp) then
            scetreweight(:) = scetreweight(:) / xmsq
        endif

    end subroutine

    subroutine lumxmsq_singletop(p,xx,z1,z2,QB,order,xmsq,central)
        use types
        use SCET
        use SCET_Jet
        use SCET_Beamfunctions
        use singletop2_nnlo_vars, only: maxbeams, beams_enabled
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'taucut.f'
        include 'beamtype.f'
    
        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(in) :: xx(2), z1, z2, QB(2)
        integer, intent(in) :: order
        real(dp), intent(out) :: xmsq
        logical, intent(in) :: central

        real(dp) :: msqlo(-nf:nf, -nf:nf)

        real(dp) :: beama0(-5:5), beamb0(-5:5)
        real(dp) :: beama1(-5:5, -1:1), beamb1(-5:5, -1:1)
        real(dp) :: beama2(-5:5, -1:3), beamb2(-5:5, -1:3)

        real(dp) :: soft1(-1:1), soft2(-1:3)
        real(dp) :: jet1q(-1:1), jet2q(-1:3)
        real(dp) :: hard(2)

        integer :: i,j,k

        real(dp) :: dot
        real(dp) :: logy16

        real(dp), allocatable :: taucuts(:), coeffs(:)

        !integer, allocatable :: xj(:)
        integer, parameter :: xj(4) = [-3,-1,2,4]

        ! on heavy and light lines separately

        xmsq = 0._dp

        call singletop2_scet_tree(p,msqlo)

        if (central .and. doMultitaucut) then
            allocate(taucuts(size(tcutarray)+1))
            allocate(coeffs(size(tcutarray)+1))
            scetreweight(:) = 0._dp
        else
            allocate(taucuts(1))
            allocate(coeffs(1))
        endif
        taucuts(1) = taucut

        if (size(taucuts) > 1) then
            taucuts(2:) = tcutarray(:)
        endif

        !allocate(xj, source=[2,4,-1,-3])

        if (any(beams_enabled(1:maxbeams) == 1)) then

        ! COMPUTE ALL SCALE DEPENDENT INGREDIENTS
        ! THIS WILL BE REPEATED BELOW FOR CROSSING
        logy16 = log(-2._dp*dot(p,1,6)/((-2._dp*p(1,4))*2._dp*p(6,4)))
        call softqq(order,soft1,soft2,logy16)

        call jetq(order,2._dp*p(6,4),jet1q,jet2q, renscale_beam1_islight_onlight, disttau=.true.)
    
        if (order >= 0) then
            call fdist(ih1,xx(1),facscale_beam1_islight_onlight,beama0,1)
            call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0,2)
        endif
    
        if (order >= 1) then
            call xbeam1bis_new(ih1,z1,xx(1),QB(1),beama1, 1, renscale_beam1_islight_onlight, &
                facscale_beam1_islight_onlight, disttau=.true.)
        endif

        if (order >= 2) then
            call xbeam2bis_new(ih1,z1,xx(1),QB(1),beama2, 1, renscale_beam1_islight_onlight, &
                facscale_beam1_islight_onlight, disttau=.true.)
        endif

        call hardqq(2._dp*dot(p,1,6),renscale_beam1_islight_onlight**2,hard)
        ! END COMPUTE ALL SCALE DEPENDENT INGREDIENTS

        ! b on second beam
        do i=1,size(xj)
            j = xj(i)
            k = 5
            call do_assemble_beam1(order, renscale_beam1_islight_onlight, as_light_beam1, &
                   taucuts, beama0(j), beamb0(k), &
                   beama1(j,:), beama2(j,:), &
                   soft1, soft2, jet1q, jet2q, hard, coeffs)

            xmsq = xmsq + msqlo(j,k)*coeffs(1)
            if (size(taucuts) > 1) then
                scetreweight(:) = scetreweight(:) + msqlo(j,k)*coeffs(2:)
            endif
        enddo

        endif

        if (any(beams_enabled(1:maxbeams) == 2)) then

        ! b on first beam, recompute ingredients
        ! notice the beama beamb exchange for the do_assemble_beam1 routine

        ! COMPUTE ALL SCALE DEPENDENT INGREDIENTS, SEE ABOVE
        logy16 = log(-2._dp*dot(p,2,6)/((-2._dp*p(2,4))*2._dp*p(6,4)))
        call softqq(order,soft1,soft2,logy16)

        call jetq(order,2._dp*p(6,4),jet1q,jet2q, renscale_beam2_islight_onlight, disttau=.true.)
    
        if (order >= 0) then
            call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0,1)
            call fdist(ih2,xx(2),facscale_beam2_islight_onlight,beamb0,2)
        endif

        if (order >= 1) then
            call xbeam1bis_new(ih2,z2,xx(2),QB(2),beamb1, 2, renscale_beam2_islight_onlight, &
                facscale_beam2_islight_onlight, disttau=.true.)
        endif

        if (order >= 2) then
            call xbeam2bis_new(ih2,z2,xx(2),QB(2),beamb2, 2, renscale_beam2_islight_onlight, &
                facscale_beam2_islight_onlight, disttau=.true.)
        endif

        call hardqq(2._dp*dot(p,2,6),renscale_beam2_islight_onlight**2,hard)
        ! END COMPUTE ALL SCALE DEPENDENT INGREDIENTS

        do i=1,size(xj)
            j = 5
            k = xj(i)
            call do_assemble_beam1(order, renscale_beam2_islight_onlight, as_light_beam2, &
                   taucuts, beamb0(k), beama0(j), &
                   beamb1(k,:), beamb2(k,:), &
                   soft1, soft2, jet1q, jet2q, hard, coeffs)
            xmsq = xmsq + msqlo(j,k)*coeffs(1)
            if (size(taucuts) > 1) then
                scetreweight(:) = scetreweight(:) + msqlo(j,k)*coeffs(2:)
            endif
        enddo
        
        endif

        if (size(taucuts) > 1 .and. xmsq /= 0._dp) then
            scetreweight(:) = scetreweight(:) / xmsq
        endif

    end subroutine

    function assemble_light_pieces(order, Qj, Qb, scale, as, taucut, h, j, s, b, b0)
        use types
        use constants
        implicit none
        include 'masses.f'
        include 'kpart.f'

        integer, intent(in) :: order
        real(dp), intent(in) :: Qj, Qb, scale, as
        real(dp), intent(in) :: taucut(:)
        real(dp), intent(in) :: h(0:2)
        real(dp), intent(in) :: j(2,0:4)
        real(dp), intent(in) :: s(2,0:4)
        real(dp), intent(in) :: b(2,0:4)
        real(dp), intent(in) :: b0
        real(dp) :: assemble_light_pieces(size(taucut(:)))

        real(dp) :: full(0:2,0:4)
        real(dp) :: LSX, LJX, LBX

        ! expressing everything in terms of LJX/LSX + log(tau) allows
        ! us to repurpose this assembly function for the production process
        LJX = log(Qj/scale**2)
        LBX = log(Qb/scale**2)
        LSX = log(1._dp/scale)

        full(:,:) = 0._dp

        if (h(0) /= 1._dp) then
            write (*,*) "WARNING: bad hard function normalization: ", h(0)
        endif

        if (coeffonly) then
            full(0,0) = 0._dp
        else
            full(0,0) = b0
        endif

        if ((order == 1) .or. (order >= 1 .and. coeffonly .eqv. .false.)) then
            full(1,2) = b(1,2) + b0*j(1,2) + b0*s(1,2)
            full(1,1) = b(1,1) + 2*LBX*b(1,2) + b0*j(1,1) + 2*b0*LJX*j(1,2) + b0*s(1,1) + &
                2*b0*LSX*s(1,2)
            full(1,0) = b(1,0) + LBX*b(1,1) + LBX**2*b(1,2) + b0*h(1) + b0*j(1,0) + &
                b0*LJX*j(1,1) + b0*LJX**2*j(1,2) + b0*s(1,0) + b0*LSX*s(1,1) + &
                b0*LSX**2*s(1,2)
        endif

        if (order >= 2) then
            full(2,4) = b(2,4) + b(1,2)*j(1,2) + b0*j(2,4) + b(1,2)*s(1,2) + b0*j(1,2)*s(1,2) + &
                b0*s(2,4)
            full(2,3) = b(2,3) + 4*LBX*b(2,4) + b(1,2)*j(1,1) + b(1,1)*j(1,2) + &
                2*LBX*b(1,2)*j(1,2) + 2*LJX*b(1,2)*j(1,2) + b0*j(2,3) + &
                4*b0*LJX*j(2,4) + b(1,2)*s(1,1) + b0*j(1,2)*s(1,1) + b(1,1)*s(1,2) + &
                2*LBX*b(1,2)*s(1,2) + 2*LSX*b(1,2)*s(1,2) + b0*j(1,1)*s(1,2) + &
                2*b0*LJX*j(1,2)*s(1,2) + 2*b0*LSX*j(1,2)*s(1,2) + b0*s(2,3) + &
                4*b0*LSX*s(2,4)
            full(2,2) = b(2,2) + 3*LBX*b(2,3) + 6*LBX**2*b(2,4) + b(1,2)*h(1) + b(1,2)*j(1,0) + &
                b(1,1)*j(1,1) + 2*LBX*b(1,2)*j(1,1) + LJX*b(1,2)*j(1,1) + &
                b(1,0)*j(1,2) + LBX*b(1,1)*j(1,2) + 2*LJX*b(1,1)*j(1,2) + &
                LBX**2*b(1,2)*j(1,2) + 4*LBX*LJX*b(1,2)*j(1,2) + LJX**2*b(1,2)*j(1,2) - &
                (2*Pi**2*b(1,2)*j(1,2))/3._dp + b0*h(1)*j(1,2) + b0*j(2,2) + &
                3*b0*LJX*j(2,3) + 6*b0*LJX**2*j(2,4) + b(1,2)*s(1,0) + &
                b0*j(1,2)*s(1,0) + b(1,1)*s(1,1) + 2*LBX*b(1,2)*s(1,1) + &
                LSX*b(1,2)*s(1,1) + b0*j(1,1)*s(1,1) + 2*b0*LJX*j(1,2)*s(1,1) + &
                b0*LSX*j(1,2)*s(1,1) + b(1,0)*s(1,2) + LBX*b(1,1)*s(1,2) + &
                2*LSX*b(1,1)*s(1,2) + LBX**2*b(1,2)*s(1,2) + 4*LBX*LSX*b(1,2)*s(1,2) + &
                LSX**2*b(1,2)*s(1,2) - (2*Pi**2*b(1,2)*s(1,2))/3._dp + b0*h(1)*s(1,2) + &
                b0*j(1,0)*s(1,2) + b0*LJX*j(1,1)*s(1,2) + 2*b0*LSX*j(1,1)*s(1,2) + &
                b0*LJX**2*j(1,2)*s(1,2) + 4*b0*LJX*LSX*j(1,2)*s(1,2) + &
                b0*LSX**2*j(1,2)*s(1,2) - (2*b0*Pi**2*j(1,2)*s(1,2))/3._dp + b0*s(2,2) + &
                3*b0*LSX*s(2,3) + 6*b0*LSX**2*s(2,4)
            full(2,1) = b(2,1) + 2*LBX*b(2,2) + 3*LBX**2*b(2,3) + 4*LBX**3*b(2,4) + &
                b(1,1)*h(1) + 2*LBX*b(1,2)*h(1) + b(1,1)*j(1,0) + 2*LBX*b(1,2)*j(1,0) + &
                b(1,0)*j(1,1) + LBX*b(1,1)*j(1,1) + LJX*b(1,1)*j(1,1) + &
                LBX**2*b(1,2)*j(1,1) + 2*LBX*LJX*b(1,2)*j(1,1) - &
                (Pi**2*b(1,2)*j(1,1))/3._dp + b0*h(1)*j(1,1) + 2*LJX*b(1,0)*j(1,2) + &
                2*LBX*LJX*b(1,1)*j(1,2) + LJX**2*b(1,1)*j(1,2) - &
                (Pi**2*b(1,1)*j(1,2))/3._dp + 2*LBX**2*LJX*b(1,2)*j(1,2) + &
                2*LBX*LJX**2*b(1,2)*j(1,2) - (2*LBX*Pi**2*b(1,2)*j(1,2))/3._dp - &
                (2*LJX*Pi**2*b(1,2)*j(1,2))/3._dp + 8*zeta3*b(1,2)*j(1,2) + &
                2*b0*LJX*h(1)*j(1,2) + b0*j(2,1) + 2*b0*LJX*j(2,2) + &
                3*b0*LJX**2*j(2,3) + 4*b0*LJX**3*j(2,4) + b(1,1)*s(1,0) + &
                2*LBX*b(1,2)*s(1,0) + b0*j(1,1)*s(1,0) + 2*b0*LJX*j(1,2)*s(1,0) + &
                b(1,0)*s(1,1) + LBX*b(1,1)*s(1,1) + LSX*b(1,1)*s(1,1) + &
                LBX**2*b(1,2)*s(1,1) + 2*LBX*LSX*b(1,2)*s(1,1) - &
                (Pi**2*b(1,2)*s(1,1))/3._dp + b0*h(1)*s(1,1) + b0*j(1,0)*s(1,1) + &
                b0*LJX*j(1,1)*s(1,1) + b0*LSX*j(1,1)*s(1,1) + b0*LJX**2*j(1,2)*s(1,1) + &
                2*b0*LJX*LSX*j(1,2)*s(1,1) - (b0*Pi**2*j(1,2)*s(1,1))/3._dp + &
                2*LSX*b(1,0)*s(1,2) + 2*LBX*LSX*b(1,1)*s(1,2) + LSX**2*b(1,1)*s(1,2) - &
                (Pi**2*b(1,1)*s(1,2))/3._dp + 2*LBX**2*LSX*b(1,2)*s(1,2) + &
                2*LBX*LSX**2*b(1,2)*s(1,2) - (2*LBX*Pi**2*b(1,2)*s(1,2))/3._dp - &
                (2*LSX*Pi**2*b(1,2)*s(1,2))/3._dp + 8*zeta3*b(1,2)*s(1,2) + &
                2*b0*LSX*h(1)*s(1,2) + 2*b0*LSX*j(1,0)*s(1,2) + &
                2*b0*LJX*LSX*j(1,1)*s(1,2) + b0*LSX**2*j(1,1)*s(1,2) - &
                (b0*Pi**2*j(1,1)*s(1,2))/3._dp + 2*b0*LJX**2*LSX*j(1,2)*s(1,2) + &
                2*b0*LJX*LSX**2*j(1,2)*s(1,2) - (2*b0*LJX*Pi**2*j(1,2)*s(1,2))/3._dp - &
                (2*b0*LSX*Pi**2*j(1,2)*s(1,2))/3._dp + 8*b0*zeta3*j(1,2)*s(1,2) + &
                b0*s(2,1) + 2*b0*LSX*s(2,2) + 3*b0*LSX**2*s(2,3) + 4*b0*LSX**3*s(2,4)
            full(2,0) = b(2,0) + LBX*b(2,1) + LBX**2*b(2,2) + LBX**3*b(2,3) + LBX**4*b(2,4) + &
                b(1,0)*h(1) + LBX*b(1,1)*h(1) + LBX**2*b(1,2)*h(1) + b0*h(2) + &
                b(1,0)*j(1,0) + LBX*b(1,1)*j(1,0) + LBX**2*b(1,2)*j(1,0) + &
                b0*h(1)*j(1,0) + LJX*b(1,0)*j(1,1) + LBX*LJX*b(1,1)*j(1,1) - &
                (Pi**2*b(1,1)*j(1,1))/6._dp + LBX**2*LJX*b(1,2)*j(1,1) - &
                (LBX*Pi**2*b(1,2)*j(1,1))/3._dp + 2*zeta3*b(1,2)*j(1,1) + &
                b0*LJX*h(1)*j(1,1) + LJX**2*b(1,0)*j(1,2) + LBX*LJX**2*b(1,1)*j(1,2) - &
                (LJX*Pi**2*b(1,1)*j(1,2))/3._dp + 2*zeta3*b(1,1)*j(1,2) + &
                LBX**2*LJX**2*b(1,2)*j(1,2) - (2*LBX*LJX*Pi**2*b(1,2)*j(1,2))/3._dp - &
                (Pi**4*b(1,2)*j(1,2))/90._dp + 4*LBX*zeta3*b(1,2)*j(1,2) + &
                4*LJX*zeta3*b(1,2)*j(1,2) + b0*LJX**2*h(1)*j(1,2) + b0*j(2,0) + &
                b0*LJX*j(2,1) + b0*LJX**2*j(2,2) + b0*LJX**3*j(2,3) + &
                b0*LJX**4*j(2,4) + b(1,0)*s(1,0) + LBX*b(1,1)*s(1,0) + &
                LBX**2*b(1,2)*s(1,0) + b0*h(1)*s(1,0) + b0*j(1,0)*s(1,0) + &
                b0*LJX*j(1,1)*s(1,0) + b0*LJX**2*j(1,2)*s(1,0) + LSX*b(1,0)*s(1,1) + &
                LBX*LSX*b(1,1)*s(1,1) - (Pi**2*b(1,1)*s(1,1))/6._dp + &
                LBX**2*LSX*b(1,2)*s(1,1) - (LBX*Pi**2*b(1,2)*s(1,1))/3._dp + &
                2*zeta3*b(1,2)*s(1,1) + b0*LSX*h(1)*s(1,1) + b0*LSX*j(1,0)*s(1,1) + &
                b0*LJX*LSX*j(1,1)*s(1,1) - (b0*Pi**2*j(1,1)*s(1,1))/6._dp + &
                b0*LJX**2*LSX*j(1,2)*s(1,1) - (b0*LJX*Pi**2*j(1,2)*s(1,1))/3._dp + &
                2*b0*zeta3*j(1,2)*s(1,1) + LSX**2*b(1,0)*s(1,2) + &
                LBX*LSX**2*b(1,1)*s(1,2) - (LSX*Pi**2*b(1,1)*s(1,2))/3._dp + &
                2*zeta3*b(1,1)*s(1,2) + LBX**2*LSX**2*b(1,2)*s(1,2) - &
                (2*LBX*LSX*Pi**2*b(1,2)*s(1,2))/3._dp - (Pi**4*b(1,2)*s(1,2))/90._dp + &
                4*LBX*zeta3*b(1,2)*s(1,2) + 4*LSX*zeta3*b(1,2)*s(1,2) + &
                b0*LSX**2*h(1)*s(1,2) + b0*LSX**2*j(1,0)*s(1,2) + &
                b0*LJX*LSX**2*j(1,1)*s(1,2) - (b0*LSX*Pi**2*j(1,1)*s(1,2))/3._dp + &
                2*b0*zeta3*j(1,1)*s(1,2) + b0*LJX**2*LSX**2*j(1,2)*s(1,2) - &
                (2*b0*LJX*LSX*Pi**2*j(1,2)*s(1,2))/3._dp - (b0*Pi**4*j(1,2)*s(1,2))/90._dp + &
                4*b0*LJX*zeta3*j(1,2)*s(1,2) + 4*b0*LSX*zeta3*j(1,2)*s(1,2) + &
                b0*s(2,0) + b0*LSX*s(2,1) + b0*LSX**2*s(2,2) + b0*LSX**3*s(2,3) + &
                b0*LSX**4*s(2,4)
        endif

        assemble_light_pieces = 0._dp

        if (order == 1 .or. (order > 1 .and. coeffonly .eqv. .false.)) then
            assemble_light_pieces = full(0,0) + &
                as/4._dp/pi*( &
                    full(1,0) + full(1,1)*log(taucut) + full(1,2)*log(taucut)**2 )
        endif


        if (order >= 2) then
            assemble_light_pieces = assemble_light_pieces + (as/4._dp/pi)**2 * ( &
                full(2,0) + full(2,1)*log(taucut) + full(2,2)*log(taucut)**2 &
                + full(2,3)*log(taucut)**3 + full(2,4)*log(taucut)**4 )
        endif

    end function


    subroutine do_assemble_beam1(order, renscale, as, taucut, beama0, beamb0, beama1, &
                beama2, soft1, soft2, jet1, jet2, hard, assemble_beam1)
        use types 
        implicit none
        include 'kpart.f'
        include 'constants.f'
        include 'nf.f'
        include 'scet_const.f'

        integer, intent(in) :: order
        real(dp), intent(in) :: renscale, as
        real(dp), intent(in) :: taucut(:)
        real(dp), intent(in) :: beama0, beamb0
        real(dp), intent(in) :: beama1(-1:1)
        real(dp), intent(in) :: beama2(-1:3)
        real(dp), intent(in) :: soft1(-1:1), soft2(-1:3)
        real(dp), intent(in) :: jet1(-1:1), jet2(-1:3)
        real(dp), intent(in) :: hard(2)
        real(dp), intent(out) :: assemble_beam1(size(taucut(:)))

        real(dp) :: lcut(size(taucut(:))), full(-1:3)

        lcut = log(taucut/renscale)

        if (coeffonly) then
            assemble_beam1 = 0._dp
        else
            assemble_beam1 = beama0*beamb0
        endif

        full(:) = 0._dp

        ! only corrections on light line for now
        if ((order == 1) .or. ((order > 1) .and. (coeffonly .eqv. .false.))) then

            full(1)= beamb0 * ( beama1(1) + beama0*soft1(1) + 0.5_dp*beama0*jet1(1) )

            full(0)= beamb0 * ( beama1(0) + beama0*soft1(0) + 0.5_dp*beama0*jet1(0) )

            full(-1)= beamb0 * ( beama1(-1) + beama0*soft1(-1) + 0.5_dp*beama0*jet1(-1) + beama0*hard(1) )

            assemble_beam1 = assemble_beam1 + as/2._dp/pi*( &
                full(-1) + full(0)*lcut + 0.5_dp*full(1)*lcut**2 )
        endif

        if (order > 1) then
            full(3)= + beamb0 * ( beama2(3) + beama1(1)*soft1(1) + 0.5_dp &
               *beama1(1)*jet1(1) + 0.5_dp*beama0*soft1(1)*jet1(1) +  &
               beama0*soft2(3) + 1._dp/4._dp*beama0*jet2(3) )

           full(2)= + beamb0 * ( beama2(2) + 3._dp/2._dp*beama1(0)*soft1(1) &
                + 3._dp/4._dp*beama1(0)*jet1(1) + 3._dp/2._dp*beama1(1)* &
               soft1(0) + 3._dp/4._dp*beama1(1)*jet1(0) + 3._dp/4._dp*beama0 &
               *soft1(0)*jet1(1) + 3._dp/4._dp*beama0*soft1(1)*jet1(0) +  &
               beama0*soft2(2) + 1._dp/4._dp*beama0*jet2(2) )

           full(1)= + beamb0 * ( beama2(1) + beama1(-1)*soft1(1) + 0.5_dp &
              *beama1(-1)*jet1(1) + 2._dp*beama1(0)*soft1(0) + beama1(0)* &
               jet1(0) + beama1(1)*soft1(-1) - 2._dp*beama1(1)*soft1(1)* &
               zeta2 + 0.5_dp*beama1(1)*jet1(-1) - beama1(1)*jet1(1)* &
               zeta2 + beama1(1)*hard(1) + 0.5_dp*beama0*soft1(-1)* &
               jet1(1) + beama0*soft1(0)*jet1(0) + 0.5_dp*beama0*soft1( &
               1)*jet1(-1) - beama0*soft1(1)*jet1(1)*zeta2 + beama0*soft1(1) &
               *hard(1) + beama0*soft2(1) + 0.5_dp*beama0*jet1(1)*hard( &
               1) + 1._dp/4._dp*beama0*jet2(1) )

           full(0)= + beamb0 * ( beama2(0) + beama1(-1)*soft1(0) + 0.5_dp &
              *beama1(-1)*jet1(0) + beama1(0)*soft1(-1) - beama1(0)*soft1( &
               1)*zeta2 + 0.5_dp*beama1(0)*jet1(-1) - 0.5_dp* &
               beama1(0)*jet1(1)*zeta2 + beama1(0)*hard(1) - beama1(1)* &
               soft1(0)*zeta2 + 2._dp*beama1(1)*soft1(1)*zeta3 - 0.5_dp &
               *beama1(1)*jet1(0)*zeta2 + beama1(1)*jet1(1)*zeta3 + 0.5_dp &
              *beama0*soft1(-1)*jet1(0) + 0.5_dp*beama0*soft1(0)* &
               jet1(-1) - 0.5_dp*beama0*soft1(0)*jet1(1)*zeta2 + beama0 &
               *soft1(0)*hard(1) - 0.5_dp*beama0*soft1(1)*jet1(0)*zeta2 &
                + beama0*soft1(1)*jet1(1)*zeta3 + beama0*soft2(0) + 0.5_dp &
              *beama0*jet1(0)*hard(1) + 1._dp/4._dp*beama0*jet2(0) )

           full(-1)= + beamb0 * ( beama2(-1) + beama1(-1)*soft1(-1) + 0.5_dp &
               *beama1(-1)*jet1(-1) + beama1(-1)*hard(1) - beama1(0)* &
               soft1(0)*zeta2 + beama1(0)*soft1(1)*zeta3 - 0.5_dp* &
               beama1(0)*jet1(0)*zeta2 + 0.5_dp*beama1(0)*jet1(1)*zeta3 &
                + beama1(1)*soft1(0)*zeta3 - 1._dp/10._dp*beama1(1)*soft1(1) &
               *zeta2**2 + 0.5_dp*beama1(1)*jet1(0)*zeta3 - 1._dp/20._dp &
               *beama1(1)*jet1(1)*zeta2**2 + 0.5_dp*beama0*soft1(-1)* &
               jet1(-1) + beama0*soft1(-1)*hard(1) - 0.5_dp*beama0* &
               soft1(0)*jet1(0)*zeta2 + 0.5_dp*beama0*soft1(0)*jet1(1)* &
               zeta3 + 0.5_dp*beama0*soft1(1)*jet1(0)*zeta3 - 1._dp/20._dp &
              *beama0*soft1(1)*jet1(1)*zeta2**2 + beama0*soft2(-1) + 0.5_dp &
              *beama0*jet1(-1)*hard(1) + 1._dp/4._dp*beama0*jet2(-1) &
                + beama0*hard(2) )

           assemble_beam1=assemble_beam1+(as/2._dp/pi)**2*( &
            full(-1) + full(0)*lcut + full(1)*lcut**2/two + full(2)*lcut**3/three + full(3)*lcut**4/four)
        endif

    end subroutine


      function passed_taucut_light(pparton,scetreweight_local,taucut_in)
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
      include 'kpart.f'

      
      logical :: passed_taucut_light
      real(dp), intent(inout), optional :: scetreweight_local(:)
      real(dp), intent(in) :: pparton(mxpart,4)
      real(dp), intent(in), optional :: taucut_in

      integer:: j
      real(dp) :: tau,tauc,p67(4),p68(4),p78(4),p678(4)
      real(dp) :: dotvec,beamsign,tauallbeam,taubeamjet,taualljet,Qsq

      logical :: bin
      common/bin/bin

      if (present(scetreweight_local)) then
          scetreweight_local(:) = 0._dp
      endif

      passed_taucut_light = .false.

        if (present(taucut_in)) then
            tauc = taucut_in
        else
            tauc = taucut
        endif

      if (corr_on_beam == 1) then
          beamsign = 1._dp
          Qsq = -dotvec(pparton(2,:)+pparton(3,:)+pparton(4,:)+pparton(5,:), &
                        pparton(2,:)+pparton(3,:)+pparton(4,:)+pparton(5,:))
      else
          beamsign = -1._dp
          Qsq = -dotvec(pparton(1,:)+pparton(3,:)+pparton(4,:)+pparton(5,:), &
                        pparton(1,:)+pparton(3,:)+pparton(4,:)+pparton(5,:))
      endif

#define E(i) pparton(i,4)
#define pz(i) (beamsign*pparton(i,3))
      p67(:) = pparton(6,:) + pparton(7,:)

      if (origKpart == ksnlo .or. (present(taucut_in) .and. kpart == kvirt)) then
          tau = min( &
              ! extra emission clustered with the beam
              min(E(7)-pz(7), E(6)-pz(6)), &
              ! extra emission clustered with the jet
              E(6)+E(7) - sqrt((p67(1)**2 + p67(2)**2 + p67(3)**2)) )
          !taubeam = min(E(7)-pz(7), E(6)-pz(6))
          !taujet = E(6)+E(7) - sqrt((p67(1)**2 + p67(2)**2 + p67(3)**2))

          !tauallbeam = huge(1._dp)
          !taubeamjet = huge(1._dp)
          !taualljet = huge(1._dp)
      elseif (origKpart == knnlo .or. (present(taucut_in) .and. kpart == kreal)) then
          p68(:) = pparton(6,:) + pparton(8,:)
          p78(:) = pparton(7,:) + pparton(8,:)
          p678(:) = pparton(6,:) + pparton(7,:) + pparton(8,:)

          tauallbeam = min(E(6)+E(7)-pz(6)-pz(7), E(6)+E(8)-pz(6)-pz(8), E(7)+E(8)-pz(7)-pz(8))
          taubeamjet = min(E(6)+E(7)+E(8) - sqrt(p67(1)**2+p67(2)**2+p67(3)**2) - pz(8), &
              E(6)+E(7)+E(8) - sqrt(p68(1)**2+p68(2)**2+p68(3)**2) - pz(7), &
              E(6)+E(7)+E(8) - sqrt(p78(1)**2+p78(2)**2+p78(3)**2) - pz(6))
          taualljet = E(6)+E(7)+E(8) - sqrt(p678(1)**2+p678(2)**2+p678(3)**2)

          tau = min(tauallbeam,taubeamjet,taualljet)
      else
          error stop "unknown kpart in maketaucut_singletop"
      endif

      if (ieee_is_nan(tau)) then
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif

#if DYNAMICTAU == 1
      tau = tau / (sqrt(Qsq) / 80._dp)
#endif


      if (bin .and. doMultitaucut .and. present(scetreweight_local)) then
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
              return
          endif
      else
         if (tau < tauc) return

      endif

      passed_taucut_light = .true.

      end function passed_taucut_light
     
end module
