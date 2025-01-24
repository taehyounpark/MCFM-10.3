!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module singletop_interf_hxd
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use LHAPDF

      private

      public :: singletop_heavy_decay_vv
      ! DONE:
      ! * take singletop_light_decay_vv (from lxd) and modify:
      ! * born function corresponds heavy line virtual corrections
      ! * hard(1) function needs to be crafted carefully since two
      !   helicity structures contribute that need to be multiplied
      !   together with a (mt + ptslash) in-between.

      public :: singletop_heavy_decay_vv_tree
      ! DONE:
      ! * take singletop_heavy_decay_vv and remove additional heavy line
      !   corrections
      ! * born function is just singletop tree
      ! * hard(1) function is just decay virtual corrections
      ! * up to the interface this function should be just like the decay below
      !   cut; we can probably just take singletop_light_decay_vv_tree
      ! * can't take light_decay_vv_tree since we identify top quark decay
      !   associated scales with "islight"

      public :: singletop_heavy_decay_rr
      ! * checked double real amplitudes

      public :: singletop_heavy_decay_rr_gs
      ! * dipole routine can be taken from heavy line real emission _gs
      !   but tree process has to be replaced with decay real corrections
      ! * make sure that 7 and 8 are treated correctly in dipole routines

      public :: singletop_heavy_decay_vr
      ! DONE:
      ! * needs calculation (from scratch, or path together?)
      public :: singletop_heavy_decay_vr_z
      ! DONE:
      ! * taken from heavy line virtual _z

      public :: singletop_heavy_decay_rv
      ! DONE:
      ! * needs calculation similar to _vr

c     public :: singletop_heavy_decay_rv_tree ! for dipoles in_gs routine
      ! DONE:
      ! * should be the same as _vv_tree, just below cut with nothing on heavy
      !   line

      public :: singletop_heavy_decay_rv_gs
      ! DONE:
      ! * dipole routine, like _gs for heavy line corrections but
      !   with tree level as in _rv_tree

      public :: singletop_decay_real_hxd

      public :: extend_trans

      contains

      ! just like extend_trans from Singletop, but with additional restoration of decay radiation
      subroutine extend_trans(pold,p,ptrans,pext)
      implicit none
      include 'types.f'
c--- take vector
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: pold(mxpart,4),p(mxpart,4),ptrans(mxpart,4),
     & pext(mxpart,4),p3(4),p4(4),p5(4),pt(4),ptt(4),p7(4),
     & p3out(4),p4out(4),p5out(4),p7out(4)
      integer:: nu

      do nu=1,4
        pt(nu)=p(3,nu)
        ptt(nu)=ptrans(3,nu)
        p3(nu)=pold(3,nu)
        p4(nu)=pold(4,nu)
        p5(nu)=pold(5,nu)

        p7(nu)=pold(7,nu)
      enddo

      call boostx(p3,pt,ptt,p3out)
      call boostx(p4,pt,ptt,p4out)
      call boostx(p5,pt,ptt,p5out)

      call boostx(p7,pt,ptt,p7out)

      do nu=1,4
        pext(1,nu)=ptrans(1,nu)
        pext(2,nu)=ptrans(2,nu)
        pext(3,nu)=p3out(nu)
        pext(4,nu)=p4out(nu)
        pext(5,nu)=p5out(nu)
        pext(6,nu)=ptrans(4,nu)
        pext(7,nu)=p7out(nu)
        pext(8:,nu)=0._dp
      enddo

      return
      end subroutine extend_trans

      function st_heavy_decay_vv(ju,jb,jn,je,jc,jd, za,zb, musq_prod, musq_decay)
          use constants
          use singletop_interf_lxd, only: coefsdk, coefs_new
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          real(dp) :: st_heavy_decay_vv
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd
          real(dp), intent(in) :: musq_prod, musq_decay

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT126, propW34
          real(dp) :: propW16

          complex(dp) :: mtsq

          real(dp) :: cvProd, cvDecay, cv0_dummy
          real(dp) :: c1Decay
          complex(dp) :: c1Prod
          complex(dp) :: vamp, amp

          mtsq = mt**2 - im*mt*twidth

          ! calculate decay factor MSbar-renormalized with musq_decay
          call coefsdk(s(jn,je),mt**2,cvDecay,c1Decay,musq_decay,0._dp,0._dp)

          ! production factor with full epinv dependence
          call coefs_new(s(ju,jd),mt**2,cv0_dummy,cvProd,c1Prod,musq_prod,epinv,epinv2)

           propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
           propW16  = 1._dp / (s(ju,jd) - wmass**2)
           propT126 = 1._dp / (s(ju,jd) + s(ju,jb) + s(jd,jb) - mtsq)

           amp = za(jn,jc)*zb(jb,ju)*(za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju))

           vamp = c1Decay*((cvProd*za(jc,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je))/2._dp -
     &     (c1Prod*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(jc,je)*
     &        (za(ju,jc)*zb(jb,ju) + za(jc,jd)*zb(jd,jb)))/(4._dp*mt**2))
     &   - (c1Prod*cvDecay*za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(je,jb))/
     &   2._dp + cvDecay*cvProd*za(jn,jc)*zb(jb,ju)*
     &   (za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju))

         st_heavy_decay_vv = real(vamp*conjg(amp)) *
     &          abs(propW34 * propW16 * propT126)**2

      end function st_heavy_decay_vv

      subroutine singletop_heavy_decay_vv(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_scet_heavy_decay, only: softfun, jetfun
          use singletop_interf_lxh, only: virtqqb_heavy
          use singletop_interf_lxd, only: assemble_decay_pieces
          implicit none
          include 'constants.f'
          include 'masses.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'ewcouple.f'
          include 'zprods_com.f'
          include 'energy.f'
          include 'scheme.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          real(dp) :: xx(2)

          real(dp) :: beama0(-5:5), beamb0(-5:5)

          real(dp) :: jet(2,0:4), soft(2,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)

          real(dp) :: x

          real(dp) :: puremass
          real(dp) :: fac

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          x = puremass(p(3,:)+p(4,:))**2 / mt**2

          soft = softfun(1)
          jet = jetfun(1)

          call spinoru(6,p,za,zb)

          ! hard(0) are just virtual corrections on light line

          scheme = 'tH-V'
          fac = (as_light_beam1/2._dp/pi*cf) * aveqq*gw**8*xn**2
          hard_ub(0) = fac*virtqqb_heavy(1,2,3,4,5,6, renscale_beam2_isheavy_onheavy**2)
          hard_ubarb(0) = fac*virtqqb_heavy(6,2,3,4,5,1, renscale_beam2_isheavy_onheavy**2)

          fac = (as_light_beam2/2._dp/pi*cf) * aveqq*gw**8*xn**2
          hard_bu(0) = fac*virtqqb_heavy(2,1,3,4,5,6, renscale_beam1_isheavy_onheavy**2)
          hard_bubar(0) = fac*virtqqb_heavy(6,1,3,4,5,2, renscale_beam1_isheavy_onheavy**2)

          ! hard(1) virtual corrections in decay (msbar renormalized: epinv = 0
          ! and pi^2/6 * cf piece below)

          ! properly normalized one-loop pieces
          fac = aveqq*gw**8*xn**2
          fac = fac * (as_light_beam1/2._dp/pi*cf) * (as_heavy_beam2/2._dp/pi*cf)
          ! proper hard function normalization
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1) = fac*st_heavy_decay_vv(1,2,3,4,5,6,za,zb,
     &                      renscale_beam2_isheavy_onheavy**2,
     &                      renscale_beam2_islight_onheavy**2)
          hard_ubarb(1) = fac*st_heavy_decay_vv(6,2,3,4,5,1,za,zb,
     &                      renscale_beam2_isheavy_onheavy**2,
     &                      renscale_beam2_islight_onheavy**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac = aveqq*gw**8*xn**2
          fac = fac * (as_light_beam2/2._dp/pi*cf) * (as_heavy_beam1/2._dp/pi*cf)
          ! proper hard function normalization
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*st_heavy_decay_vv(2,1,3,4,5,6,za,zb,
     &                      renscale_beam1_isheavy_onheavy**2,
     &                      renscale_beam1_islight_onheavy**2)

          hard_bubar(1) = fac*st_heavy_decay_vv(6,1,3,4,5,2,za,zb,
     &                      renscale_beam1_isheavy_onheavy**2,
     &                      renscale_beam1_islight_onheavy**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)

          call fdist(ih1,xx(1),facscale_beam1_islight_onheavy,beama0,1)
          call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy,beamb0,2)

          msq([2,4],5) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)*[beama0(2),beama0(4)]*beamb0(5)

          msq([-3,-1],5) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)*[beama0(-3),beama0(-1)]*beamb0(5)

          call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy,beama0,1)
          call fdist(ih2,xx(2),facscale_beam2_islight_onheavy,beamb0,2)

          msq(5,[2,4]) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)*beama0(5)*[beamb0(2),beamb0(4)]

          msq(5,[-3,-1]) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)*beama0(5)*[beamb0(-3),beamb0(-1)]

      end subroutine singletop_heavy_decay_vv

      function st_heavy_decay_rv_MMMM_M(ju,jb,jn,je,jc,jd,jg, za,zb, musq)
          use constants
          use singletop_interf_lxd, only: coefsdk
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          real(dp) :: st_heavy_decay_rv_MMMM_M
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: musq

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT345, propW34
          real(dp) :: propW16

          complex(dp) :: mtsq

          real(dp) :: cv
          real(dp) :: c1
          complex(dp) :: treeamp, c1amp

          mtsq = mt**2 - im*mt*twidth

          ! calculate decay factor MSbar-renormalized with musq_decay
          call coefsdk(s(jn,je),mt**2,cv,c1,musq,0._dp,0._dp)

           propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
           propW16  = 1._dp / (s(ju,jd) - wmass**2)
           propT345 = 1._dp / (s(jn,je) + s(jn,jc) + s(je,jc) - mtsq)

          treeamp =  (za(jn,jc)*zb(jb,ju)*
     &     (za(jc,jg)*za(ju,jd)*zb(jb,ju)*zb(jc,je) -
     &       mt**2*za(jd,jg)*zb(je,jb) -
     &       za(jn,jg)*za(ju,jd)*zb(jb,ju)*zb(je,jn)))/
     &   ((mt**2 - za(ju,jb)*zb(jb,ju) - za(jb,jd)*zb(jd,jb) -
     &       za(ju,jd)*zb(jd,ju))*zb(jg,jb))

          c1amp = -(c1*za(jn,jc)*zb(jb,ju)*zb(jc,je)*
     &     (-(za(jc,jg)*za(ju,jd)*zb(jb,ju)) +
     &       za(jd,jg)*(za(je,jc)*zb(je,jb) + za(jn,jc)*zb(jn,jb))))/
     &   (2._dp*(mt**2 - za(ju,jb)*zb(jb,ju) - za(jb,jd)*zb(jd,jb) -
     &       za(ju,jd)*zb(jd,ju))*zb(jg,jb))


         st_heavy_decay_rv_MMMM_M = real((c1amp + cv*treeamp)*conjg(treeamp)) *
     &          abs(propW34 * propW16 * propT345)**2

      end function st_heavy_decay_rv_MMMM_M

      function st_heavy_decay_rv_MMMM_P(ju,jb,jn,je,jc,jd,jg, za,zb, musq)
          use constants
          use singletop_interf_lxd, only: coefsdk
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          real(dp) :: st_heavy_decay_rv_MMMM_P
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: musq

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT345, propW34
          real(dp) :: propW16

          complex(dp) :: mtsq

          real(dp) :: cv
          real(dp) :: c1
          complex(dp) :: treeamp, c1amp

          mtsq = mt**2 - im*mt*twidth

          ! calculate decay factor MSbar-renormalized with musq_decay
          call coefsdk(s(jn,je),mt**2,cv,c1,musq,0._dp,0._dp)

          propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
          propW16  = 1._dp / (s(ju,jd) - wmass**2)
          propT345 = 1._dp / (s(jn,je) + s(jn,jc) + s(je,jc) - mtsq)

          c1amp = -(c1*za(jc,jd)*za(jn,jc)*zb(jc,je)*
     &      (-(za(jb,jc)*zb(jb,ju)*
     &           (-mt**2 + za(ju,jb)*zb(jb,ju) + za(jb,jd)*zb(jd,jb) +
     &             za(ju,jd)*zb(jd,ju))) +
     &        za(jb,jg)*zb(jb,ju)*
     &         (za(je,jc)*zb(jg,je) + za(jn,jc)*zb(jg,jn)) +
     &        za(jc,jg)*(-mt**2 + za(ju,jb)*zb(jb,ju) +
     &           za(jb,jd)*zb(jd,jb) + za(ju,jd)*zb(jd,ju))*zb(jg,ju)))/
     &   (2._dp*za(jb,jg)*za(jc,jg)*
     &     (mt**2 - za(ju,jb)*zb(jb,ju) - za(jb,jd)*zb(jd,jb) -
     &       za(ju,jd)*zb(jd,ju)))

          treeamp = za(jn,jc)*(-((za(jb,jc)*zb(jb,ju)*
     &          (za(jc,jd)*zb(jc,je) - za(jn,jd)*zb(je,jn)))/
     &        (za(jb,jg)*za(jc,jg))) +
     &     (za(jb,jd)*za(jn,jc)*zb(jb,ju)*zb(je,jn)*zb(jg,jb))/
     &      (za(jc,jg)*(-mt**2 + za(ju,jb)*zb(jb,ju) +
     &          za(jb,jd)*zb(jd,jb) + za(ju,jd)*zb(jd,ju))) +
     &     (mt**2*za(jc,jd)*zb(jb,ju)*zb(jg,je))/
     &      (za(jc,jg)*(-mt**2 + za(ju,jb)*zb(jb,ju) +
     &          za(jb,jd)*zb(jd,jb) + za(ju,jd)*zb(jd,ju))) +
     &     (za(jc,jd)*zb(jc,je)*zb(jg,ju))/za(jb,jg) -
     &     (za(jn,jd)*zb(je,jn)*zb(jg,ju))/za(jb,jg) +
     &     (za(jn,jc)*za(ju,jd)*zb(jb,ju)*zb(je,jn)*zb(jg,ju))/
     &      (za(jc,jg)*(-mt**2 + za(ju,jb)*zb(jb,ju) +
     &          za(jb,jd)*zb(jd,jb) + za(ju,jd)*zb(jd,ju))))

         st_heavy_decay_rv_MMMM_P = real((c1amp + cv*treeamp)*conjg(treeamp)) *
     &          abs(propW34 * propW16 * propT345)**2

      end function st_heavy_decay_rv_MMMM_P

      subroutine singletop_heavy_decay_rv(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_nnlo_vars
          use singletop2_scet_heavy_decay, only: softfun, jetfun
          use singletop_interf_lxd, only: assemble_decay_pieces
          implicit none
          include 'constants.f'
          include 'masses.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'ewcouple.f'
          include 'zprods_com.f'
          include 'energy.f'
          include 'scheme.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib)

          real(dp) :: xx(2)

          real(dp) :: jet(2,0:4), soft(2,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)
          real(dp) :: hard_ug(0:1), hard_ubarg(0:1), hard_gu(0:1), hard_gubar(0:1)

          real(dp) :: x

          real(dp) :: puremass
          real(dp) :: fac

          real(dp) :: ubtdg_h

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          x = puremass(p(3,:)+p(4,:))**2 / mt**2

          soft = softfun(1)
          jet = jetfun(1)

          call spinoru(7,p,za,zb)

          ! hard(0) are just real emission corrections in heavy production

          scheme = 'tH-V'

          fac=2._dp*(4*pi*as_heavy_beam2)*cf*gw**8*xn**2
          hard_ub(0) = aveqq*fac*(ubtdg_h(1,2,3,4,5,6,7,p))
          hard_ubarb(0) = aveqq*fac*(ubtdg_h(6,2,3,4,5,1,7,p))
          hard_ug(0) =  aveqg*fac*(ubtdg_h(1,7,3,4,5,6,2,p))
          hard_ubarg(0) = aveqg*fac*(ubtdg_h(6,7,3,4,5,1,2,p))

          fac=2._dp*(4*pi*as_heavy_beam1)*cf*gw**8*xn**2
          hard_bu(0) = aveqq*fac*(ubtdg_h(2,1,3,4,5,6,7,p))
          hard_bubar(0) = aveqq*fac*(ubtdg_h(6,1,3,4,5,2,7,p))
          hard_gu(0) = aveqg*fac*(ubtdg_h(2,7,3,4,5,6,1,p))
          hard_gubar(0) = aveqg*fac*(ubtdg_h(6,7,3,4,5,2,1,p))

          ! hard(1) virtual corrections in decay (msbar renormalized: epinv = 0
          ! and pi^2/6 * cf piece below)

      fac=2._dp*(4*pi*as_heavy_beam2)*cf*gw**8*xn**2
      fac = fac * (as_light_beam2/2._dp/pi * cf) ! decay treated as "light"

      fac = fac / (as_light_beam2/4._dp/pi) ! for proper hard function normalization

      hard_ub(1) = aveqq*fac*(
     &      st_heavy_decay_rv_MMMM_P(1,2,3,4,5,6,7,za,zb,renscale_beam2_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(1,2,3,4,5,6,7,za,zb,renscale_beam2_islight_onheavy**2))

      hard_ubarb(1) = aveqq*fac*(
     &      st_heavy_decay_rv_MMMM_P(6,2,3,4,5,1,7,za,zb,renscale_beam2_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(6,2,3,4,5,1,7,za,zb,renscale_beam2_islight_onheavy**2))

      hard_ug(1) = aveqg*fac*(
     &      st_heavy_decay_rv_MMMM_P(1,7,3,4,5,6,2,za,zb,renscale_beam2_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(1,7,3,4,5,6,2,za,zb,renscale_beam2_islight_onheavy**2))

      hard_ubarg(1) = aveqg*fac*(
     &      st_heavy_decay_rv_MMMM_P(6,7,3,4,5,1,2,za,zb,renscale_beam2_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(6,7,3,4,5,1,2,za,zb,renscale_beam2_islight_onheavy**2))

      hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
      hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)
      hard_ug(1) = hard_ug(1) + pi**2/6._dp * cf * hard_ug(0)
      hard_ubarg(1) = hard_ubarg(1) + pi**2/6._dp * cf * hard_ubarg(0)


      fac=2._dp*(4*pi*as_heavy_beam1)*cf*gw**8*xn**2
      fac = fac * (as_light_beam1/2._dp/pi * cf)

      fac = fac / (as_light_beam1/4._dp/pi) ! for proper hard function normalization

      hard_bu(1) = aveqq*fac*(
     &      st_heavy_decay_rv_MMMM_P(2,1,3,4,5,6,7,za,zb,renscale_beam1_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(2,1,3,4,5,6,7,za,zb,renscale_beam1_islight_onheavy**2))

      hard_bubar(1) = aveqq*fac*(
     &      st_heavy_decay_rv_MMMM_P(6,1,3,4,5,2,7,za,zb,renscale_beam1_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(6,1,3,4,5,2,7,za,zb,renscale_beam1_islight_onheavy**2))

      hard_gu(1) = aveqg*fac*(
     &      st_heavy_decay_rv_MMMM_P(2,7,3,4,5,6,1,za,zb,renscale_beam1_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(2,7,3,4,5,6,1,za,zb,renscale_beam1_islight_onheavy**2))

      hard_gubar(1) = aveqg*fac*(
     &      st_heavy_decay_rv_MMMM_P(6,7,3,4,5,2,1,za,zb,renscale_beam1_islight_onheavy**2) +
     &      st_heavy_decay_rv_MMMM_M(6,7,3,4,5,2,1,za,zb,renscale_beam1_islight_onheavy**2))

      hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
      hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)
      hard_gu(1) = hard_gu(1) + pi**2/6._dp * cf * hard_gu(0)
      hard_gubar(1) = hard_gubar(1) + pi**2/6._dp * cf * hard_gubar(0)


          msq([2,4],5, 1) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onheavy, as_light_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)

          msq([-3,-1],5, 1) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onheavy, as_light_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)

          msq([2,4],0, 3) = hard_ug(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onheavy, as_light_beam2, taucut,
     &        hard_ug/hard_ug(0), jet, soft)

          msq([-1,-3],0, 3) = hard_ubarg(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onheavy, as_light_beam2, taucut,
     &        hard_ubarg/hard_ubarg(0), jet, soft)



          msq(5,[2,4], 1) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onheavy, as_light_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)

          msq(5,[-3,-1], 1) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onheavy, as_light_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)

          msq(0,[2,4], 3) = hard_gu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onheavy, as_light_beam1, taucut,
     &        hard_gu/hard_gu(0), jet, soft)

          msq(0,[-3,-1], 3) = hard_gubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onheavy, as_light_beam1, taucut,
     &        hard_gubar/hard_gubar(0), jet, soft)

      end subroutine singletop_heavy_decay_rv

      function st_heavy_decay_vr_MMMM_M(ju,jb,jn,je,jc,jd,jg, za,zb, musq)
          use constants
          use singletop_interf_lxd, only: coefs_new
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          real(dp) :: st_heavy_decay_vr_MMMM_M
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: musq

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT126, propW34
          real(dp) :: propW16

          complex(dp) :: mtsq

          real(dp) :: cv, cv0_dummy
          complex(dp) :: c1
          complex(dp) :: treeamp, c1amp

          mtsq = mt**2 - im*mt*twidth

          call coefs_new(s(ju,jd),mt**2,cv0_dummy,cv,c1,musq,epinv,epinv2)

           propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
           propW16  = 1._dp / (s(ju,jd) - wmass**2)
           propT126 = 1._dp / (s(ju,jb) + s(ju,jd) + s(jd,jb) - mtsq)

           !tested to reproduce real emission without virt for
           !cv = 1._dp
           !c1 = 0._dp

          c1amp = (c1*za(jb,jd)*zb(jb,ju)*
     &     (za(jn,jg)*zb(je,jb)*
     &        (mt**2 - za(je,jc)*zb(jc,je) - za(jn,je)*zb(je,jn)) +
     &       za(jn,jc)*(za(ju,jg)*zb(jb,ju)*zb(jc,je) -
     &          za(jc,jg)*zb(jc,jb)*zb(jc,je) -
     &          za(jd,jg)*zb(jc,je)*zb(jd,jb) -
     &          za(jn,jg)*zb(jc,jn)*zb(je,jb) +
     &          za(jn,jg)*zb(jc,jb)*zb(je,jn))))/
     &   (2._dp*(mt**2 - za(je,jc)*zb(jc,je) - za(jn,jc)*zb(jc,jn) -
     &       za(jn,je)*zb(je,jn))*zb(jg,jc))

           treeamp =
     &  (zb(jb,ju)*(mt**2*za(jd,jg)*za(jn,jc)*zb(jc,je) +
     &       za(jb,jd)*(za(jc,jg)*za(jn,jc)*zb(jc,jb)*zb(jc,je) +
     &          za(jn,jg)*(-(za(jn,jc)*zb(jc,jb)*zb(je,jn)) +
     &             zb(je,jb)*
     &              (-mt**2 + za(je,jc)*zb(jc,je) +
     &                za(jn,jc)*zb(jc,jn) + za(jn,je)*zb(je,jn)))) +
     &       za(ju,jd)*(za(jc,jg)*za(jn,jc)*zb(jc,je)*zb(jc,ju) +
     &          za(jn,jg)*((-mt**2 + za(je,jc)*zb(jc,je) +
     &                za(jn,je)*zb(je,jn))*zb(je,ju) +
     &             za(jn,jc)*
     &              (-(zb(jc,ju)*zb(je,jn)) + zb(jc,jn)*zb(je,ju))))))/
     &   ((mt**2 - za(je,jc)*zb(jc,je) - za(jn,jc)*zb(jc,jn) -
     &       za(jn,je)*zb(je,jn))*zb(jg,jc))

         st_heavy_decay_vr_MMMM_M = real((c1amp+cv*treeamp)*conjg(treeamp)) *
     &          abs(propW34 * propW16 * propT126)**2

      end function st_heavy_decay_vr_MMMM_M

      function st_heavy_decay_vr_MMMM_P(ju,jb,jn,je,jc,jd,jg, za,zb, musq)
          use constants
          use singletop_interf_lxd, only: coefs_new
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          real(dp) :: st_heavy_decay_vr_MMMM_P
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: musq

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT126, propW34
          real(dp) :: propW16

          complex(dp) :: mtsq

          real(dp) :: cv, cv0_dummy
          complex(dp) :: c1
          complex(dp) :: treeamp, c1amp

          mtsq = mt**2 - im*mt*twidth

          call coefs_new(s(ju,jd),mt**2,cv0_dummy,cv,c1,musq,epinv,epinv2)

           !tested to reproduce real emission without virt for
           !cv = 1._dp
           !c1 = 0._dp

           propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
           propW16  = 1._dp / (s(ju,jd) - wmass**2)
           propT126 = 1._dp / (s(ju,jb) + s(ju,jd) + s(jd,jb) - mtsq)

          c1amp = (c1*za(jb,jd)*za(jn,jc)*zb(jb,ju)*
     &     (za(jn,jc)*zb(je,jn)*zb(jg,jb) +
     &       (za(ju,jc)*zb(jb,ju) + za(jc,jd)*zb(jd,jb))*zb(jg,je)))/
     &   (2._dp*za(jc,jg)*(mt**2 - za(je,jc)*zb(jc,je) -
     &       za(jn,jc)*zb(jc,jn) - za(jn,je)*zb(je,jn)))

          treeamp = -(za(jn,jc)*zb(jb,ju)*
     &     (za(jb,jd)*za(jn,jc)*zb(je,jn)*zb(jg,jb) +
     &       mt**2*za(jc,jd)*zb(jg,je) +
     &       za(jn,jc)*za(ju,jd)*zb(je,jn)*zb(jg,ju)))/
     &   (za(jc,jg)*(mt**2 - za(je,jc)*zb(jc,je) -
     &       za(jn,jc)*zb(jc,jn) - za(jn,je)*zb(je,jn)))


         st_heavy_decay_vr_MMMM_P = real((c1amp+cv*treeamp)*conjg(treeamp)) *
     &          abs(propW34 * propW16 * propT126)**2

      end function st_heavy_decay_vr_MMMM_P

      ! decay scale is treated as "islight", this is a special version for that
      subroutine singletop_decay_real_hxd(p,msq)
          use singletop2_nnlo_vars
          use singletop2_scet_heavy_decay, only: qqbtbbargd
      implicit none
      include 'types.f'

c       Matrix element for t-bbar production
c        b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
c       averaged(summed) over initial(final) colours and spins
c--N  B average over spins only -- colour factors cancel
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

      fac=2._dp*(4*pi*as_light_beam2)*cf*aveqq*gw**8*xn**2
      ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
      ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
      msq([2,4],5, 1,2) = ub
      msq([-1,-3],5, 1,2) = ubarb

      fac=2._dp*(4*pi*as_light_beam1)*cf*aveqq*gw**8*xn**2
      bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
      bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
      msq(5,[2,4], 1,1) = bu
      msq(5,[-1,-3], 1,1) = bubar

      end subroutine singletop_decay_real_hxd

      subroutine singletop_heavy_decay_vr(p,msq)
          use singletop2_nnlo_vars
          use singletop2_scet_heavy_prod, only: virtqqb_heavy
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
      include 'scheme.f'

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp):: fac
      integer:: ib


      call spinoru(7,p,za,zb)
      ib=5*nwz

      scheme = 'tH-V'

      msq = 0._dp

      fac = 2._dp*(4*pi*as_light_beam2)*cf*aveqq*gw**8*xn**2 ! decay real
      fac = fac * (as_heavy_beam2/2._dp/pi*cf) ! virt

      msq([2,4],5, 1,2) = fac*(st_heavy_decay_vr_MMMM_P(1,2,3,4,5,6,7,za,zb, renscale_beam2_isheavy_onheavy**2) +
     &     st_heavy_decay_vr_MMMM_M(1,2,3,4,5,6,7,za,zb, renscale_beam2_isheavy_onheavy**2))

      msq([-1,-3],5, 1,2) = fac*(st_heavy_decay_vr_MMMM_P(6,2,3,4,5,1,7,za,zb, renscale_beam2_isheavy_onheavy**2) +
     &     st_heavy_decay_vr_MMMM_M(6,2,3,4,5,1,7,za,zb, renscale_beam2_isheavy_onheavy**2))

      fac = 2._dp*(4*pi*as_light_beam1)*cf*aveqq*gw**8*xn**2 ! decay real
      fac = fac * (as_heavy_beam1/2._dp/pi*cf) ! virt

      msq(5,[2,4], 1,1) = fac*(st_heavy_decay_vr_MMMM_P(2,1,3,4,5,6,7,za,zb, renscale_beam1_isheavy_onheavy**2) +
     &     st_heavy_decay_vr_MMMM_M(2,1,3,4,5,6,7,za,zb, renscale_beam1_isheavy_onheavy**2))

      msq(5,[-1,-3], 1,1) = fac*(st_heavy_decay_vr_MMMM_P(6,1,3,4,5,2,7,za,zb, renscale_beam1_isheavy_onheavy**2) +
     &     st_heavy_decay_vr_MMMM_M(6,1,3,4,5,2,7,za,zb, renscale_beam1_isheavy_onheavy**2))

      end subroutine singletop_heavy_decay_vr

      subroutine singletop_heavy_decay_vr_z(p,z)
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
        real(dp) :: dot
        real(dp) :: if_mqq, fi_mqq, if_mqg

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
        do is=1,3
               B2(b,b,q,is) = as_heavy_beam2/2/pi * cf*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
               !B2(q,g,q,is) = as_heavy_beam2/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam2_isheavy_onheavy**2),is)

               ! massive spec needed for DDIS
               B2(q,g,q,is) = as_heavy_beam2/2/pi * tr * if_mqg(z, xl25, mbar25, is)


               B1(b,b,q,is) = as_heavy_beam1/2/pi * cf*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
               !B1(q,g,q,is) = as_heavy_beam1/2/pi * tr * ii_qg(z, log(2*dot(p,1,2)/renscale_beam1_isheavy_onheavy**2),is)

               ! massive spec for DDIS
               B1(q,g,q,is) = as_heavy_beam1/2/pi * tr * if_mqg(z, xl15, mbar15, is)
        enddo

      end subroutine

      ! like light_decay_vv_tree, but with scale setting adjusted
      ! we handle decay corrections as "light"
      subroutine singletop_heavy_decay_vv_tree(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_scet_heavy_decay, only: softfun, jetfun
          use singletop_interf_lxh, only: qqbtbbar
          use singletop_interf_lxd, only: assemble_decay_pieces, virtqqbdk
          implicit none
          include 'constants.f'
          include 'masses.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'ewcouple.f'
          include 'zprods_com.f'
          include 'energy.f'
          include 'scheme.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          real(dp) :: xx(2)

          real(dp) :: jet(2,0:4), soft(2,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)

          real(dp) :: x

          real(dp) :: puremass
          real(dp) :: fac

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          x = puremass(p(3,:)+p(4,:))**2 / mt**2

          soft = softfun(1)
          jet = jetfun(1)

          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2

          scheme = 'tH-V'


          ! hard(0) just Born amplitudes

          hard_ub(0) = fac*qqbtbbar(1,2,3,4,5,6)
          hard_ubarb(0) = fac*qqbtbbar(6,2,3,4,5,1)

          hard_bu(0) = fac*qqbtbbar(2,1,3,4,5,6)
          hard_bubar(0) =fac*qqbtbbar(6,1,3,4,5,2)

          ! hard(1) corresponds to virtual corrections in decay
          ! but MSbar renormalized epinv = 0 set in virtqqbdk, pi^2/6*cf added
          ! below

          ! properly normalized one-loop pieces
          fac=(as_light_beam2/2._dp/pi*cf)*aveqq*gw**8*xn**2
          ! proper hard function normalization
          fac = fac / (as_light_beam2/4._dp/pi)
          hard_ub(1) = fac*virtqqbdk(p, 1,2,3,4,5,6,renscale_beam2_islight_onlight**2)
          hard_ubarb(1) = fac*virtqqbdk(p, 6,2,3,4,5,1,renscale_beam2_islight_onlight**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac=(as_light_beam1/2._dp/pi*cf)*aveqq*gw**8*xn**2
          ! proper hard function normalization
          fac = fac / (as_light_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqbdk(p, 2,1,3,4,5,6,renscale_beam1_islight_onlight**2)
          hard_bubar(1) = fac*virtqqbdk(p, 6,1,3,4,5,2,renscale_beam1_islight_onlight**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)


          msq([2,4],5) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onlight, as_light_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)

          msq([-3,-1],5) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_islight_onlight, as_light_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)

          msq(5,[2,4]) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onlight, as_light_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)

          msq(5,[-3,-1]) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_islight_onlight, as_light_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)

      end subroutine singletop_heavy_decay_vv_tree

      subroutine singletop_heavy_decay_rv_gs(p,ndmx,msq)
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
          real(dp), intent(inout) :: msq(ndmx,-nf:nf,-nf:nf,max_bcontrib)

          real(dp) :: dummyv(-nf:nf,-nf:nf),dsubv
          real(dp) :: sub(4), msqx(-nf:nf,-nf:nf)

          integer :: noglue(4) = [-1,-3,2,4]

          external donothing_gvec

          ndmax = 6

          msq = 0._dp
          sub = 0._dp

        qqproc = .true.

        corr_islight = .false.
        corr_beam1 = .false.
        corr_on_beam = 2

        ! for dips_mass:
        ! 3 -> top (3+4+5)
        ! 4 -> 6
        ! 5 -> 7

        call dips_mass(1,p,2,5,3,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(1,noglue,5, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

        call dips_mass(2,p,3,5,2,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(2,noglue,5, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

        ! binned as final state b~ contrib
        call dips_mass(5,p,2,5,3,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(5,noglue,0, 3) = 2._dp*tr*sub(qg)*msqx(noglue,5)

        corr_islight = .false.
        corr_beam1 = .true.
        corr_on_beam = 1

        call dips_mass(3,p,1,5,3,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(3,5,noglue, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

        call dips_mass(4,p,3,5,1,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(4,5,noglue, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

        ! binned as final state b~ contrib
        call dips_mass(6,p,1,5,3,sub,dsubv,msqx,dummyv,singletop_heavy_decay_vv_tree,donothing_gvec)
        msq(6,0,noglue, 3) = 2._dp*tr*sub(qg)*msqx(5,noglue)

      end subroutine singletop_heavy_decay_rv_gs


      ! copy of routine without _simple
      ! removed b_contrib and corr_on_beam here for use in
      ! rr_gs
      subroutine singletop_decay_real_hxd_simple(p,msq)
          use singletop2_nnlo_vars
          use singletop2_scet_heavy_decay, only: qqbtbbargd
      implicit none
      include 'types.f'

c       Matrix element for t-bbar production
c        b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
c       averaged(summed) over initial(final) colours and spins
c--N  B average over spins only -- colour factors cancel
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,ub,bu,bubar,ubarb
      integer:: ib

      call spinoru(7,p,za,zb)
      ib=5*nwz

      msq(:,:) = 0._dp

      fac=2._dp*(4*pi*as_light_beam2)*cf*aveqq*gw**8*xn**2
      ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
      ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
      msq([2,4],5) = ub
      msq([-1,-3],5) = ubarb

      fac=2._dp*(4*pi*as_light_beam1)*cf*aveqq*gw**8*xn**2
      bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
      bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
      msq(5,[2,4]) = bu
      msq(5,[-1,-3]) = bubar

      end subroutine singletop_decay_real_hxd_simple

#define WITH_BBAR 3
      subroutine singletop_heavy_decay_rr_gs(p,ndmx,msq)
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

          real(dp) :: dummyv(-nf:nf,-nf:nf), dsubv
          real(dp) :: sub(4), msqx(-nf:nf,-nf:nf)

          integer :: noglue(4) = [-1,-3,2,4]

          external donothing_gvec

          ndmax = 6

          msq = 0._dp
          sub = 0._dp

        qqproc = .true.

        corr_islight = .false.
        corr_beam1 = .false.
        corr_on_beam = 2

        ! for dips_mass:
        ! 3 -> top (3+4+5+7)
        ! 4 -> 6
        ! 5 -> 8

        call dips_mass(1,p,2,5,3,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(1,noglue,5, 1, 2) = 2._dp*cf*sub(qq)*msqx(noglue,5)

        call dips_mass(2,p,3,5,2,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(2,noglue,5, 1, 2) = 2._dp*cf*sub(qq)*msqx(noglue,5)

        ! binned as final state b~ contrib
        call dips_mass(5,p,2,5,3,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(5,noglue,0, WITH_BBAR, 2) = 2._dp*tr*sub(qg)*msqx(noglue,5)

        corr_islight = .false.
        corr_beam1 = .true.
        corr_on_beam = 1

        call dips_mass(3,p,1,5,3,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(3,5,noglue, 1, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

        call dips_mass(4,p,3,5,1,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(4,5,noglue, 1, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

        ! binned as final state b~ contrib
        call dips_mass(6,p,1,5,3,sub,dsubv,msqx,dummyv,singletop_decay_real_hxd_simple,donothing_gvec)
        msq(6,0,noglue, WITH_BBAR, 1) = 2._dp*tr*sub(qg)*msqx(5,noglue)

      end subroutine singletop_heavy_decay_rr_gs

      subroutine singletop_heavy_decay_rr(p,msqall)
        use singletop2_scale_m
        use singletop2_nnlo_vars
        implicit none
        include 'types.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'nwz.f'
        include 'sprods_com.f'
        include 'zprods_com.f'
        include 'constants.f'
        include 'ewcouple.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

        real(dp) :: fac,facqg,facqq,tmp(-nf:nf,max_bcontrib)

        integer :: iperm,i1,i2

        call spinoru(8,p,za,zb)

        msqall(:,:,:,:) = 0._dp

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1
        do iperm=1,maxbeams
            corr_on_beam = beams_enabled(iperm)

            tmp(:,:)=0._dp

            if (corr_on_beam == 2) then
              i1=1
              i2=2
              fac=(fourpi*as_light_beam2)*(fourpi*as_heavy_beam2)
            else
              i1=2
              i2=1
              fac=(fourpi*as_light_beam1)*(fourpi*as_heavy_beam1)
            endif
            facqq=fac*aveqq
            facqg=fac*aveqg

c Note: as set up here, gluon 7 is radiated from heavy line and gluon 8 in decay
c Update Tobias: reverted to 7 radiated from decay and 8 from heavy line

            if (iand(partons_enabled, quarkChannel) > 0) then
                tmp([2,4], 1) = facqq*msqheavyxdecay(i1,i2,3,4,5,6,7,8)   ! u b -> t d g g
                tmp([-1,-3], 1) = facqq*msqheavyxdecay(6,i2,3,4,5,i1,7,8) ! d~ b -> t u~ g g
            endif

            if (iand(partons_enabled, gluonChannel) > 0) then
                tmp([2,4], WITH_BBAR)=facqg*msqheavyxdecay(i1,8,3,4,5,6,7,i2)         ! u b -> t d b~ g
                tmp([-1,-3], WITH_BBAR)=facqg*msqheavyxdecay(6,8,3,4,5,i1,7,i2)       ! d~ b -> t u~ b~ g
            endif

c update     msq array,
            if (corr_on_beam == 1) then
                msqall(5,:, 1, 1) = tmp(:,1)
                msqall(0,:, WITH_BBAR, 1) = tmp(:, WITH_BBAR)
            else
                msqall(:,5, 1, 2) = tmp(:,1)
                msqall(:,0, WITH_BBAR, 2) = tmp(:, WITH_BBAR)
            endif

        enddo

      end subroutine singletop_heavy_decay_rr


      function msqheavyxdecay(p1,p2,p3,p4,p5,p6,p7,p8)

c Matrix element squared for (four) diagrams of the form:

c     p1 ---->------------->------ p6
c                    $
c                    $          o p7
c                    $         o
c                    $   t    o
c     p2 ---->---------||-->------ p5     (on-shell top, p3457^2 = mt^2)
c               o          $
c                o          $ -->-- p3
c                 o          |
c                   p8       ^
c                            | p4

c with p7 and p8 attached to the lines shown, and no overall factor of 1/2 for Bose symmetry

c All overall factors are included except for averaging over the initial state
c and a strong coupling factor of gs^4

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: msqheavyxdecay,s16,s34,s345,s34578,prop34,prop16,prop345,prop34578
      complex(dp):: app,apm,amp,amm,zab2,zab3

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)

      s34=s(p3,p4)
      prop34=(s34-wmass**2)**2
      if (s34 > 0._dp) then
        prop34=prop34+(wmass*wwidth)**2
      endif

      msqheavyxdecay=0._dp

      prop16=(s(p1,p6)-wmass**2)**2
      if (s16 > 0._dp) then
        prop16=prop16+(wmass*wwidth)**2
      endif

c No width necessary for these top propagators since we assume s3457 = mt^2
      s34578=s(p1,p2)+s(p1,p6)+s(p2,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      prop34578=s34578-mt**2
      prop345=s345-mt**2

      app =  + prop345**(-1)*prop34578**(-1) * (
     &     + 1/(za(p5,p7))/(za(p6,p8))*za(p3,p5)*za(p5,p6)*zb(p1,p2)*
     &    zb(p4,p7)*zab2(p6,p1,p2,p8)*mt**2
     &     + 1/(za(p5,p7))/(za(p6,p8))*za(p3,p5)*zb(p1,p2)*zab2(p5,p3,
     &    p5,p4)*zab2(p6,p1,p2,p8)*zab3(p6,p3,p4,p5,p7)
     &     )
      app = app + prop345**(-1) * (
     &     + 1/(za(p2,p8))/(za(p5,p7))/(za(p6,p8))*za(p3,p5)*za(p5,p6)*
     &    zb(p4,p7)*zab2(p6,p2,p8,p1)*mt**2
     &     + 1/(za(p2,p8))/(za(p5,p7))/(za(p6,p8))*za(p3,p5)*zab2(p5,p3
     &    ,p5,p4)*zab2(p6,p2,p8,p1)*zab3(p6,p3,p4,p5,p7)
     &     )

      apm =  + prop345**(-1)*prop34578**(-1) * (
     &     + 1/(za(p5,p7))/(zb(p2,p8))*za(p1,p6)*za(p3,p5)*za(p5,p8)*
     &    zb(p1,p2)**2*zb(p4,p7)*mt**2
     &     + 1/(za(p5,p7))/(zb(p2,p8))*za(p1,p6)*za(p3,p5)*zb(p1,p2)**2
     &    *zab2(p5,p3,p5,p4)*zab3(p8,p3,p4,p5,p7)
     &     + 1/(za(p5,p7))/(zb(p2,p8))*za(p3,p5)*za(p6,p8)*zb(p1,p2)*
     &    zb(p2,p7)*zab2(p5,p3,p5,p4)*mt**2
     &     - 1/(za(p5,p7))/(zb(p2,p8))*za(p3,p5)*za(p6,p8)*zb(p1,p2)*
     &    zb(p4,p7)*zab3(p5,p3,p4,p7,p2)*mt**2
     &     )

      amp =  + prop345**(-1)*prop34578**(-1) * (
     &     - 1/(za(p6,p8))/(zb(p4,p7))*za(p3,p5)*zb(p1,p2)*zab2(p6,p1,
     &    p2,p8)*zab2(p7,p3,p5,p4)*zab3(p6,p3,p5,p7,p4)
     &     )
      amp = amp + prop345**(-1) * (
     &     - 1/(za(p2,p8))/(za(p6,p8))/(zb(p4,p7))*za(p3,p5)*zab2(p6,p2
     &    ,p8,p1)*zab2(p7,p3,p5,p4)*zab3(p6,p3,p5,p7,p4)
     &     )
      amp = amp + prop34578**(-1) * (
     &     - 1/(za(p6,p8))/(zb(p4,p7))/(zb(p5,p7))*zb(p1,p2)*zab2(p3,p5
     &    ,p7,p4)*zab2(p6,p1,p2,p8)*zab3(p6,p3,p5,p7,p4)
     &     )
      amp = amp - 1/(za(p2,p8))/(za(p6,p8))/(zb(p4,p7))/(zb(p5,p7))*
     & zab2(p3,p5,p7,p4)*zab2(p6,p2,p8,p1)*zab3(p6,p3,p5,p7,p4)

      amm =  + prop345**(-1)*prop34578**(-1) * (
     &     - 1/(zb(p2,p8))/(zb(p4,p7))*za(p1,p6)*za(p3,p5)*zb(p1,p2)**2
     &    *zab2(p7,p3,p5,p4)*zab3(p8,p3,p5,p7,p4)
     &     - 1/(zb(p2,p8))/(zb(p4,p7))*za(p3,p5)*za(p6,p8)*zb(p1,p2)*
     &    zb(p2,p4)*zab2(p7,p3,p5,p4)*mt**2
     &     )
      amm = amm + prop34578**(-1) * (
     &     - 1/(zb(p2,p8))/(zb(p4,p7))/(zb(p5,p7))*za(p1,p6)*zb(p1,p2)
     &    **2*zab2(p3,p5,p7,p4)*zab3(p8,p3,p5,p7,p4)
     &     - 1/(zb(p2,p8))/(zb(p4,p7))/(zb(p5,p7))*za(p6,p8)*zb(p1,p2)*
     &    zb(p2,p4)*zab2(p3,p5,p7,p4)*mt**2
     &     )
      msqheavyxdecay=msqheavyxdecay
     & +(abs(app)**2+abs(apm)**2+abs(amp)**2+abs(amm)**2)/(prop34*prop16)

c Overall factors
      msqheavyxdecay=msqheavyxdecay*V**2*gwsq**4/(mt*twidth)**2

      return
      end

      end module
