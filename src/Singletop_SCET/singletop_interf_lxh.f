!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module singletop_interf_lxh
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use LHAPDF

      private

      public :: singletop_jet_light_heavy_rr_all
      public :: singletop_jet_light_heavy_rr_gs_all

      public :: singletop_jet_light_heavy_vr_all

      public :: singletop_jet_light_heavy_vr_z
      public :: qqb_tbb_g_heavy_all_swap

      public :: singletop_jet_light_heavy_rv
      public :: singletop_jet_light_heavy_rv_tree ! tree piece for _gs dipoles
      public :: singletop_jet_light_heavy_rv_gs


      public :: singletop_jet_light_heavy_vv
      ! for splitting pieces
      public :: singletop_jet_light_heavy_vv_tree

      public :: passed_taucut_lxh

      real(dp), public, save :: z1int, z2int
!$omp threadprivate(z1int, z2int)


      ! additional routines used in other places

      public :: virtqqb_light, virtqqb_heavy, qqbtbbar

      contains

      function ubtdg_l_virt(ju,jb,jn,je,jc,jd,jg,p,musq)
          use types
      implicit none
      real(dp):: ubtdg_l_virt

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'

      integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(in) :: musq

      complex(dp):: amp1(2), amp0(2)

      real(dp) :: eta
      real(dp) :: massvec,s12

      scheme = 'tH-V'

      if (scheme =='dred') then
c------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
c------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

      s12 = massvec(-p(ju,:)-p(jd,:)-p(jg,:))

      amp1(1) = streal_lightResonant_MMMM_P_virt(ju,jb,jn,je,jc,jd,jg, za,zb,s12,musq)
      amp1(2) = streal_lightResonant_MMMM_M_virt(ju,jb,jn,je,jc,jd,jg, za,zb,s12,musq)

      amp0(1) = streal_lightResonant_MMMM_P_tree(ju,jb,jn,je,jc,jd,jg, za,zb)
      amp0(2) = streal_lightResonant_MMMM_M_tree(ju,jb,jn,je,jc,jd,jg, za,zb)

      ubtdg_l_virt = sum(real(amp1*conjg(amp0)))

      end

      ! modified to include vertex correction on light line
      function ubtdg_h_virt(ju,jb,jn,je,jc,jd,jg,p,musq)
          use types
      implicit none
      real(dp):: ubtdg_h_virt
c     Matrix element squared for single top production with gluon
c     radiation in production (radiation from final line)
c      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc))+d(jd)+g(jg)
c     masses of b quarks c.c=b.b=0

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'

      integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(in) :: musq

      integer :: nu
      real(dp):: pt(4),ptDpt
      real(dp):: sne,sdu,prop,twoptg
      complex(dp):: ampf(2)

      real(dp) :: cv0
      real(dp) :: eta
      complex(dp) :: xl12, lnrat

      scheme = 'tH-V'

      if (scheme =='dred') then
c------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
c------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

      xl12=lnrat(-s(ju,jd), musq)

      ! copied from coefs_new
c----- correction to the vertex in units of as/4/pi*cf
      cv0=-2._dp*(epinv2*epinv-epinv*real(xl12))-real(xl12**2)
     &           -3._dp*(epinv-real(xl12))-7._dp-eta

      ! we just multiply one of the amps below with cv0


      do nu=1,4
      pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      twoptg=
     & 2._dp*(pt(4)*p(jg,4)-pt(1)*p(jg,1)-pt(2)*p(jg,2)-pt(3)*p(jg,3))

      sne=s(jn,je)
      sdu=s(jd,ju)
      if (sdu < 0._dp) then
      prop=(sdu-wmass**2)**2
      else
      prop=(sdu-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     &    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop

c  -Lefthanded gluon tb-line
      ampf(1)=
     & -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))*zb(jb,ju)*za(ju,jd)
     & +mt**2*zb(je,jb)*za(jg,jd)
      ampf(1)=za(jc,jn)*zb(ju,jb)/(twoptg*zb(jg,jb))*ampf(1)

c---eikonal form
c      ampf(1)=
c     & -(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
c     & *(za(jg,ju)*zb(ju,jb)+za(jg,jd)*zb(jd,jb))
c     & /zb(jg,jb)
c     & -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))
c     & *za(jg,jd)
c      ampf(1)=za(jc,jn)*zb(ju,jb)/twoptg*ampf(1)

c  -Righthanded gluon tb-line
c---eikonal form
      ampf(2)=
     & +(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))*zb(ju,jg)/zb(ju,jb)
     & -(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & *(zb(jg,ju)*za(ju,jb)+zb(jg,jd)*za(jd,jb))/twoptg
      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jb,jg)*ampf(2)

c---alternative form
c      ampf(2)=
c     & +(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))*zb(ju,jg)/zb(ju,jb)
c     & +(-(zb(je,jc)*za(jc,jb)+zb(je,jn)*za(jn,jb))
c     &   *(zb(jg,ju)*za(ju,jd)+zb(jg,jb)*za(jb,jd))
c     & +mt**2*zb(je,jg)*za(jb,jd))/twoptg
c      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jb,jg)*ampf(2)


      !ubtdg_h_virt=(abs(ampf(1))**2+abs(ampf(2))**2)/prop

      ! where is the factor of two?
      ubtdg_h_virt = real(ampf(1)*conjg(ampf(1))*cv0)/prop +
     &      real(ampf(2)*conjg(ampf(2))*cv0)/prop

      return
      end

      subroutine singletop_jet_light_heavy_vr_all(p,msqvall)
        use singletop2_nnlo_vars
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
      real(dp), intent(out) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp):: fac

      call spinoru(7,p,za,zb)

      msqvall = 0._dp

      if (nwz == +1) then
                fac = 2._dp*(4*pi*as_heavy_beam2)*cf*gw**8*xn**2
                fac = fac * cf*(as_light_beam1/2._dp/pi)
                msqvall([2,4],5, 1,1) =
     &              aveqq*fac*(ubtdg_h_virt(1,2,3,4,5,6,7,p,renscale_beam1_islight_onlight**2))
                msqvall([-1,-3],5, 1,1) =
     &              aveqq*fac*(ubtdg_h_virt(6,2,3,4,5,1,7,p,renscale_beam1_islight_onlight**2))
                msqvall([2,4],0, 3,1) =
     &              aveqg*fac*(ubtdg_h_virt(1,7,3,4,5,6,2,p,renscale_beam1_islight_onlight**2))
                msqvall([-1,-3],0, 3,1) =
     &              aveqg*fac*(ubtdg_h_virt(6,7,3,4,5,1,2,p,renscale_beam1_islight_onlight**2))

                fac = 2._dp*(4*pi*as_heavy_beam1)*cf*gw**8*xn**2
                fac = fac * cf*(as_light_beam2/2._dp/pi)
                msqvall(5,[2,4], 1,2) =
     &              aveqq*fac*(ubtdg_h_virt(2,1,3,4,5,6,7,p,renscale_beam2_islight_onlight**2))
                msqvall(5,[-1,-3], 1,2) =
     &              aveqq*fac*(ubtdg_h_virt(6,1,3,4,5,2,7,p,renscale_beam2_islight_onlight**2))
                msqvall(0,[2,4], 3,2) =
     &              aveqg*fac*(ubtdg_h_virt(2,7,3,4,5,6,1,p,renscale_beam2_islight_onlight**2))
                msqvall(0,[-1,-3], 3,2) =
     &              aveqg*fac*(ubtdg_h_virt(6,7,3,4,5,2,1,p,renscale_beam2_islight_onlight**2))
      elseif (nwz == -1) then
          error stop "nwz = -1 not implemented in qqb_tbb_g_heavy"
      endif

      end subroutine

      ! like qqb_tbb_g_heavy_all but with corr_on_beam reidentified w.r.t. light
      ! line
      subroutine qqb_tbb_g_heavy_all_swap(p,msqall)
        use singletop2_nnlo_vars
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
                fac=2._dp*(4*pi*as_heavy_beam2)*cf*gw**8*xn**2
                msqall([2,4],5, 1,1) = aveqq*fac*(ubtdg_h(1,2,3,4,5,6,7,p))
                msqall([-1,-3],5, 1,1) = aveqq*fac*(ubtdg_h(6,2,3,4,5,1,7,p))
                msqall([2,4],0, 3,1) = aveqg*fac*(ubtdg_h(1,7,3,4,5,6,2,p))
                msqall([-1,-3],0, 3,1) = aveqg*fac*(ubtdg_h(6,7,3,4,5,1,2,p))

                fac=2._dp*(4*pi*as_heavy_beam1)*cf*gw**8*xn**2
                msqall(5,[2,4], 1,2) = aveqq*fac*(ubtdg_h(2,1,3,4,5,6,7,p))
                msqall(5,[-1,-3], 1,2) = aveqq*fac*(ubtdg_h(6,1,3,4,5,2,7,p))
                msqall(0,[2,4], 3,2) = aveqg*fac*(ubtdg_h(2,7,3,4,5,6,1,p))
                msqall(0,[-1,-3], 3,2) = aveqg*fac*(ubtdg_h(6,7,3,4,5,2,1,p))
      elseif (nwz == -1) then
          error stop "nwz = -1 not implemented in qqb_tbb_g_heavy"
      endif

      end subroutine

      ! copy of lumxmsq_singletop.f90, singletop2_scet_z
      subroutine singletop_jet_light_heavy_vr_z(p,z)
          use singletop2_nnlo_vars
        implicit none

        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'PR_new.f'
        include 'PR_stop.f'
        include 'agq.f'
        include 'nwz.f'
        include 'scheme.f'

        real(dp), intent(in) :: p(mxpart,4), z
        real(dp) :: dot, if_qq, fi_qq, if_qg

        integer :: is
        real(dp) :: xl16,xl26

        xl16 = log(-2*dot(p,1,6)/renscale_beam1_islight_onlight**2)
        xl26 = log(-2*dot(p,2,6)/renscale_beam2_islight_onlight**2)

        Q1 = zip
        Q2 = zip
        B1 = zip
        B2 = zip

        scheme = 'tH-V'

        do is=1,3
            ! corr_on_beam = 1
            B1(q,q,b,is) = as_light_beam1/2/pi * cf*(if_qq(z,xl16,is) + fi_qq(z,xl16,is))
            Q1(q,g,q,is) = as_light_beam1/2/pi * tr * if_qg(z, xl16,is)

            ! corr_on_beam = 2
            B2(q,q,b,is) = as_light_beam2/2/pi * cf*(if_qq(z,xl26,is) + fi_qq(z,xl26,is))
            Q2(q,g,q,is) = as_light_beam2/2/pi * tr * if_qg(z, xl26,is)
        enddo

      end subroutine singletop_jet_light_heavy_vr_z

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
          prop=prop *((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)
     &        *((st-mt**2)**2+(mt*twidth)**2)
          qqbtbbar=s(jc,jn)*s(ju,jb)
     &     *(-(s(ju,jd)+s(jb,jd))*(s(jn,je)+s(jc,je))-st*s(je,jd))/prop
      end

      ! * like lumxmsq_singletop_prod.f90, lumxmsq_singletop_prod
      ! * this one doesn't have radiation on the light line, for dipole subtractions
      ! * only one single taucut

      ! right.. this is actually identical to singletop_jet_light_heavy_vv_tree
      subroutine singletop_jet_light_heavy_rv_tree(p,msq)
          use types
          use SCET
          use SCET_Beamfunctions
          use SCET_Jet
          use singletop2_scet_heavy_decay, only : softfun
          use singletop2_nnlo_vars
          use singletop2_scale_m
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'scet_const.f'
          include 'vegas_common.f'
          include 'energy.f'
          include 'scheme.f'
          include 'ewcouple.f'
          include 'ckm.f'
          include 'zprods_com.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          real(dp) :: xx(2), QB(2)

          real(dp) :: beama0_heavy(-5:5)
          real(dp) :: beamb0_heavy(-5:5)
          real(dp) :: beama1(-5:5, -1:1)
          real(dp) :: beamb1(-5:5, -1:1)

          real(dp) :: soft2(2,0:4), soft(1,0:4), jet(1,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)

          real(dp) :: epinv_sav, epinv2_sav

          ! functions
          real(dp) :: massvec

          integer :: k
          real(dp) :: x

          real(dp) :: fac

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          QB(1) = -2*p(1,4)
          QB(2) = -2*p(2,4)

          scheme = 'tH-V'

          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2
          hard_ub(0) = fac*qqbtbbar(1,2,3,4,5,6)
          hard_ubarb(0) = fac*qqbtbbar(6,2,3,4,5,1)

          hard_bu(0) = fac*qqbtbbar(2,1,3,4,5,6)
          hard_bubar(0) =fac*qqbtbbar(6,1,3,4,5,2)

          ! NEEDS MSBAR RENORMALIZATION
          ! error stop "needs msbar renormalization"
          ! there we go:
          epinv_sav = epinv
          epinv2_sav = epinv2

          epinv = 0._dp
          epinv2 = 0._dp

          fac=aveqq*gw**8*xn**2*(as_heavy_beam2/2._dp/pi*cf)
          ! HARD FUNCTION WITHOUT as/4/pi !!
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1)  = fac*virtqqb_heavy(1,2,3,4,5,6, renscale_beam2_isheavy_onlight**2)
          hard_ubarb(1) = fac*virtqqb_heavy(6,2,3,4,5,1, renscale_beam2_isheavy_onlight**2)

          ! for c1 = 0, cv = 1
          !write (*,*) "COMPARISON TREE", virtqqb_heavy(1,2,3,4,5,6,172d0**2) / qqbtbbar(1,2,3,4,5,6)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac=aveqq*gw**8*xn**2*(as_heavy_beam1/2._dp/pi*cf)
          ! HARD FUNCTION WITHOUT as/4/pi !!
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqb_heavy(2,1,3,4,5,6, renscale_beam1_isheavy_onlight**2)
          hard_bubar(1) = fac*virtqqb_heavy(6,1,3,4,5,2, renscale_beam1_isheavy_onlight**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)

          epinv = epinv_sav
          epinv2 = epinv2_sav

          soft2 = softfun(1)
          soft(1,0:4) = soft2(1,0:4)

          !error stop "make sure currentNd is properly set"

          !call fdist(ih1,xx(1),singletop2_dipscale(currentNd,1),beama0_light,1)

          call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0_heavy,2)

          call xbeam1bis_new(ih2,z2int,xx(2),QB(2),beamb1, 2, renscale_beam2_isheavy_onlight,
     &            facscale_beam2_isheavy_onlight,disttau=.false.)


          call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0_heavy,1)

          !call fdist(ih2,xx(2),singletop2_dipscale(currentNd,2),beamb0_light,2)

          call xbeam1bis_new(ih1,z1int,xx(1),QB(1),beama1, 1, renscale_beam1_isheavy_onlight,
     &        facscale_beam1_isheavy_onlight,disttau=.false.)


c ub   channel
          x = massvec(-p(1,:)-p(6,:)) / mt**2

          msq = 0._dp

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beamb1(k,-1)
          jet(1,1) = beamb1(k,0)
          jet(1,2) = beamb1(k,1)/2._dp

          msq([2,4],5) = hard_ub(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

          msq([-1,-3],5) = hard_ubarb(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft, beamb0_heavy(k))


          x = massvec(-p(2,:)-p(6,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beama1(k,-1)
          jet(1,1) = beama1(k,0)
          jet(1,2) = beama1(k,1)/2._dp

          msq(5,[2,4]) = hard_bu(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft, beama0_heavy(5))

          msq(5,[-1,-3]) = hard_bubar(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft, beama0_heavy(5))


      end subroutine

      ! so this routine is used for VV splitting pieces
      subroutine singletop_jet_light_heavy_vv_tree(p,msq)
          use types
          use SCET
          use SCET_Beamfunctions
          use SCET_Jet
          use singletop2_scet_heavy_decay, only : softfun, hardfun
          use singletop2_nnlo_vars
          use singletop2_scale_m
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'scet_const.f'
          include 'vegas_common.f'
          include 'energy.f'
          include 'scheme.f'
          include 'ewcouple.f'
          include 'ckm.f'
          include 'zprods_com.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          real(dp) :: xx(2), QB(2)
          real(dp) :: epinv_sav, epinv2_sav

          real(dp) :: beama0_heavy(-5:5)
          real(dp) :: beamb0_heavy(-5:5)
          real(dp) :: beama1(-5:5, -1:1)
          real(dp) :: beamb1(-5:5, -1:1)

          real(dp) :: soft2(2,0:4), soft(1,0:4), jet(1,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)

          ! functions
          real(dp) :: massvec

          integer :: k
          real(dp) :: x

          real(dp) :: fac

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          QB(1) = -2*p(1,4)
          QB(2) = -2*p(2,4)

          scheme = 'tH-V'

          call spinoru(6,p,za,zb)

          fac=aveqq*gw**8*xn**2
          hard_ub(0) = fac*qqbtbbar(1,2,3,4,5,6)
          hard_ubarb(0) = fac*qqbtbbar(6,2,3,4,5,1)

          hard_bu(0) = fac*qqbtbbar(2,1,3,4,5,6)
          hard_bubar(0) =fac*qqbtbbar(6,1,3,4,5,2)

          ! NEEDS MSBAR RENORMALIZATION
          ! error stop "needs msbar renormalization"
          ! there we go:

          epinv_sav = epinv
          epinv2_sav = epinv2

          epinv = 0._dp
          epinv2 = 0._dp

          fac=aveqq*gw**8*xn**2*(as_heavy_beam2/2._dp/pi*cf)

          ! hard function without as/4/pi
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1)  = fac*virtqqb_heavy(1,2,3,4,5,6, renscale_beam2_isheavy_onlight**2)
          hard_ubarb(1) = fac*virtqqb_heavy(6,2,3,4,5,1, renscale_beam2_isheavy_onlight**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

c         x = massvec(-p(1,:)-p(6,:)) / mt**2
c         write (*,*) hard_ub(0), hard_ub(1) +
c    &          pi**2/6._dp*cf*hard_ub(0)
c         write (*,*) real(hardfun(1, renscale_beam2_isheavy_onheavy, x, p, 1,2,3,4,5,6, production=.true.),dp)
c         write (*,*) "UBARB"
c         write (*,*) hard_ubarb(0), hard_ubarb(1) +
c    &          pi**2/6._dp*cf*hard_ubarb(0)
c         write (*,*) real(hardfun(1, renscale_beam2_isheavy_onheavy, x, p, 6,2,3,4,5,1, production=.true.),dp)
c         pause

          fac=aveqq*gw**8*xn**2*(as_heavy_beam1/2._dp/pi*cf)
          ! hard function without as/4/pi
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqb_heavy(2,1,3,4,5,6, renscale_beam1_isheavy_onlight**2)
          hard_bubar(1) = fac*virtqqb_heavy(6,1,3,4,5,2, renscale_beam1_isheavy_onlight**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)

          epinv = epinv_sav
          epinv2 = epinv2_sav

          soft2 = softfun(1)
          soft(1,0:4) = soft2(1,0:4)

          ! z1int and z2int are module level variables, set in
          ! singletop_vvint.f90

          call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0_heavy,2)

          call xbeam1bis_new(ih2,z2int,xx(2),QB(2),beamb1, 2, renscale_beam2_isheavy_onlight,
     &            facscale_beam2_isheavy_onlight,disttau=.false.)


          call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0_heavy,1)

          call xbeam1bis_new(ih1,z1int,xx(1),QB(1),beama1, 1, renscale_beam1_isheavy_onlight,
     &        facscale_beam1_isheavy_onlight,disttau=.false.)


c ub   channel
          x = massvec(-p(1,:)-p(6,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beamb1(k,-1)
          jet(1,1) = beamb1(k,0)
          jet(1,2) = beamb1(k,1)/2._dp

          msq([2,4],5) = hard_ub(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

          msq([-1,-3],5) = hard_ubarb(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft, beamb0_heavy(k))

          x = massvec(-p(2,:)-p(6,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beama1(k,-1)
          jet(1,1) = beama1(k,0)
          jet(1,2) = beama1(k,1)/2._dp

          msq(5,[2,4]) = hard_bu(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft, beama0_heavy(5))

          msq(5,[-1,-3]) = hard_bubar(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft, beama0_heavy(5))


      end subroutine

      ! two loop below cut, also needs tree version for the
      ! splitting pieces
      subroutine singletop_jet_light_heavy_vv(p,msq)
          use types
          use SCET
          use SCET_Beamfunctions
          use SCET_Jet
          use singletop2_scet_heavy_decay, only : softfun
          use singletop2_nnlo_vars
          use singletop2_scale_m
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'scet_const.f'
          include 'vegas_common.f'
          include 'energy.f'
          include 'scheme.f'
          include 'ewcouple.f'
          include 'ckm.f'
          include 'zprods_com.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          real(dp) :: xx(2), QB(2)

          real(dp) :: beama0_light(-5:5), beama0_heavy(-5:5)
          real(dp) :: beamb0_light(-5:5), beamb0_heavy(-5:5)
          real(dp) :: beama1(-5:5, -1:1)
          real(dp) :: beamb1(-5:5, -1:1)

          real(dp) :: soft2(2,0:4), soft(1,0:4), jet(1,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)

          ! functions
          real(dp) :: massvec

          integer :: k
          real(dp) :: x

          real(dp) :: fac

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          QB(1) = -2*p(1,4)
          QB(2) = -2*p(2,4)

          scheme = 'tH-V'

          call spinoru(6,p,za,zb)

          fac = (as_light_beam1/2._dp/pi*cf) * aveqq*gw**8*xn**2
          hard_ub(0) = fac*virtqqb_light(1,2,3,4,5,6, renscale_beam1_islight_onlight**2)
          hard_ubarb(0) = fac*virtqqb_light(6,2,3,4,5,1, renscale_beam1_islight_onlight**2)

          fac = (as_light_beam2/2._dp/pi*cf) * aveqq*gw**8*xn**2
          hard_bu(0) = fac*virtqqb_light(2,1,3,4,5,6, renscale_beam2_islight_onlight**2)
          hard_bubar(0) = fac*virtqqb_light(6,1,3,4,5,2, renscale_beam2_islight_onlight**2)


          ! light line with epinv
          ! heavy line MSbar renormalized (epinv = epinv2 = 0)

          fac=aveqq*gw**8*xn**2*(as_heavy_beam2/2._dp/pi*cf)
          fac = fac * (as_light_beam1/2._dp/pi*cf)
          ! hard function without as/4/pi
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1)  = fac*virtqqb_heavy_light(1,2,3,4,5,6,
     &          renscale_beam2_isheavy_onlight**2,
     &          renscale_beam1_islight_onlight**2)
          hard_ubarb(1) = fac*virtqqb_heavy_light(6,2,3,4,5,1,
     &          renscale_beam2_isheavy_onlight**2,
     &          renscale_beam1_islight_onlight**2)

          ! part of scheme conversion
          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac=aveqq*gw**8*xn**2*(as_heavy_beam1/2._dp/pi*cf)
          fac = fac * (as_light_beam2/2._dp/pi*cf)
          ! hard function without as/4/pi
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqb_heavy_light(2,1,3,4,5,6,
     &          renscale_beam1_isheavy_onlight**2,
     &          renscale_beam2_islight_onlight**2)
          hard_bubar(1) = fac*virtqqb_heavy_light(6,1,3,4,5,2,
     &          renscale_beam1_isheavy_onlight**2,
     &          renscale_beam2_islight_onlight**2)

          ! part of scheme conversion
          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)


          soft2 = softfun(1)
          soft(1,0:4) = soft2(1,0:4)

          call fdist(ih1,xx(1),facscale_beam1_islight_onlight,beama0_light,1)
          call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0_heavy,2)

          call xbeam1bis_new(ih2,z2int,xx(2),QB(2),beamb1, 2, renscale_beam2_isheavy_onlight,
     &            facscale_beam2_isheavy_onlight,disttau=.false.)


          call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0_heavy,1)
          call fdist(ih2,xx(2),facscale_beam2_islight_onlight,beamb0_light,2)

          call xbeam1bis_new(ih1,z1int,xx(1),QB(1),beama1, 1, renscale_beam1_isheavy_onlight,
     &        facscale_beam1_isheavy_onlight,disttau=.false.)

c ub   channel
          x = massvec(-p(1,:)-p(6,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beamb1(k,-1)
          jet(1,1) = beamb1(k,0)
          jet(1,2) = beamb1(k,1)/2._dp

          msq([2,4],5) = beama0_light([2,4])*hard_ub(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

          msq([-1,-3],5) = beama0_light([-1,-3])*hard_ubarb(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft, beamb0_heavy(k))


          x = massvec(-p(2,:)-p(6,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beama1(k,-1)
          jet(1,1) = beama1(k,0)
          jet(1,2) = beama1(k,1)/2._dp

          msq(5,[2,4]) = beamb0_light([2,4])*hard_bu(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft, beama0_heavy(5))

          msq(5,[-1,-3]) = beamb0_light([-1,-3])*hard_bubar(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft, beama0_heavy(5))


      end subroutine singletop_jet_light_heavy_vv

      subroutine singletop_jet_light_heavy_rv(p,msq)
          use types
          use SCET
          use SCET_Beamfunctions
          use SCET_Jet
          use singletop2_scet_heavy_decay, only : softfun, hardfun
          use singletop2_nnlo_vars
          use singletop2_scale_m
          use singletop_jet, only: ampsq_ugd_tdkb
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'masses.f'
          include 'epinv.f'
          include 'epinv2.f'
          include 'scet_const.f'
          include 'vegas_common.f'
          include 'energy.f'
          include 'scheme.f'
          include 'ewcouple.f'
          include 'ckm.f'
          include 'zprods_com.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib)

          real(dp) :: xx(2), QB(2)
          real(dp) :: epinv_sav, epinv2_sav

          real(dp) :: beama0_heavy(-5:5)
          real(dp) :: beamb0_heavy(-5:5)
          real(dp) :: beama1(-5:5, -1:1)
          real(dp) :: beamb1(-5:5, -1:1)

          real(dp) :: soft2(2,0:4), soft(1,0:4), jet(1,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)
          real(dp) :: hard_gb(0:1), hard_bg(0:1)

          ! functions
          real(dp) :: massvec

          integer :: k
          real(dp) :: x

          real(dp) :: fac
          real(dp), parameter :: genfac = 2._dp

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          QB(1) = -2*p(1,4)
          QB(2) = -2*p(2,4)

          scheme = 'tH-V'

          call spinoru(7,p,za,zb)

          ! BORN PIECES
          fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
          hard_ub(0) = aveqq*fac*ampsq_ugd_tdkb(1,2,3,4,5,6,7,p)
          hard_ubarb(0) = aveqq*fac*ampsq_ugd_tdkb(6,2,3,4,5,1,7,p)
          hard_gb(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,2,3,4,5,6,1,p)

          fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
          hard_bu(0) = aveqq*fac*ampsq_ugd_tdkb(2,1,3,4,5,6,7,p)
          hard_bubar(0) = aveqq*fac*ampsq_ugd_tdkb(6,1,3,4,5,2,7,p)
          hard_bg(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,1,3,4,5,6,2,p)


          ! VIRT PIECES

          epinv_sav = epinv
          epinv2_sav = epinv2

          epinv = 0._dp
          epinv2 = 0._dp

          fac = 2._dp*(4*pi*as_light_beam1)*cf*gw**8*xn**2
          fac = fac * cf*(as_heavy_beam2/2._dp/pi)
          ! HARD FUNCTION WITHOUT as/4/pi !!
          fac = fac / (as_heavy_beam2/4._dp/pi)

          ! for c1 = 0, cv = 1
          !write (*,*) "COMPARISON", ubtdg_l_virt(1,2,3,4,5,6,7,p,170d0**2) / ampsq_ugd_tdkb(1,2,3,4,5,6,7,p)

          hard_ub(1) = aveqq*fac*ubtdg_l_virt(1,2,3,4,5,6,7,p,renscale_beam2_isheavy_onheavy**2)
          hard_ubarb(1) = aveqq*fac*ubtdg_l_virt(6,2,3,4,5,1,7,p,renscale_beam2_isheavy_onheavy**2)
          hard_gb(1) = genfac*aveqg*fac*ubtdg_l_virt(7,2,3,4,5,6,1,p,renscale_beam2_isheavy_onheavy**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)
          hard_gb(1) = hard_gb(1) + pi**2/6._dp * cf * hard_gb(0)

          fac = 2._dp*(4*pi*as_light_beam2)*cf*gw**8*xn**2
          fac = fac * cf*(as_heavy_beam1/2._dp/pi)
          ! HARD FUNCTION WITHOUT as/4/pi !!
          fac = fac / (as_heavy_beam1/4._dp/pi)

          hard_bu(1) = aveqq*fac*ubtdg_l_virt(2,1,3,4,5,6,7,p,renscale_beam1_isheavy_onheavy**2)
          hard_bubar(1) = aveqq*fac*ubtdg_l_virt(6,1,3,4,5,2,7,p,renscale_beam1_isheavy_onheavy**2)
          hard_bg(1) = genfac*aveqg*fac*ubtdg_l_virt(7,1,3,4,5,6,2,p,renscale_beam1_isheavy_onheavy**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)
          hard_bg(1) = hard_bg(1) + pi**2/6._dp * cf * hard_bg(0)

          epinv = epinv_sav
          epinv2 = epinv2_sav

          soft2 = softfun(1)
          soft(1,0:4) = soft2(1,0:4)

          ! z1int and z2int are module level variables, set in
          ! singletop_vvint.f90

          call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0_heavy,2)

          call xbeam1bis_new(ih2,z2int,xx(2),QB(2),beamb1, 2, renscale_beam2_isheavy_onlight,
     &            facscale_beam2_isheavy_onlight,disttau=.false.)


          call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0_heavy,1)

          call xbeam1bis_new(ih1,z1int,xx(1),QB(1),beama1, 1, renscale_beam1_isheavy_onlight,
     &        facscale_beam1_isheavy_onlight,disttau=.false.)


c ub   channel
          x = massvec(-p(1,:)-p(6,:)-p(7,:)) / mt**2

          msq = 0._dp

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beamb1(k,-1)
          jet(1,1) = beamb1(k,0)
          jet(1,2) = beamb1(k,1)/2._dp

          msq([2,4],5, 1) = hard_ub(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft, beamb0_heavy(k))

          msq([-1,-3],5, 1) = hard_ubarb(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft, beamb0_heavy(k))

          msq(0,5, 1) = hard_gb(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_gb/hard_gb(0), jet, soft, beamb0_heavy(k))

          x = massvec(-p(2,:)-p(6,:)-p(7,:)) / mt**2

          k = 5
          ! convert to Laplace-space log-notation
          jet(1,0) = beama1(k,-1)
          jet(1,1) = beama1(k,0)
          jet(1,2) = beama1(k,1)/2._dp

          msq(5,[2,4], 1) = hard_bu(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft, beama0_heavy(5))

          msq(5,[-1,-3], 1) = hard_bubar(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft, beama0_heavy(5))

          msq(5,0, 1) = hard_bg(0)*
     &        assemble_production_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bg/hard_bg(0), jet, soft, beama0_heavy(5))


      end subroutine

       function virtqqb_heavy_light(ju,jb,jn,je,jc,jd,musq_heavy,musq_light)
       implicit none
       real(dp):: virtqqb_heavy_light

       integer:: ju,jd,jn,je,jc,jb
       include 'constants.f'
       include 'nf.f'
       include 'mxpart.f'
       include 'cplx.h'
       include 'masses.f'
       include 'sprods_com.f'
       include 'zprods_com.f'
       include 'scheme.f'
       include 'epinv.f'
       include 'epinv2.f'
       real(dp):: snec,prop,cv,cv0,mtsq
       complex(dp):: c1,amp,ampho
       real(dp) :: cv0_dummy

       real(dp), intent(in) :: musq_heavy, musq_light


       mtsq=mt**2
       snec=+s(jn,je)+s(je,jc)+s(jc,jn)

       ! calculate light factor with epinv dependence and musq_light
       call coefs_new(s(ju,jd),mtsq,cv0,cv,c1,musq_light,epinv,epinv2)
       ! calculate heavy line factors MSbar-renormalized with musq_heavy
       call coefs_new(s(ju,jd),mtsq,cv0_dummy,cv,c1,musq_heavy,0d0,0d0)

       if (s(ju,jd) < 0._dp) then
       prop=(s(ju,jd)-wmass**2)**2
       else
       prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
       endif
       prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
       prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

       amp=za(jc,jn)*zb(ju,jb)
     &  *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
       ! no cv0 here
       ampho=za(jc,jn)*zb(ju,jb)
     &  *(cplx1(cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     &  +c1*chalf*zb(je,jb)*za(jb,jd))

       ! here multiplied with cv0
       virtqqb_heavy_light=real(cv0*amp*conjg(ampho))/prop

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

       ! only heavy line corrections
       cv0 = 0._dp

       amp=za(jc,jn)*zb(ju,jb)
     &  *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
       ampho=za(jc,jn)*zb(ju,jb)
     &  *(cplx1(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     &  +c1*chalf*zb(je,jb)*za(jb,jd))

       virtqqb_heavy=real(amp*conjg(ampho))/prop

       end



      ! like passed_taucut_heavyprod, but without p(8,:) input
      function passed_taucut_lxh(pparton,scetreweight_local,taucut_in)
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
          include 'kpart.f'
          include 'ipsgen.f'

          logical :: passed_taucut_lxh
          real(dp), intent(in) :: pparton(mxpart,4)
          real(dp), intent(inout), optional :: scetreweight_local(:)
          real(dp), intent(in), optional :: taucut_in

          integer j
          real(dp) :: tau,tauc
          real(dp) :: qsq

          real(dp) :: dotvec, massvec, nn(4)

          logical :: bin
          common/bin/bin

          passed_taucut_lxh = .false.

          if (present(taucut_in)) then
              tauc = taucut_in
          else
              tauc = taucut
          endif

          ! remember: corr_on_beam and kpart is w.r.t. dipole corrections
          if (corr_on_beam == 2) then
              nn(4) = 1._dp
              nn(1:3) = -pparton(1,1:3)/pparton(1,4)

              if (kpart == kreal) then
                  qsq = massvec(-pparton(2,:)-pparton(6,:)-pparton(8,:))
                  nn(:) = -dotvec(pparton(7,:)+pparton(1,:),nn(:))*pparton(1,:)/pparton(1,4) / two
              else
                  qsq = massvec(-pparton(2,:)-pparton(6,:))
                  nn(:) = -dotvec(pparton(7,:)+pparton(1,:),nn(:))*pparton(1,:)/pparton(1,4) / two
              endif
          else
              nn(4) = 1._dp
              nn(1:3) = -pparton(2,1:3)/pparton(2,4)

              if (kpart == kreal) then
                  qsq = massvec(-pparton(1,:)-pparton(6,:)-pparton(8,:))
                  nn(:) = -dotvec(pparton(7,:)+pparton(2,:),nn(:))*pparton(2,:)/pparton(2,4) / two
              else
                  qsq = massvec(-pparton(1,:)-pparton(6,:))
                  nn(:) = -dotvec(pparton(7,:)+pparton(2,:),nn(:))*pparton(2,:)/pparton(2,4) / two
              endif
          endif

          if (kpart == kreal) then
              ! double real
              ! 8 is on the light line
              tau = 2._dp * dotvec(pparton(7,:), nn) / (mt**2 - qsq)
          elseif (kpart == kvirt) then
              tau = 2._dp * dotvec(pparton(7,:), nn) / (mt**2 - qsq)
          endif

          !write (*,*) "EVENT TAU = ", tau


c---      heck to make sure no NaN
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

          passed_taucut_lxh = .true.

      end function passed_taucut_lxh

      subroutine singletop_jet_light_heavy_rr_all(p,msqall)
#define WITH_BBAR 3
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

        real(dp) :: fac,facqg,facqq,facgg,tmp(-nf:nf,max_bcontrib)
        real(dp) :: genfac

        integer :: iperm,i1,i2

        call spinoru(8,p,za,zb)

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Rcola results for a specific final state
        genfac=2._dp

        msqall(:,:,:,:) = 0._dp

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1
        do iperm=1,maxbeams
            corr_on_beam = beams_enabled(iperm)

            tmp(:,:)=0._dp

            if (corr_on_beam == 1) then
              i1=1
              i2=2
              fac=(fourpi*as_light_beam1)*(fourpi*as_heavy_beam2)
            else
              i1=2
              i2=1
              fac=(fourpi*as_light_beam2)*(fourpi*as_heavy_beam1)
            endif
            facqq=fac*aveqq
            facqg=fac*aveqg
            facgg=fac*avegg

            if (iand(partons_enabled, quarkChannel) > 0) then
                tmp([2,4], 1) = facqq*msqlightxheavy(i1,i2,3,4,5,6,7,8)   ! u b -> t d g g
                tmp([-1,-3], 1) = facqq*msqlightxheavy(6,i2,3,4,5,i1,7,8) ! d~ b -> t u~ g g
            endif

            if (iand(partons_enabled, gluonChannel) > 0) then
                tmp(0, 1)=facqg*genfac*msqlightxheavy(8,i2,3,4,5,6,7,i1)         ! g b -> t d g u~
                tmp([2,4], WITH_BBAR)=facqg*msqlightxheavy(i1,7,3,4,5,6,i2,8)    ! u g -> t d b~ g
                tmp([-1,-3], WITH_BBAR)=facqg*msqlightxheavy(6,7,3,4,5,i1,i2,8)  ! d~ g -> t u~ b~ g
                tmp(0, WITH_BBAR)=facgg*genfac*msqlightxheavy(8,7,3,4,5,6,i2,i1) ! g g -> t d b~ u~
            endif

c update     msq array,
            if (corr_on_beam == 1) then
                msqall(:,5, 1, 1) = tmp(:,1)
                msqall(:,0, WITH_BBAR, 1) = tmp(:, WITH_BBAR)
            else
                msqall(5,:, 1, 2) = tmp(:,1)
                msqall(0,:, WITH_BBAR, 2) = tmp(:, WITH_BBAR)
            endif

        enddo

      end subroutine

      function assemble_production_pieces(x, scale, as, taucut, h, j, s, j0)
          use types
          use constants
          implicit none
          include 'masses.f'
          include 'kpart.f'

          real(dp), intent(in) :: x, scale, as
          real(dp), intent(in) :: taucut
          real(dp), intent(in) :: h(0:1)
          real(dp), intent(in) :: j(1,0:4)
          real(dp), intent(in) :: s(1,0:4)
          real(dp), intent(in) :: j0
          real(dp) :: assemble_production_pieces

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

          !if (coeffonly) then
              full(0,0) = 0._dp
          !else
              !full(0,0) = j0
          !endif


          full(1,2) = j(1,2) + j0*s(1,2)
          full(1,1) = j(1,1) + 2*LJX*j(1,2) + j0*(s(1,1) + 2*LSX*s(1,2))
          full(1,0) = j(1,0) + LJX*(j(1,1) + LJX*j(1,2)) +
     &        j0*(h(1) + s(1,0) + LSX*(s(1,1) + LSX*s(1,2)))

          assemble_production_pieces = full(0,0) +
     &        as/4._dp/pi*(
     &            full(1,0) + full(1,1)*log(taucut) + full(1,2)*log(taucut)**2 )

      end function

      subroutine singletop_virt_light(p,msqv)
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


      subroutine singletop_jet_light_heavy_rv_gs(p,ndmx,msq)
        use singletop2_scale_m
        use singletop2_nnlo_vars
        use singletop_jet2
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'nwz.f'
        include 'ptilde.f'
        include 'qqgg.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ndmx
        real(dp), intent(inout) :: msq(ndmx,-nf:nf,-nf:nf)

        real(dp) :: dummyv(-nf:nf,-nf:nf), dsubv
        real(dp) :: msq18_6(-nf:nf,-nf:nf), msq68_1(-nf:nf,-nf:nf)
        real(dp) :: msq28_6(-nf:nf,-nf:nf), msq68_2(-nf:nf,-nf:nf)
        real(dp) :: msq18_2(-nf:nf,-nf:nf), msq16_2(-nf:nf,-nf:nf)
        real(dp) :: msq28_1(-nf:nf,-nf:nf), msq26_1(-nf:nf,-nf:nf)
        real(dp) :: sub18_6(4), sub68_1(4)
        real(dp) :: sub28_6(4), sub68_2(4)
        real(dp) :: sub18_2(4), sub16_2(4)
        real(dp) :: sub28_1(4), sub26_1(4)

c        real(dp) :: msq17_6(-nf:nf,-nf:nf), msq67_1(-nf:nf,-nf:nf)
c        real(dp) :: msq27_6(-nf:nf,-nf:nf), msq67_2(-nf:nf,-nf:nf)
c        real(dp) :: msq17_2(-nf:nf,-nf:nf)
c        real(dp) :: msq27_1(-nf:nf,-nf:nf)
c        real(dp) :: sub17_6(4), sub67_1(4)
c        real(dp) :: sub27_6(4), sub67_2(4)
c        real(dp) :: sub17_2(4)
c        real(dp) :: sub27_1(4)

        integer :: noglue(4) = [-1,-3,2,4]

        external donothing_gvec

        ndmax = 8

        msq = 0._dp

        corr_islight = .true.

        corr_beam1 = .true.
        call dips(1,p, 1,7,6,sub18_6,dsubv,msq18_6,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(2,p, 6,7,1,sub68_1,dsubv,msq68_1,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(3,p, 1,7,6,sub18_2,dsubv,msq18_2,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(4,p, 1,6,7,sub16_2,dsubv,msq16_2,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)

        msq(1,noglue,5) = 2._dp*cf*sub18_6(qq)*msq18_6(noglue,5)
        msq(2,noglue,5) = 2._dp*cf*sub68_1(qq)*msq68_1(noglue,5)
        msq(3,0,5) = 2._dp*tr*sub18_2(qg) * sum(msq18_2(1:5,5))
        msq(4,0,5) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1, 5))

        corr_beam1 = .false.

        call dips(5,p, 2,7,6,sub28_6,dsubv,msq28_6,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(6,p, 6,7,2,sub68_2,dsubv,msq68_2,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(7,p, 2,7,6,sub28_1,dsubv,msq28_1,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)
        call dips(8,p, 2,6,7,sub26_1,dsubv,msq26_1,dummyv,singletop_jet_light_heavy_rv_tree,donothing_gvec)

        msq(5,5,noglue) = 2._dp*cf*sub28_6(qq)*msq28_6(5,noglue)
        msq(6,5,noglue) = 2._dp*cf*sub68_2(qq)*msq68_2(5,noglue)
        msq(7,5,0) = 2._dp*tr*sub28_1(qg) * sum(msq28_1(5,1:5))
        msq(8,5,0) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(5,-5:-1))

      end subroutine singletop_jet_light_heavy_rv_gs

      subroutine singletop_jet_light_heavy_rr_gs_all(p,ndmx,msqall)
        use singletop2_scale_m
        use singletop2_nnlo_vars
        use singletop_jet2
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'nwz.f'
        include 'ptilde.f'
        include 'qqgg.f'

        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ndmx
        real(dp), intent(inout) :: msqall(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

        real(dp) :: dummyv(-nf:nf,-nf:nf), dsubv
        real(dp) :: msq18_6(-nf:nf,-nf:nf), msq68_1(-nf:nf,-nf:nf)
        real(dp) :: msq28_6(-nf:nf,-nf:nf), msq68_2(-nf:nf,-nf:nf)
        real(dp) :: msq18_2(-nf:nf,-nf:nf), msq16_2(-nf:nf,-nf:nf)
        real(dp) :: msq28_1(-nf:nf,-nf:nf), msq26_1(-nf:nf,-nf:nf)
        real(dp) :: sub18_6(4), sub68_1(4)
        real(dp) :: sub28_6(4), sub68_2(4)
        real(dp) :: sub18_2(4), sub16_2(4)
        real(dp) :: sub28_1(4), sub26_1(4)

c        real(dp) :: sub17_6(4), sub67_1(4)
c        real(dp) :: sub27_6(4), sub67_2(4)
c        real(dp) :: sub17_2(4)
c        real(dp) :: sub27_1(4)

        integer :: noglue(4) = [-1,-3,2,4]

        external donothing_gvec

        ndmax = 8

        msqall = 0._dp

        corr_islight = .true.

        ! the cobswitch routine treats corr_on_beam w.r.t. light line internally

        corr_on_beam = 1
        corr_beam1 = .true.
        call dips(1,p, 1,8,6,sub18_6,dsubv,msq18_6,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)
        call dips(2,p, 6,8,1,sub68_1,dsubv,msq68_1,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)
        call dips(3,p, 1,8,6,sub18_2,dsubv,msq18_2,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)

        ! _swap routine is needed to get 6 and 7 into correct order after dipole transformation
        call dips(4,p, 1,6,8,sub16_2,dsubv,msq16_2,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)

        msqall(1,noglue,5, 1,1) = 2._dp*cf*sub18_6(qq)*msq18_6(noglue,5)
        msqall(1,noglue,0, WITH_BBAR,1) = 2._dp*cf*sub18_6(qq)*msq18_6(noglue,0)

        msqall(2,noglue,5, 1,1) = 2._dp*cf*sub68_1(qq)*msq68_1(noglue,5)
        msqall(2,noglue,0, WITH_BBAR,1) = 2._dp*cf*sub68_1(qq)*msq68_1(noglue,0)

        msqall(3,0,5, 1,1) = 2._dp*tr*sub18_2(qg) * sum(msq18_2(1:5,5))
        msqall(4,0,5, 1,1) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1, 5))

        msqall(3,0,0, WITH_BBAR,1) = 2._dp*tr*sub18_2(qg) * sum(msq18_2(1:5,0))
        msqall(4,0,0, WITH_BBAR,1) = 2._dp*tr*sub16_2(qg) * sum(msq16_2(-5:-1, 0))

        corr_on_beam = 2
        corr_beam1 = .false.

        call dips(5,p, 2,8,6,sub28_6,dsubv,msq28_6,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)
        call dips(6,p, 6,8,2,sub68_2,dsubv,msq68_2,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)
        call dips(7,p, 2,8,6,sub28_1,dsubv,msq28_1,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)

        ! _swap routine is needed to get 6 and 7 into correct order after dipole transformation
        call dips(8,p, 2,6,8,sub26_1,dsubv,msq26_1,dummyv,singletop_jet_heavy_cobswitch,donothing_gvec)

        msqall(5,5,noglue, 1,2) = 2._dp*cf*sub28_6(qq)*msq28_6(5,noglue)
        msqall(5,0,noglue, WITH_BBAR,2) = 2._dp*cf*sub28_6(qq)*msq28_6(0,noglue)

        msqall(6,5,noglue, 1,2) = 2._dp*cf*sub68_2(qq)*msq68_2(5,noglue)
        msqall(6,0,noglue, WITH_BBAR,2) = 2._dp*cf*sub68_2(qq)*msq68_2(0,noglue)

        msqall(7,5,0, 1,2) = 2._dp*tr*sub28_1(qg) * sum(msq28_1(5,1:5))
        msqall(7,0,0, WITH_BBAR,2) = 2._dp*tr*sub28_1(qg) * sum(msq28_1(0,1:5))

        msqall(8,5,0, 1,2) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(5,-5:-1))
        msqall(8,0,0, WITH_BBAR,2) = 2._dp*tr*sub26_1(qg) * sum(msq26_1(0,-5:-1))

      end subroutine singletop_jet_light_heavy_rr_gs_all



      function msqlightxheavy(p1,p2,p3,p4,p5,p6,p7,p8)

c Matrix element squared for (four) diagrams of the form:

c                           o p8
c                          o
c                         o
c     p1 ---->------------->------ p6
c                    $
c                    $
c                    $
c                    $   t
c     p2 ---->---------||-->------ p5     (on-shell top, p345^2 = mt^2)
c               o          $
c                o          $ -->-- p3
c                 o          |
c                   p7       ^
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
      integer:: p1,p2,p3,p4,p5,p6,p7,p8,iperm,j7,j8
      real(dp):: msqlightxheavy,s168,s34,s3457,prop34,prop168,prop3457
      complex(dp):: app,apm,amp,amm,zab2,zab3

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)

      s34=s(p3,p4)
      prop34=(s34-wmass**2)**2
      if (s34 > 0._dp) then
        prop34=prop34+(wmass*wwidth)**2
      endif

      msqlightxheavy=0._dp

c Could loop to 2 here, but do this externally in order to make crossings easier
      do iperm=1,1

      if (iperm == 1) then
        j7=p7
        j8=p8
      else
        j7=p8
        j8=p7
      endif

      s168=s(p1,p6)+s(p1,j8)+s(p6,j8)
      prop168=(s168-wmass**2)**2
      if (s168 > 0._dp) then
        prop168=prop168+(wmass*wwidth)**2
      endif

c No width necessary for this top propagator since we assume s345 = mt^2
      s3457=s(p3,p4)+s(p3,p5)+s(p3,j7)+s(p4,p5)+s(p4,j7)+s(p5,j7)
      prop3457=s3457-mt**2

      app =  + prop3457**(-1) * (
     &     + 1./(za(p1,j8))/(za(p6,j7))*za(p3,p5)*zb(p2,j8)*zab2(p6,p3,
     &    p5,p4)*zab3(p6,p3,p4,p5,j7)
     &     + 1./(za(p1,j8))/(za(p6,j7))/(za(p6,j8))*za(p1,p6)*za(p3,p5)*
     &    zb(p1,p2)*zab2(p6,p3,p5,p4)*zab3(p6,p3,p4,p5,j7)
     &     )
      app = app + 1./(za(p1,j8))/(za(p2,j7))/(za(p6,j7))*za(p3,p5)*zab2(
     & p6,p2,j7,j8)*zab2(p6,p3,p5,p4)
     &     - 1./(za(p1,j8))/(za(p2,j7))/(za(p6,j7))/(za(p6,j8))*za(p1,p6
     &    )*za(p3,p5)*zab2(p6,p2,j7,p1)*zab2(p6,p3,p5,p4)

      apm =  + prop3457**(-1) * (
     &     - 1./(za(p6,j7))/(zb(p1,j8))/(zb(p6,j8))*za(p3,p5)*zb(p1,p2)*
     &    zb(p1,p6)*zab2(p6,p3,p5,p4)*zab3(p6,p3,p4,p5,j7)
     &     - 1./(za(p6,j7))/(zb(p6,j8))*za(p3,p5)*za(p6,j8)*zb(p1,p2)*
     &    zb(p4,j7)*mt**2
     &     - 1./(za(p6,j7))/(zb(p6,j8))*za(p3,p5)*zb(p1,p2)*zab2(p6,p3,
     &    p5,p4)*zab3(j8,p3,p4,p5,j7)
     &     )
      apm = apm + 1./(za(p2,j7))/(za(p6,j7))/(zb(p1,j8))/(zb(p6,j8))*za(
     & p3,p5)*zb(p1,p6)*zab2(p6,p2,j7,p1)*zab2(p6,p3,p5,p4)
     &     + 1./(za(p2,j7))/(za(p6,j7))/(zb(p6,j8))*za(p3,p5)*zab2(p6,p2
     &    ,j7,p1)*zab2(j8,p3,p5,p4)

      amp =  + prop3457**(-1) * (
     &     - 1./(za(p1,j8))/(za(p6,j8))/(zb(p2,j7))*za(p1,p6)*za(p3,p5)*
     &    za(p6,j7)*zb(p1,p2)*zb(p2,p4)*mt**2
     &     + 1./(za(p1,j8))/(za(p6,j8))/(zb(p2,j7))*za(p1,p6)*za(p3,p5)*
     &    zb(p1,p2)*zab2(p6,p1,j8,p2)*zab2(j7,p3,p5,p4)
     &     - 1./(za(p1,j8))/(zb(p2,j7))*za(p3,p5)*za(p6,j7)*zb(p2,p4)*
     &    zb(p2,j8)*mt**2
     &     + 1./(za(p1,j8))/(zb(p2,j7))*za(p3,p5)*zb(p2,j8)*zab2(p6,p1,
     &    j8,p2)*zab2(j7,p3,p5,p4)
     &     )

      amm =  + prop3457**(-1) * (
     &     + 1./(zb(p1,j8))/(zb(p2,j7))/(zb(p6,j8))*za(p3,p5)*za(p6,j7)*
     &    zb(p1,p2)*zb(p1,p6)*zb(p2,p4)*mt**2
     &     - 1./(zb(p1,j8))/(zb(p2,j7))/(zb(p6,j8))*za(p3,p5)*zb(p1,p2)*
     &    zb(p1,p6)*zab2(p6,p1,j8,p2)*zab2(j7,p3,p5,p4)
     &     - 1./(zb(p2,j7))/(zb(p6,j8))*za(p3,p5)*za(j7,j8)*zb(p1,p2)*
     &    zb(p2,p4)*mt**2
     &     - 1./(zb(p2,j7))/(zb(p6,j8))*za(p3,p5)*zb(p1,p2)*zab2(j7,p3,
     &    p5,p4)*zab2(j8,p1,p6,p2)
     &     )

      msqlightxheavy=msqlightxheavy
     & +(abs(app)**2+abs(apm)**2+abs(amp)**2+abs(amm)**2)/(prop34*prop168)

      enddo

c Overall factors
      msqlightxheavy=msqlightxheavy*V**2*gwsq**4/(mt*twidth)**2

      return
      end

      function msqheavy(p1,p2,p3,p4,p5,p6,p7)

c Matrix element squared for (two) diagrams of the form:

c     p1 ---->------------->------ p6
c                    $
c                    $
c                    $
c                    $   t
c     p2 ---->---------||-->------ p5     (on-shell top, p345^2 = mt^2)
c               o          $
c                o          $ -->-- p3
c                 o          |
c                   p7       ^
c                            | p4

c with p7 attached to the line shown

c All overall factors are included except for averaging over the initial state
c and a strong coupling factor of gs^2

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      integer:: p1,p2,p3,p4,p5,p6,p7,j7
      real(dp):: msqheavy,s16,s34,s3457,prop34,prop16,prop3457
      complex(dp):: ap,am,zab2,zab3

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)

      s34=s(p3,p4)
      prop34=(s34-wmass**2)**2
      if (s34 > 0._dp) then
        prop34=prop34+(wmass*wwidth)**2
      endif

      j7=p7

      s16=s(p1,p6)
      prop16=(s16-wmass**2)**2
      if (s16 > 0._dp) then
        prop16=prop16+(wmass*wwidth)**2
      endif

c No width necessary for this top propagator since we assume s345 = mt^2
      s3457=s(p3,p4)+s(p3,p5)+s(p3,j7)+s(p4,p5)+s(p4,j7)+s(p5,j7)
      prop3457=s3457-mt**2

      ap =  + prop3457**(-1) * (
     &     - 1./(za(p6,j7))*za(p3,p5)*zb(p1,p2)*zab2(p6,p3,p5,p4)*zab3(
     &    p6,p3,p4,p5,j7)
     &     )
      ap = ap + 1./(za(p2,j7))/(za(p6,j7))*za(p3,p5)*zab2(p6,p2,j7,
     & p1)*zab2(p6,p3,p5,p4)

      am =  + prop3457**(-1) * (
     &     - 1./(zb(p2,j7))*za(p3,p5)*za(p6,p1)*zb(p1,p2)**2*zab2(j7,p3,
     &    p5,p4)
     &     + 1./(zb(p2,j7))*za(p3,p5)*za(p6,j7)*zb(p1,p2)*zb(p2,p4)*
     &    mt**2
     &     )

      msqheavy=(abs(ap)**2+abs(am)**2)/(prop34*prop16)

c Overall factors
      msqheavy=msqheavy*xn*V*gwsq**4/(mt*twidth)**2

      return
      end

      subroutine coefs_new(s12,mtsq,cv0,cv,c1,musq, epinv_in, epinv2_in)
      implicit none
      include 'types.f'
c-----In this routine:-
c-----cv0 is the results for all massless vertex function
c-----cv and c1 are  is the results for one-mass vertex function

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

      real(dp), optional, intent(in) :: epinv_in, epinv2_in
      real(dp), intent(in) :: musq

      real(dp) :: epinv_use, epinv2_use

      if (present(epinv_in)) then
          epinv_use = epinv_in
      else
          epinv_use = epinv
      endif

      if (present(epinv2_in)) then
          epinv2_use = epinv2_in
      else
          epinv2_use = epinv2
      endif

      if (scheme =='dred') then
c------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
c------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

c**********************************************************************
c   Massless case
c   Taken from
c   %\cite{Altarelli:1979ub}
c   \bibitem{Altarelli:1979ub}
c   G.~Altarelli, R.~K.~Ellis and G.~Martinelli,
c   %``Large Perturbative Corrections To The Drell-Yan Process In QCD,''
c   Nucl.\ Phys.\ B {\bf 157}, 461 (1979).
c   %%CITATION = NUPHA,B157,461;%%
c   Using Eqn(58) with normalization changed to
c   as/2/pi*cf*(4*pi)^ep/Gamma(1-ep)
c   Taking account that Gamma(1-ep)^2/Gamma(1-2*ep)=1-ep^2*pi^2/6
c**********************************************************************
      xl12=lnrat(-s12,musq)

c-----2/22/2012
c-----This appears to be the correction to the vertex in units of as/4/pi*cf
c-----despite the above comment
      cv0=-2._dp*(epinv2_use*epinv_use-epinv_use*real(xl12))-real(xl12**2)
     &           -3._dp*(epinv_use-real(xl12))-7._dp-eta

c---- this routine has been constructed following closely
c---- the notation of
c---- %\cite{Gottschalk:1980rv}
c---- \bibitem{Gottschalk:1980rv}
c---- T.~Gottschalk,
c---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
c---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
c---- %%CITATION = PHRVA,D23,56;%%
c----- Adapted from Eqs.(A8,A9)

c NB  s12=-Q^2, -taucs=mtsq+Q^2
      taucs=s12-mtsq
      la=-s12/(mtsq-s12)
      oml=1._dp-la
c-----oml=mtsq/(mtsq-s12)
      logoml=-lnrat(-taucs,mtsq)
      logsca=lnrat(-taucs,musq)
      Kfun=cplx1(oml/la)*logoml

c--- Minus sign relative to Gottschalk since incoming b has momentum
c--- vector reversed for the t-channel process
c--- s-channel process follows by crossing
      c1=-ctwo*Kfun

      if (la < 1._dp) then
      Li2la=ddilog(la)
      else
      logla=lnrat(-s12,-taucs)
      Li2la=pisqo6-ddilog(oml)-real(logla*logoml)
      endif
c-----Again from A8 and A9 these are in units of alpha_s/4/pi*CF
      cv=-epinv2_use*epinv_use
     & -epinv_use*(2.5_dp+real(logoml-logsca))
     & -0.5_dp*(11._dp+eta)-pisqo6+2._dp*Li2la-real(Kfun)
     &  -0.5_dp*real(logoml*(cone-logoml))
     &  +2.5_dp*real(logsca)+real(logsca*logoml)-0.5_dp*real(logsca**2)

c answer from gotts.mac
c   ans:
c   -1/ep^2;
c   -2.5/ep -(log(oml)/ep-+log(sca)/ep)
c  -6-pisqo6+2*Li2-Ka
c  -0.5*log(oml)*(1-log(oml))
c  +2.5*log(sca)+log(oml)*log(sca)-log(sca)^2/2
      return

      end subroutine coefs_new

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

      ! only light line corrections
      cv = 0._dp
      c1 = 0._dp

      amp=za(jc,jn)*zb(ju,jb)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

      ampho=za(jc,jn)*zb(ju,jb)
     & *(cplx1(cv0+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & +c1*chalf*zb(je,jb)*za(jb,jd))

      virtqqb_light=real(amp*conjg(ampho))/prop

      end

       function streal_lightResonant_MMMM_P_virt(ju,jb,jn,je,jc,jd,jg,za,zb,s12,musq)
           use constants
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_P_virt
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
           real(dp), intent(in) :: s12,musq

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           real(dp) :: cv,cv0
           complex(dp) :: c1

           mtsq = mt**2 - im*mt*twidth

           call coefs_new(s12,mt**2,cv0,cv,c1,musq)

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_P_virt =  -(c1*za(jb,jd)*za(jn,jc)*zb(je,jb)*
     &      (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)))/
     &   (2._dp*za(jd,jg)*za(ju,jg)) +
     &  (cv*za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*
     &     (za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) +
     &       za(jd,jg)*zb(jg,je)))/(za(jd,jg)*za(ju,jg))

            streal_lightResonant_MMMM_P_virt = streal_lightResonant_MMMM_P_virt *
     &          propW34 * propW167 * propT1267


       end function streal_lightResonant_MMMM_P_virt

       function streal_lightResonant_MMMM_P_tree(ju,jb,jn,je,jc,jd,jg, za,zb)
           use constants
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_P_tree
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_P_tree =
     &  (za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*
     &     (za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) +
     &       za(jd,jg)*zb(jg,je)))/(za(jd,jg)*za(ju,jg))

            streal_lightResonant_MMMM_P_tree = streal_lightResonant_MMMM_P_tree *
     &          propW34 * propW167 * propT1267

       end function streal_lightResonant_MMMM_P_tree

       function streal_lightResonant_MMMM_M_virt(ju,jb,jn,je,jc,jd,jg, za,zb,s12,musq)
           use constants
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_M_virt
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
           real(dp), intent(in) :: s12,musq

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           real(dp) :: cv,cv0
           complex(dp) :: c1

           mtsq = mt**2 - im*mt*twidth

           call coefs_new(s12,mt**2,cv0,cv,c1,musq)

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_M_virt = -(c1*za(jn,jc)*zb(jb,ju)*zb(je,jb)*
     &      (za(jb,jg)*zb(jg,je)*zb(jg,ju) +
     &        za(jb,jd)*(zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju))))/
     &   (2._dp*zb(jg,jd)*zb(jg,je)*zb(jg,ju)) +
     &  (cv*za(jn,jc)*zb(jb,ju)*
     &     (za(jb,jd)*zb(je,jb)*
     &        (zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) +
     &       za(ju,jd)*zb(je,ju)*
     &        (zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) +
     &       zb(jg,je)*(za(jd,jg)*zb(je,ju)*zb(jg,jd) +
     &          (za(jb,jg)*zb(je,jb) + za(ju,jg)*zb(je,ju))*zb(jg,ju))))
     &    /(zb(jg,jd)*zb(jg,je)*zb(jg,ju))

            streal_lightResonant_MMMM_M_virt = streal_lightResonant_MMMM_M_virt *
     &          propW34 * propW167 * propT1267

       end function streal_lightResonant_MMMM_M_virt

       function streal_lightResonant_MMMM_M_tree(ju,jb,jn,je,jc,jd,jg, za,zb)
           use constants
           implicit none
           include 'nf.f'
           include 'mxpart.f'
           include 'masses.f'
           complex(dp) :: streal_lightResonant_MMMM_M_tree
           complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
           integer, intent(in) :: ju,jb,jn,je,jc,jd,jg

           integer :: j,k
           real(dp) :: s
           s(j,k) = real(za(j,k)*zb(k,j))

           complex(dp) :: propT1267, propW34
           real(dp) :: propW167

           complex(dp) :: mtsq

           mtsq = mt**2 - im*mt*twidth

            propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
            propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
            propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

           streal_lightResonant_MMMM_M_tree =
     &  (za(jn,jc)*zb(jb,ju)*
     &     (za(jb,jd)*zb(je,jb)*
     &        (zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) +
     &       za(ju,jd)*zb(je,ju)*
     &        (zb(je,ju)*zb(jg,jd) + zb(jd,je)*zb(jg,ju)) +
     &       zb(jg,je)*(za(jd,jg)*zb(je,ju)*zb(jg,jd) +
     &          (za(jb,jg)*zb(je,jb) + za(ju,jg)*zb(je,ju))*zb(jg,ju))))
     &    /(zb(jg,jd)*zb(jg,je)*zb(jg,ju))

            streal_lightResonant_MMMM_M_tree = streal_lightResonant_MMMM_M_tree *
     &          propW34 * propW167 * propT1267

       end function streal_lightResonant_MMMM_M_tree

      end module
