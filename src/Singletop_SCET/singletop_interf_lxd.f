!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module singletop_interf_lxd
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use LHAPDF

      private

      public :: singletop_light_decay_vv
      ! DONE:
      ! * take lumxmsq_singletop_decay_jetmass, replace hard function
      !   with nlo vertex corrections from singletop2_scet_virt_heavy_decay_all.
      ! * Check that it agrees with previous implementation.
      ! * Attach an additional light correction vertex factor of cv0.

c     ! for splitting pieces
      public :: singletop_light_decay_vv_tree
      ! DONE:
      ! * singletop_light_decay_vv without additional factors of cv0
      ! * up to the interface this should be just equivalent to the decay below
      !   cut

      public :: singletop_light_decay_rr
      ! * checked double real amplitudes

      public :: singletop_light_decay_rr_gs
      ! * dipole routine can be taken from light line real emission _gs
      !   but tree process has to be replaced with decay real corrections
      ! * make sure that 7 and 8 are treated correctly in dipole routines

      public :: singletop_light_decay_vr
      ! DONE:
      ! * attach additional factor cv0 to decay real emission routines

      public :: singletop_light_decay_vr_z
      ! DONE:
      ! * taken from light line virtual _z

      public :: singletop_light_decay_rv
      ! DONE:
      ! * requires new calculation where light line current with real emission
      !   is contracted with the rest.

      !public :: singletop_light_decay_rv_tree ! for dipoles in_gs routine
      ! DONE:
      ! * should be the same as _vv_tree, just below cut with nothing on light
      !   line

      public :: singletop_light_decay_rv_gs
      ! DONE:
      ! * dipole routine, like _gs for light line corrections but
      !   with tree level as in _rv_tree

c     public :: passed_taucut_lxd
      ! DONE:
      ! * should correspond to taucut routine from decay

      ! like singletop2_heavy_decay_g_all, but
      ! swaps "correction on beam labels" in msqall
      public :: singletop_decay_real_swap

      public :: coefs_new, coefsdk, virtqqbdk

      public :: assemble_decay_pieces

      contains

      subroutine singletop_decay_real_swap(p,msq)
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

      fac=2._dp*(4*pi*as_heavy_beam2)*cf*aveqq*gw**8*xn**2
      ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
      ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
      msq([2,4],5, 1,1) = ub
      msq([-1,-3],5, 1,1) = ubarb

      fac=2._dp*(4*pi*as_heavy_beam1)*cf*aveqq*gw**8*xn**2
      bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
      bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
      msq(5,[2,4], 1,2) = bu
      msq(5,[-1,-3], 1,2) = bubar

      end subroutine singletop_decay_real_swap

      subroutine singletop_light_decay_rv_gs(p,ndmx,msq)
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

          real(dp) :: dummyv(-nf:nf,-nf:nf), dsubv
          real(dp) :: sub(4), msqx(-nf:nf,-nf:nf)

          integer :: noglue(4) = [-1,-3,2,4]

          external donothing_gvec

          ndmax = 8

          msq = 0._dp
          sub = 0._dp

          corr_islight = .true.

          corr_on_beam = 1
          corr_beam1 = .true.

          call dips(1,p, 1,7,6,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(1,noglue,5, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

          call dips(2,p, 6,7,1,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(2,noglue,5, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

          call dips(3,p, 1,7,6,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(3,0,5, 1) = 2._dp*tr*sub(qg) * sum(msqx(1:5,5))

          call dips(4,p, 1,6,7,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(4,0,5, 1) = 2._dp*tr*sub(qg) * sum(msqx(-5:-1, 5))

          corr_on_beam = 2
          corr_beam1 = .false.

          call dips(5,p, 2,7,6,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(5,5,noglue, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

          call dips(6,p, 6,7,2,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(6,5,noglue, 1) = 2._dp*cf*sub(qq)*msqx(5,noglue)

          call dips(7,p, 2,7,6,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(7,5,0, 1) = 2._dp*tr*sub(qg) * sum(msqx(5,1:5))

          call dips(8,p, 2,6,7,sub,dsubv,msqx,dummyv,singletop_light_decay_vv_tree,donothing_gvec)
          msq(8,5,0, 1) = 2._dp*tr*sub(qg) * sum(msqx(5,-5:-1))

      end subroutine singletop_light_decay_rv_gs


      function ubtdg_l(ju,jb,jn,je,jc,jd,jg,p)
      implicit none
      include 'types.f'
      real(dp):: ubtdg_l
c     Matrix element squared for single top production with gluon
c     radiation in production
c      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc))+d(jd)+g(jg)
c     masses of b quarks c.c=b.b=0

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: ju,jb,jn,je,jc,jd,jg,nu
      real(dp):: p(mxpart,4),pt(4),ptDpt
      real(dp):: sne,sdug,prop
      complex(dp):: ampi(2)

      do nu=1,4
      pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2

      sne=s(jn,je)
      sdug=s(jd,ju)+s(jd,jg)+s(ju,jg)

      if (sdug < 0._dp) then
      prop=(sdug-wmass**2)**2
      else
      prop=(sdug-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     &    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop

c  -Lefthanded gluon
      ampi(1)=(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & *zb(ju,jd)/(zb(jg,ju)*zb(jg,jd))
     & -(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))/zb(jg,jd)
      ampi(1)=ampi(1)*za(jc,jn)*zb(ju,jb)

c  -Righthanded gluon
      ampi(2)=za(ju,jd)/(za(ju,jg)*za(jd,jg))
     & -zb(jg,jb)/(zb(ju,jb)*za(ju,jg))
      ampi(2)=-ampi(2)*za(jc,jn)*zb(ju,jb)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))


      ubtdg_l=(abs(ampi(1))**2+abs(ampi(2))**2)/prop
      return
      end

      function st_light_decay_rv_MMMM_M(ju,jb,jn,je,jc,jd,jg, za,zb, cv,c1)
          use constants
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          real(dp) :: st_light_decay_rv_MMMM_M
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: cv,c1

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT1267, propW34
          real(dp) :: propW167

          complex(dp) :: mtsq

          complex(dp) :: treeamp, c1amp

          mtsq = mt**2 - im*mt*twidth

          propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
          propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
          propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

          c1amp = -(c1*za(jn,jc)*zb(jb,ju)*zb(jc,je)*
     &      (za(jc,jd)*zb(jd,ju) + za(jc,jg)*zb(jg,ju)))/
     &   (2._dp*zb(jg,jd)*zb(jg,ju))

          treeamp = -(za(jn,jc)*zb(jb,ju)*
     &     (za(jb,jd)*zb(jd,ju)*zb(je,jb) +
     &       za(ju,jd)*zb(jd,ju)*zb(je,ju) +
     &       za(jd,jg)*zb(jd,ju)*zb(jg,je) -
     &       za(jd,jg)*zb(jd,je)*zb(jg,ju) +
     &       za(jb,jg)*zb(je,jb)*zb(jg,ju) +
     &       za(ju,jg)*zb(je,ju)*zb(jg,ju)))/(zb(jg,jd)*zb(jg,ju))

         st_light_decay_rv_MMMM_M = real((c1amp+cv*treeamp)*conjg(treeamp)) *
     &          abs(propW34 * propW167 * propT1267)**2

      end function st_light_decay_rv_MMMM_M

      function st_light_decay_rv_MMMM_P(ju,jb,jn,je,jc,jd,jg, za,zb, cv,c1)
          use constants
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'masses.f'
          real(dp) :: st_light_decay_rv_MMMM_P
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: ju,jb,jn,je,jc,jd,jg
          real(dp), intent(in) :: cv, c1

          integer :: j,k
          real(dp) :: s
          s(j,k) = real(za(j,k)*zb(k,j))

          complex(dp) :: propT1267, propW34
          real(dp) :: propW167

          complex(dp) :: mtsq

          complex(dp) :: c1amp, treeamp

          mtsq = mt**2 - im*mt*twidth

          propW34  = 1._dp / (s(jn,je) - wmass**2 + im*wmass*wwidth)
          propW167  = 1._dp / (s(ju,jd)+s(ju,jg)+s(jd,jg) - wmass**2)
          propT1267 = 1._dp / (s(jc,je)+s(jc,jn)+s(je,jn) - mtsq)

          c1amp = -(c1*za(jc,jd)*za(jn,jc)*zb(jc,je)*
     &      (za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb)))/
     &   (2._dp*za(jd,jg)*za(ju,jg))
          treeamp = -(za(jn,jc)*(za(ju,jd)*zb(jb,ju) + za(jd,jg)*zb(jg,jb))*
     &     (za(jb,jd)*zb(je,jb) + za(ju,jd)*zb(je,ju) +
     &       za(jd,jg)*zb(jg,je)))/(za(jd,jg)*za(ju,jg))

         st_light_decay_rv_MMMM_P = real((c1amp+cv*treeamp)*conjg(treeamp))*
     &          abs(propW34 * propW167 * propT1267)**2

      end function st_light_decay_rv_MMMM_P

      subroutine singletop_light_decay_rv(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_scet_heavy_decay, only: softfun, jetfun
          use singletop_jet, only: ampsq_ugd_tdkb
          use singletop2_nnlo_vars
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
          include 'sprods_com.f'
          include 'energy.f'
          include 'scheme.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib)

          real(dp) :: xx(2)

          real(dp) :: jet(2,0:4), soft(2,0:4)
          real(dp) :: hard_ub(0:1), hard_bu(0:1), hard_ubarb(0:1), hard_bubar(0:1)
          real(dp) :: hard_gb(0:1), hard_bg(0:1)

          real(dp) :: x
          real(dp) :: puremass
          real(dp) :: fac
          real(dp), parameter :: genfac = 2._dp

          real(dp) :: cv, c1

          xx(1) = -2._dp*p(1,4)/sqrts
          xx(2) = -2._dp*p(2,4)/sqrts

          x = puremass(p(3,:)+p(4,:))**2 / mt**2

          soft = softfun(1)
          jet = jetfun(1)

          call spinoru(7,p,za,zb)

          scheme = 'tH-V'


          ! tree amplitudes: just real emission on light line
          fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
          hard_ub(0) = aveqq*fac*ampsq_ugd_tdkb(1,2,3,4,5,6,7,p)
          hard_ubarb(0) = aveqq*fac*ampsq_ugd_tdkb(6,2,3,4,5,1,7,p)
          hard_gb(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,2,3,4,5,6,1,p)

          fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
          hard_bu(0) = aveqq*fac*ampsq_ugd_tdkb(2,1,3,4,5,6,7,p)
          hard_bubar(0) = aveqq*fac*ampsq_ugd_tdkb(6,1,3,4,5,2,7,p)
          hard_bg(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,1,3,4,5,6,2,p)

          fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
          fac = fac*(as_heavy_beam2/2._dp/pi*cf)
c         ! proper hard function normalization
          fac = fac / (as_heavy_beam2/4._dp/pi)

          call coefsdk(s(3,4),mt**2,cv,c1,renscale_beam2_isheavy_onlight**2,0._dp,0._dp)

          hard_ub(1) = aveqq*fac*(st_light_decay_rv_MMMM_M(1,2,3,4,5,6,7,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(1,2,3,4,5,6,7,za,zb,cv,c1))

          hard_ubarb(1) = aveqq*fac*(st_light_decay_rv_MMMM_M(6,2,3,4,5,1,7,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(6,2,3,4,5,1,7,za,zb,cv,c1))

          hard_gb(1) = genfac*aveqg*fac*(st_light_decay_rv_MMMM_M(7,2,3,4,5,6,1,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(7,2,3,4,5,6,1,za,zb,cv,c1))

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)
          hard_gb(1) = hard_gb(1) + pi**2/6._dp * cf * hard_gb(0)


          fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
          fac = fac*(as_heavy_beam1/2._dp/pi*cf)
c         ! proper hard function normalization
          fac = fac / (as_heavy_beam1/4._dp/pi)

          call coefsdk(s(3,4),mt**2,cv,c1,renscale_beam1_isheavy_onlight**2,0._dp,0._dp)

          hard_bu(1) = aveqq*fac*(st_light_decay_rv_MMMM_M(2,1,3,4,5,6,7,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(2,1,3,4,5,6,7,za,zb,cv,c1))

          hard_bubar(1) = aveqq*fac*(st_light_decay_rv_MMMM_M(6,1,3,4,5,2,7,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(6,1,3,4,5,2,7,za,zb,cv,c1))

          hard_bg(1) = genfac*aveqg*fac*(st_light_decay_rv_MMMM_M(7,1,3,4,5,6,2,za,zb,cv,c1) +
     &                      st_light_decay_rv_MMMM_P(7,1,3,4,5,6,2,za,zb,cv,c1))

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)
          hard_bg(1) = hard_bg(1) + pi**2/6._dp * cf * hard_bg(0)



c         msq([2,4],5, 1) = hard_ub(0)
c         msq([-3,-1],5, 1) = hard_ubarb(0)
c         msq(0,5, 1) = hard_gb(0)

c         msq(5,[2,4], 1) = hard_bu(0)
c         msq(5,[-3,-1], 1) = hard_bubar(0)
c         msq(5,0, 1) = hard_bg(0)


          msq([2,4],5, 1) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)

          msq([-3,-1],5, 1) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)

          msq(0,5, 1) = hard_gb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onheavy, as_heavy_beam2, taucut,
     &        hard_gb/hard_gb(0), jet, soft)


          msq(5,[2,4], 1) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)

          msq(5,[-3,-1], 1) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)

          msq(5,0, 1) = hard_bg(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onheavy, as_heavy_beam1, taucut,
     &        hard_bg/hard_bg(0), jet, soft)

      end subroutine singletop_light_decay_rv

      subroutine singletop_decay_real(p,msq)
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

      fac=2._dp*(4*pi*as_heavy_beam2)*cf*aveqq*gw**8*xn**2
      ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
      ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
      msq([2,4],5) = ub
      msq([-1,-3],5) = ubarb

      fac=2._dp*(4*pi*as_heavy_beam1)*cf*aveqq*gw**8*xn**2
      bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
      bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
      msq(5,[2,4]) = bu
      msq(5,[-1,-3]) = bubar

      end subroutine singletop_decay_real

      subroutine singletop_light_decay_rr_gs(p,ndmx,msq)
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

          ndmax = 8

          msq = 0._dp
          sub = 0._dp

          corr_islight = .true.

          corr_on_beam = 1
          corr_beam1 = .true.

          call dips(1,p, 1,8,6,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(1,noglue,5, 1, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

          call dips(2,p, 6,8,1,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(2,noglue,5, 1, 1) = 2._dp*cf*sub(qq)*msqx(noglue,5)

          call dips(3,p, 1,8,6,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(3,0,5, 1, 1) = 2._dp*tr*sub(qg) * sum(msqx(1:5,5))

          call dips(4,p, 1,6,8,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(4,0,5, 1, 1) = 2._dp*tr*sub(qg) * sum(msqx(-5:-1, 5))

          corr_on_beam = 2
          corr_beam1 = .false.

          call dips(5,p, 2,8,6,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(5,5,noglue, 1, 2) = 2._dp*cf*sub(qq)*msqx(5,noglue)

          call dips(6,p, 6,8,2,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(6,5,noglue, 1, 2) = 2._dp*cf*sub(qq)*msqx(5,noglue)

          call dips(7,p, 2,8,6,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(7,5,0, 1, 2) = 2._dp*tr*sub(qg) * sum(msqx(5,1:5))

          call dips(8,p, 2,6,8,sub,dsubv,msqx,dummyv,singletop_decay_real,donothing_gvec)
          msq(8,5,0, 1, 2) = 2._dp*tr*sub(qg) * sum(msqx(5,-5:-1))

      end subroutine singletop_light_decay_rr_gs

      subroutine singletop_light_decay_vv_tree(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_scet_heavy_decay, only: softfun, jetfun
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
          fac=(as_heavy_beam2/2._dp/pi*cf)*aveqq*gw**8*xn**2
          ! proper hard function normalization
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1) = fac*virtqqbdk(p, 1,2,3,4,5,6,renscale_beam2_isheavy_onlight**2)
          hard_ubarb(1) = fac*virtqqbdk(p, 6,2,3,4,5,1,renscale_beam2_isheavy_onlight**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac=(as_heavy_beam1/2._dp/pi*cf)*aveqq*gw**8*xn**2
          ! proper hard function normalization
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqbdk(p, 2,1,3,4,5,6,renscale_beam1_isheavy_onlight**2)
          hard_bubar(1) = fac*virtqqbdk(p, 6,1,3,4,5,2,renscale_beam1_isheavy_onlight**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)


c         msq([2,4],5) = hard_ub(0)
c         msq([-3,-1],5) = hard_ubarb(0)
c         msq(5,[2,4]) = hard_bu(0)
c         msq(5,[-3,-1]) = hard_bubar(0)

c         msq([2,4],5) = hard_ub(0)
c         msq([-3,-1],5) = hard_ubarb(0)
c         msq(5,[2,4]) = hard_bu(0)
c         msq(5,[-3,-1]) = hard_bubar(0)


          msq([2,4],5) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)

          msq([-3,-1],5) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)

          msq(5,[2,4]) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)

          msq(5,[-3,-1]) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)

      end subroutine singletop_light_decay_vv_tree

      subroutine singletop_light_decay_vv(p,msq)
          use types
          use SCET
          use SCET_Jet
          use singletop2_scet_heavy_decay, only: softfun, jetfun
          use singletop_interf_lxh, only: virtqqb_light
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
          hard_ub(0) = fac*virtqqb_light(1,2,3,4,5,6, renscale_beam1_islight_onlight**2)
          hard_ubarb(0) = fac*virtqqb_light(6,2,3,4,5,1, renscale_beam1_islight_onlight**2)

          fac = (as_light_beam2/2._dp/pi*cf) * aveqq*gw**8*xn**2
          hard_bu(0) = fac*virtqqb_light(2,1,3,4,5,6, renscale_beam2_islight_onlight**2)
          hard_bubar(0) = fac*virtqqb_light(6,1,3,4,5,2, renscale_beam2_islight_onlight**2)

          ! hard(1) virtual corrections in decay (msbar renormalized: epinv = 0
          ! and pi^2/6 * cf piece below) with cv0 factor multiplied (for virtual
          ! corrections on light line)

          ! properly normalized one-loop pieces
          fac = aveqq*gw**8*xn**2
          fac = fac * (as_light_beam1/2._dp/pi*cf) * (as_heavy_beam2/2._dp/pi*cf)
          ! proper hard function normalization
          fac = fac / (as_heavy_beam2/4._dp/pi)
          hard_ub(1) = fac*virtqqbdk_wlight(p, 1,2,3,4,5,6,
     &                      renscale_beam1_islight_onlight**2,
     &                      renscale_beam2_isheavy_onlight**2)
          hard_ubarb(1) = fac*virtqqbdk_wlight(p, 6,2,3,4,5,1,
     &                      renscale_beam1_islight_onlight**2,
     &                      renscale_beam2_isheavy_onlight**2)

          hard_ub(1) = hard_ub(1) + pi**2/6._dp * cf * hard_ub(0)
          hard_ubarb(1) = hard_ubarb(1) + pi**2/6._dp * cf * hard_ubarb(0)

          fac = aveqq*gw**8*xn**2
          fac = fac * (as_light_beam2/2._dp/pi*cf) * (as_heavy_beam1/2._dp/pi*cf)
          ! proper hard function normalization
          fac = fac / (as_heavy_beam1/4._dp/pi)
          hard_bu(1) = fac*virtqqbdk_wlight(p, 2,1,3,4,5,6,
     &                      renscale_beam2_islight_onlight**2,
     &                      renscale_beam1_isheavy_onlight**2)

          hard_bubar(1) = fac*virtqqbdk_wlight(p, 6,1,3,4,5,2,
     &                      renscale_beam2_islight_onlight**2,
     &                      renscale_beam1_isheavy_onlight**2)

          hard_bu(1) = hard_bu(1) + pi**2/6._dp * cf * hard_bu(0)
          hard_bubar(1) = hard_bubar(1) + pi**2/6._dp * cf * hard_bubar(0)

          call fdist(ih1,xx(1),facscale_beam1_islight_onlight,beama0,1)
          call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight,beamb0,2)

          msq([2,4],5) = hard_ub(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ub/hard_ub(0), jet, soft)*[beama0(2),beama0(4)]*beamb0(5)

          msq([-3,-1],5) = hard_ubarb(0)*assemble_decay_pieces(x,
     &        renscale_beam2_isheavy_onlight, as_heavy_beam2, taucut,
     &        hard_ubarb/hard_ubarb(0), jet, soft)*[beama0(-3),beama0(-1)]*beamb0(5)

          call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight,beama0,1)
          call fdist(ih2,xx(2),facscale_beam2_islight_onlight,beamb0,2)

          msq(5,[2,4]) = hard_bu(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bu/hard_bu(0), jet, soft)*beama0(5)*[beamb0(2),beamb0(4)]

          msq(5,[-3,-1]) = hard_bubar(0)*assemble_decay_pieces(x,
     &        renscale_beam1_isheavy_onlight, as_heavy_beam1, taucut,
     &        hard_bubar/hard_bubar(0), jet, soft)*beama0(5)*[beamb0(-3),beamb0(-1)]

      end subroutine singletop_light_decay_vv

      subroutine singletop_light_decay_vr(p,msq)
          use singletop2_nnlo_vars
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
      include 'scheme.f'

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp):: fac,ub,bu,bubar,ubarb
      integer:: ib


      call spinoru(7,p,za,zb)
      ib=5*nwz

      scheme = 'tH-V'

      msq = 0._dp

      fac=2._dp*(4*pi*as_heavy_beam2)*cf*aveqq*gw**8*xn**2
      fac = fac * (as_light_beam1/2._dp/pi*cf)
      ub=fac*qqbtbbargd_cv0(1,2,3,4,5,6,7,p, renscale_beam1_islight_onlight**2)
      ubarb=fac*qqbtbbargd_cv0(6,2,3,4,5,1,7,p, renscale_beam1_islight_onlight**2)
      msq([2,4],5, 1,1) = ub
      msq([-1,-3],5, 1,1) = ubarb

      fac=2._dp*(4*pi*as_heavy_beam1)*cf*aveqq*gw**8*xn**2
      fac = fac * (as_light_beam2/2._dp/pi*cf)
      bu=fac*qqbtbbargd_cv0(2,1,3,4,5,6,7,p, renscale_beam2_islight_onlight**2)
      bubar=fac*qqbtbbargd_cv0(6,1,3,4,5,2,7,p, renscale_beam2_islight_onlight**2)
      msq(5,[2,4], 1,2) = bu
      msq(5,[-1,-3], 1,2) = bubar

      end subroutine singletop_light_decay_vr

      subroutine singletop_light_decay_vr_z(p,z)
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
        real(dp) :: dot, if_qq, fi_qq, if_qg

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
            Q1(q,g,q,is) = as_light_beam1/2/pi * tr * if_qg(z, xl16,is)

            ! corr_on_beam = 2
            B2(q,q,b,is) = as_light_beam2/2/pi * cf*(if_qq(z,xl26,is) + fi_qq(z,xl26,is))
            Q2(q,g,q,is) = as_light_beam2/2/pi * tr * if_qg(z, xl26,is)
        enddo

      end subroutine singletop_light_decay_vr_z

      function assemble_decay_pieces(x, scale, as, taucut, h, j, s)
          use types
          use constants
          implicit none
          include 'masses.f'
          include 'kpart.f'

          real(dp), intent(in) :: x, scale, as
          real(dp), intent(in) :: taucut
          real(dp), intent(in) :: h(0:1)
          real(dp), intent(in) :: j(2,0:4)
          real(dp), intent(in) :: s(2,0:4)
          real(dp) :: assemble_decay_pieces

          real(dp) :: LogMtMu
          real(dp) :: full(0:2,0:4)
          real(dp) :: LSX, LJX

          LogMtMu = log(mt/scale)

          LJX = 2*LogMtMu
          LSX = LogMtMu - log(1-x)

          full(:,:) = 0._dp

          if (h(0) /= 1._dp) then
              write (*,*) "WARNING: bad hard function normalization!"
          endif

          full(1,2) = j(1,2) + s(1,2)
          full(1,1) = j(1,1) + 2*LJX*j(1,2) + s(1,1) + 2*LSX*s(1,2)
          full(1,0) = h(1) + j(1,0) + LJX*j(1,1) + LJX**2*j(1,2) + s(1,0) + LSX*s(1,1) +
     &       LSX**2*s(1,2)

          assemble_decay_pieces =
     &            as/4._dp/pi*(
     &                full(1,0) + full(1,1)*log(taucut) + full(1,2)*log(taucut)**2 )

      end function

      function qqbtbbargd_cv0(ju,jb,jn,je,jc,jd,jg,p,musq_light)
      implicit none
      include 'types.f'
      real(dp):: qqbtbbargd_cv0
c       Matrix element squared for single top production with gluon
c       radiation in decay (radiation from final line)
c        u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc)+g(jg))+d(jd)
c       masses of b quarks c.c=b.b=0

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: ju,jb,jn,je,jc,jd,jg,nu
      real(dp):: p(mxpart,4),pt(4),ptDpt
      real(dp):: sne,sdu,prop,twoptg
      complex(dp):: ampf(2),amp0
      real(dp) :: musq_light

      real(dp) :: cv0,cv
      complex(dp) :: c1_prod

      do nu=1,4
          pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)+p(jg,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      twoptg=
     &  2._dp*(pt(4)*p(jg,4)-pt(1)*p(jg,1)-pt(2)*p(jg,2)-pt(3)*p(jg,3))

      sne=s(jn,je)
      sdu=s(jd,ju)
      if (sdu < 0._dp) then
          prop=(sdu-wmass**2)**2
      else
          prop=(sdu-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     &*((ptDpt-mt**2)**2+(mt*twidth)**2)*prop


      ! calculate light factor with epinv dependence and musq_light
      call coefs_new(s(ju,jd),mt**2,cv0,cv,c1_prod,musq_light,epinv,epinv2)


      amp0=(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))
     &*za(jc,jn)*zb(ju,jb)
c---  eikonal form
      ampf(1)=-amp0
     &*(za(jg,je)*zb(je,jc)+za(jg,jn)*zb(jn,jc))/zb(jg,jc)/twoptg
      ampf(1)=ampf(1)-za(jg,jn)*zb(ju,jb)/zb(jg,jc)
     &*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))

c---  eikonal form
      ampf(2)=-amp0
     &*(zb(jg,je)*za(je,jc)+zb(jg,jn)*za(jn,jc))/za(jc,jg)/twoptg
     &-za(jc,jn)*zb(ju,jb)/twoptg*zb(je,jg)
     &*(zb(jg,jc)*za(jc,jd)+zb(jg,je)*za(je,jd)+zb(jg,jn)*za(jn,jd))


c     qqbtbbargd_cv0=(abs(ampf(1))**2+abs(ampf(2))**2)/prop

      qqbtbbargd_cv0 = real(ampf(1)*conjg(ampf(1))*cv0 +
     &                      ampf(2)*conjg(ampf(2))*cv0 )/prop

      return
      end function qqbtbbargd_cv0

      function virtqqbdk_wlight(p, ju,jb,jn,je,jc,jd,musq_light,musq_decay)
      implicit none
      include 'types.f'
      real(dp):: virtqqbdk_wlight


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      real(dp):: snec,prop,mtsq,cv,c1,cv0
      complex(dp):: amp,ampho

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: ju,jd,jn,je,jc,jb
      real(dp), intent(in) :: musq_light, musq_decay

      complex(dp) :: c1_prod

      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      ! calculate light factor with epinv dependence and musq_light
      call coefs_new(s(ju,jd),mtsq,cv0,cv,c1_prod,musq_light,epinv,epinv2)
      ! calculate decay factor MSbar-renormalized with musq_decay
      call coefsdk(s(jn,je),mtsq,cv,c1,musq_decay,0._dp,0._dp)

      if (s(ju,jd) < 0._dp) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      ! note cv0 multiplied here
      amp=za(jc,jn)*zb(ju,jb)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)) * cv0

      ! coefficients are in terms of alphas/4/pi
      ampho=za(jc,jn)*zb(ju,jb)
     & *(cplx1(cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & +cplx1(0.5_dp*c1)*zb(je,jc)*za(jc,jd))

      ! but factor of two is left out here
      virtqqbdk_wlight=real(amp*conjg(ampho))/prop
      return
      end function virtqqbdk_wlight


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
      include 'epinv.f'
      include 'epinv2.f'
      real(dp):: snec,prop,mtsq,cv,c1
      complex(dp):: amp,ampho

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: ju,jd,jn,je,jc,jb
      real(dp), intent(in) :: musq

      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      call coefsdk(s(jn,je),mtsq,cv,c1,musq,0._dp,0._dp)
      if (s(ju,jd) < 0._dp) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      amp=za(jc,jn)*zb(ju,jb)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

      ! coefficients are in terms of alphas/4/pi
      ampho=za(jc,jn)*zb(ju,jb)
     & *(cplx1(cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & +cplx1(0.5_dp*c1)*zb(je,jc)*za(jc,jd))

      ! but factor of two is left out here
      virtqqbdk=real(amp*conjg(ampho))/prop
      return
      end function virtqqbdk


      subroutine coefsdk(s12,mtsq,cv,c1,musq,epinv,epinv2)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scheme.f'
      include 'alfacut.f'
      !include 'includect.f'
      real(dp):: cv,c1,s12,mtsq,ddilog,rsq,omrsq,eta,Kfun,
     & wlog,rlog,mulog,epinv,epinv2

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

c--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967
      mulog=log(musq/mtsq)

c---- this routine has been constructed from
c---- %\cite{Gottschalk:1980rv}
c---- \bibitem{Gottschalk:1980rv}
c---- T.~Gottschalk,
c---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
c---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
c---- %%CITATION = PHRVA,D23,56;%%
c----- Adapted from Eqs.(A8,A9)

      Kfun=1._dp/rsq*wlog
      c1=2._dp*Kfun
      cv=-(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2)
     & -(epinv+mulog)*(2.5_dp-2._dp*wlog)
     & -0.5_dp*(11._dp+eta)-pisqo6-2._dp*ddilog(rsq)
     &  +3._dp*wlog-2._dp*wlog**2-Kfun

      return
      end subroutine coefsdk

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
      end function qqbtbbar


      subroutine singletop_light_decay_rr(p,msqall)
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

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Rcola results for a specific final state
        integer, parameter :: genfac = 2._dp

        call spinoru(8,p,za,zb)

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

c Note: as set up here, gluon 7 is radiated from light line and gluon 8 in decay
c Update Tobias: I switched 7 and 8, so 7 is now in decay and 8 radiated from light line

            if (iand(partons_enabled, quarkChannel) > 0) then
                tmp([2,4], 1) = facqq*msqlightxdecay(i1,i2,3,4,5,6,7,8)   ! u b -> t d g g
                tmp([-1,-3], 1) = facqq*msqlightxdecay(6,i2,3,4,5,i1,7,8) ! d~ b -> t u~ g g
            endif

            if (iand(partons_enabled, gluonChannel) > 0) then
                tmp(0, 1)=facqg*genfac*msqlightxdecay(8,i2,3,4,5,6,7,i1)         ! g b -> t d u~ g
            endif

c update     msq array,
            if (corr_on_beam == 1) then
                msqall(:,5, 1, 1) = tmp(:,1)
            else
                msqall(5,:, 1, 2) = tmp(:,1)
            endif

        enddo

      end subroutine singletop_light_decay_rr


      function msqlightxdecay(p1,p2,p3,p4,p5,p6,p7,p8)

c Matrix element squared for (four) diagrams of the form:

c                           o p8
c                          o
c                         o
c     p1 ---->------------->------ p6
c                    $
c                    $           o p7
c                    $          o
c                    $   t     o
c     p2 ---->---------||-->------ p5     (on-shell top, p3457^2 = mt^2)
c                          $
c                           $ -->-- p3
c                            |
c                            ^
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
      real(dp):: msqlightxdecay,s168,s34,s345,prop34,prop168,prop345
      complex(dp):: app,apm,amp,amm,zab2,zab3

c Statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zab3(p1,p2,p3,p4,p5)=za(p1,p2)*zb(p2,p5)+za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5)

      s34=s(p3,p4)
      prop34=(s34-wmass**2)**2
      if (s34 > 0._dp) then
        prop34=prop34+(wmass*wwidth)**2
      endif

      msqlightxdecay=0._dp

      s168=s(p1,p6)+s(p1,p8)+s(p6,p8)
      prop168=(s168-wmass**2)**2
      if (s168 > 0._dp) then
        prop168=prop168+(wmass*wwidth)**2
      endif

c No width necessary for this top propagator since we assume s3457 = mt^2
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      prop345=s345-mt**2

      app =  + prop345**(-1) * (
     &     + 1/(za(p1,p8))/(za(p5,p7))*za(p3,p5)*za(p5,p6)*zb(p2,p8)*
     &    zb(p4,p7)*mt**2
     &     + 1/(za(p1,p8))/(za(p5,p7))*za(p3,p5)*zb(p2,p8)*zab2(p5,p3,
     &    p5,p4)*zab3(p6,p3,p4,p5,p7)
     &     + 1/(za(p1,p8))/(za(p5,p7))/(za(p6,p8))*za(p1,p6)*za(p3,p5)*
     &    za(p5,p6)*zb(p1,p2)*zb(p4,p7)*mt**2
     &     + 1/(za(p1,p8))/(za(p5,p7))/(za(p6,p8))*za(p1,p6)*za(p3,p5)*
     &    zb(p1,p2)*zab2(p5,p3,p5,p4)*zab3(p6,p3,p4,p5,p7)
     &     )

      apm =  + prop345**(-1) * (
     &     - 1/(za(p5,p7))/(zb(p1,p8))/(zb(p6,p8))*za(p3,p5)*za(p5,p6)*
     &    zb(p1,p2)*zb(p1,p6)*zb(p4,p7)*mt**2
     &     - 1/(za(p5,p7))/(zb(p1,p8))/(zb(p6,p8))*za(p3,p5)*zb(p1,p2)*
     &    zb(p1,p6)*zab2(p5,p3,p5,p4)*zab3(p6,p3,p4,p5,p7)
     &     - 1/(za(p5,p7))/(zb(p6,p8))*za(p3,p5)*za(p5,p8)*zb(p1,p2)*
     &    zb(p4,p7)*mt**2
     &     - 1/(za(p5,p7))/(zb(p6,p8))*za(p3,p5)*zb(p1,p2)*zab2(p5,p3,
     &    p5,p4)*zab3(p8,p3,p4,p5,p7)
     &     )

      amp =  + prop345**(-1) * (
     &     - 1/(za(p1,p8))/(za(p6,p8))/(zb(p4,p7))*za(p1,p6)*za(p3,p5)*
     &    zb(p1,p2)*zab2(p7,p3,p5,p4)*zab3(p6,p3,p5,p7,p4)
     &     - 1/(za(p1,p8))/(zb(p4,p7))*za(p3,p5)*zb(p2,p8)*zab2(p7,p3,
     &    p5,p4)*zab3(p6,p3,p5,p7,p4))
      amp = amp - 1/(za(p1,p8))/(za(p6,p8))/(zb(p4,p7))/(zb(p5,p7))*za(
     & p1,p6)*zb(p1,p2)*zab2(p3,p5,p7,p4)*zab3(p6,p3,p5,p7,p4)
     &     - 1/(za(p1,p8))/(zb(p4,p7))/(zb(p5,p7))*zb(p2,p8)*zab2(p3,p5
     &    ,p7,p4)*zab3(p6,p3,p5,p7,p4)

      amm =  + prop345**(-1) * (
     &     + 1/(zb(p1,p8))/(zb(p4,p7))/(zb(p6,p8))*za(p3,p5)*zb(p1,p2)*
     &    zb(p1,p6)*zab2(p7,p3,p5,p4)*zab3(p6,p3,p5,p7,p4)
     &     + 1/(zb(p4,p7))/(zb(p6,p8))*za(p3,p5)*zb(p1,p2)*zab2(p7,p3,
     &    p5,p4)*zab3(p8,p3,p5,p7,p4)
     &     )
      amm = amm + 1/(zb(p1,p8))/(zb(p4,p7))/(zb(p5,p7))/(zb(p6,p8))*zb(
     & p1,p2)*zb(p1,p6)*zab2(p3,p5,p7,p4)*zab3(p6,p3,p5,p7,p4)
     &     + 1/(zb(p4,p7))/(zb(p5,p7))/(zb(p6,p8))*zb(p1,p2)*zab2(p3,p5
     &    ,p7,p4)*zab3(p8,p3,p5,p7,p4)

      msqlightxdecay=msqlightxdecay
     & +(abs(app)**2+abs(apm)**2+abs(amp)**2+abs(amm)**2)/(prop34*prop168)

c Overall factors
      msqlightxdecay=msqlightxdecay*V**2*gwsq**4/(mt*twidth)**2

      return
      end


      end module
