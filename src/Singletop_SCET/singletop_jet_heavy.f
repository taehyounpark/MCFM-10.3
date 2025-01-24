!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module singletop_jet2
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use singletop_jet2_deps
      use singletop_jetdeps

      public :: singletop_jet_heavy
      public :: singletop_jet_heavy_cobswitch
      public :: singletop_jet_heavy_cobswitch_swap
      public :: singletop_jet_heavy_all
      public :: singletop_jet_heavy_virt_all
      public :: singletop_jet_heavy_real_all
      public :: singletop_jet_heavy_gvec
      public :: singletop_jet_heavy_gs_all
      public :: singletop_jet_heavy_z

      public :: Wtoponshell_gen
      public :: Watoponshell_gen

      private

      contains

#define ONLY_MAIN_B 1
#define WITH_B 2
#define WITH_BBAR 3
#define WITH_B_BBAR 4
#define WITH_BBAR_BBAR 5

      subroutine singletop_jet_heavy_all(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (qqb_w_twdk) by K. Ellis

c note that ordering is:
c      u(p1) + b(p2) -> t(W(p34)+b(p5)) + d(p6) + g(p7)
c  and u(p1) + g(p2) -> t(W(p34)+b(p5)) + d(p6) + bbar(p7)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      integer:: ht,hb,h2,hc,i1,i2
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
      real(dp):: fac,msq_qb,msq_qbb,msq_qg,msq_qbg
      complex(dp):: prop
      complex(dp):: mtop(2,2),
     & mqb(2,2,2),mqbb(2,2,2),mqg(2,2,2),mqbg(2,2,2),
     & mtotqb(2,2,2),mtotqbb(2,2,2),mtotqg(2,2,2),mtotqbg(2,2,2)

      integer :: m

c---initialize
      prop=cplx2(zip,mt*twidth)

      if (nwz == +1) then
          msqall = 0._dp

          call tdecay(p,3,4,5,mtop)

          do m=1,maxbeams
              corr_on_beam = beams_enabled(m)

              mtotqb(:,:,:)=czip
              mtotqbb(:,:,:)=czip
              mtotqg(:,:,:)=czip
              mtotqbg(:,:,:)=czip

              if (corr_on_beam == 1) then
                i1=2
                i2=1
              elseif (corr_on_beam == 2) then
                i1=1
                i2=2
              endif

              call Wtoponshell_gen(i2,7,6,i1,p,0,mqb)
              call Wtoponshell_gen(i2,7,i1,6,p,0,mqbb)
              call Wtoponshell_gen(7,i2,6,i1,p,0,mqg)
              call Wtoponshell_gen(7,i2,i1,6,p,0,mqbg)
              !call tdecay(p,3,4,5,mtop)
c              mtop(:,:)=mt*twidth/sqrt(gwsq) ! DEBUG: remove decay for Recola check
              do hb=1,1; do h2=1,2; do hc=1,2; do ht=1,1
                  mtotqb(hb,h2,hc)=mtotqb(hb,h2,hc)+mtop(hb,ht)*mqb(ht,h2,hc)
                  mtotqbb(hb,h2,hc)=mtotqbb(hb,h2,hc)+mtop(hb,ht)*mqbb(ht,h2,hc)
                  mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)+mtop(hb,ht)*mqg(ht,h2,hc)
                  mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)+mtop(hb,ht)*mqbg(ht,h2,hc)
              enddo; enddo; enddo; enddo

              msq_qb = sum(abs(mtotqb)**2)
              msq_qbb = sum(abs(mtotqbb)**2)
              msq_qg = sum(abs(mtotqg)**2)
              msq_qbg = sum(abs(mtotqbg)**2)

              if (corr_on_beam == 2) then
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam2)/abs(prop)**2

                  ! only main b from top decay
                  msqall([2,4],5, 1, corr_on_beam)=aveqq*fac*msq_qb
                  msqall([-1,-3],5, 1, corr_on_beam)=aveqq*fac*msq_qbb

                  ! additional b~ from g splitting
                  msqall([2,4],0, WITH_BBAR, corr_on_beam)=aveqg*fac*msq_qg
                  msqall([-1,-3],0, WITH_BBAR, corr_on_beam)=aveqg*fac*msq_qbg
              else
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam1)/abs(prop)**2

                  msqall(5,[2,4], 1, corr_on_beam)=aveqq*fac*msq_qb
                  msqall(5,[-1,-3], 1, corr_on_beam)=aveqq*fac*msq_qbb

                  msqall(0,[2,4], WITH_BBAR, corr_on_beam)=aveqg*fac*msq_qg
                  msqall(0,[-1,-3], WITH_BBAR, corr_on_beam)=aveqg*fac*msq_qbg
              endif

          enddo
      else
         write(6,*) 'Abort in singletop_jet_heavy_all'
         stop
      endif

      end subroutine singletop_jet_heavy_all

      subroutine singletop_jet_heavy_cobswitch_swap(p,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
      real(dp) :: pswap(mxpart,4)

      pswap(1:5,:)=p(1:5,:)
      pswap(6,:)=p(7,:)
      pswap(7,:)=p(6,:)

      call singletop_jet_heavy_cobswitch(pswap,msq)

      end subroutine singletop_jet_heavy_cobswitch_swap


      subroutine singletop_jet_heavy_cobswitch(p,msq)

      use singletop2_nnlo_vars
      use singletop2_scet_heavy_prod, only: qqb_tbb_g_heavy_all
c      use singletop_interf_lxh, only: msqheavy
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (qqb_w_twdk) by K. Ellis

c note that ordering is:
c      u(p1) + b(p2) -> t(W(p34)+b(p5)) + d(p6) + g(p7)
c  and u(p1) + g(p2) -> t(W(p34)+b(p5)) + d(p6) + bbar(p7)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
c      include 'zprods_com.f'
      integer:: ht,hb,h2,hc,i1,i2
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
      real(dp):: fac,msq_qb,msq_qbb,msq_qg,msq_qbg
      complex(dp):: prop
      complex(dp):: mtop(2,2),
     & mqb(2,2,2),mqbb(2,2,2),mqg(2,2,2),mqbg(2,2,2),
     & mtotqb(2,2,2),mtotqbb(2,2,2),mtotqg(2,2,2),mtotqbg(2,2,2)

c     real(dp) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      integer :: cob

c---initialize
      prop=cplx2(zip,mt*twidth)

      if (nwz == +1) then
          msq = 0._dp

          call tdecay(p,3,4,5,mtop)

          cob = corr_on_beam

          !do cob=1,2
              mtotqb(:,:,:)=czip
              mtotqbb(:,:,:)=czip
              mtotqg(:,:,:)=czip
              mtotqbg(:,:,:)=czip

              if (cob == 1) then
                i1=1
                i2=2
              elseif (cob == 2) then
                i1=2
                i2=1
              endif

              call Wtoponshell_gen(i2,7,6,i1,p,0,mqb)
              call Wtoponshell_gen(i2,7,i1,6,p,0,mqbb)
              call Wtoponshell_gen(7,i2,6,i1,p,0,mqg)
              call Wtoponshell_gen(7,i2,i1,6,p,0,mqbg)
              !call tdecay(p,3,4,5,mtop)
c              mtop(:,:)=mt*twidth/sqrt(gwsq) ! DEBUG: remove decay for Recola check
              do hb=1,1; do h2=1,2; do hc=1,2; do ht=1,1
                  mtotqb(hb,h2,hc)=mtotqb(hb,h2,hc)+mtop(hb,ht)*mqb(ht,h2,hc)
                  mtotqbb(hb,h2,hc)=mtotqbb(hb,h2,hc)+mtop(hb,ht)*mqbb(ht,h2,hc)
                  mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)+mtop(hb,ht)*mqg(ht,h2,hc)
                  mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)+mtop(hb,ht)*mqbg(ht,h2,hc)
              enddo; enddo; enddo; enddo

              msq_qb = sum(abs(mtotqb)**2)
              msq_qbb = sum(abs(mtotqbb)**2)
              msq_qg = sum(abs(mtotqg)**2)
              msq_qbg = sum(abs(mtotqbg)**2)

              if (cob == 1) then
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam2)/abs(prop)**2

                  ! only main b from top decay
                  msq([2,4],5)=aveqq*fac*msq_qb
                  msq([-1,-3],5)=aveqq*fac*msq_qbb

                  ! additional b~ from g splitting
                  msq([2,4],0)=aveqg*fac*msq_qg
                  msq([-1,-3],0)=aveqg*fac*msq_qbg
              else
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam1)/abs(prop)**2

                  msq(5,[2,4])=aveqq*fac*msq_qb
                  msq(5,[-1,-3])=aveqq*fac*msq_qbb

                  msq(0,[2,4])=aveqg*fac*msq_qg
                  msq(0,[-1,-3])=aveqg*fac*msq_qbg
              endif

              ! I am a little worried about the numerical stability of these ...

c             if (cob == 2) then
c                 call qqb_tbb_g_heavy_all(p,msqall)
c                 if ( abs(msq(2,5) / msqall(2,5,1,2) - 1d0) > 1d-3 ) then
c                     write (*,*) "25", msq(2,5) / msqall(2,5, 1,2), msq(2,5), msqall(2,5,1,2)
c                 endif
c                 if ( abs(msq(2,0) / msqall(2,0,3,2) - 1d0) > 1d-3 ) then
c                     write (*,*) "20", msq(2,0) / msqall(2,0, 3,2), msq(2,0), msqall(2,0,3,2)
c                 endif
c             elseif (cob == 1) then
c                 call qqb_tbb_g_heavy_all(p,msqall)
c                 if ( abs(msq(5,2) / msqall(5,2,1,1) - 1d0) > 1d-3 ) then
c                     write (*,*) "52", msq(5,2) / msqall(5,2, 1,1), msq(5,2), msqall(5,2,1,1)
c                 endif
c                 if ( abs(msq(0,2) / msqall(0,2,3,1) - 1d0) > 1d-3 ) then
c                     write (*,*) "02", msq(0,2) / msqall(0,2, 3,1), msq(0,2), msqall(0,2,3,1)
c                 endif
c             endif


          !enddo
      else
         write(6,*) 'Abort in singletop_jet_heavy_cobswitch'
      stop
      endif

      end subroutine singletop_jet_heavy_cobswitch

      subroutine singletop_jet_heavy(p,msq)
      use singletop2_nnlo_vars
      use singletop2_scet_heavy_prod, only: qqb_tbb_g_heavy_all
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (qqb_w_twdk) by K. Ellis

c note that ordering is:
c      u(p1) + b(p2) -> t(W(p34)+b(p5)) + d(p6) + g(p7)
c  and u(p1) + g(p2) -> t(W(p34)+b(p5)) + d(p6) + bbar(p7)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      integer:: ht,hb,h2,hc,i1,i2
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
      real(dp):: fac,msq_qb,msq_qbb,msq_qg,msq_qbg
      complex(dp):: prop
      complex(dp):: mtop(2,2),
     & mqb(2,2,2),mqbb(2,2,2),mqg(2,2,2),mqbg(2,2,2),
     & mtotqb(2,2,2),mtotqbb(2,2,2),mtotqg(2,2,2),mtotqbg(2,2,2)

c     real(dp) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      integer :: cob

c---initialize
      prop=cplx2(zip,mt*twidth)

      if (nwz == +1) then
          msq = 0._dp

          call tdecay(p,3,4,5,mtop)

          cob = corr_on_beam

          !do cob=1,max_corr_on_beam
              mtotqb(:,:,:)=czip
              mtotqbb(:,:,:)=czip
              mtotqg(:,:,:)=czip
              mtotqbg(:,:,:)=czip

              if (cob == 1) then
                i1=2
                i2=1
              elseif (cob == 2) then
                i1=1
                i2=2
              endif

              call Wtoponshell_gen(i2,7,6,i1,p,0,mqb)
              call Wtoponshell_gen(i2,7,i1,6,p,0,mqbb)
              call Wtoponshell_gen(7,i2,6,i1,p,0,mqg)
              call Wtoponshell_gen(7,i2,i1,6,p,0,mqbg)
              !call tdecay(p,3,4,5,mtop)
c              mtop(:,:)=mt*twidth/sqrt(gwsq) ! DEBUG: remove decay for Recola check
              do hb=1,1; do h2=1,2; do hc=1,2; do ht=1,1
                  mtotqb(hb,h2,hc)=mtotqb(hb,h2,hc)+mtop(hb,ht)*mqb(ht,h2,hc)
                  mtotqbb(hb,h2,hc)=mtotqbb(hb,h2,hc)+mtop(hb,ht)*mqbb(ht,h2,hc)
                  mtotqg(hb,h2,hc)=mtotqg(hb,h2,hc)+mtop(hb,ht)*mqg(ht,h2,hc)
                  mtotqbg(hb,h2,hc)=mtotqbg(hb,h2,hc)+mtop(hb,ht)*mqbg(ht,h2,hc)
              enddo; enddo; enddo; enddo

              msq_qb = sum(abs(mtotqb)**2)
              msq_qbb = sum(abs(mtotqbb)**2)
              msq_qg = sum(abs(mtotqg)**2)
              msq_qbg = sum(abs(mtotqbg)**2)

              if (cob == 2) then
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam2)/abs(prop)**2

                  ! only main b from top decay
                  msq([2,4],5)=aveqq*fac*msq_qb
                  msq([-1,-3],5)=aveqq*fac*msq_qbb

                  ! additional b~ from g splitting
                  msq([2,4],0)=aveqg*fac*msq_qg
                  msq([-1,-3],0)=aveqg*fac*msq_qbg
              else
                  fac=V*xn*gwsq**4*(4._dp*pi*as_heavy_beam1)/abs(prop)**2

                  msq(5,[2,4])=aveqq*fac*msq_qb
                  msq(5,[-1,-3])=aveqq*fac*msq_qbb

                  msq(0,[2,4])=aveqg*fac*msq_qg
                  msq(0,[-1,-3])=aveqg*fac*msq_qbg
              endif

              ! I am a little worried about the numerical stability of these ...

c             if (cob == 2) then
c                 call qqb_tbb_g_heavy_all(p,msqall)
c                 if ( abs(msq(2,5) / msqall(2,5,1,2) - 1d0) > 1d-3 ) then
c                     write (*,*) "25", msq(2,5) / msqall(2,5, 1,2), msq(2,5), msqall(2,5,1,2)
c                 endif
c                 if ( abs(msq(2,0) / msqall(2,0,3,2) - 1d0) > 1d-3 ) then
c                     write (*,*) "20", msq(2,0) / msqall(2,0, 3,2), msq(2,0), msqall(2,0,3,2)
c                 endif
c             elseif (cob == 1) then
c                 call qqb_tbb_g_heavy_all(p,msqall)
c                 if ( abs(msq(5,2) / msqall(5,2,1,1) - 1d0) > 1d-3 ) then
c                     write (*,*) "52", msq(5,2) / msqall(5,2, 1,1), msq(5,2), msqall(5,2,1,1)
c                 endif
c                 if ( abs(msq(0,2) / msqall(0,2,3,1) - 1d0) > 1d-3 ) then
c                     write (*,*) "02", msq(0,2) / msqall(0,2, 3,1), msq(0,2), msqall(0,2,3,1)
c                 endif
c             endif


          !enddo
      else
         stop 'Abort in singletop_jet_heavy'
      endif

      end subroutine singletop_jet_heavy

      subroutine singletop_jet_heavy_virt_all(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (qqb_w_twdk_v) by F. Tramontano
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'scheme.f'
c nodecay = .true. returns the squared matrix elements for an undecayed top quark

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      logical, parameter :: nodecay = .false.
      integer:: j,m,ibeam1,ibeam2,i6,i1,i3,i4,iq,ig,ib,k1,k6,iproc
      real(dp) :: q(mxpart,4),fac,r(mxpart,4),wprop,msqbit,gsq
      complex(dp):: ampl0(2),amplv(2),amplp(2,2),ampld(2)
      real(dp):: twotDg,dot
      complex(dp):: wpropfix,spp_ft,spm_ft,
     & smm_ft,smp_ft
      real(dp) :: musq

      scheme = 'tH-V'

      msqall(:,:,:,:)=zero

      do m=1,maxbeams
          corr_on_beam = beams_enabled(m)

      if (corr_on_beam == 1) then
        ibeam1=2
        ibeam2=1
      elseif (corr_on_beam == 2) then
        ibeam1=1
        ibeam2=2
      endif

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        k6=ibeam1
        k1=6
        i3=4
        i4=3
        iq=-1 ! antitop quark
      elseif (nwz == +1) then
        k6=6
        k1=ibeam1
        i3=3
        i4=4
        iq=1 ! top quark
      else
        write(6,*) 'Error in singletop_jet_heavy_virt, nwz is not +1 or -1 :   ',nwz
        stop
      endif

c--- overall factor contained in diag.prc
      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
          fac=(as_heavy_beam1/2._dp/pi)*(V/two)*xn*two*gsq*gwsq**4
      else
          gsq = 4*pi*as_heavy_beam2
          fac=(as_heavy_beam2/2._dp/pi)*(V/two)*xn*two*gsq*gwsq**4
      endif

c--- loop over all crossings
      do iproc=1,4

c--- u + b -> t + d + g
      if     (iproc == 1) then
        ig=7
        ib=ibeam2
        i1=k1
        i6=k6
c--- u + g -> t + d + b~
      elseif (iproc == 2) then
        ig=ibeam2
        ib=7
        i1=k1
        i6=k6
c--- d~ + b -> t + u~ + g
      elseif (iproc == 3) then
        ig=7
        ib=ibeam2
        i1=k6
        i6=k1
c--- d~ + g -> t + u~ + b~
      elseif (iproc == 4) then
        ig=ibeam2
        ib=7
        i1=k6
        i6=k1
      endif

c-- calculate auxiliary momentum array - gq case
      twotDg=two*(dot(p,3,ig)+dot(p,4,ig)+dot(p,5,ig))
      q(1:7,:)=p(1:7,:)
      q(8,:)=p(3,:)+p(4,:)+p(5,:)-p(ig,:)*mt**2/twotDg

c---fill matrices of spinor products
      call spinoru(8,q,za,zb)

c--- Note: call to tree now passes top mass as a parameter
      call tree(mt,ig,ib,i6,i1,8,amplp)
      call wampd(mt,twidth,ig,i3,i4,5,8,ampld)

      r(1,:)=p(ig,:)
      r(2,:)=p(ib,:)
      r(3,:)=p(i6,:)
      r(4,:)=p(i1,:)
      r(5,:)=p(3,:)+p(4,:)+p(5,:)

      if (corr_on_beam == 1) then
          musq = renscale_beam1_isheavy_onheavy**2
      else
          musq = renscale_beam2_isheavy_onheavy**2
      endif

      spp_ft=virt_pp(mt,1,2,3,4,5,r,musq)
      spm_ft=virt_pm(mt,1,2,3,4,5,r,musq)
      smm_ft=virt_mm(mt,1,2,3,4,5,r,musq)
      smp_ft=virt_mp(mt,1,2,3,4,5,r,musq)

      wprop=s(i6,i1)-wmass**2

c This removes the width (not needed for s < 0) that has been applied in tree
      wpropfix=sqrt((s(i6,i1)-wmass**2)**2+(wmass*wwidth)**2)/wprop
      amplp(:,:)=amplp(:,:)*wpropfix

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark, for Born and virtual

      if (nodecay) then
        msqbit=
     &  +real(amplp(1,1)*conjg(smm_ft/wprop))
     &  +real(amplp(2,1)*conjg(spm_ft/wprop))
     &  +real(amplp(1,2)*conjg(smp_ft/wprop))
     &  +real(amplp(2,2)*conjg(spp_ft/wprop))
c        msqbit=(
c     &  +real(amplp(1,1)*conjg(amplp(1,1)))
c     &  +real(amplp(2,1)*conjg(amplp(2,1)))
c     &  +real(amplp(1,2)*conjg(amplp(1,2)))
c     &  +real(amplp(2,2)*conjg(amplp(2,2))))/ason2pi ! To check LO
        msqbit=msqbit*fac/gwsq**2
      else
c--- lowest order amplitudes
        ampl0(1)=amplp(1,1)*ampld(1)
     &          +amplp(2,1)*ampld(2)

        ampl0(2)=amplp(1,2)*ampld(1)
     &          +amplp(2,2)*ampld(2)

c--- virtual amplitudes
        amplv(1)=(smm_ft/wprop)*ampld(1)
     &          +(spm_ft/wprop)*ampld(2)

        amplv(2)=(smp_ft/wprop)*ampld(1)
     &          +(spp_ft/wprop)*ampld(2)

        msqbit=zip
        do j=1,2
        msqbit=msqbit+real(amplv(j)*conjg(ampl0(j)))
c        msqbit=msqbit+real(ampl0(j)*conjg(ampl0(j)))/ason2pi ! To check LO
        enddo
        msqbit=msqbit*fac
      endif

      if (corr_on_beam == 2) then
          select case (iproc)
            case (1)
                msqall([2,4],5, 1,corr_on_beam) = aveqq*msqbit
            case (2)
                ! b~ in final state gets binned into '3'
                msqall([2,4],0, WITH_BBAR,corr_on_beam) = aveqg*msqbit
            case (3)
                msqall([-1,-3],5, 1,corr_on_beam) = aveqq*msqbit
            case (4)
                ! b~ in final state gets binned into '3'
                msqall([-1,-3],0, WITH_BBAR,corr_on_beam) = aveqg*msqbit
            case DEFAULT
                write(6,*) 'Unexpected iproc in singletop_jet_heavy_virt', iproc
                stop
          end select
      else
          select case (iproc)
            case (1)
                msqall(5,[2,4], 1,corr_on_beam) = aveqq*msqbit
            case (2)
                ! b~ in final state gets binned into '3'
                msqall(0,[2,4], WITH_BBAR,corr_on_beam) = aveqg*msqbit
            case (3)
                msqall(5,[-1,-3], 1,corr_on_beam) = aveqq*msqbit
            case (4)
                ! b~ in final state gets binned into '3'
                msqall(0,[-1,-3], WITH_BBAR,corr_on_beam) = aveqg*msqbit
            case DEFAULT
                write(6,*) 'Unexpected iproc in singletop_jet_heavy_virt', iproc
                stop
          end select
      endif

      enddo ! end of loop over crossings

      !if (corr_on_beam == 1) then
        !msqall(:,:,1,corr_on_beam) = transpose(msqall(:,:,1,corr_on_beam))
        !msqall(:,:,3,corr_on_beam) = transpose(msqall(:,:,3,corr_on_beam))
      !endif

      enddo

      return
      end subroutine singletop_jet_heavy_virt_all

      subroutine singletop_jet_heavy_real_all(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'types.f'
c Heavily based on routine for t-channel single top + b production (qg_tbqdk_g)
c***********************************************************************
c     Real MEs for t-channel single top, with explicit b-quark         *
c                                                                      *
c     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5) + g(p6)                  *
c                                                                      *
c      (and related crossings and MEs)                                 *
c                                                                      *
c     Author: J. Campbell, March 18, 2008                              *
c                         (added decay May 2011)                       *
c                                                                      *
c***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'noglue.f'
      include 'nwz.f'
      include 'ewcouple.f'
c nodecay = .true. returns the squared matrix elements for an undecayed top quark

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      logical, parameter :: nodecay = .false.
      real(dp):: q(mxpart,4),msq_qbarg,
     & msq_qg,
     & msq_qb,msq_qbarb,msq_qb_qq,msq_qbarb_qq,
     & msq_qq,msq_qqbar,msq_qbarq,msq_qbarqbar,
     & msq_qb_bb,msq_qbarb_bb,msq_qbbar_bb,msq_qbarbbar_bb
      integer:: k,i1,i2,i6,i7,i8,m


c--- initialize
      msqall(:,:,:,:)=0._dp

      if (nodecay) then
        q(1:2,:)=p(1:2,:)
        q(3,:)=p(3,:)+p(4,:)+p(5,:)
        q(4,:)=p(6,:)
        q(5,:)=p(7,:)
        q(6,:)=p(8,:)
        i6=4
        i7=5
        i8=6
      else
        i6=6
        i7=7
        i8=8
      endif

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1
      do m=1,maxbeams
          corr_on_beam = beams_enabled(m)

      if (corr_on_beam == 2) then
        i1=1
        i2=2
      else
        i1=2
        i2=1
      endif

c     msq_qg = 0._dp
c     msq_qbarg = 0._dp
c     msq_qb = 0._dp
c     msq_qbarb = 0._dp

c     msq_qq = 0._dp
c     msq_qbarq = 0._dp
c     msq_qqbar = 0._dp
c     msq_qbarqbar = 0._dp
c     msq_qb_qq = 0._dp
c     msq_qbarb_qq = 0._dp
c     msq_qb_bb = 0._dp
c     msq_qbarb_bb = 0._dp
c     msq_qbbar_bb = 0._dp
c     msq_qbarbbar_bb = 0._dp

c set mass of quark and antiquark according to nwz
      if (nwz == +1) then
               if (nodecay) then
                  call inter_gen(q,i1,i2,3,i7,i6,i8,mt,msq_qg)      ! u  g -> t d  b~ g
                  call inter_gen(q,i6,i2,3,i7,i1,i8,mt,msq_qbarg)   ! d~ g -> t u~ b~ g
                  call inter_gen(q,i1,i7,3,i2,i6,i8,mt,msq_qb)      ! u  b -> t d  g  g
                  call inter_gen(q,i6,i7,3,i2,i1,i8,mt,msq_qbarb)   ! d~ b -> t u~ g  g
               else
                  call interdk_gen(p,i1,i2,3,4,5,i7,i6,i8,mt,msq_qg)
                  call interdk_gen(p,i6,i2,3,4,5,i7,i1,i8,mt,msq_qbarg)
                  call interdk_gen(p,i1,i7,3,4,5,i2,i6,i8,mt,msq_qb)
                  call interdk_gen(p,i6,i7,3,4,5,i2,i1,i8,mt,msq_qbarb)
               endif
      else
            write(6,*) 'Real radiation not written for nwz = -1'
            stop
      endif

c--- the matrix elements with an extra quark line:
      if (nwz == +1) then
            if (nodecay) then
                  call inter_qq_gen(q,i1,i2,3,i7,i6,i8,mt,msq_qq)         ! u  q  -> t d  b~ q
                  call inter_qq_gen(q,i6,i2,3,i7,i1,i8,mt,msq_qbarq)      ! d~ q  -> t u~ b~ q
                  call inter_qq_gen(q,i1,i8,3,i7,i6,i2,mt,msq_qqbar)      ! u  q~ -> t d  b~ q~
                  call inter_qq_gen(q,i6,i8,3,i7,i1,i2,mt,msq_qbarqbar)   ! d~ q~ -> t u~ b~ q~

                  call inter_qq_gen(q,i1,i7,3,i2,i6,i8,mt,msq_qb_qq)      ! u  b  -> t d  q~ q
                  call inter_qq_gen(q,i6,i7,3,i2,i1,i8,mt,msq_qbarb_qq)   ! d~ b  -> t u~ q~ q

                  call inter_qqid_gen(q,i1,i7,3,i2,i6,i8,mt,msq_qb_bb)      ! u  b  -> t d  b~ b
                  call inter_qqid_gen(q,i6,i7,3,i2,i1,i8,mt,msq_qbarb_bb)   ! d~ b  -> t u~ b~ b
                  call inter_qqid_gen(q,i1,i7,3,i8,i6,i2,mt,msq_qbbar_bb)   ! u  b~ -> t d  b~ b~
                  call inter_qqid_gen(q,i6,i7,3,i8,i1,i2,mt,msq_qbarbbar_bb)! d~ b~ -> t u~ b~ b~

            else
                  call interdk_qq_gen(p,i1,i2,3,4,5,i7,i6,i8,mt,msq_qq)         ! u  q  -> t d  b~ q
                  call interdk_qq_gen(p,i6,i2,3,4,5,i7,i1,i8,mt,msq_qbarq)      ! d~ q  -> t u~ b~ q
                  call interdk_qq_gen(p,i1,i8,3,4,5,i7,i6,i2,mt,msq_qqbar)      ! u  q~ -> t d  b~ q~
                  call interdk_qq_gen(p,i6,i8,3,4,5,i7,i1,i2,mt,msq_qbarqbar)   ! d~ q~ -> t u~ b~ q~

                  call interdk_qq_gen(p,i1,i7,3,4,5,i2,i6,i8,mt,msq_qb_qq)      ! u  b  -> t d  q~ q
                  call interdk_qq_gen(p,i6,i7,3,4,5,i2,i1,i8,mt,msq_qbarb_qq)   ! d~ b  -> t u~ q~ q

                  call interdk_qqid_gen(p,i1,i7,3,4,5,i2,i6,i8,mt,msq_qb_bb)      ! u  b  -> t d  b~ b
                  call interdk_qqid_gen(p,i6,i7,3,4,5,i2,i1,i8,mt,msq_qbarb_bb)   ! d~ b  -> t u~ b~ b
                  call interdk_qqid_gen(p,i1,i7,3,4,5,i8,i6,i2,mt,msq_qbbar_bb)   ! u  b~ -> t d  b~ b~
                  call interdk_qqid_gen(p,i6,i7,3,4,5,i8,i1,i2,mt,msq_qbarbbar_bb)! d~ b~ -> t u~ b~ b~
            endif
      else
            write(6,*) 'Real radiation not written for nwz = -1'
            stop
      endif

      if (corr_on_beam == 2) then

          if (iand(partons_enabled, quarkChannel) > 0) then
              msqall([2,4],5, 1,corr_on_beam) = aveqq*(half*msq_qb + (nf-1)*msq_qb_qq)
              msqall([2,4],5, WITH_B_BBAR,corr_on_beam) = aveqq*msq_qb_bb

              msqall([-1,-3],5, 1,corr_on_beam) = aveqq*(half*msq_qbarb + (nf-1)*msq_qbarb_qq)
              msqall([-1,-3],5, WITH_B_BBAR,corr_on_beam) = aveqq*msq_qbarb_bb

              do k=1,4
                  msqall([2,4],k, WITH_BBAR,corr_on_beam) = aveqq*msq_qq
                  msqall([2,4],-k, WITH_BBAR,corr_on_beam) = aveqq*msq_qqbar
                  msqall([-1,-3],k, WITH_BBAR,corr_on_beam) = aveqq*msq_qbarq
                  msqall([-1,-3],-k, WITH_BBAR,corr_on_beam) = aveqq*msq_qbarqbar
              enddo

              ! we reserve category 5 for the additional b~ b~ pieces
              msqall([2,4],-5, WITH_BBAR_BBAR,corr_on_beam) = aveqq*half*msq_qbbar_bb
              msqall([-1,-3],-5, WITH_BBAR_BBAR,corr_on_beam) = aveqq*half*msq_qbarbbar_bb
          endif

          if (iand(partons_enabled, gluonChannel) > 0) then
              msqall([2,4],0, WITH_BBAR,corr_on_beam) = aveqg*msq_qg
              msqall([-1,-3],0, WITH_BBAR,corr_on_beam) = aveqg*msq_qbarg
          endif

      else
          if (iand(partons_enabled, quarkChannel) > 0) then
              msqall(5,[2,4], 1,corr_on_beam) = aveqq*(half*msq_qb + (nf-1)*msq_qb_qq)
              msqall(5,[2,4], WITH_B_BBAR,corr_on_beam) = aveqq*msq_qb_bb

              msqall(5,[-1,-3], 1,corr_on_beam) = aveqq*(half*msq_qbarb + (nf-1)*msq_qbarb_qq)
              msqall(5,[-1,-3], WITH_B_BBAR,corr_on_beam) = aveqq*msq_qbarb_bb

              do k=1,4
                  msqall(k,[2,4], WITH_BBAR,corr_on_beam) = aveqq*msq_qq
                  msqall(-k,[2,4], WITH_BBAR,corr_on_beam) = aveqq*msq_qqbar
                  msqall(k,[-1,-3], WITH_BBAR,corr_on_beam) = aveqq*msq_qbarq
                  msqall(-k,[-1,-3], WITH_BBAR,corr_on_beam) = aveqq*msq_qbarqbar
              enddo

              ! we reserve category 5 for the additional b~ b~ pieces
              msqall(-5,[2,4], WITH_BBAR_BBAR,corr_on_beam) = aveqq*half*msq_qbbar_bb
              msqall(-5,[-1,-3], WITH_BBAR_BBAR,corr_on_beam) = aveqq*half*msq_qbarbbar_bb
          endif

          if (iand(partons_enabled, gluonChannel) > 0) then
              msqall(0,[2,4], WITH_BBAR,corr_on_beam) = aveqg*msq_qg
              msqall(0,[-1,-3], WITH_BBAR,corr_on_beam) = aveqg*msq_qbarg
          endif
      endif

c     if (corr_on_beam == 1) then ! i.e. corr_on_beam = 1
c       msqall(:,:, 1, corr_on_beam) = transpose(msqall(:,:, 1,corr_on_beam))
c       !msqall(:,:, 2, corr_on_beam) = transpose(msqall(:,:, 2,corr_on_beam))
c       msqall(:,:, 3, corr_on_beam) = transpose(msqall(:,:, 3,corr_on_beam))
c       msqall(:,:, 4, corr_on_beam) = transpose(msqall(:,:, 4,corr_on_beam))
c       msqall(:,:, 5, corr_on_beam) = transpose(msqall(:,:, 5,corr_on_beam))
c     endif

      enddo

      return
      end subroutine singletop_jet_heavy_real_all

      subroutine singletop_jet_heavy_gvec(p,n,in,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'nwz.f'
      include 'masses.f'
      include 'mxpart.f'
      integer:: in,i1,i2
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),n(4),fac,s16,dot

      msq(:,:)=0._dp

      if (corr_on_beam == 1) then
        i1=2
        i2=1
      elseif (corr_on_beam == 2) then
        i1=1
        i2=2
      endif

c--- fac adds additional color factor and corrects for (unnecessary)
c--- W width in propagator
      s16=two*dot(p,i1,6)
      fac=xn*((s16-wmass**2)**2+(wmass*wwidth)**2)/(s16-wmass**2)**2

      if (nwz == +1) then
            if     (in == 7) then
                  msq( 2,5) = aveqq*fac*wtgvecn(mt,twidth,7,i2,6,i1,3,4,5,p,n)
                  msq(-1,5) = aveqq*fac*wtgvecn(mt,twidth,7,i2,i1,6,3,4,5,p,n)
            elseif (in == i2) then
                  msq( 2,0) = aveqg*fac*wtgvecn(mt,twidth,i2,7,6,i1,3,4,5,p,n)
                  msq(-1,0) = aveqg*fac*wtgvecn(mt,twidth,i2,7,i1,6,3,4,5,p,n)
            endif
      else
            write(6,*) 'singletop_jet_heavy_gvec not written for nwz=-1'
            stop
      endif

c duplicate to 2nd generation
      msq(4,:)=msq(2,:)
      msq(-3,:)=msq(-1,:)

      if (corr_on_beam == 1) then
c--- interchange flavor labels
            msq(:,2)=msq(2,:); msq(2,:)=0._dp
            msq(:,4)=msq(4,:); msq(4,:)=0._dp
            msq(:,-1)=msq(-1,:); msq(-1,:)=0._dp
            msq(:,-3)=msq(-3,:); msq(-3,:)=0._dp
      endif

      return
      end subroutine singletop_jet_heavy_gvec

      subroutine singletop_jet_heavy_gs_all(p,ndmx,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qqgg.f'
      include 'ckm.f'
      include 'breit.f'
      include 'masses.f'
      integer:: j,k

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: ndmx
      real(dp), intent(out) :: msq(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp)::
     & sub28_7(4),sub28_7v,msq28_7(-nf:nf,-nf:nf),msq28_7v(-nf:nf,-nf:nf),
     & sub27_8(4),sub27_8v,msq27_8(-nf:nf,-nf:nf),msq27_8v(-nf:nf,-nf:nf),
     & sub78_2(4),sub78_2v,msq78_2(-nf:nf,-nf:nf),msq78_2v(-nf:nf,-nf:nf),
     & sub26_3(4),sub26_3v,msq26_3(-nf:nf,-nf:nf),msq26_3v(-nf:nf,-nf:nf),
     & sub36_2(4),msq36_2(-nf:nf,-nf:nf),
     & sub56_3(4),sub56_3v,msq56_3(-nf:nf,-nf:nf),msq56_3v(-nf:nf,-nf:nf),
     & sub36_5(4),msq36_5(-nf:nf,-nf:nf),
     & sub35_6(4),msq35_6(-nf:nf,-nf:nf),
     & sub25_3(4),msq25_3(-nf:nf,-nf:nf),
     & sub35_2(4),msq35_2(-nf:nf,-nf:nf),
     & sub18_7(4),sub18_7v,msq18_7(-nf:nf,-nf:nf),msq18_7v(-nf:nf,-nf:nf),
     & sub17_8(4),sub17_8v,msq17_8(-nf:nf,-nf:nf),msq17_8v(-nf:nf,-nf:nf),
     & sub78_1(4),sub78_1v,msq78_1(-nf:nf,-nf:nf),msq78_1v(-nf:nf,-nf:nf),
     & sub16_3(4),sub16_3v,msq16_3(-nf:nf,-nf:nf),msq16_3v(-nf:nf,-nf:nf),
     & sub36_1(4),msq36_1(-nf:nf,-nf:nf),
     & sub15_3(4),msq15_3(-nf:nf,-nf:nf),
     & sub35_1(4),msq35_1(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),dsubv

      logical :: beam1, beam2
      logical :: quarks, gluon

      external donothing_gvec

      corr_islight = .false.
      beam1 = any(beams_enabled(1:maxbeams) == 1)
      beam2 = any(beams_enabled(1:maxbeams) == 2)

      gluon = iand(partons_enabled, gluonChannel) > 0
      quarks = iand(partons_enabled, quarkChannel) > 0

c--- massless dipoles
      if (beam2) then
          corr_on_beam = 2
          corr_beam1 = .false.
          call dips(1,p,2,8,7,sub28_7,sub28_7v,msq28_7,msq28_7v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips(16,p,7,8,2,sub78_2,sub78_2v,msq78_2,msq78_2v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips(2,p,2,7,8,sub27_8,sub27_8v,msq27_8,msq27_8v,singletop_jet_heavy,singletop_jet_heavy_gvec)
      endif

      if (beam1) then
          corr_on_beam = 1
          corr_beam1 = .true.
          call dips(3,p,1,8,7,sub18_7,sub18_7v,msq18_7,msq18_7v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips(17,p,7,8,1,sub78_1,sub78_1v,msq78_1,msq78_1v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips(4,p,1,7,8,sub17_8,sub17_8v,msq17_8,msq17_8v,singletop_jet_heavy,singletop_jet_heavy_gvec)
      endif


c--- massive dipoles
      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.
      if (beam2) then
          corr_on_beam = 2
          corr_beam1 = .false.
          call dips_mass(5,p,2,6,3,sub26_3,sub26_3v,msq26_3,msq26_3v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips_mass(6,p,3,6,2,sub36_2,dsubv,msq36_2,dummyv,singletop_jet_heavy,donothing_gvec)
          call dips_mass(9,p,2,5,3,sub25_3,dsubv,msq25_3,dummyv,singletop_jet_heavy,donothing_gvec)
          call dips_mass(10,p,3,5,2,sub35_2,dsubv,msq35_2,dummyv,singletop_jet_heavy,donothing_gvec)
      endif

      if (beam1) then
          corr_on_beam = 1
          corr_beam1 = .true.
          call dips_mass(12,p,1,6,3,sub16_3,sub16_3v,msq16_3,msq16_3v,singletop_jet_heavy,singletop_jet_heavy_gvec)
          call dips_mass(13,p,3,6,1,sub36_1,dsubv,msq36_1,dummyv,singletop_jet_heavy,donothing_gvec)
          call dips_mass(14,p,1,5,3,sub15_3,dsubv,msq15_3,dummyv,singletop_jet_heavy,donothing_gvec)
          call dips_mass(15,p,3,5,1,sub35_1,dsubv,msq35_1,dummyv,singletop_jet_heavy,donothing_gvec)
      endif


      ! dips 7, 8 and 11 are called below for each corr_on_beam

      msq(:,:,:,:,:)=0._dp

c Corrections with heavy beam 2
      if (beam2) then
      corr_on_beam = 2
      corr_beam1=.false.

      !!! DIP BLOCK
      call dips_mass(7,p,3,6,5,sub36_5,dsubv,msq36_5,dummyv,singletop_jet_heavy,donothing_gvec)
      call dips_mass(8,p,3,5,6,sub35_6,dsubv,msq35_6,dummyv,singletop_jet_heavy,donothing_gvec)
      mass2=zip ! emission from a massless quark
      ggproc=.true.     ! only need this flag for final-final dipoles
      gqproc=.true.     ! only need this flag for final-final dipoles
      call dips_mass(11,p,5,6,3,sub56_3,sub56_3v,msq56_3,msq56_3v,singletop_jet_heavy,singletop_jet_heavy_gvec)
      ggproc = .false.
      gqproc = .false.
      mass2=mt  ! restore correct value
      !!! END DIP BLOCK

c These are singularities from 4-quark contributions to (q, q) and (q, q~)
      if (quarks) then
      do k=1,4
          j=2
          msq(1, j, k, 3,2) = two*Cf*(msq28_7(j,0)*sub28_7(gq)+msq28_7v(j,0)*sub28_7v)
          msq(1, j,-k, 3,2) = msq(1, j, k, 3,2)
          j=-1
          msq(1, j, k, 3,2) = two*Cf*(msq28_7(j,0)*sub28_7(gq)+msq28_7v(j,0)*sub28_7v)
          msq(1, j,-k, 3,2) = msq(1, j, k, 3,2)
      enddo

      j=2
      msq(1, j,-5, 3,2) = Cf*(msq28_7(j,0)*sub28_7(gq)+msq28_7v(j,0)*sub28_7v)
      msq(2, j,-5, 3,2) = Cf*(msq27_8(j,0)*sub27_8(gq)+msq27_8v(j,0)*sub27_8v)
      j=-1
      msq(1, j,-5, 3,2) = Cf*(msq28_7(j,0)*sub28_7(gq)+msq28_7v(j,0)*sub28_7v)
      msq(2, j,-5, 3,2) = Cf*(msq27_8(j,0)*sub27_8(gq)+msq27_8v(j,0)*sub27_8v)

c These are singularities from 4-quark contributions to (q,b)
      msq(1, 2, 5, 3,2) = msq(1, 2, 3, 3,2) ! right side: light quarks
      msq(1,-1, 5, 3,2) = msq(1,-1, 3, 3,2) ! right side: light quarks
      k=5
      do j=2,-1,-3
c--- standard massless dipoles
          msq(16,j,k, 1,2)=msq(16,j,k, 1,2)+tr*dfloat(nf)*(msq78_2(j,k)*sub78_2(gq)-msq78_2v(j,k)*sub78_2v)
c--- massive dipoles
          msq(11,j,k, 1,2)=msq(11,j,k, 1,2)+tr*dfloat(nf)*(msq56_3(j,k)*sub56_3(gq)-msq56_3v(j,k)*sub56_3v)
      enddo
      endif

c These are singularities from 2-quark, 2-gluon contributions to (q,g)
      if (gluon) then
      k=0
      do j=2,-1,-3
c--- standard massless dipoles
      msq(2,j,k, 1,2)=2._dp*tr*msq27_8(j,5)*sub27_8(qg)

      msq(1,j,k, 3,2)=xn*(msq28_7(j,k)*(sub28_7(gg))+msq28_7v(j,k)*sub28_7v)
      msq(16,j,k, 3,2)=xn*msq78_2(j,k)*sub78_2(qq)
c--- massive dipoles in relabelled scheme (345 -> 3, 6 -> 4, 7 -> 5, 8 -> 6)
      msq(5,j,k, 3,2)=xn*(msq26_3(j,k)*sub26_3(gg)+msq26_3v(j,k)*sub26_3v)          ! 28 singularity
      msq(6,j,k, 3,2)=xn*(msq36_2(j,k)*sub36_2(qq))                                 ! (345,8)
      msq(7,j,k, 3,2)= -(msq36_5(j,k)*sub36_5(qq))/xn                               ! (345,8)
      msq(11,j,k, 3,2)=-(msq56_3(j,k)*sub56_3(qq))/xn                               ! 78 singularity
      enddo
      endif

c These are singularities from 2-quark, 2-gluon contributions to (q,b)
      if (quarks) then
      k=5
      do j=2,-1,-3
c--- standard massless dipoles
      msq(1,j,k, 1,2)=msq(1,j,k, 1,2)+half*xn*msq28_7(j,k)*sub28_7(qq)
      msq(2,j,k, 1,2)=msq(2,j,k, 1,2)+half*xn*msq27_8(j,k)*sub27_8(qq)

      msq(16,j,k, 1,2)=msq(16,j,k, 1,2)+half*xn*(msq78_2(j,k)*sub78_2(gg)+msq78_2v(j,k)*sub78_2v)
c--- massive dipoles in relabelled scheme (345 -> 3, 6 -> 4, 7 -> 5, 8 -> 6)
      msq(7,j,k, 1,2)=half*xn*(msq36_5(j,k)*sub36_5(qq))                              ! (345,8)
      msq(8,j,k, 1,2)=half*xn*(msq35_6(j,k)*sub35_6(qq))                              ! (345,7)
      msq(11,j,k, 1,2)=msq(11,j,k, 1,2)
     &          +half*xn*(msq56_3(j,k)*sub56_3(gg)+msq56_3v(j,k)*sub56_3v)       ! 78 singularity
      msq(5,j,k, 1,2)=-half*(msq26_3(j,k)*sub26_3(qq))/xn                             ! 28 singularity
      msq(9,j,k, 1,2)=-half*(msq25_3(j,k)*sub25_3(qq))/xn                             ! 27 singularity
      msq(6,j,k, 1,2)= -half*(msq36_2(j,k)*sub36_2(qq))/xn                            ! (345,8)
      msq(10,j,k, 1,2)=-half*(msq35_2(j,k)*sub35_2(qq))/xn                            ! (345,7)
      enddo

      endif

c duplicate results to 2nd gen
      msq(:, 4, :, :,2) = msq(:, 2, :, :,2)
      msq(:,-3, :, :,2) = msq(:,-1, :, :,2)

      endif
c End of corrections with heavy beam 2


c Corrections with heavy beam 1
      if (beam1) then
      corr_on_beam = 1
      corr_beam1=.true.

      !!! DIP BLOCK
      call dips_mass(18,p,3,6,5,sub36_5,dsubv,msq36_5,dummyv,singletop_jet_heavy,donothing_gvec)
      call dips_mass(19,p,3,5,6,sub35_6,dsubv,msq35_6,dummyv,singletop_jet_heavy,donothing_gvec)

      mass2=zip ! emission from a massless quark
      ggproc=.true.     ! only need this flag for final-final dipoles
      gqproc=.true.     ! only need this flag for final-final dipoles
      call dips_mass(20,p,5,6,3,sub56_3,sub56_3v,msq56_3,msq56_3v,singletop_jet_heavy,singletop_jet_heavy_gvec)
      ggproc = .false.
      gqproc = .false.
      mass2=mt  ! restore correct value
      !!! END DIP BLOCK

c These are singularities from 4-quark contributions to (q, q) and (q, q~)
      if (quarks) then
      do j=1,4
      k=2
      msq(3, j, k, 3,1) = two*Cf*(msq18_7(0,k)*sub18_7(gq)+msq18_7v(0,k)*sub18_7v)
      msq(3,-j, k, 3,1) = msq(3, j, k, 3,1)
      k=-1
      msq(3, j, k, 3,1) = two*Cf*(msq18_7(0,k)*sub18_7(gq)+msq18_7v(0,k)*sub18_7v)
      msq(3,-j, k, 3,1) = msq(3, j, k, 3,1)
      enddo

      k=2
      msq(3,-5, k, 3,1) = Cf*(msq18_7(0,k)*sub18_7(gq)+msq18_7v(0,k)*sub18_7v)
      msq(4,-5, k, 3,1) = Cf*(msq17_8(0,k)*sub17_8(gq)+msq17_8v(0,k)*sub17_8v)
      k=-1
      msq(3,-5, k, 3,1) = Cf*(msq18_7(0,k)*sub18_7(gq)+msq18_7v(0,k)*sub18_7v)
      msq(4,-5, k, 3,1) = Cf*(msq17_8(0,k)*sub17_8(gq)+msq17_8v(0,k)*sub17_8v)

c These are singularities from 4-quark contributions to (q,b)
      msq(3, 5, 2, 3,1) = msq(3, 3, 2, 3,1)
      msq(3, 5,-1, 3,1) = msq(3, 3,-1, 3,1)
      j=5
      do k=2,-1,-3
c--- standard massless dipoles
          msq(17,j,k, 1,1)=msq(17,j,k, 1,1)+tr*dfloat(nf)*(msq78_1(j,k)*sub78_1(gq)-msq78_1v(j,k)*sub78_1v)
c--- massive dipoles
          msq(20,j,k, 1,1)=msq(20,j,k, 1,1)+tr*dfloat(nf)*(msq56_3(j,k)*sub56_3(gq)-msq56_3v(j,k)*sub56_3v)
      enddo
      endif

c These are singularities from 2-quark, 2-gluon contributions to (q,g)
      if (gluon) then
      j=0
      do k=2,-1,-3
c--- standard massless dipoles
      msq(4,j,k, 1,1)=2._dp*tr*msq17_8(5,k)*sub17_8(qg)
      msq(3,j,k, 3,1)=msq(3,j,k, 3,1)+xn*(msq18_7(j,k)*sub18_7(gg)+msq18_7v(j,k)*sub18_7v)
      msq(17,j,k, 3,1)=msq(17,j,k, 3,1)+xn*msq78_1(j,k)*sub78_1(qq)
c--- massive dipoles in relabelled scheme (345 -> 3, 6 -> 4, 7 -> 5, 8 -> 6)
      msq(12,j,k, 3,1)=xn*(msq16_3(j,k)*sub16_3(gg)+msq16_3v(j,k)*sub16_3v)          ! 18 singularity
      msq(13,j,k, 3,1)=xn*(msq36_1(j,k)*sub36_1(qq))                                 ! (345,8)
      msq(18,j,k, 3,1)= -(msq36_5(j,k)*sub36_5(qq))/xn                               ! (345,8)
      msq(20,j,k, 3,1)=-(msq56_3(j,k)*sub56_3(qq))/xn                               ! 78 singularity
      enddo
      endif

c These are singularities from 2-quark, 2-gluon contributions to (q,b)
      if (quarks) then
      j=5
      do k=2,-1,-3
c--- standard massless dipoles
      msq(3,j,k, 1,1)=msq(3,j,k, 1,1)+half*xn*msq18_7(j,k)*sub18_7(qq)
      msq(3,j,k, 1,1)=msq(3,j,k, 1,1)+half*xn*msq17_8(j,k)*sub17_8(qq)

      msq(17,j,k, 1,1)=msq(17,j,k, 1,1)+half*xn*(msq78_1(j,k)*sub78_1(gg)+msq78_1v(j,k)*sub78_1v)
c--- massive dipoles in relabelled scheme (345 -> 3, 6 -> 4, 7 -> 5, 8 -> 6)
      msq(18,j,k, 1,1)=half*xn*(msq36_5(j,k)*sub36_5(qq))                              ! (345,8)
      msq(19,j,k, 1,1)=half*xn*(msq35_6(j,k)*sub35_6(qq))                              ! (345,7)
      msq(20,j,k, 1,1)=msq(20,j,k, 1,1)
     &          +half*xn*(msq56_3(j,k)*sub56_3(gg)+msq56_3v(j,k)*sub56_3v)       ! 78 singularity
      msq(12,j,k, 1,1)=-half*(msq16_3(j,k)*sub16_3(qq))/xn                             ! 18 singularity
      msq(14,j,k, 1,1)=-half*(msq15_3(j,k)*sub15_3(qq))/xn                             ! 17 singularity
      msq(13,j,k, 1,1)= -half*(msq36_1(j,k)*sub36_1(qq))/xn                            ! (345,8)
      msq(15,j,k, 1,1)=-half*(msq35_1(j,k)*sub35_1(qq))/xn                            ! (345,7)
      enddo
      endif

c duplicate results to 2nd gen
      msq(:, :, 4, :,1) = msq(:, :, 2, :,1)
      msq(:, :,-3, :,1) = msq(:, :,-1, :,1)

      endif
c End of corrections with heavy beam 1

      return
      end subroutine singletop_jet_heavy_gs_all

      subroutine singletop_jet_heavy_z(p,z)
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      integer is
      real(dp), intent(in) :: p(mxpart,4),z
      real(dp):: dot,p1Dp3,p2Dp3,p35sq,xl13,xl23,xl17,xl27,
     & mbar13,mbar23,mbar35,
     & if_qq,if_gg,if_gq,if_qg,fi_qq,fi_gg,if_mqq,if_mgg,fi_mqq,
c     & fi_mgg,
     & ff_mqq0,ff_mgg,ff_1mqq
      real(dp) :: xl35_beam1, xl35_beam2

c--- calculate the two dot products involving the top momentum ("p3")
c--- f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6)+g(p7)
      p35sq=(p(3,4)+p(4,4)+p(5,4)+p(7,4))**2-(p(3,1)+p(4,1)+p(5,1)+p(7,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(7,2))**2-(p(3,3)+p(4,3)+p(5,3)+p(7,3))**2
      mbar35=mt/sqrt(p35sq)

      if (any(beams_enabled(1:maxbeams) == 1)) then
c CDTS (5.45,5.77)
          p1Dp3=dot(p,1,3)+dot(p,1,4)+dot(p,1,5)
          mbar13=mt/sqrt(-two*p1Dp3)
          xl13=log(-two*p1Dp3/renscale_beam1_isheavy_onheavy**2)
          xl17=log(-two*dot(p,1,7)/renscale_beam1_isheavy_onheavy**2)
          xl35_beam1=log(    p35sq/renscale_beam1_isheavy_onheavy**2)

          do is=1,3
              B1(b,b,q,is) = as_heavy_beam1/4._dp/pi*xn*(if_qq(z,xl17,is)+half*fi_gg(z,xl17,is)
     &                                         +ff_mqq0(z,xl35_beam1,mbar35,is)+half*ff_mgg(z,xl35_beam1,mbar35,is))
     &                      -as_heavy_beam1/4._dp/pi/xn*(if_mqq(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is))

              B1(g,g,q,is) = as_heavy_beam1/4._dp/pi*xn*(if_gg(z,xl17,is)+fi_qq(z,xl17,is)
     &                                         +if_mgg(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is))
     &                      -as_heavy_beam1/4._dp/pi/xn*(two*ff_1mqq(z,xl35_beam1,mbar35,is))

              B1(g,q,q,is) = as_heavy_beam1/4._dp/pi*two*Cf*if_gq(z,xl17,is)
              B1(q,g,q,is) = as_heavy_beam1/4._dp/pi*two*Tr*if_qg(z,xl17,is)
          enddo
      endif

      if (any(beams_enabled(1:maxbeams) == 2)) then
          p2Dp3=dot(p,2,3)+dot(p,2,4)+dot(p,2,5)
          mbar23=mt/sqrt(-two*p2Dp3)
          xl23=log(-two*p2Dp3/renscale_beam2_isheavy_onheavy**2)
          xl27=log(-two*dot(p,2,7)/renscale_beam2_isheavy_onheavy**2)
          xl35_beam2=log(    p35sq/renscale_beam2_isheavy_onheavy**2)

          do is=1,3
c 2-quark, 2-gluon contribution and final-state g->qqb splittings from 4-quark contribution
              B2(b,b,q,is) = as_heavy_beam2/4._dp/pi*xn*(if_qq(z,xl27,is)+half*fi_gg(z,xl27,is)
     &                                         +ff_mqq0(z,xl35_beam2,mbar35,is)+half*ff_mgg(z,xl35_beam2,mbar35,is))
     &                      -as_heavy_beam2/4._dp/pi/xn*(if_mqq(z,xl23,mbar23,is)+fi_mqq(z,xl23,mbar23,is))

              B2(g,g,q,is) = as_heavy_beam2/4._dp/pi*xn*(if_gg(z,xl27,is)+fi_qq(z,xl27,is)
     &                                         +if_mgg(z,xl23,mbar23,is)+fi_mqq(z,xl23,mbar23,is))
     &                      -as_heavy_beam2/4._dp/pi/xn*(two*ff_1mqq(z,xl35_beam2,mbar35,is))

c off diagonal terms from 4-quark contribution
              B2(g,q,q,is) = as_heavy_beam2/4._dp/pi*two*Cf*if_gq(z,xl27,is)
              B2(q,g,q,is) = as_heavy_beam2/4._dp/pi*two*Tr*if_qg(z,xl27,is)
          enddo
      endif

      end subroutine singletop_jet_heavy_z


      subroutine Wtoponshell_gen(q1,q2,q6,q7,p,iswitch,m)
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (Wtoponshell) by K. Ellis,
c but here extended to work for more generic momentum labels

c***********************************************************************
c     Author: R.K. Ellis, May 2012                                     *
c                                                                      *
c     b(q1)+g(q2)--> W^-(e-(q6)+nu~(q7))+t(p3,p4,p5)                   *
c                                                                      *
c     keeping polarization information for t,pc and gluon              *
c     iswitch= 0 for no gluon emission                                 *
c     iswitch=+1 for gluon emission in top decay                       *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & alt,dot,mt2,twoptDp2,senb
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p2,t,c,e,nb,eb,si,aa,bb,iswitch,q1,q2,q6,q7
      parameter(c=1,p2=2,e=3,nb=4,t=5,eb=6)
c-----matrix element for b(p1)+p2(g) -> t+W where t is on shell
c-----t rendered massless wrt eb, and c massless
c--- statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
c---zero all arrays
      m(:,:,:)=czip
      do si=1,4
      q(c,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(e,si)=p(q6,si)
      q(nb,si)=p(q7,si)
c      if (iswitch == 0) then
c      q(t,si)=p(3,si)+p(4,si)+p(5,si)
c      elseif (iswitch == 1) then
c      q(t,si)=p(3,si)+p(4,si)+p(5,si)+p(8,si)
c      endif
      q(t,si)=-p(q1,si)-p(q2,si)-p(q6,si)-p(q7,si)
      q(eb,si)=p(4,si)
      enddo
      mt2=mt**2
c---- now render "t" massless wrt to vector eb
      alt=mt2/(two*dot(q,t,eb))
      senb=two*dot(q,e,nb)
      if (senb < zip) then
      cprop=cplx2(senb-wmass**2,zip)
      else
      cprop=cplx2(senb-wmass**2,wmass*wwidth)
      endif
      twoptDp2=two*dot(q,t,p2)
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      enddo
      call spinoru(6,q,za,zb)

c----order of indices is polt,polg,polc
      m(1,1,2)=czip

      m(1,2,2)=czip

      m(1,1,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(e,t)*za(t,p2)*zb(c,
     &    t)*zb(c,nb)*izb(c,p2) + za(e,eb)*za(t,p2)*zb(c,nb)*zb(c,eb)*
     &    izb(c,p2)*alt + za(e,p2)*za(t,p2)*zb(c,nb) + za(e,p2)*zb(c,nb
     &    )*zb(c,eb)*izb(c,p2)*izb(t,eb)*mt2 )

      m(1,2,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(e,c)*zb(c,nb)*zb(eb
     &    ,p2)*iza(c,p2)*izb(t,eb)*mt2 - za(e,t)*za(c,t)*zb(c,nb)*zb(t,
     &    p2)*iza(c,p2) - za(e,eb)*za(c,t)*zb(c,nb)*zb(eb,p2)*iza(c,p2)
     &    *alt )
      m(1,2,1) = m(1,2,1) + cprop**(-1) * ( za(e,t)*zb(nb,p2)*iza(c,p2)
     &     )


      m(2,1,2)=czip

      m(2,2,2)=czip

      m(2,1,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(e,t)*za(eb,p2)*
     &    zb(c,t)*zb(c,nb)*iza(t,eb)*izb(c,p2)*mt - za(e,eb)*za(eb,p2)*
     &    zb(c,nb)*zb(c,eb)*iza(t,eb)*izb(c,p2)*alt*mt - za(e,p2)*za(eb
     &    ,p2)*zb(c,nb)*iza(t,eb)*mt - za(e,p2)*zb(c,t)*zb(c,nb)*izb(c,
     &    p2)*mt )

      m(2,2,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(e,c)*zb(c,nb)*
     &    zb(t,p2)*iza(c,p2)*mt + za(e,t)*za(c,eb)*zb(c,nb)*zb(t,p2)*
     &    iza(c,p2)*iza(t,eb)*mt + za(e,eb)*za(c,eb)*zb(c,nb)*zb(eb,p2)
     &    *iza(c,p2)*iza(t,eb)*alt*mt )
      m(2,2,1) = m(2,2,1) + cprop**(-1) * (  - za(e,eb)*zb(nb,p2)*iza(c
     &    ,p2)*iza(t,eb)*mt )

      return
      end subroutine Wtoponshell_gen
      end module


      subroutine Watoponshell_gen(q1,q2,q6,q7,p,iswitch,m)
      implicit none
      include 'types.f'
c Heavily based on routine for W+t production (Watoponshell) by K. Ellis,
c but here extended to work for more generic momentum labels

c***********************************************************************
c     Author: R.K. Ellis, January 2012                                 *
c                                                                      *
c     b~(q1)+g(q2)--> W^+(nu(q6)+e+(q7))+t~(p3,p4,p5)                 *
c                                                                      *
c     keeping polarization information for t~,pb and gluon             *
c     iswitch= 0 for no gluon emission                                 *
c     iswitch=-1 for gluon emission in atop decay                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & ala,dot,s34,mt2,twopaDp2
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p2,a,p3,p4,b,e,si,aa,bb,iswitch,q1,q2,q6,q7
      parameter(b=1,p2=2,p3=3,p4=4,a=5,e=6)
c-----matrix element for p1+p2(g) -> t+c where both t and c are on shell
c-----t rendered massless wrt eb, and c rendered massless wrt p2
c--- statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
c---zero all arrays
      m(:,:,:)=czip
      do si=1,4
      q(b,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(p3,si)=p(q6,si)
      q(p4,si)=p(q7,si)
c      if (iswitch == 0) then
c      q(a,si)=p(3,si)+p(4,si)+p(5,si)
c      elseif (iswitch == -1) then
c      q(a,si)=p(3,si)+p(4,si)+p(5,si)+p(8,si)
c      endif
      q(a,si)=-p(q1,si)-p(q2,si)-p(q6,si)-p(q7,si)
      q(e,si)=p(3,si)
      enddo
      mt2=mt**2
c---- now render "a" massless wrt to vector e
      ala=mt2/(two*dot(q,a,e))
      s34=two*dot(q,p3,p4)
      if (s34 < zip) then
      cprop=cplx2(s34-wmass**2,zip)
      else
      cprop=cplx2(s34-wmass**2,wmass*wwidth)
      endif
      twopaDp2=two*dot(q,a,p2)
      do si=1,4
      q(a,si)=q(a,si)-ala*q(e,si)
      enddo
      call spinoru(6,q,za,zb)

c----order of indices is polb,polg,polatop
      m(1,1,2)= + cprop**(-1)*twopaDp2**(-1) * ( za(e,p2)*za(b,p3)*zb(e
     &    ,p4)*zb(b,a)*izb(b,p2)*ala + za(b,p3)*za(a,p2)*zb(b,a)*zb(a,
     &    p4)*izb(b,p2) )
      m(1,1,2) = m(1,1,2) + cprop**(-1) * ( za(p2,p3)*zb(a,p4)*izb(b,p2
     &    ) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*mt2*twopaDp2**(-1) * ( za(e,p2)
     &    *za(b,p3)*zb(b,p4)*iza(e,a)*izb(b,p2) )

      m(1,2,2)= + cprop**(-1)*twopaDp2**(-1) * ( za(e,b)*za(b,p3)*zb(e,
     &    p4)*zb(a,p2)*iza(b,p2)*ala - za(b,a)*za(b,p3)*zb(a,p2)*zb(a,
     &    p4)*iza(b,p2) - za(b,p3)*zb(a,p2)*zb(p2,p4) )
      m(1,2,2) = m(1,2,2) + cprop**(-1)*mt2*twopaDp2**(-1) * (  - za(e,
     &    b)*za(b,p3)*zb(p2,p4)*iza(e,a)*iza(b,p2) )

      m(1,1,1)= + cprop**(-1)*mt*twopaDp2**(-1) * (  - za(e,p2)*za(b,p3
     &    )*zb(e,b)*zb(e,p4)*izb(e,a)*izb(b,p2)*ala - za(b,p3)*za(a,p2)
     &    *zb(e,b)*zb(a,p4)*izb(e,a)*izb(b,p2) + za(b,p3)*za(a,p2)*zb(b
     &    ,p4)*izb(b,p2) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt * ( za(p2,p3)*zb(e,p4)*izb(e
     &    ,a)*izb(b,p2) )

      m(1,2,1)= + cprop**(-1)*mt*twopaDp2**(-1) * ( za(e,b)*za(b,p3)*
     &    zb(e,p2)*zb(e,p4)*iza(b,p2)*izb(e,a)*ala - za(b,a)*za(b,p3)*
     &    zb(e,p2)*zb(a,p4)*iza(b,p2)*izb(e,a) + za(b,a)*za(b,p3)*zb(p2
     &    ,p4)*iza(b,p2) - za(b,p3)*zb(e,p2)*zb(p2,p4)*izb(e,a) )


      m(2,1,2)=czip

      m(2,2,2)=czip

      m(2,1,1)=czip

      m(2,2,1)=czip

      return
      end subroutine Watoponshell_gen
