!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_int
      use types
      use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
      use MCFMStorage
      use Scalevar
      use PDFerrors
      use LHAPDF
      use SCET
      use singletop2_nnlo
      use singletop2_nnlo_vars
      use singletop2_scale_m
!!!#ifdef HAVE_RECOLA
!!!      use singletop2_recola
!!!      use singletop_heavy_recola
!!!#endif
      use m_gencuts, only : reweight_user, enable_reweight_user

      implicit none

      public :: lowint_assemble

      public :: virtint_assemble
      public :: virtint_assemble_kbq_tpq

      public :: realint_assemble

      public :: lowint_fillmsq
      public :: virtint_fillmsq
      public :: realint_fillmsq

      public :: calc_singletop_pdfs_virt
      public :: calc_singletop_pdfs_lord
      public :: calc_singletop_pdfs_real

      include 'nf.f'

      ! last entry: beam1, beam2
      real(dp), save, private :: singletop_pdfs(-nf:nf, max_corr_on_beam, 2)
!$omp threadprivate(singletop_pdfs)

#define IPSDECAY 3
#define IPSHEAVY 2
#define IPSLIGHT 1

#define ONLY_MAIN_B 1
#define WITH_B 2
#define WITH_BBAR 3
#define WITH_B_BBAR 4
#define WITH_BBAR_BBAR 5

      contains

      function lowint_assemble(msqall,bcontribs,pswt,xx)
          implicit none
          include 'nf.f'
          include 'constants.f'
          include 'energy.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: lowint_assemble(max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: pswt, xx(2)
          logical, intent(in) :: bcontribs(max_bcontrib)
          real(dp), intent(in) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          integer :: j,k,m
          real(dp) :: flux

          lowint_assemble = 0._dp

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)

          loop_corr_on_beam: do m=1,maxbeams
              corr_on_beam = beams_enabled(m)
              loop_b_contrib: do b_contrib=1,max_bcontrib
                  if (bcontribs(b_contrib) .eqv. .false.) cycle

                  do j=-nf,nf; do k=-nf,nf
                      if ((selectpdfs(1,k) .eqv. .false.) &
                          .or. (selectpdfs(2,j) .eqv. .false.)) cycle

                      lowint_assemble(b_contrib, corr_on_beam) = &
                          lowint_assemble(b_contrib, corr_on_beam) + &
                          singletop_pdfs(k,corr_on_beam,1) * &
                          singletop_pdfs(j,corr_on_beam,2) * &
                          msqall(k,j,b_contrib,corr_on_beam)
                  enddo; enddo

                  lowint_assemble(b_contrib, corr_on_beam) = &
                      lowint_assemble(b_contrib, corr_on_beam)*flux*pswt/BrnRat
              enddo loop_b_contrib
          enddo loop_corr_on_beam
          
      end function

      function scetint_fillxmsq(p,xx,z1,z2,QB,pswt,central)
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'kpart.f'
          include 'energy.f'
          include 'constants.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: scetint_fillxmsq

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: xx(2),z1,z2,QB(2)
          real(dp), intent(in) :: pswt
          logical, intent(in) :: central

          real(dp) :: flux
          integer :: iorder


          scetint_fillxmsq = 0._dp

          if (kpart == ksnlo .or. kpart == ksnloV) then
              iorder = 1
          else
              iorder = 2
          endif

          if (currentContrib == 1) then
!             call lumxmsq_singletop(p,xx,z1,z2,QB,iorder, &
!                 scetint_fillxmsq,central)
              call lumxmsq_singletop_new(p,xx,z1,z2,QB,iorder, &
                  scetint_fillxmsq,central)
          elseif (currentContrib == 2) then
              call lumxmsq_singletop_prod(p,xx,z1,z2,QB,iorder, &
                  scetint_fillxmsq,central)

          elseif (currentContrib == 3) then
              call lumxmsq_singletop_decay_jetmass(p,xx,z1,z2,iorder, &
                  scetint_fillxmsq,central)
          endif

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)

          scetint_fillxmsq = scetint_fillxmsq * &
              (4*sqrt(z1*z2)) * flux * pswt / BrnRat

      end function 

      subroutine virtint_fillmsq(p, z, bcontribs, msqvall)
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'kprocess.f'
          include 'ipsgen.f'

          real(dp), intent(in) :: p(mxpart,4), z
          logical, intent(in) :: bcontribs(max_bcontrib)
          real(dp), intent(out) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          if (kcase == kbq_tpq) then
              if (currentContrib == 1) then
                  call singletop2_scet_virt_light_all(p,msqvall)
                  call singletop2_scet_z(p,z)
              elseif (currentContrib == 2) then
                  call singletop2_scet_virt_heavy_prod_all(p,msqvall)
                  call singletop2_scet_heavy_prod_z(p,z)
              elseif (currentContrib == 3) then
                  call singletop2_scet_virt_heavy_decay_all(p,msqvall)
              elseif (currentContrib == 4) then
                  call singletop_jet_light_heavy_vr_all(p,msqvall)
                  call singletop_jet_light_heavy_vr_z(p,z)
              elseif (currentContrib == 5) then
                  call singletop_light_decay_vr(p,msqvall)
                  call singletop_light_decay_vr_z(p,z)
              elseif (currentContrib == 6) then
                  call singletop_heavy_decay_vr(p,msqvall)
                  call singletop_heavy_decay_vr_z(p,z)
              endif
          elseif (kcase == kbq_tpq_jet) then
              if (currentContrib == 1) then
                  call singletop_jet_light_virt_all(p,msqvall)
                  call singletop_jet_light_z(p,z)
              elseif (currentContrib == 2) then
                  call singletop_jet_heavy_virt_all(p,msqvall)
                  call singletop_jet_heavy_z(p,z)
              elseif (currentContrib == 3) then
                  call singletop_jet_decay_virt_all(p,msqvall)
              endif
          endif

          return

      end subroutine

      subroutine realint_fillmsq(p, ndmax, bcontribs, msqall, msqcall)
          use Multichannel
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'kprocess.f'
          include 'maxd.f'
          include 'incldip.f'
          include 'ipsgen.f'

          real(dp), intent(in) :: p(mxpart,4)
          integer, intent(in) :: ndmax
          logical, intent(in) :: bcontribs(:,:)
          real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(inout) :: msqcall(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          incldip(:) = .false.

          if (kcase == kbq_tpq) then
              if (currentContrib == 1) then
                  msqall = 0._dp
                  call singletop_jet_light_msqall(p,msqall)
                  call singletop2_scet_gsall(p,ndmax,msqcall)

                  !call singcheck(singletop_jet_light_msqall, &
                      !singletop2_scet_gsall, p)
              elseif (currentContrib == 2) then
                  call qqb_tbb_g_heavy_all(p,msqall)
                  call singletop2_scet_heavy_prod_gs_all(p,ndmax,msqcall)
              elseif (currentContrib == 3) then
                  call singletop2_heavy_decay_g_all(p,msqall)
                  call singletop2_heavy_decay_gs_all_new(p,ndmax,msqcall)
              elseif (currentContrib == 4) then
                  call singletop_jet_light_heavy_rr_all(p,msqall)
                  call singletop_jet_light_heavy_rr_gs_all(p,ndmax,msqcall)
                  
                  !call singcheck(singletop_jet_light_heavy_rr_all, &
                      !singletop_jet_light_heavy_rr_gs_all, p)
              elseif (currentContrib == 5) then
                  call singletop_light_decay_rr(p,msqall)
                  call singletop_light_decay_rr_gs(p,ndmax,msqcall)

                  ! LABORDAY ADDED
                  !call singcheck(singletop_light_decay_rr, &
                      !singletop_light_decay_rr_gs, p)
              elseif (currentContrib == 6) then
                  call singletop_heavy_decay_rr(p,msqall)
                  call singletop_heavy_decay_rr_gs(p,ndmax,msqcall) 

                  !call singcheck(singletop_heavy_decay_rr, &
                      !singletop_heavy_decay_rr_gs, p)
              endif
          elseif (kcase == kbq_tpq_jet) then
              if (currentContrib == 1) then
                  call singletop_jet_light_real_all(p,msqall)
                  call singletop_jet_light_gs_all(p,ndmax,msqcall)
              elseif (currentContrib == 2) then
                  call singletop_jet_heavy_real_all(p,msqall)
                  call singletop_jet_heavy_gs_all(p,ndmax,msqcall)
              elseif (currentContrib == 3) then
                  call singletop_jet_decay_real_all(p,msqall)
                  call singletop_jet_decay_gs(p,ndmax,msqcall)

                  !write (*,*) "SINGCHECK CALL"
                  !call singcheck(singletop_jet_decay_real_all,singletop_jet_decay_gs,p)
              endif
          endif

!       if (currentContrib == 1) then
!           if (includereal) call singletop_jet_light(p,msq)
!           call singletop2_scet_gs(p,msqc)
!           !call singcheck(singletop_jet_light,singletop2_scet_gs,p)
!       elseif (currentContrib == 2) then
!           !if (includereal) call singletop_jet_heavy(p,msq)
!           !if (includereal) call qqb_tbb_g_heavy(p,msq)
!           call singletop2_scet_heavy_prod_gs(p,msqc)
!       elseif (currentContrib == 3) then
!           if (includereal) call singletop2_heavy_decay_g(p,msq)
!           call singletop2_heavy_decay_gs(p,msqc)
!           !if (includereal) call bq_tpq_gdk(p,msq)
!           !call bq_tpq_gsdk(p,msqc)
!       endif
      end subroutine

      subroutine lowint_fillmsq(p, bcontribs, msqall)
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'kprocess.f'
          include 'ipsgen.f'

          real(dp), intent(in) :: p(mxpart,4)
          logical, intent(in) :: bcontribs(max_bcontrib)
          real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          if (kcase == kbq_tpq) then
              msqall = 0._dp
              if (currentContrib == 1) then
                  call singletop2_scet_tree_ub(p,msqall(:,:,1,1))
                  call singletop2_scet_tree_bu(p,msqall(:,:,1,2))
              elseif (any(currentContrib == [2,3])) then
                  call singletop2_scet_tree_ub(p,msqall(:,:,1,2))
                  call singletop2_scet_tree_bu(p,msqall(:,:,1,1))
              elseif (currentContrib == 4) then
                  if (ipsgen == 4) then
                      ! this is only for call within singletop_virtint
                      ! when virt corrections are on the light line
                      
                      ! return real emission on heavy line
                      call qqb_tbb_g_heavy_all_swap(p,msqall)
                  endif
              elseif (currentContrib == 5) then
                  !if (ipsgen == 6) then ! VR
                      ! this is only for call within singletop_virtint
                      ! when virt corrections are on the light or heavy line

                      ! real emission in decay
                      call singletop_decay_real_swap(p,msqall)
                  !endif
              elseif (currentContrib == 6) then
                  call singletop_decay_real_hxd(p,msqall)
                  !call singletop2_scet_tree_ub(p,msqall(:,:,1,2))
                  !call singletop2_scet_tree_bu(p,msqall(:,:,1,1))
              endif
          elseif (kcase == kbq_tpq_jet) then
              if (currentContrib == 1) then
                  call singletop_jet_light_msqall(p,msqall)
              elseif (currentContrib == 2) then
                  ! any of these works, but I suspect the second is faster
                  ! (and simpler)
                  !call singletop_jet_heavy_all(p,msqall)
                  call qqb_tbb_g_heavy_all(p,msqall)
              elseif (currentContrib == 3) then
                  !call singletop2_heavy_decay_g_all(p,msqall)
                  call singletop_jet_decay_all(p,msqall)
              endif
          endif

          return

      end subroutine

      subroutine calc_singletop_pdfs_real(xx, ndmax)
          implicit none
          include 'maxd.f'
          include 'incldip.f'
          include 'beamtype.f'

          real(dp), intent(in) :: xx(2)
          integer, intent(in) :: ndmax

          integer :: nd

          !singletop2_pdfs(:,:) = 0._dp
          singletop_pdfs = 0._dp
          singletop2_dipole_pdfs(:,:,:) = 0._dp

          if (any(currentContrib == [1,4,5])) then
              if (usemask) then
                  if (any(beams_enabled(1:maxbeams) == 1)) then
                      ! corr_on_beam = 1
                      call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop_pdfs(:,1,1),1,maskb1)
                      ! only b for beam 2
                       singletop_pdfs(5,1,2) = fdist_one(ih2,xx(2),facscale_beam2_isheavy_onlight,5, 2)
                  endif

                  if (any(beams_enabled(1:maxbeams) == 2)) then
                      ! corr_on_beam = 2
                       singletop_pdfs(5,2,1) = fdist_one(ih1,xx(1),facscale_beam1_isheavy_onlight,5, 1)
                      call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop_pdfs(:,2,2),2,maskb2)
                  endif

              else
                  ! corr_on_beam == 1
                  call fdist(ih1, xx(1), facscale_beam1_islight_onlight, singletop_pdfs(:,1,1),1)
                  call fdist(ih2, xx(2), facscale_beam2_isheavy_onlight, singletop_pdfs(:,1,2),2)

                  ! corr_on_beam == 2
                  call fdist(ih1, xx(1), facscale_beam1_isheavy_onlight, singletop_pdfs(:,2,1),1)
                  call fdist(ih2, xx(2), facscale_beam2_islight_onlight, singletop_pdfs(:,2,2),2)
              endif
          elseif (any(currentContrib == [2,3,6])) then
              ! corr_on_beam == 2
              call fdist(ih1, xx(1), facscale_beam1_islight_onheavy, singletop_pdfs(:,2,1),1)
              call fdist(ih2, xx(2), facscale_beam2_isheavy_onheavy, singletop_pdfs(:,2,2),2)

              ! corr_on_beam == 1
              call fdist(ih1, xx(1), facscale_beam1_isheavy_onheavy, singletop_pdfs(:,1,1),1)
              call fdist(ih2, xx(2), facscale_beam2_islight_onheavy, singletop_pdfs(:,1,2),2)
          endif

          if (usemask) then
              do nd=1,ndmax
                  if (incldip(nd) .eqv. .false.) cycle
                  call fdist(ih1,xx(1),singletop2_dipscale(nd,1),singletop2_dipole_pdfs(nd,1,:),1,maskb1)
                  call fdist(ih2,xx(2),singletop2_dipscale(nd,2),singletop2_dipole_pdfs(nd,2,:),2,maskb2)
              enddo
          else
              do nd=1,ndmax
                  if (incldip(nd) .eqv. .false.) cycle
                  call fdist(ih1,xx(1),singletop2_dipscale(nd,1),singletop2_dipole_pdfs(nd,1,:),1)
                  call fdist(ih2,xx(2),singletop2_dipscale(nd,2),singletop2_dipole_pdfs(nd,2,:),2)
              enddo
          endif

      end subroutine

      subroutine calc_singletop_pdfs_lord(xx)
          implicit none
          include 'beamtype.f'

          real(dp), intent(in) :: xx(2)

          singletop_pdfs(:,:,:) = 0._dp

          if (currentContrib == 1) then
              if (usemask) then
                  if (any(beams_enabled(1:maxbeams) == 1)) then
                      ! corr_on_beam = 1
                      call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop_pdfs(:,1,1),1,maskb1)
                      ! only b for beam 2
                      singletop_pdfs(5,1,2) = fdist_one(ih2,xx(2),facscale_beam2_isheavy_onlight,5, 2)
                  endif

                  if (any(beams_enabled(1:maxbeams) == 2)) then
                      ! corr_on_beam = 2
                      singletop_pdfs(5,2,1) = fdist_one(ih1,xx(1),facscale_beam1_isheavy_onlight,5, 1)
                      call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop_pdfs(:,2,2),2,maskb2)
                  endif
              else
                  ! corr_on_beam = 1
                  call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop_pdfs(:,1,1),1)
                  call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight, singletop_pdfs(:,1,2),2)

                  ! corr_on_beam = 2
                  call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight, singletop_pdfs(:,2,1),1)
                  call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop_pdfs(:,2,2),2)
              endif
          elseif (any(currentContrib == [2,3])) then
              call fdist(ih1, xx(1), facscale_beam1_islight_onheavy, singletop_pdfs(:,2,1),1)
              call fdist(ih2, xx(2), facscale_beam2_isheavy_onheavy, singletop_pdfs(:,2,2),2)

              ! corr_on_beam == 1
              call fdist(ih1, xx(1), facscale_beam1_isheavy_onheavy, singletop_pdfs(:,1,1),1)
              call fdist(ih2, xx(2), facscale_beam2_islight_onheavy, singletop_pdfs(:,1,2),2)
          endif
      end subroutine

      subroutine calc_singletop_pdfs_virt(xx, z)
          implicit none
          include 'beamtype.f'

          real(dp), intent(in) :: xx(2)
          real(dp), intent(in), optional :: z

          singletop2_pdfs(:,:) = 0._dp

          if (any(currentContrib == [1,4,5])) then
              if (usemask) then
                  if (any(beams_enabled(1:maxbeams) == 1)) then
                      ! corr_on_beam = 1
                      call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop_pdfs(:,1,1),1,maskb1)
                      ! only b for beam 2
                      singletop_pdfs(5,1,2) = fdist_one(ih2,xx(2),facscale_beam2_isheavy_onlight,5, 2)
                  endif

                  if (any(beams_enabled(1:maxbeams) == 2)) then
                      ! corr_on_beam = 2
                      singletop_pdfs(5,2,1) = fdist_one(ih1,xx(1),facscale_beam1_isheavy_onlight,5, 1)
                      call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop_pdfs(:,2,2),2,maskb2)
                  endif
              else
                  call fdist(ih1,xx(1),facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_light,:),1)
                  call fdist(ih2,xx(2),facscale_beam2_isheavy_onlight, singletop2_pdfs(i_beam2_heavy,:),2)

                  call fdist(ih1,xx(1),facscale_beam1_isheavy_onlight, singletop2_pdfs(i_beam1_heavy,:),1)
                  call fdist(ih2,xx(2),facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_light,:),2)
              endif
          elseif (any(currentContrib == [2,3,6])) then
              call fdist(ih1,xx(1),facscale_beam1_islight_onheavy, singletop2_pdfs(i_beam1_light,:),1)
              call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_heavy,:),2)

              call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_heavy,:),1)
              call fdist(ih2,xx(2),facscale_beam2_islight_onheavy, singletop2_pdfs(i_beam2_light,:),2)
          endif

          if (present(z)) then
              if (any(currentContrib == [1,4,5])) then
                  if (usemask) then
                      if (any(beams_enabled(1:maxbeams) == 1)) then
                          if (z > xx(1)) then
                              call fdist(ih1,xx(1)/z,facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_light_z,:),1,maskb1)
                          endif
!                         if (z > xx(2)) then
!                             call fdist(ih2,xx(2)/z,facscale_beam2_isheavy_onlight, singletop2_pdfs(i_beam2_heavy_z,:),2,maskb2)
!                         endif
                      endif

                      if (any(beams_enabled(1:maxbeams) == 2)) then
!                         if (z > xx(1)) then
!                             call fdist(ih1,xx(1)/z,facscale_beam1_isheavy_onlight, singletop2_pdfs(i_beam1_heavy_z,:),1,maskb1)
!                         endif
                          if (z > xx(2)) then
                              call fdist(ih2,xx(2)/z,facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_light_z,:),2,maskb2)
                          endif
                      endif
                  else
                      if (z > xx(1)) then
                          call fdist(ih1,xx(1)/z,facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_light_z,:),1)
                          call fdist(ih1,xx(1)/z,facscale_beam1_isheavy_onlight, singletop2_pdfs(i_beam1_heavy_z,:),1)
                      endif

                      if (z > xx(2)) then
                          call fdist(ih2,xx(2)/z,facscale_beam2_isheavy_onlight, singletop2_pdfs(i_beam2_heavy_z,:),2)
                          call fdist(ih2,xx(2)/z,facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_light_z,:),2)
                      endif
                  endif
              elseif (any(currentContrib == [2,3,6])) then
                  if (z > xx(1)) then
                      call fdist(ih1,xx(1)/z,facscale_beam1_islight_onheavy, singletop2_pdfs(i_beam1_light_z,:),1)
                      call fdist(ih1,xx(1)/z,facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_heavy_z,:),1)
                  endif

                  if (z > xx(2)) then
                      call fdist(ih2,xx(2)/z,facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_heavy_z,:),2)
                      call fdist(ih2,xx(2)/z,facscale_beam2_islight_onheavy, singletop2_pdfs(i_beam2_light_z,:),2)
                  endif

              endif
          endif

      end subroutine

      function realint_assemble(ndmax, msqall,msqcall,bcontribs,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kprocess.f'
          include 'constants.f'
          include 'energy.f'
          include 'maxd.f'
          include 'incldip.f'

          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          integer, intent(in) :: ndmax
          real(dp), intent(in) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqcall(ndmax,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          logical, intent(in) :: bcontribs(0:ndmax,max_bcontrib)
          real(dp), intent(in) :: pswt, xx(2)

          real(dp) :: realint_assemble(0:ndmax,max_bcontrib,max_corr_on_beam)

          real(dp) :: flux
          integer :: j,k,m,nd


          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)

          realint_assemble = 0._dp

          do nd=0,ndmax
              do m=1,maxbeams
                  corr_on_beam = beams_enabled(m)
                  do b_contrib=1,max_bcontrib
                      ! bcontribs is filled by includedipole
                      ! incldip is filled by the dipole routines alpha cuts


                      ! We assume that each incldip is associated with a unique
                      ! bcontrib. I think this holds. 
                      if ((bcontribs(nd,b_contrib) .eqv. .false.) .or. &
                          ((nd > 0) .and. (incldip(nd) .eqv. .false.)) ) cycle

                      if (nd == 0) then
                          do k=-nf,nf; do j=-nf,nf
                              realint_assemble(nd,b_contrib,corr_on_beam) = &
                                  realint_assemble(nd,b_contrib,corr_on_beam) + &
                                    msqall(j,k,b_contrib,corr_on_beam) * &
                                    singletop_pdfs(j,corr_on_beam,1) * &
                                    singletop_pdfs(k,corr_on_beam,2)
                          enddo; enddo
                      else
                          do k=-nf,nf; do j=-nf,nf
                              realint_assemble(nd,b_contrib,corr_on_beam) = &
                                  realint_assemble(nd,b_contrib,corr_on_beam) - &
                                  msqcall(nd,j,k,b_contrib,corr_on_beam) * &
                                  singletop2_dipole_pdfs(nd,1,j) * &
                                  singletop2_dipole_pdfs(nd,2,k)
                          enddo; enddo
                      endif

!                     do j=-nf,nf; do k=-nf,nf
!                         if ((selectpdfs(1,j) .eqv. .false.) &
!                             .or. (selectpdfs(2,k) .eqv. .false.)) cycle

!                         if ((nd == 0) .and. (msqall(j,k,b_contrib,corr_on_beam) == 0._dp)) cycle
!                         if ((nd > 0) .and. (msqcall(nd,j,k,b_contrib,corr_on_beam) == 0._dp)) cycle

!                         if (nd == 0) then
!                             realint_assemble(nd,b_contrib,corr_on_beam) = &
!                                 realint_assemble(nd,b_contrib,corr_on_beam) + &
!                                   msqall(j,k,b_contrib,corr_on_beam) * &
!                                   gingletop_pdfs(j,corr_on_beam,1) * &
!                                   singletop_pdfs(k,corr_on_beam,2)
!                         else
!                             realint_assemble(nd,b_contrib,corr_on_beam) = &
!                                 realint_assemble(nd,b_contrib,corr_on_beam) - &
!                                 msqcall(nd,j,k,b_contrib,corr_on_beam) * &
!                                 singletop2_dipole_pdfs(nd,1,j) * &
!                                 singletop2_dipole_pdfs(nd,2,k)
!                         endif
!                     enddo; enddo

                      realint_assemble(nd,b_contrib, corr_on_beam) = &
                          realint_assemble(nd,b_contrib, corr_on_beam) * &
                            flux*pswt/BrnRat

                  enddo
              enddo
          enddo


      end function

      function virtint_assemble(msqall,msqvall,z,bcontribs,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kprocess.f'

          real(dp) :: virtint_assemble(max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: z, pswt, xx(2)
          logical, intent(in) :: bcontribs(max_bcontrib)

          if (kcase == kbq_tpq) then
              virtint_assemble = virtint_assemble_kbq_tpq(msqall,msqvall,z,bcontribs,pswt,xx)
          elseif (kcase == kbq_tpq_jet) then
              virtint_assemble = virtint_assemble_kbq_tpq_jet(msqall,msqvall,z,bcontribs,pswt,xx)
          endif

      end function

      function virtint_assemble_kbq_tpq(msqall,msqvall,z,bcontribs,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'PR_new.f'
          include 'PR_stop.f'
          include 'agq.f'
          include 'epinv.f'
          include 'ipsgen.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: virtint_assemble_kbq_tpq(max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: z, pswt, xx(2)
          logical, intent(in) :: bcontribs(max_bcontrib)

          real(dp) :: xmsq(max_bcontrib,max_corr_on_beam)
          real(dp) :: AP(-1:1, -1:1, 3)
          integer :: j,k
          real(dp) :: flux

          logical :: bin
          common/bin/bin

          virtint_assemble_kbq_tpq = 0._dp

          ! "deliberately" not compatible with t~ (no abs(5))
          ! and explicitly constrained to the available partonic channels
          ! for maximum testability and readability

          xmsq = 0._dp

          if (currentContrib == 1 .or. currentContrib == 5) then
              do j=-nf,nf; do k=-nf,nf
                  if ( any(j == [2,4,-1,-3]) .and. k == 5 ) then
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqvall(j,k, 1,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      endif
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)
                  elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqvall(j,k, 1,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      endif
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
                  elseif (j==5 .and. k==0) then
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + sum(msqall(j,[2,4,-1,-3], 1, 2)) * (AP(q,g,2) + Q2(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
                  elseif (j==0 .and. k==5) then
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + sum(msqall([2,4,-1,-3],k, 1,1)) * (AP(q,g,2) + Q1(q,g,q,2)) * &
                          singletop2_pdfs(i_beam1_light_z,j)/z * singletop2_pdfs(i_beam2_heavy,k)
                  endif
              enddo; enddo
          elseif (currentContrib == 2 .or. currentContrib == 6) then
              do j=-nf,nf; do k=-nf,nf
                  if (any(j == [2,4,-1,-3]) .and. k == 5) then
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqvall(j,k, 1,2) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      endif

                      call singletop2_fillAP(z, i_beam2_heavy, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)) * &
                          singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy_z,k)/z

                  elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqvall(j,k, 1,1) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      endif

                      call singletop2_fillAP(z, i_beam1_heavy, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z * singletop2_pdfs(i_beam2_light,k)
                  elseif (any(j == [2,4,-1,-3]) .and. k == 0) then
                      call singletop2_fillAP(z, i_beam2_heavy, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,+5, 1,2)*(AP(q,g,2) + B2(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy_z,k)/z
                  elseif (j == 0 .and. any(k == [2,4,-1,-3])) then
                      call singletop2_fillAP(z, i_beam1_heavy, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(+5,k, 1,1)*(AP(q,g,2) + B1(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z * singletop2_pdfs(i_beam2_light,k)
                  endif
              enddo; enddo
          elseif (currentContrib == 3) then
              do j=-nf,nf; do k=-nf,nf
                  if (any(j == [2,4,-1,-3]) .and. k == 5) then
                      xmsq(ONLY_MAIN_B, 2) = xmsq(ONLY_MAIN_B, 2) + msqvall(j,k,1,2) * &
                          singletop2_pdfs(i_beam1_light,j) * &
                          singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B, 2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2) * &
                              singletop2_pdfs(i_beam1_light,j) * &
                              singletop2_pdfs(i_beam2_heavy,k)
                      endif
                  elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqvall(j,k,1,1) * &
                          singletop2_pdfs(i_beam1_heavy,j) * &
                          singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k,1,1) * &
                              singletop2_pdfs(i_beam1_heavy,j) * &
                              singletop2_pdfs(i_beam2_light,k)
                      endif
                  endif
              enddo; enddo
          elseif (currentContrib == 4) then
              if (ipsgen == 4) then
              do j=-nf,nf; do k=-nf,nf
                  if ( any(j == [2,4,-1,-3]) .and. k == 5 ) then
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqvall(j,k, 1,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      !if (coeffonly .eqv. .false.) then
                          !xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      !endif
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + msqall(j,k, 1,1)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)
                  elseif ( any(j == [2,4,-1,-3]) .and. k == 0 ) then
                      xmsq(WITH_BBAR,1) = xmsq(WITH_BBAR,1) + msqvall(j,k, WITH_BBAR,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      !if (coeffonly .eqv. .false.) then
                          !xmsq(WITH_BBAR,1) = xmsq(WITH_BBAR,1) + msqall(j,k, 1,1) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      !endif
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(WITH_BBAR,1) = xmsq(WITH_BBAR,1) + msqall(j,k, WITH_BBAR,1)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      xmsq(WITH_BBAR,1) = xmsq(WITH_BBAR,1) + msqall(j,k, WITH_BBAR,1)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)
                  elseif ( j == 5  .and. any(k == [2,4,-1,-3])) then
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqvall(j,k, 1,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      !if (coeffonly .eqv. .false.) then
                          !xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      !endif
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + msqall(j,k, 1,2)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
                  elseif ( j == 0  .and. any(k == [2,4,-1,-3])) then
                      xmsq(WITH_BBAR,2) = xmsq(WITH_BBAR,2) + msqvall(j,k, WITH_BBAR,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      !if (coeffonly .eqv. .false.) then
                          !xmsq(WITH_BBAR,2) = xmsq(WITH_BBAR,2) + msqall(j,k, 1,2) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      !endif
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(WITH_BBAR,2) = xmsq(WITH_BBAR,2) + msqall(j,k, WITH_BBAR,2)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      xmsq(WITH_BBAR,2) = xmsq(WITH_BBAR,2) + msqall(j,k, WITH_BBAR,2)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
                  elseif (j==5 .and. k==0) then
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(ONLY_MAIN_B,2) = xmsq(ONLY_MAIN_B,2) + sum(msqall(j,[2,4,-1,-3], 1, 2)) * (AP(q,g,2) + Q2(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z
                  elseif (j==0 .and. k==5) then
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(ONLY_MAIN_B,1) = xmsq(ONLY_MAIN_B,1) + sum(msqall([2,4,-1,-3],k, 1,1)) * (AP(q,g,2) + Q1(q,g,q,2)) * &
                          singletop2_pdfs(i_beam1_light_z,j)/z * singletop2_pdfs(i_beam2_heavy,k)
                  elseif (j==0 .and. k==0) then
                      call singletop2_fillAP(z, i_beam2_light, AP)
                      xmsq(WITH_BBAR,2) = xmsq(WITH_BBAR,2) + sum(msqall(j,[2,4,-1,-3], WITH_BBAR, 2)) * (AP(q,g,2) + Q2(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light_z,k)/z

                      call singletop2_fillAP(z, i_beam1_light, AP)
                      xmsq(WITH_BBAR,1) = xmsq(WITH_BBAR,1) + sum(msqall([2,4,-1,-3],k, WITH_BBAR,1)) * (AP(q,g,2) + Q1(q,g,q,2)) * &
                          singletop2_pdfs(i_beam1_light_z,j)/z * singletop2_pdfs(i_beam2_heavy,k)
                  endif
              enddo; enddo

              endif
          endif
                  
          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)

          virtint_assemble_kbq_tpq = 0._dp
          where (xmsq /= 0._dp) virtint_assemble_kbq_tpq(:,:) = &
                  xmsq(:,:)*(2*sqrt(z))*flux*pswt/BrnRat

          return
          
      end function

      function virtint_assemble_kbq_tpq_jet(msqall,msqvall,z,bcontribs,pswt,xx)
          implicit none
          include 'nf.f'
          include 'kpart.f'! coeffonly
          include 'constants.f'
          include 'energy.f'
          include 'PR_new.f'
          include 'PR_stop.f'
          include 'agq.f'
          include 'epinv.f'
          real(dp) :: BrnRat
          common/BrnRat/BrnRat

          real(dp) :: virtint_assemble_kbq_tpq_jet(max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: msqvall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)
          real(dp), intent(in) :: z, pswt, xx(2)
          logical, intent(in) :: bcontribs(max_bcontrib)

          real(dp) :: xmsq(max_corr_on_beam)
          real(dp) :: AP(-1:1, -1:1, 3)
          integer :: j,k
          real(dp) :: flux

          logical :: bin
          common/bin/bin

          virtint_assemble_kbq_tpq_jet = 0._dp

          ! "deliberately" not compatible with t~ (no abs(5))
          ! and explicitly constrained to the available partonic channels
          ! for maximum testability and readability

          xmsq = 0._dp

          if (currentContrib == 1) then
              do j=-nf,nf; do k=-nf,nf
                  if (j == 5 .and. k == 5 ) then
                      call singletop2_fillAP(z, i_beam1_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqall(g,k,1,1)*(AP(g,q,2)+B1(g,q,b,2)) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)

                      call singletop2_fillAP(z, i_beam2_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          msqall(j,g,1,2)*(AP(g,q,2)+B2(g,q,b,2)) * &
                              singletop2_pdfs(i_beam1_heavy,j)*singletop2_pdfs(i_beam2_light_z,k)/z

                  elseif (j /= 0 .and. k == 5 ) then
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          (msqvall(j,k,1,1)) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                              (msqall(j,k,1,1)) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      endif

                      call singletop2_fillAP(z, i_beam1_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqall(j,k,1,1)*(AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          (msqall(j,k,1,1)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3)) &
                                   + msqall(g,k,1,1)*(AP(g,q,2)+B1(g,q,b,2))) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)

                  elseif (j == 0 .and. k == 5 ) then
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          (msqvall(j,k,1,1)) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                              (msqall(j,k,1,1)) * singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      endif

                      call singletop2_fillAP(z, i_beam1_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqall(j,k,1,1) * (AP(g,g,1)-AP(g,g,3)+B1(g,g,b,1)-B1(g,g,b,3)) * &
                              singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          (msqall(j,k,1,1) * (AP(g,g,2)+AP(g,g,3)+B1(g,g,b,2)+B1(g,g,b,3)) &
                                   + sum(msqall([-1,-3,2,4],k, 1,1)) * (AP(q,g,2) + B1(q,g,b,2))) * &
                              singletop2_pdfs(i_beam1_light_z,j)/z*singletop2_pdfs(i_beam2_heavy,k)

                  elseif (j == 5 .and. k /= 0 ) then
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          (msqvall(j,k,1,2)) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                              (msqall(j,k,1,2)) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      endif

                      call singletop2_fillAP(z, i_beam2_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          msqall(j,k,1,2)*(AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          (msqall(j,k,1,2)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3)) &
                                   + msqall(j,g,1,2)*(AP(g,q,2)+B2(g,q,b,2))) *&
                              singletop2_pdfs(i_beam1_heavy,j)*singletop2_pdfs(i_beam2_light_z,k)/z

                  elseif (j == 5 .and. k == 0  ) then
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          (msqvall(j,k,1,2)) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                              (msqall(j,k,1,2)) * singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      endif

                      call singletop2_fillAP(z, i_beam2_light, AP)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          msqall(j,k,1,2) * (AP(g,g,1)-AP(g,g,3)+B2(g,g,b,1)-B2(g,g,b,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          (msqall(j,k,1,2) * (AP(g,g,2)+AP(g,g,3)+B2(g,g,b,2)+B2(g,g,b,3)) &
                                   + sum(msqall(j,[-1,-3,2,4], 1,2)) * (AP(q,g,2) + B2(q,g,b,2))) * &
                              singletop2_pdfs(i_beam1_heavy,j)*singletop2_pdfs(i_beam2_light_z,k)/z
                  endif

              enddo; enddo
          elseif (currentContrib == 2) then
              do j=-nf,nf; do k=-nf,nf
                  if (any(j == [2,4,-1,-3]) .and. k == 5) then
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          msqvall(j,k, 1,2) * singletop2_pdfs(i_beam1_light,j) * &
                          singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                              msqall(j,k, 1,2) * singletop2_pdfs(i_beam1_light,j) * &
                              singletop2_pdfs(i_beam2_heavy,k)
                      endif

                      call singletop2_fillAP(z, i_beam2_heavy, AP)
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                          msqall(j,k, 1,2)*(AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)) * &
                            singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)

                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                           msqall(j,k, 1,2)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3)) * &
                         singletop2_pdfs(i_beam1_light,j)*singletop2_pdfs(i_beam2_heavy_z,k)/z

                      virtint_assemble_kbq_tpq_jet(3,2) = virtint_assemble_kbq_tpq_jet(3,2) + &
                           msqall(j,0, 3,2)*(AP(g,q,2)+B2(g,q,q,2)) * &
                         singletop2_pdfs(i_beam1_light,j)*singletop2_pdfs(i_beam2_heavy_z,k)/z

                  elseif (any(j == [2,4,-1,-3]) .and. k == 0) then
                      virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) + &
                          msqvall(j,k, WITH_BBAR,2) * singletop2_pdfs(i_beam1_light,j) * &
                          singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) + &
                              msqall(j,k, WITH_BBAR,2) * singletop2_pdfs(i_beam1_light,j) * &
                              singletop2_pdfs(i_beam2_heavy,k)
                      endif

                      call singletop2_fillAP(z, i_beam2_heavy, AP)
                       virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) + &
                           msqall(j,k, WITH_BBAR,2) * (AP(g,g,1)-AP(g,g,3)+B2(g,g,q,1)-B2(g,g,q,3)) * &
                               singletop2_pdfs(i_beam1_light,j) * singletop2_pdfs(i_beam2_heavy,k)

                       virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,2) + &
                            msqall(j,k, WITH_BBAR,2) * (AP(g,g,2)+AP(g,g,3)+B2(g,g,q,2)+B2(g,g,q,3)) * &
                               singletop2_pdfs(i_beam1_light,j)*singletop2_pdfs(i_beam2_heavy_z,k)/z

                       virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + &
                            msqall(j,5, 1,2) * (AP(q,g,2) + B2(q,g,q,2)) * &
                               singletop2_pdfs(i_beam1_light,j)*singletop2_pdfs(i_beam2_heavy_z,k)/z

                  elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqvall(j,k, 1,1) * singletop2_pdfs(i_beam1_heavy,j) * &
                          singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                              msqall(j,k, 1,1) * singletop2_pdfs(i_beam1_heavy,j) * &
                              singletop2_pdfs(i_beam2_light,k)
                      endif

                      call singletop2_fillAP(z, i_beam1_heavy, AP)
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqall(j,k, 1,1)*(AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)

                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                          msqall(j,k, 1,1)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z*singletop2_pdfs(i_beam2_light,k)

                      virtint_assemble_kbq_tpq_jet(3,1) = virtint_assemble_kbq_tpq_jet(3,1) + &
                                     msqall(0,k, 3,1)*(AP(g,q,2)+B1(g,q,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z*singletop2_pdfs(i_beam2_light,k)

                  elseif (j == 0 .and. any(k == [2,4,-1,-3])) then
                      virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) + &
                          msqvall(j,k, WITH_BBAR,1) * singletop2_pdfs(i_beam1_heavy,j) * &
                          singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) + &
                              msqall(j,k, WITH_BBAR,1) * singletop2_pdfs(i_beam1_heavy,j) * &
                              singletop2_pdfs(i_beam2_light,k)
                      endif

                      call singletop2_fillAP(z, i_beam1_heavy, AP)
                      virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) + &
                          msqall(j,k, WITH_BBAR,1) * (AP(g,g,1)-AP(g,g,3)+B1(g,g,q,1)-B1(g,g,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy,j) * singletop2_pdfs(i_beam2_light,k)

                      virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) = virtint_assemble_kbq_tpq_jet(WITH_BBAR,1) + &
                          msqall(j,k, WITH_BBAR,1) * (AP(g,g,2)+AP(g,g,3)+B1(g,g,q,2)+B1(g,g,q,3)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z*singletop2_pdfs(i_beam2_light,k)

                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + &
                                     msqall(5,k, 1,1) * (AP(q,g,2) + B1(q,g,q,2)) * &
                              singletop2_pdfs(i_beam1_heavy_z,j)/z*singletop2_pdfs(i_beam2_light,k)
                  endif

                  if (any(j == [2,4,-1,-3]) .and. any(k == [4,3,2,1,-1,-2,-3,-4,-5])) then
                      call singletop2_fillAP(z, i_beam2_heavy, AP)
                      virtint_assemble_kbq_tpq_jet(3,2) = virtint_assemble_kbq_tpq_jet(3,2) + &
                          msqall(j,0, 3,2)*(AP(g,q,2)+B2(g,q,q,2)) * &
                          singletop2_pdfs(i_beam1_light,j)*singletop2_pdfs(i_beam2_heavy_z,k)/z
                  endif
                  if (any(j == [4,3,2,1,-1,-2,-3,-4,-5]) .and. any(k == [2,4,-1,-3])) then
                      call singletop2_fillAP(z, i_beam1_heavy, AP)
                      virtint_assemble_kbq_tpq_jet(3,1) = virtint_assemble_kbq_tpq_jet(3,1) + &
                          msqall(0,k,3,1)*(AP(g,q,2)+B1(g,q,q,2)) * &
                          singletop2_pdfs(i_beam1_heavy_z,j)/z*singletop2_pdfs(i_beam2_light,k)
                  endif
              enddo; enddo
          elseif (currentContrib == 3) then
              do j=-nf,nf; do k=-nf,nf
                  if (any(j == [2,4,-1,-3]) .and. k == 5) then
                      virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + msqvall(j,k,1,2) * &
                          singletop2_pdfs(i_beam1_light,j) * &
                          singletop2_pdfs(i_beam2_heavy,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,2) = virtint_assemble_kbq_tpq_jet(1,2) + msqall(j,k, 1,2) * &
                              singletop2_pdfs(i_beam1_light,j) * &
                             singletop2_pdfs(i_beam2_heavy,k)
                      endif
                  elseif (j == 5 .and. any(k == [2,4,-1,-3])) then
                      virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + msqvall(j,k,1,1) * &
                          singletop2_pdfs(i_beam1_heavy,j) * &
                          singletop2_pdfs(i_beam2_light,k)
                      if (coeffonly .eqv. .false.) then
                          virtint_assemble_kbq_tpq_jet(1,1) = virtint_assemble_kbq_tpq_jet(1,1) + msqall(j,k,1,1) * &
                              singletop2_pdfs(i_beam1_heavy,j) * &
                              singletop2_pdfs(i_beam2_light,k)
                      endif
                  endif
              enddo; enddo
          endif

          flux = fbGev2/(2._dp*xx(1)*xx(2)*sqrts**2)

          virtint_assemble_kbq_tpq_jet = virtint_assemble_kbq_tpq_jet * &
              (2*sqrt(z))*flux*pswt/BrnRat

      end function virtint_assemble_kbq_tpq_jet

end module
