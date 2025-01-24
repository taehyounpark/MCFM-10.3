!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop2_nnlo
      use singletop2_nnlo_vars
      use singletop2_scet_light
      use singletop2_scet_heavy_prod
      use singletop2_scet_heavy_decay
      use singletop_phase
      use singletop_jet
      use singletop_jet2
      use singletop_jet3
      use singletop_interf_lxh
      use singletop_interf_lxd
      use singletop_interf_hxd
      implicit none
      private

      public :: singletop2_scet_tree
      public :: singletop2_scet_tree_ub
      public :: singletop2_scet_tree_bu

      ! light line
      public :: singletop2_scet_virt_light
      public :: singletop2_scet_virt_light_all
      public :: singletop2_scet_z
      public :: singletop2_scet_gsall
      ! heavy line decay
      public :: singletop2_scet_virt_heavy_decay_all
      public :: singletop2_heavy_decay_g_all
      public :: singletop2_heavy_decay_gs_all_new
      public :: singletop_jet_decay_all
      public :: singletop_jet_decay_virt_all
      public :: singletop_jet_decay_real_all
      public :: singletop_jet_decay_gs
      ! heavy line production
      public :: singletop2_scet_virt_heavy_prod_all
      public :: singletop2_scet_heavy_prod_z
      public :: singletop2_scet_heavy_prod_gs_all
      public :: qqb_tbb_g_heavy_all
      public :: singletop_jet_heavy_virt_all
      public :: singletop_jet_heavy
      public :: singletop_jet_heavy_cobswitch
      public :: singletop_jet_heavy_all
      public :: singletop_jet_heavy_z
      public :: singletop_jet_heavy_real_all
      public :: singletop_jet_heavy_gs_all

      ! light line jet pieces
      !public :: singletop_jet_light
      public :: singletop_jet_light_msqall
      public :: singletop_jet_light_virt_all
      public :: singletop_jet_light_z
      public :: singletop_jet_light_real_all
      public :: singletop_jet_light_gvec
      public :: singletop_jet_light_gs_all

      ! light line, SCET pieces
      public :: passed_taucut_light
      public :: lumxmsq_singletop
      public :: lumxmsq_singletop_new

      ! heavy line decay, SCET pieces
      public :: passed_taucut_decay
      public :: lumxmsq_singletop_decay_jetmass

      ! heavy line production, SCET pieces
      public :: passed_taucut_heavyprod
      public :: lumxmsq_singletop_prod

      public :: gen_singletop

      ! interference lxh
      public :: singletop_jet_light_heavy_rr_all
      public :: singletop_jet_light_heavy_rr_gs_all

      public :: singletop_jet_light_heavy_vr_all
      public :: singletop_jet_light_heavy_vr_z
      public :: qqb_tbb_g_heavy_all_swap

      public :: singletop_jet_light_heavy_vv
      public :: singletop_jet_light_heavy_vv_tree

      public :: passed_taucut_lxh

      ! interference lxd

      public :: singletop_light_decay_vv
      public :: singletop_light_decay_vv_tree

      public :: singletop_light_decay_vr
      public :: singletop_light_decay_vr_z

      public :: singletop_light_decay_rr
      public :: singletop_light_decay_rr_gs

      public :: singletop_light_decay_rv
      public :: singletop_light_decay_rv_gs

      public :: singletop_decay_real_swap

      ! interference hxd

      public :: singletop_heavy_decay_vv
      public :: singletop_heavy_decay_vv_tree

      public :: singletop_heavy_decay_vr
      public :: singletop_heavy_decay_vr_z

      public :: singletop_decay_real_hxd

      public :: singletop_heavy_decay_rv
      public :: singletop_heavy_decay_rv_gs
    
      public :: singletop_heavy_decay_rr
      public :: singletop_heavy_decay_rr_gs


      ! global variables, see below

      public :: singletop2_nnlo_enable_light
      public :: singletop2_nnlo_enable_heavy_prod
      public :: singletop2_nnlo_enable_heavy_decay
      public :: singletop2_nnlo_enable_lxh
      public :: singletop2_nnlo_enable_lxd
      public :: singletop2_nnlo_enable_hxd
      public :: singletop2_nnlo_fully_inclusive
      public :: corr_on_beam
      public :: max_corr_on_beam, max_bcontrib

end module

