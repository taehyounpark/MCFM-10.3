!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter(p,wt,wt2,nd)
          use types
          use MCFMStorage
          use SCET
          use omp_lib
          use parseinput, only: cfg
          use m_config
          use ieee_arithmetic
      implicit none
c--- Variable passed in to this routine:

c---      p:  4-momenta of leptons and jets in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)

c---     wt:  weight of this event

c---    wt2:  weight^2 of this event

c---     nd:  an integer:: specifying the dipole number of this contribution
c---          (if applicable), otherwise equal to zero

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'nplot.f'
      include 'nproc.f'
      include 'first.f'
      include 'taucut.f'
      include 'lhcb.f'


      real(dp), intent(in) :: p(mxpart,4),wt,wt2
      integer, intent(in) :: nd
      integer :: switch
      integer, save :: plotindex
!$omp threadprivate(plotindex)

      integer :: nplotmax
      common/nplotmax/nplotmax

      logical,save :: disableHistograms = .false.

c--- This routine simply picks out a process-specific plotting routine
c---  (if available) and falls back to the generic routine otherwise.
c---  So far available: W_only, Z_only, Wbbbar, Wbbmas, WpWp2j, WpWp3j
c---  For the convenience of the user who wants to bail out and do their
c---  own plotting we provide the dummy routine userplotter

c---  the index of the plot - stored in the nplot.f common and used for
c---  keeping track of the plot index across different files (GPS)

      ! do NOT remove, see SCET module
      ! used to transfer current dipole contribution information
      ! from nplotter to bookplot without changing the bookplot interface
      currentNd = nd

      if (ieee_is_nan(wt)) then
          write (*,*) "NaN in nplotter, this should never happen!"
          return
      endif

      ! WARNING: do NOT remove this first total cross section histogram
      ! it is used for PDF and scale uncertainties
      if (first) then
           call bookplot(1,1,'total cross', 0.5d0,0d0,0d0,0d0,1d0,1d0,'lin')
      else
           call bookplot(1,2,'total cross', 0.5d0,wt,wt**2,0d0,1d0,1d0,'lin')
      endif

      if (first) then
          call cfg_get(cfg, "extra%nohistograms", disableHistograms)
      endif

      ! for debugging, only consider total cross histogram:
      if (disableHistograms) then
          nplotmax = 1
          goto 888
      endif

      nextnplot = 2

c---  first allow for user plots
      call nplotter_user(p,wt,wt2,nd)

c LHCb plots if in that mode
      if (cut_mode > 0) then
        call lhcb_plots(p,wt,wt2,nd)
        goto 888
      endif

c--- switch:  an integer:: equal to either 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      if (nd>0) then
         switch=1
      else
         switch=0
      endif

      if (kcase==ktopanom .or. kcase==kbq_tpq_jet .or. nproc == 161 .or. nproc == 162) then
        call nplotter_ktopanom(p,wt,wt2,nd)
        goto 888
      endif

      if (nprocbelow == 300 .or. nprocbelow == 305) then
        call nplotter_Vgamma(p,wt,wt2,switch,nd)
        goto 888
      endif

c--- Special plotting routine for WW -> leptons
      if((nproc==61).or.(nproc==66).or.(nproc==126).or.(nproc==127)) then
         call nplotter_VV(p,wt,wt2,switch,0)
         goto 888
      endif

c----- special plotting routine for VHgaga
      if((nproc==93).or.(nproc==98).or.(nproc==104)
     &   .or.(nproc==612).or.(nproc==617).or.(nproc==622)
     &   .or.(nproc==370).or.(nproc==371).or.(nproc==301)) then
         call nplotter_VHgaga(p,wt,wt2,switch,nd)
         goto 888
      endif

c====== special plotting routine for VHWW
      if((nproc==94).or.(nproc==99).or.(nproc==106)
     &   .or.(nproc==613).or.(nproc==618).or.(nproc==623)) then
         call nplotter_VHWW(p,wt,wt2,switch,nd)
         goto 888
      endif

c--- work out which plotting routine to use when first called
c-----> saves string comparison in general and important for combining SCET
      if (first) then
        first = .false.
        if     ((kcase==kW_only) .or. (kcase==kWln_ew)) then
          plotindex=1
c          plotindex=1000 ! revert to default plotting routine
        elseif (kcase==kZ_only .or. kcase==kgg2lep) then
          plotindex=2
        elseif (kcase==kW_cjet) then
          plotindex=3
        elseif (kcase==kWbbbar) then
          plotindex=4
        elseif ((kcase==kWbbmas) .or. (kcase==kW_bjet)) then
          plotindex=5
        elseif ((kcase==kWpWp2j) .or. (kcase==kWpWp3j))then
          plotindex=6
        elseif ((kcase==kW_1jet) .or. (kcase==kW_2jet)
     &     .or. (kcase==kW_3jet)) then
          plotindex=7
        elseif ((kcase==kHWW_4l) .or. (kcase==kHWW_tb)
     &     .or. (kcase==kHWW2lq) .or. (kcase==kHWWint)
     &     .or. (kcase==kHWWHpi) .or. (kcase==kggWW4l)
     &     .or. (kcase==kWWqqbr) .or. (kcase==kggWWbx)
     &     .or. (kcase==kWW_jet) .or. (kcase==kWZbbar)
     &     .or. (kcase==kWZ_jet)) then
          plotindex=8
        elseif (kcase==kWW2jet) then
          plotindex=9
c--- photon processes also need to know the dipole number
        elseif ((kcase==kWgamma) .or. (kcase==kZgamma)
     &     .or. (kcase==kWgajet) .or. (kcase==kWgaj_a)
     &     .or. (kcase==kWga_ew) .or. (kcase==kWgajew)
     &     .or. (kcase==kWgajja)) then
          plotindex=10
        elseif ((kcase==kgamgam) .or. (kcase==kgg2gam)) then
          plotindex=11
        elseif (kcase==kgmgmjt) then
          if (usescet) then
            plotindex=11 ! same routine below and above cut for SCET
          else
            plotindex=12
          endif
        elseif (kcase==kdirgam) then
          plotindex=13
        elseif (kcase==ktrigam) then
          plotindex=14
        elseif (kcase==kW_2gam)  then
          plotindex=15
        elseif (kcase==kZ_2gam)  then
          plotindex=16
        elseif (kcase==kZgajet)  then
          plotindex=17
        elseif ((kcase==ktt_bbl) .or. (kcase==ktt_ldk)
     &     .or. (kcase==ktt_bbh) .or. (kcase==ktt_bbu)
     &     .or. (kcase==ktt_hdk) .or. (kcase==ktthWdk)
     &     .or. (kcase==ktt_udk)) then
          plotindex=18
        elseif ((kcase==k4ftwdk) .or. (kcase==kdk_4ft)) then
          plotindex=19
        elseif ((kcase==kt_bbar) .or. (kcase==ktdecay)) then
          plotindex=20
        elseif (kcase==kqq_ttw) then
          plotindex=21
        elseif ((kcase==kH_tjet) .or. (kcase==kZ_tjet)) then
          plotindex=22
        elseif ((kcase==kH_tdkj) .or. (kcase==kZ_tdkj)) then
          plotindex=23
        elseif (kcase==kqqtthz) then
          plotindex=24
        elseif ((kcase==kdm_jet).or.(kcase==kdm2jet)) then
          plotindex=25
        elseif ((kcase==kdm_gam).or.(kcase==kdm_gaj)) then
          plotindex=26
        elseif ((kcase==kqqZZqq).or.(kcase==kqqWWqq)
     &     .or. (kcase==kqqVVqq).or.(kcase==kqqWWss)
     &     .or. (kcase==kqqWZqq).or.(kcase==kWpmZjj)
     &     .or.(kcase==kqq_ttg)) then
          plotindex=27
        elseif ((kcase==kHZZ_4l)
     &    .or.  (kcase==kHZZ_tb)
     &    .or.  (kcase==kHZZint)
     &    .or.  (kcase==kHZZHpi)
     &    .or.  (kcase==kggZZ4l)
     &    .or.  (kcase==kggZZbx)
     &    .or.  (kcase==kHZZqgI)
     &    .or.  (kcase==kZZlept)
     &    .or.  (kcase==kVVlept)
     &    .or.  (kcase==kggVV4l)
     &    .or.  (kcase==kggVVbx)
     &    .or.  (kcase==kHVV_tb)
     &    .or.  (kcase==kZZ_jet)) then
           plotindex=28
        elseif((kcase == kWHbbar).or.(kcase == kZHbbar)
     &     .or.(kcase == kWH1jet).or.(kcase == kZH1jet)
     &     .or.(kcase == kWHbbdk).or.(kcase == kZHbbdk)
     &          ) then
c           plotindex=29
           plotindex=1000
        elseif (kcase==kZbbbar) then
           plotindex=30
        elseif ((kcase==ktt_tot) .or. (kcase==kbb_tot)
     &     .or. (kcase==kcc_tot) .or. (kcase==ktt_mix)) then
           plotindex=31
        elseif ((kcase==ktwojet).or.(kcase==ktwo_ew)) then
           plotindex=32
        elseif ((kcase==kZ_1jet) .or. (kcase==kZ_2jet)
     &     .or. (kcase==kZ_3jet)) then
           plotindex=33
        elseif ((kcase==kggfus0) .or. (kcase==kggfus1)
     &     .or. (kcase==kggfus2) .or. (kcase==kggfus3)
     &     .or. (kcase==kggscl0) .or. (kcase==khjetma)
     &     .or. (kcase==khttjet) .or. (kcase==khttscl)
     &     .or. (kcase==kh2jmas) .or. (kcase==kh2jscl)) then
          plotindex=34
       else
          plotindex=1000
        endif
      endif

      select case (plotindex)
      case (1)
        call nplotter_W_only(p,wt,wt2,switch)
      case (2)
        call nplotter_Z_only(p,wt,wt2,switch)
      case (3)
        call nplotter_Wbbmas(p,wt,wt2,switch)
      case (4)
        call nplotter_Wbbmas(p,wt,wt2,switch)
      case (5)
        call nplotter_Wbbmas(p,wt,wt2,switch)
      case (6)
        call nplotter_WpWp(p,wt,wt2,switch)
      case (7)
        call nplotter_Wjets(p,wt,wt2,switch)
      case (8)
        call nplotter_VV(p,wt,wt2,switch,0)
      case (9)
        call nplotter_WW_jet(p,wt,wt2,switch)
      case (10)
        call nplotter_Vgamma(p,wt,wt2,switch,nd)
      case (11)
        call nplotter_gamgam(p,wt,wt2,switch,nd)
      case (12)
        call nplotter_gmgmjt(p,wt,wt2,switch)
      case (13)
        call nplotter_dirgam(p,wt,wt2,switch,nd)
      case (14)
        call nplotter_trigam(p,wt,wt2,switch)
      case (15)
        call nplotter_wgamgam(p,wt,wt2,switch,nd)
      case (16)
        call nplotter_zgamgam(p,wt,wt2,switch,nd)
      case (17)
        call nplotter_zgamjet(p,wt,wt2,switch,nd)
      case (18)
        call nplotter_ttbar(p,wt,wt2,switch)
      case (19)
        call nplotter_4ftwdk(p,wt,wt2,switch)
      case (20)
        call nplotter_tbbar(p,wt,wt2,switch)
      case (21)
        call nplotter_ttw(p,wt,wt2,switch)
      case (22)
        call nplotter_Ztj(p,wt,wt2,switch)
      case (23)
        call nplotter_Ztjdk(p,wt,wt2,switch)
      case (24)
        call nplotter_ttZ(p,wt,wt2,switch)
      case (25)
         call nplotter_dm_monj(p,wt,wt2,switch)
      case (26)
         call nplotter_dm_mongam(p,wt,wt2,switch,nd)
      case (27)
         call nplotter_qqZZqq(p,wt,wt2,switch)
      case (28)
         call nplotter_ZZlept(p,wt,wt2,switch)
      case (29)
c         call nplotter_VHbbar(p,wt,wt2,switch,nd)
         call nplotter_VHbbarHXSWG(p,wt,wt2,switch,nd)
      case (30)
         call nplotter_Zbbbar(p,wt,wt2,switch,nd)
c         call nplotter_VHgaga(p,wt,wt2,switch,nd)
      case (31)
         call nplotter_tt_tot(p,wt,wt2,switch,nd)
      case (32)
         call nplotter_twojet(p,wt,wt2,switch,nd)
      case (33)
         call nplotter_Zjets(p,wt,wt2,switch)
      case (34)
         call nplotter_Hjets(p,wt,wt2,switch)
      case (1000)
         call nplotter_auto(p,wt,wt2)
      case DEFAULT
        write(6,*) 'unexpected plotindex in nplotter =',plotindex
        stop
      end select

 888  continue

      first = .false.

      return
      end

