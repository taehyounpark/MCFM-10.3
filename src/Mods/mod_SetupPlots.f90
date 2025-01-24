!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module MCFMSetupPlots
      implicit none

      public :: setup_plots
      public :: nplotter_new

      private

      interface 
          subroutine ibook(p,wt,ids,vals,wts)
              use types
              implicit none
              include 'mxpart.f'
              real(dp), intent(in) :: p(mxpart,4), wt
              integer, allocatable, intent(out) :: ids(:)
              real(dp), allocatable, intent(out) :: vals(:)
              real(dp), allocatable, intent(out) :: wts(:)
          end subroutine
      end interface

      procedure (ibook), pointer :: pbook => null()
          
      logical, save :: disableHistograms = .false.

      contains

      subroutine setup_plots()
          use nplotter_Z, only: setup_Z => setup, book_Z => book
          use nplotter_W, only: setup_W => setup, book_W => book
          use nplotter_singletop, only: setup_singletop => setup, book_singletop => book
          use nplotter_Higgs, only: setup_Higgs => setup, book_Higgs => book
          use nplotter_Higgs_to_ZGa, only: setup_Higgs_to_ZGa => setup, book_Higgs_to_Zga => book
          use nplotter_ZGamma, only: setup_ZGamma => setup, book_ZGamma => book
          use nplotter_ZZ, only: setup_ZZ => setup, book_ZZ => book
          use nplotter_WZ, only: setup_WZ => setup, book_WZ => book
          use nplotter_WW, only: setup_WW => setup, book_WW => book
          use nplotter_Diphoton, only: setup_Diphoton => setup, book_Diphoton => book
          use nplotter_WH_bbar_tautau, only: setup_WH_bbar_tautau => setup, &
                                             book_WH_bbar_tautau => book
          use nplotter_WH_gamgam, only: setup_WH_gamgam => setup, &
                                        book_WH_gamgam => book
          use nplotter_ZH_bbar_tautau_gamgam, only: setup_ZH_bbar_tautau_gamgam => setup, &
                                                    book_ZH_bbar_tautau_gamgam => book
          use nplotter_ZH_WW, only: setup_ZH_WW => setup, book_ZH_WW => book
          use MCFMPlotting
          use parseinput, only: cfg, cfg_get
          use types
          implicit none
          include 'nproc.f'
          include 'kprocess.f'

          integer :: one

          call cfg_get(cfg, "extra%nohistograms", disableHistograms)

          one = plot_setup_uniform(0._dp,1._dp,1._dp,'total cross')

          if (disableHistograms) then
              return
          endif

          if (kcase==kW_only .or. kcase == kW_2jet .or. kcase == kW_1jet) then
              call setup_W()
              pbook => book_W
          elseif (kcase==kZ_only .or. kcase==kgg2lep .or. kcase == kZ_2jet .or. kcase == kZ_1jet) then
              call setup_Z()
              pbook => book_Z
          elseif (kcase == kbq_tpq .or. kcase == kbq_tpq_jet) then
              call setup_singletop()
              pbook => book_singletop
          elseif ((kcase == kggfus0) .or. (kcase == kggfus1)) then
              call setup_Higgs()
              pbook => book_Higgs
          elseif (kcase == kHi_Zga .or. kcase == kHi_Zaj) then
              call setup_Higgs_to_Zga()
              pbook => book_Higgs_to_Zga
          elseif ((kcase==kZGamma) .or. (kcase == kZgajet) .or. (kcase==kWgamma) .or. (kcase==kWgajet)) then
              call setup_ZGamma()
              pbook => book_ZGamma
          elseif ((kcase == kgamgam) .or. (kcase == kgmgmjt) .or. (kcase==kgg2gam)) then
              call setup_Diphoton()
              pbook => book_Diphoton
          elseif ((kcase==kZZlept) .or. (kcase==kHZZ_tb) .or. (kcase==kHZZint) .or.(kcase==kHZZHpi) .or. (kcase==kggZZ4l) .or. (kcase==kggZZbx)) then
              call setup_ZZ()
              pbook => book_ZZ
          elseif (kcase==kWZbbar) then
              call setup_WZ()
              pbook => book_WZ
          elseif ((kcase==kWWqqbr) .or. (kcase==kWW2jet)) then
              call setup_WW()
              pbook => book_WW
          elseif (kcase == kWHbbar) then
              call setup_WH_bbar_tautau()
              pbook => book_WH_bbar_tautau
          elseif (kcase == kWHgaga) then
              call setup_WH_gamgam()
              pbook => book_WH_gamgam
          elseif (kcase == kZH__WW .or. nproc == 623) then
              ! above cut of this process is also kZH1jet, we therefore have to
              ! select by nproc
              call setup_ZH_WW()
              pbook => book_ZH_WW
          elseif (kcase == kZHbbar .or. kcase == kZHgaga .or. kcase == kZH1jet) then
              call setup_ZH_bbar_tautau_gamgam()
              pbook => book_ZH_bbar_tautau_gamgam
          endif

          ! could also add auto-histograms here

      end subroutine

      subroutine nplotter_new(p,wt)
          use types
          use MCFMPlotting
          implicit none
          include 'mxpart.f'
          real(dp), intent(in) :: p(mxpart,4), wt

          integer, allocatable :: ids(:)
          real(dp), allocatable :: vals(:)
          real(dp), allocatable :: wts(:)

          ! could call plotting for auto histograms here

          ! call process specific plotting
          if (disableHistograms .or. (.not. associated(pbook))) then
              call plot_book([1],[0.5_dp],wt,[wt])
          else
              call pbook(p,wt,ids,vals,wts)
              ! equal weights for all observables if not passed out

              if (.not. allocated(wts)) then
                  allocate(wts(1+size(ids)), source=wt)
              endif
              call plot_book([1,ids],[0.5_dp,vals],wt,[wt,wts])
          endif

      end subroutine

end module
