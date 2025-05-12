!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
module parseinput
  use m_config
  use types
  implicit none

  type(CFG_T), public, save :: cfg

  public

contains

  subroutine default_config()
    implicit none

    real(dp), parameter :: sqrts = 14000d0
    integer :: i

    call cfg_add(cfg, "mcfm_version", "10.3", "MCFM version number")
    call cfg_add(cfg, "writerefs", .true., "write references")

    ! [histogram]
    call cfg_add(cfg, "histogram%writetop", .false., "write top-drawer histograms")
    call cfg_add(cfg, "histogram%writetxt", .false., "write raw table file for each histogram")
    call cfg_add(cfg, "histogram%newstyle", .false., "use new histogram routines")

    ! [general]
    call cfg_add(cfg, "general%nproc", 1, "process number")
    call cfg_add(cfg, "general%part", "nlo", "part: lo, nlo, nlocoeff, nnlocoeff")
    call cfg_add(cfg, "general%runstring", "run", "string identifying the run")
    call cfg_add(cfg, "general%rundir", "./", "directory for output and snapshot files")
    call cfg_add(cfg, "general%sqrts", sqrts, "center of mass energy")
    call cfg_add(cfg, "general%ih1", +1, "ih1: +1 for proton, -1 for antiproton")
    call cfg_add(cfg, "general%ih2", +1, "ih2: +1 for proton, -1 for antiproton")
    call cfg_add(cfg, "general%zerowidth", .false., "use zero width approx. (proc. dependent)")
    call cfg_add(cfg, "general%removebr", .false., "remove decay branching ratio (proc. dependent)")
    call cfg_add(cfg, "general%ewcorr", "none", "electroweak corrections: none, sudakov or exact")

    ! [nnlo]
    call cfg_add(cfg, "nnlo%dynamictau", .true., "event-by-event tau definition")
    call cfg_add(cfg, "nnlo%tcutarray", [real(dp) ::], "optional array &
             &of taucut values that should be sampled on the fly in addition", dynamic_size=.true.)
    call cfg_add(cfg, "nnlo%useqt", .false., "use qt subtractions (implementation of 2207.07056)")
    call cfg_add(cfg, "nnlo%useqt_nnlo", .false., "use qt subtractions (implementation of 2202.07738)")
    call cfg_add(cfg, "nnlo%useqt_nnlo_useGLY", .false., "useGLY")
    call cfg_add(cfg, "nnlo%useptveto", .false., "use ptveto subtractions")
    call cfg_add(cfg, "nnlo%useptveto_useBNR", .true., "useBNR")

    ! [resummation]
    call cfg_add(cfg, "resummation%usegrid", .false., "use pre-generated grid files for beam functions")
    call cfg_add(cfg, "resummation%makegrid", .false., "create beam function grids and then exit")
    call cfg_add(cfg, "resummation%gridoutpath", "/tmp", "location for writing grids")
    call cfg_add(cfg, "resummation%gridinpath", "/usr/local/share/LHAPDF", "location for LHAPDF grids")

    call cfg_add(cfg, "resummation%res_range", [0._dp, 80._dp], "integration range of purely resummed contribution")
    call cfg_add(cfg, "resummation%resexp_range", [1._dp, 80._dp], "integration range of fixed-order expansion of resummed contribution")
    call cfg_add(cfg, "resummation%fo_cutoff", 1d0, "minimum qT for fixed-order calculation")
    call cfg_add(cfg, "resummation%transitionswitch", 0.4d0, "transition switch to be used in plotting routines")
    call cfg_add(cfg, "resummation%scalevar_rapidity", .false., "rapidity scale variation")

    ! additional optional settings in [general]
    !call cfg_add(cfg, "general%nevtrequested", 0, "") ! JC: removed 3/8/22 (no longer operational)
    call cfg_add(cfg, "general%vdecayid", .false., "")
    call cfg_add(cfg, "general%v34id", "", "")
    call cfg_add(cfg, "general%v56id", "", "")

    ! [masses]
    call cfg_add(cfg, "masses%hmass", 125d0, "Higgs mass")
    call cfg_add(cfg, "masses%mt", 173d0, "Top-quark mass")
    call cfg_add(cfg, "masses%mb", 4.66d0, "Bottom-quark mass")
    call cfg_add(cfg, "masses%mc", 1.275d0, "Charm-quark mass")


    call cfg_add(cfg, "masses%CKMdiag", .false., "Set CKM to diagonal and ignore CKMrotate")
    call cfg_add(cfg, "masses%CKMrotate", .false., "Implement VCKM by rotating pdfs")

    ! [scales]
    call cfg_add(cfg, "scales%renscale", 1d0, "Renormalization scale")
    call cfg_add(cfg, "scales%facscale", 1d0, "Factorization scale")
    call cfg_add(cfg, "scales%dynamicscale", "none", "Dynamic scale")
    call cfg_add(cfg, "scales%doscalevar", .false., "perform scale variation?")
    call cfg_add(cfg, "scales%maxscalevar", 6, "can be 2, 6 or 8 for 2, 6 or 8-point variation")
    call cfg_add(cfg, "scales%vetoscalevar", .false., "veto scale variation?")
    call cfg_add(cfg, "scales%timelikemusq", .false., "Time-like musq")

    ! [integration]
    call cfg_add(cfg, "integration%usesobol", .true., "Low Sobol low discrepancy sequence")
    call cfg_add(cfg, "integration%seed", 0, "Seed to use for pseudo random number integration")
    call cfg_add(cfg, "integration%precisiongoal", 0.002d0, "Relative precision goal")
    call cfg_add(cfg, "integration%readin", .false., "Resume from previous integration snapshot")
    call cfg_add(cfg, "integration%writeintermediate", .true., &
                   "write histograms and results after each vegas iteration")
    call cfg_add(cfg, "integration%warmupprecisiongoal", 0.25d0, &
                   "Relative precision goal for each contribution during warmup")
    call cfg_add(cfg, "integration%warmupchisqgoal", 2.5d0, &
                   "Chisq/it goal for each contribution during warmup")
    call cfg_add(cfg, "integration%scalarselect", 1, &
                   "Selection of library for one-loop integrals")

    ! [pdf]
    call cfg_add(cfg, "pdf%pdlabel", "CT14.NN", "PDF label for internal routines")

    ! [lhapdf]
    call cfg_add(cfg, "lhapdf%lhapdfset", ["CT14nnlo"], "LHAPDF PDF label", dynamic_size=.true.)
    call cfg_add(cfg, "lhapdf%lhapdfmember", [0], "LHAPDF PDF member number, -1 for PDF uncertainties", &
                               dynamic_size=.true.)
    call cfg_add(cfg, "lhapdf%dopdferrors", .false., "calculate PDF uncertainties.")

    ! [basicjets]
    call cfg_add(cfg, "basicjets%inclusive", .true., "jet-inclusive cross-section")
    call cfg_add(cfg, "basicjets%algorithm", "ankt", "jet-algorithm: ankt, ktal, cone, hqrk, none")
    call cfg_add(cfg, "basicjets%ptjetmin", 30d0, "minimum jet pT")
    call cfg_add(cfg, "basicjets%etajetmax", 2.4d0, "maximum jet rapidity, absolute value")
    call cfg_add(cfg, "basicjets%Rcutjet", 0.5d0, "minimum jet separation in R")
    call cfg_add(cfg, "basicjets%userap", .true. , "use rapidity not pseudorapidity for jets")

    ! [masscuts]
    call cfg_add(cfg, "masscuts%m34min", 0d0, "minimum mass of 3-4 system")
    call cfg_add(cfg, "masscuts%m56min", 0d0, "minimum mass of 5-6 system")
    call cfg_add(cfg, "masscuts%m3456min", 0d0, "minimum mass of 3-4-5-6 system")
    call cfg_add(cfg, "masscuts%pt34min", 0d0, "minimum pT of 3-4 system")
    ! max values are handled below

    call cfg_add(cfg, "masscuts%closestZ", .false., "identfy Z pairs using mass closest to Z")
    call cfg_add(cfg, "masscuts%moppmin", 0d0, "minimum invariant mass for opposite-sign leptons")

    ! [cuts]
    call cfg_add(cfg, "cuts%makecuts", .true., "enable these additional cuts")

    call cfg_add(cfg, "cuts%ptleptmin", 20d0, "minimum lepton pT")
    call cfg_add(cfg, "cuts%etaleptmax", 2.4d0, "maximum lepton rapidity, absolute value")
    call cfg_add(cfg, "cuts%etaleptveto", [0d0, 0d0], "lepton rapidity veto")
    call cfg_add(cfg, "cuts%ptminmiss", 30d0, "minimum missing pt")

    call cfg_add(cfg, "cuts%y34min", 0d0, "minimum 34 rapidity")
    call cfg_add(cfg, "cuts%y34max", 99d0, "maximum 34 rapidity")

    call cfg_add(cfg, "cuts%ptlept2min", 20d0, "additional leptons minimum pt")
    call cfg_add(cfg, "cuts%etalept2max", 2.4d0, "maximum rapidity of additional leptons, absolute value")
    call cfg_add(cfg, "cuts%etalept2veto", [0d0, 0d0], "additional lepton rapidity veto")

    call cfg_add(cfg, "cuts%m34transmin", 0d0, "minimum (3,4) transverse mass")
    call cfg_add(cfg, "cuts%Rjlmin", 0d0, "R(jet,lept)_min")
    call cfg_add(cfg, "cuts%Rllmin", 0d0, "R(lept,lept)_min")
    call cfg_add(cfg, "cuts%delyjjmin", 0d0, "Delta eta(jet,jet)_min")
    call cfg_add(cfg, "cuts%jetsopphem", .false., "force jets to be in opposite hemispheres")
    call cfg_add(cfg, "cuts%lbjscheme", 0, "lepbtwnjets scheme")
    call cfg_add(cfg, "cuts%ptbjetmin", 0d0, "b-jet minimum pT; can also add ptbjetmax")
    call cfg_add(cfg, "cuts%etabjetmax", 100d0, "b-jet maximum rapidity, absolute value; can also add ptbjetmin")

    call cfg_add(cfg, "cuts%elptmin", 0d0, "minimum electron pT")
    call cfg_add(cfg, "cuts%muptmin", 0d0, "minimum muon pT")
    call cfg_add(cfg, "cuts%Relelmin", 0d0, "R(electron,electron)_min")
    call cfg_add(cfg, "cuts%Relmumin", 0d0, "R(muon,muon)_min")
    call cfg_add(cfg, "cuts%Rmumumin", 0d0, "R(muon,muon)_min")

    call cfg_add(cfg, "cuts%Rlepiso", 0d0, "R(lept)_isolation")
    call cfg_add(cfg, "cuts%fraclepiso", 0d0, "lepton isolation fraction")

    call cfg_add(cfg, "cuts%elrapmax", 99d0, "maximum electron rapidity, absolute value")
    call cfg_add(cfg, "cuts%murapmax", 99d0, "maximum muon rapidity, absolute value")
    call cfg_add(cfg, "cuts%etaelveto", [0d0, 0d0], "electron rapidity veto")
    call cfg_add(cfg, "cuts%etamuveto", [0d0, 0d0], "muon rapidity veto")
    call cfg_add(cfg, "cuts%ptminmissrel", -1d0, "minimum missing ptrel")
    call cfg_add(cfg, "cuts%mllmin", 0d0, "minimum dilepton invariant mass")
    call cfg_add(cfg, "cuts%ptllmin", 0d0, "minimum dilepton transverse momentum")

    call cfg_add(cfg, "cuts%m3lmin", 0d0, "minimum mass of 3-lepton system")

    ! [photon]
    call cfg_add(cfg, "photon%fragmentation", .false., "fragmentation included")
    call cfg_add(cfg, "photon%fragmentation_set", "GdRG__LO", "fragmentation set")
    call cfg_add(cfg, "photon%fragmentation_scale", 1d0, "fragmentation scale")

    call cfg_add(cfg, "photon%gammptmin", 40d0, "minimum photon pT; can also add gammptmax")
    call cfg_add(cfg, "photon%gammrapmax", 2.5d0, "maximum photon rapidity; can also add gammrapmin")
    call cfg_add(cfg, "photon%gammrapveto", [0d0, 0d0], "photon rapidity veto")
    call cfg_add(cfg, "photon%gammpt2", 25d0, "second photon minimum pT")
    call cfg_add(cfg, "photon%gammpt3", 25d0, "third photon minimum pT")
    call cfg_add(cfg, "photon%Rgalmin", 0d0, "R(photon,lepton)_min")
    call cfg_add(cfg, "photon%Rgagamin", 0.4d0, "R(photon,photon)_min")
    call cfg_add(cfg, "photon%Rgajetmin", 0d0, "R(photon,jet)_min")
    call cfg_add(cfg, "photon%cone_ang", 0.4d0, "cone size for isolation")
    call cfg_add(cfg, "photon%epsilon_h", 0.5d0, "epsilon_h, energy fraction for isolation")
    call cfg_add(cfg, "photon%n_pow", 1d0, "n_pow, exponent for smooth-cone isolation")
    call cfg_add(cfg, "photon%fixed_coneenergy", .false., "use fixed cone energy epsilon_h in isolation")
    call cfg_add(cfg, "photon%hybrid", .false., "hybrid scheme with inner cone radius R_inner")
    call cfg_add(cfg, "photon%R_inner", 0.1d0, "inner cone radius for hybrid isolation")

! Removed: JC implementation of hybrid isolation
!    call cfg_add(cfg, "photon%hybridiso", .false., "hybrid isolation")
!    call cfg_add(cfg, "photon%eps_dyn", 0.5d0, "")
!    call cfg_add(cfg, "photon%R_dyn", 0.1d0, "")
!    call cfg_add(cfg, "photon%n_dyn", 1d0, "")
!    call cfg_add(cfg, "photon%eps_iso", 0d0, "")
!    call cfg_add(cfg, "photon%Etfix_iso", 4d0, "")
!    call cfg_add(cfg, "photon%R_iso", 0.4d0, "")

    ! [wz2jet]
    call cfg_add(cfg, "wz2jet%qflag", .true., "")
    call cfg_add(cfg, "wz2jet%gflag", .true., "")

    ! [hjetmass]
    call cfg_add(cfg, "hjetmass%mtex", 0, "controls approximation for 2-loop virtual corrections")

    ! [bsm_higgs]
    call cfg_add(cfg, "bsm_higgs%bsm", "eft", "Choose whether to vary EFT or DM parameters as BSM scenario")
    call cfg_add(cfg, "bsm_higgs%t1", 1d0, "t1")
    call cfg_add(cfg, "bsm_higgs%t2", 1d0, "t2")
    call cfg_add(cfg, "bsm_higgs%t3", 1d0, "t3")
    call cfg_add(cfg, "bsm_higgs%t4", 1d0, "t4")
    call cfg_add(cfg, "bsm_higgs%t5", 1d0, "t5")
    call cfg_add(cfg, "bsm_higgs%t6", 1d0, "t6")
    call cfg_add(cfg, "bsm_higgs%w1", 1d0, "w1")
    call cfg_add(cfg, "bsm_higgs%w2", 1d0, "w2")
    call cfg_add(cfg, "bsm_higgs%w3", 1d0, "w3")
    call cfg_add(cfg, "bsm_higgs%w4", 1d0, "w4")
    call cfg_add(cfg, "bsm_higgs%w5", 1d0, "w5")
    call cfg_add(cfg, "bsm_higgs%c6", [-20d0, 20d0, 10d0], "c6 BSM values")
    call cfg_add(cfg, "bsm_higgs%ct", [-1d0, 1d0, 1d0], "ct BSM values")
    call cfg_add(cfg, "bsm_higgs%cg", [-0d01, 0d01, 0d01], "cg SM value")
    call cfg_add(cfg, "bsm_higgs%c6_sm", 0d0, "c6 SM value")
    call cfg_add(cfg, "bsm_higgs%ct_sm", 0d0, "ct SM value")
    call cfg_add(cfg, "bsm_higgs%cg_sm", 0d0, "cg SM value")
    call cfg_add(cfg, "bsm_higgs%cx_sm", 0d0, "cx SM value")
    call cfg_add(cfg, "bsm_higgs%cx",    1d0, "cx BSM value")
    call cfg_add(cfg, "bsm_higgs%mx", 1000d0, "mx BSM value")

    ! [output]
    call cfg_add(cfg, "output%metadata", "metadata.csv", "Metadata output file") 
    call cfg_add(cfg, "output%events",   "events.csv",   "CSV output file") 

    ! [anom_higgs]
    call cfg_add(cfg, "anom_higgs%hwidth_ratio", 1d0, "Gamma_H / Gamma_H(SM)")
    call cfg_add(cfg, "anom_higgs%cttH", 1d0, "cttH")
    call cfg_add(cfg, "anom_higgs%cWWH", 1d0, "cWWH")

    ! [anom_wz]
    call cfg_add(cfg, "anom_wz%enable", .false., "enable anomalous W/Z couplings")
    call cfg_add(cfg, "anom_wz%delg1_z", 0d0, "Delta g1(Z)")
    call cfg_add(cfg, "anom_wz%delk_z", 0d0, "Delta K(Z)")
    call cfg_add(cfg, "anom_wz%delk_g", 0d0, "Delta K(gamma)")
    call cfg_add(cfg, "anom_wz%lambda_z", 0d0, "Lambda(Z)")
    call cfg_add(cfg, "anom_wz%lambda_g", 0d0, "Lambda(gamma)")
    call cfg_add(cfg, "anom_wz%h1Z", 0d0, "h1(Z)")
    call cfg_add(cfg, "anom_wz%h1gam", 0d0, "h1(gamma)")
    call cfg_add(cfg, "anom_wz%h2Z", 0d0, "h2(Z)")
    call cfg_add(cfg, "anom_wz%h2gam", 0d0, "h2(gamma)")
    call cfg_add(cfg, "anom_wz%h3Z", 0d0, "h3(Z)")
    call cfg_add(cfg, "anom_wz%h3gam", 0d0, "h3(gamma)")
    call cfg_add(cfg, "anom_wz%h4Z", 0d0, "h4(Z)")
    call cfg_add(cfg, "anom_wz%h4gam", 0d0, "h4(gamma)")
    call cfg_add(cfg, "anom_wz%tevscale", 2d0, "Form-factor scale, in TeV")

    ! [singletop]

    call cfg_add(cfg, "singletop%c_phiq", 0d0, "C_phiq (O1), real-valued")
    call cfg_add(cfg, "singletop%c_phiphi", [0d0,0d0], "C_phiphi (O2), real and imaginary part")
    call cfg_add(cfg, "singletop%c_tw", [0d0,0d0], "C_tW (O3), real and imaginary part")
    call cfg_add(cfg, "singletop%c_bw", [0d0,0d0], "C_bW (O4), real and imaginary part")
    call cfg_add(cfg, "singletop%c_tg", [0d0,0d0], "C_tG (O6), real and imaginary part")
    call cfg_add(cfg, "singletop%c_bg", [0d0,0d0], "C_bG (O7), real and imaginary part")
    call cfg_add(cfg, "singletop%lambda", 1000d0, "Lambda, scale of EFT breakdown in GeV")
    call cfg_add(cfg, "singletop%enable_lambda4", .false., "enable 1/Lambda^4 contributions")
    call cfg_add(cfg, "singletop%disable_sm", .false., "disable Standard Model contributions")
    call cfg_add(cfg, "singletop%mode_anomcoup", .false., "anomalous couplings mode (only LO)")

    ! [extra]
    call cfg_add(cfg, "extra%verbose", .false., "verbose output")
    call cfg_add(cfg, "extra%new_pspace", .true., "use new_pspace")
    call cfg_add(cfg, "extra%nohistograms", .false., "generate no additional histograms")
    call cfg_add(cfg, "extra%benchmark", "none", "benchmark to run")
    call cfg_add(cfg, "extra%pdfchannels", "", "constraints for PDF channels")
    call cfg_add(cfg, "extra%griddebug", .false., "print grid")
    call cfg_add(cfg, "extra%dummypdf", .false., "dummy pdfs x^0.1*(1-x)")
  end subroutine default_config

  subroutine read_config()
    use, intrinsic :: iso_c_binding
    use anomcoup_tbW
    use MCFMStorage, only: selectpdfs, griddebug
    use PDFerrors, only : doPDFerrors, PDFnames, PDFmembers, numPDFsets
    use Scalevar, only : doScalevar
    use m_gencuts, only : enable_reweight_user
    use MCFMSettings, only: newStyleHistograms, resexp_linPC
    use qtResummation_params, only: qtcutoff, qtminRes, qtmaxRes, qtminResexp, &
                                    qtmaxResexp, transitionSwitch, &
!                                    scalevar_rapidity, &
                                    scalevar_rapidity, &
            enable_fixed_y, fixed_y, enable_dsigma_dQ, generations, &
            fix_alphas_nf5
    use singletop2_nnlo, only: singletop2_nnlo_enable_light, &
        singletop2_nnlo_enable_heavy_prod, singletop2_nnlo_enable_heavy_decay, &
        singletop2_nnlo_enable_lxh, singletop2_nnlo_enable_lxd, &
        singletop2_nnlo_enable_hxd, singletop2_nnlo_fully_inclusive
    use ggHwilson, only: expansionorder, Wilsonorder
    use SCET, only: useQT
    use ptveto

    implicit none
    include 'cplx.h'
    include 'constants.f'
    include 'maxwt.f'! for nevtrequested
    include 'outputoptions.f'! writetop, etc.
    include 'vdecayid.f'! vdecayid, v34id, v56id
    include 'nproc.f'
    include 'masses.f'
    include 'runstring.f'
    include 'energy.f'! sqrts
    include 'zerowidth.f'
    include 'removebr.f'
    include 'pdlabel.f'
    include 'lhapdf.f'
    include 'clustering.f'
    include 'jetcuts.f'
    include 'Rcut.f'
    include 'makecuts.f'
    include 'lhcb.f'
    include 'limits.f'
    include 'leptcuts.f'
    include 'frag.f'
    include 'flags.f'! qflag, gflag
    include 'noglue.f'
    include 'realwt.f'
    include 'lc.f'
    include 'verbose.f'
    include 'debug.f'
    include 'new_pspace.f'
    include 'asymptotic.f'
    include 'alfacut.f'
    include 'betacut.f'
    include 'anomHiggs.f'
    include 'anom_higgs.f'
    include 'anomcoup.f'
    include 'scalevar.f'
    include 'ewcorr.f'
    include 'cutoff.f'
    include 'rtsmin.f'
    include 'kpart.f'
    include 'kprocess.f'
    include 'taucut.f'
    include 'ewinput.f'
    include 'squark.f'
    include 'Cabibbo.f'
!    include 'hybridiso.f'
    include 'beamtype.f'
    include 'yukawas.f'
    include 'mpicommon.f'
    include 'scalarselect.f'
    include 'interference.f'
    include 'userap.f'
    include 'bsm_higgs.f'
    include 'csvfile.f'
    character(len=cfg_string_len) :: mcfm_version, part, dynstring
    character(len=cfg_string_len) :: pdfchannels

    character(len=255) :: rundir
    common/rundir/rundir
    integer :: err_mkdir

    logical :: writerefs
    common/writerefs/writerefs

    logical :: spira
    common/spira/spira

    integer :: nmin, nmax
    common/nmin/nmin
    common/nmax/nmax

    real(dp) :: m34min, m34max, m56min, m56max
    real(dp) :: leptveto(2), lept2veto(2)

    real(dp) :: c_phiq_in, c_phiphi_in(2), c_tw_in(2), c_bg_in(2)
    real(dp) :: c_bw_in(2), c_tg_in(2), lambda_in

    integer :: pdfchan1, pdfchan2
    integer :: i

    logical :: dummypdf
    common/dummypdf/dummypdf

    real(dp) :: res_range(2), resexp_range(2), res_fo_cutoff

    interface
       ! posix mkdir
       function c_mkdir(pathname, mode) bind(C,name="mkdir")
         use iso_c_binding
         implicit none
         integer(c_int) :: c_mkdir
         character(kind=c_char), intent(in) :: pathname(*)
         integer(c_int16_t), value, intent(in) :: mode
       end function c_mkdir
    end interface

    call cfg_get(cfg, "mcfm_version", mcfm_version)
    if (mcfm_version(1:4) /= "10.3") then
       write (*,*) "Unsupported input file format version ", mcfm_version(1:4)
       stop
    endif

    ! [general]
    call cfg_get(cfg, "writerefs", writerefs)
    !call cfg_get(cfg, "general%nevtrequested", nevtrequested)    ! JC: removed 3/8/22 (no longer operational)
    call cfg_get(cfg, "general%vdecayid", vdecayid)
    call cfg_get(cfg, "general%v34id", v34id)
    call cfg_get(cfg, "general%v56id", v56id)

    call cfg_get(cfg, "general%nproc", nproc)
    call cfg_get(cfg, "general%part", part)
    call parse_part(part)
    call cfg_get(cfg, "general%runstring", runstring)

    call cfg_get(cfg, "general%rundir", rundir)
    if (trim(rundir) /= "./") then
       err_mkdir = c_mkdir(trim(rundir)//C_NULL_CHAR, int(o'755',c_int16_t)) 
       if (err_mkdir > 0) then
          write (*,*) "Error creating "//trim(rundir)// ": ", err_mkdir
          error stop
       endif
    endif

    call cfg_get(cfg, "general%sqrts", sqrts)
    call cfg_get(cfg, "general%ih1", ih1)
    call cfg_get(cfg, "general%ih2", ih2)
    call cfg_get(cfg, "general%zerowidth", zerowidth)
    call cfg_get(cfg, "general%removebr", removebr)
    call cfg_get_add(cfg, "general%ewscheme",ewscheme,ewscheme,"ew scheme")

    ! [nnlo]
    call cfg_get(cfg, "nnlo%dynamictau", dynamictau)
    call cfg_get(cfg, "nnlo%useqt", useQT)
    call cfg_get(cfg, "nnlo%useqt_nnlo", useqt_nnlo)
    call cfg_get(cfg, "nnlo%useqt_nnlo_useGLY", useGLY)
    call cfg_get(cfg, "nnlo%useptveto", usept)
    call cfg_get(cfg, "nnlo%useptveto_useBNR", useBNR)
    call parse_taucut()
    call parse_ewcorr()

    ! [histogram]
    call cfg_get(cfg, "histogram%writetop", writetop)
    call cfg_get(cfg, "histogram%writetxt", writeroot)
    call cfg_get(cfg, "histogram%newstyle", newStyleHistograms)

    if (nproc == 1610 .or. nproc == 1650) then
        newStyleHistograms = .true.
    endif

    ! [masses]
    call cfg_get(cfg, "masses%hmass", hmass)
    call cfg_get_add(cfg, "masses%wmass", wmass_inp, wmass_inp, "W mass")
    call cfg_get_add(cfg, "masses%zmass", zmass_inp, zmass_inp, "Z mass")
    call cfg_get_add(cfg, "masses%wwidth", wwidth, wwidth, "W width")
    call cfg_get_add(cfg, "masses%zwidth", zwidth, zwidth, "Z width")
    call cfg_get(cfg, "masses%mt", mt)
    call cfg_get(cfg, "masses%mb", mb)
    call cfg_get(cfg, "masses%mc", mc)
    call cfg_get(cfg, "masses%CKMdiag", CKMdiag)
    call cfg_get(cfg, "masses%CKMrotate", CKMrotate)

    ! Yukawas
    mt_yuk=mt
    if (abs(mb) > 1.e-8_dp) then
      mb_yuk=mb
    else
      mb_yuk=4.75_dp
    endif
    if (abs(mc) > 1.e-8_dp) then
      mc_yuk=mc
    else
      mc_yuk=1.5_dp
    endif
    mbsq=mb_yuk**2
    mcsq=mc_yuk**2

    ! [pdf]
    call cfg_get(cfg, "pdf%pdlabel", pdlabel)

    ! [lhapdf]
    if (cfg_var_size(cfg, "lhapdf%lhapdfset") > 1) then
       numpdfsets = cfg_var_size(cfg, "lhapdf%lhapdfset") 
       if (numpdfsets /= cfg_var_size(cfg, "lhapdf%lhapdfmember")) then
          error stop "Please specify a pdf member for each pdf set"
       endif
       write (*,*) "Running with multiple PDF sets!"
    else
       numpdfsets = 1
    endif

    allocate(PDFnames(numpdfsets))
    allocate(PDFmembers(numpdfsets))
    call cfg_get(cfg, "lhapdf%lhapdfset", PDFnames)
    call cfg_get(cfg, "lhapdf%lhapdfmember", PDFmembers)
    call cfg_get(cfg, "lhapdf%dopdferrors", doPDFerrors)

    ! [scales]
    ! note: setup ewcorr before
    call parse_scales()

    ! [basicjets]
    call cfg_get(cfg, "basicjets%inclusive", inclusive)
    call parse_jetalgo()
    call cfg_get(cfg, "basicjets%ptjetmin", ptjetmin)
    call cfg_get_add(cfg, "basicjets%ptjetmax", ptjetmax, sqrts, "")
    call cfg_get(cfg, "basicjets%etajetmax", etajetmax)
    call cfg_get_add(cfg, "basicjets%etajetmin", etajetmin, 0d0, "")
    call cfg_get(cfg, "basicjets%Rcutjet", Rcut)
    call cfg_get(cfg, "basicjets%userap", userap)

    if ((etajetmin < 0._dp) .or. (etajetmax < 0._dp)) then
       write(6,*) 'etajetmin and etajetmax are absolute values,'
       write(6,*) ' please reset to a positive value.'
       stop
    endif

    ! [masscuts]
    call cfg_get(cfg, "masscuts%m34min", m34min)
    wsqmin=m34min**2
    call cfg_get_add(cfg, "masscuts%m34max", m34max, sqrts, "")
    if (m34max > sqrts*0.9999_dp) m34max = sqrts*0.9999_dp
    wsqmax=m34max**2

    call cfg_get(cfg, "masscuts%m56min", m56min)
    bbsqmin=m56min**2
    call cfg_get_add(cfg, "masscuts%m56max", m56max, sqrts, "")
    if (m56max > sqrts*0.9999_dp) m56max = sqrts*0.9999_dp
    bbsqmax=m56max**2

    call cfg_get(cfg, "masscuts%m3456min", m3456min)
    call cfg_get_add(cfg, "masscuts%m3456max", m3456max, sqrts, "")
    if (m3456max > sqrts) m3456max = sqrts

    call cfg_get(cfg, "masscuts%pt34min", pt34min)
    call cfg_get_add(cfg, "masscuts%pt34max", pt34max, sqrts, "")
    if (pt34max > sqrts) pt34max = sqrts

    call cfg_get(cfg, "masscuts%closestZ", closestZ)
    call cfg_get(cfg, "masscuts%moppmin", moppmin)

    ! [cuts]
    call cfg_get(cfg, "cuts%makecuts", makecuts)

    call cfg_get(cfg, "cuts%ptleptmin", leptptmin)
    call cfg_get_add(cfg, "cuts%ptleptmax", leptptmax, sqrts, "")
    call cfg_get(cfg, "cuts%etaleptmax", leptrapmax)
    call cfg_get_add(cfg, "cuts%etaleptmin", leptrapmin, 0d0, "")
    call cfg_get(cfg, "cuts%etaleptveto", leptveto)
    leptveto1min = leptveto(1)
    leptveto1max = leptveto(2)
    call cfg_get(cfg, "cuts%ptminmiss", misspt)

    call cfg_get(cfg, "cuts%ptlept2min", leptpt2min)
    call cfg_get_add(cfg, "cuts%ptlept2max", leptpt2max, sqrts, "")
    call cfg_get(cfg, "cuts%etalept2max", leptrap2max)
    call cfg_get_add(cfg, "cuts%etalept2min", leptrap2min, 0d0, "")

    call cfg_get_add(cfg, "cuts%ptlept3min", leptpt3min, leptpt2min, "")
    call cfg_get_add(cfg, "cuts%ptlept3max", leptpt3max, leptpt2max, "")
    call cfg_get_add(cfg, "cuts%etalept3max", leptrap3max, leptrap2max, "")
    call cfg_get_add(cfg, "cuts%etalept3min", leptrap3min, leptrap2min, "")

    call cfg_get(cfg, "cuts%etalept2veto", lept2veto)
    leptveto2min = lept2veto(1)
    leptveto2max = lept2veto(2)

    call cfg_get(cfg, "cuts%m34transmin", mtrans34cut)
    call cfg_get(cfg, "cuts%Rjlmin", Rjlmin)
    call cfg_get(cfg, "cuts%Rllmin", Rllmin)
    call cfg_get(cfg, "cuts%delyjjmin", delyjjmin)
    call cfg_get(cfg, "cuts%jetsopphem", jetsopphem)
    call cfg_get(cfg, "cuts%lbjscheme", lbjscheme)
    call cfg_get(cfg, "cuts%ptbjetmin", ptbjetmin)
    call cfg_get_add(cfg, "cuts%ptbjetmax", ptbjetmax, sqrts, "")
    call cfg_get(cfg, "cuts%etabjetmax", etabjetmax)
    call cfg_get_add(cfg, "cuts%etabjetmin", etabjetmin, 0d0, "")

    call cfg_get(cfg, "cuts%elptmin", elptmin)
    call cfg_get(cfg, "cuts%muptmin", muptmin)
    call cfg_get(cfg, "cuts%Relelmin", Relelmin)
    call cfg_get(cfg, "cuts%Relmumin", Relmumin)
    call cfg_get(cfg, "cuts%Rmumumin", Rmumumin)

    call cfg_get(cfg, "cuts%Rlepiso", Rlepiso)
    call cfg_get(cfg, "cuts%fraclepiso", fraclepiso)

    call cfg_get(cfg, "cuts%elrapmax", elrapmax)
    call cfg_get(cfg, "cuts%murapmax", murapmax)
    call cfg_get(cfg, "cuts%etaelveto", leptveto)
    elvetomin = leptveto(1)
    elvetomax = leptveto(2)
    call cfg_get(cfg, "cuts%etamuveto", leptveto)
    muvetomin = leptveto(1)
    muvetomax = leptveto(2)
    call cfg_get(cfg, "cuts%ptminmissrel", missrelpt)
    call cfg_get(cfg, "cuts%mllmin", mllmin)
    call cfg_get(cfg, "cuts%ptllmin", ptllmin)

    call cfg_get_add(cfg, "cuts%mllmax", mllmax, sqrts, "")
    if (mllmax > sqrts) mllmax = sqrts

    call cfg_get(cfg, "cuts%m3lmin", m3lmin)

    call cfg_get(cfg, "cuts%y34min", y34min)
    call cfg_get(cfg, "cuts%y34max", y34max)

    ! [photon]
    call cfg_get(cfg, "photon%fragmentation", frag)
    call cfg_get(cfg, "photon%fragmentation_set", fragset)
    call cfg_get(cfg, "photon%fragmentation_scale", frag_scale)
    frag_scalestart = frag_scale

    call cfg_get(cfg, "photon%gammptmin", gammptmin)
    call cfg_get_add(cfg, "photon%gammptprod", gammptprod, 0d0, "")
    call cfg_get_add(cfg, "photon%gammptmax", gammptmax, sqrts, "")

    call cfg_get(cfg, "photon%gammrapmax", gammrapmax)
    call cfg_get_add(cfg, "photon%gammrapmin", gammrapmin, 0d0, "")
    call cfg_get(cfg, "photon%gammrapveto", leptveto)
    gammvetomin = leptveto(1)
    gammvetomax = leptveto(2)

    call cfg_get(cfg, "photon%gammpt2", gammpt2)
    call cfg_get(cfg, "photon%gammpt3", gammpt3)
    call cfg_get(cfg, "photon%Rgalmin", Rgalmin)
    call cfg_get(cfg, "photon%Rgagamin", Rgagamin)
    call cfg_get(cfg, "photon%Rgajetmin", Rgajetmin)
    call cfg_get(cfg, "photon%cone_ang", cone_ang)
    call cfg_get(cfg, "photon%epsilon_h", epsilon_h)
    call cfg_get(cfg, "photon%n_pow", n_pow)
    call cfg_get(cfg, "photon%fixed_coneenergy", fixed_coneenergy)
    call cfg_get(cfg, "photon%hybrid", photiso_hybrid)
    call cfg_get(cfg, "photon%R_inner", R_inner)

! Removed: JC implementation of hybrid isolation
!    call cfg_get(cfg, "photon%hybridiso", hybridiso)
!    call cfg_get(cfg, "photon%eps_dyn", eps_dyn)
!    call cfg_get(cfg, "photon%R_dyn", R_dyn)
!    call cfg_get(cfg, "photon%n_dyn", n_dyn)
!    call cfg_get(cfg, "photon%eps_iso", eps_iso)
!    call cfg_get(cfg, "photon%Etfix_iso", Etfix_iso)
!    call cfg_get(cfg, "photon%R_iso", R_iso)

    ! [lhcb]
    if (cfg_var_configadded(cfg, "lhcb%cut_mode")) then
       ! assume that block with lhcb config settings exists
       call cfg_get(cfg, "lhcb%cut_mode", cut_mode)
       call cfg_get(cfg, "lhcb%dir_mode", dir_mode)
       call cfg_get(cfg, "lhcb%nl_min", nl_min)
       call cfg_get(cfg, "lhcb%nj_min", nj_min)
       call cfg_get(cfg, "lhcb%nb_min", nb_min)

       if (cut_mode > 0) then
          call lhcb_config()
       endif
    endif

    ! [integration]
    call parse_seed()
    call cfg_get(cfg, "integration%scalarselect", scalarselect)

    ! [wz2jet]
    call cfg_get(cfg, "wz2jet%qflag", qflag)
    call cfg_get(cfg, "wz2jet%gflag", gflag)

    ! [hjetmass]
    call cfg_get(cfg, "hjetmass%mtex", mtex)

    ! [anom_higgs]
    call cfg_get(cfg, "anom_higgs%hwidth_ratio", hwidth_ratio)
    call cfg_get(cfg, "anom_higgs%cttH", cttH)
    call cfg_get(cfg, "anom_higgs%cWWH", cWWH)

    ! [bsm_higgs]
    call cfg_get(cfg, "bsm_higgs%bsm", bsm_higgs_scenario)
    call cfg_get(cfg, "bsm_higgs%t1", t1)
    call cfg_get(cfg, "bsm_higgs%t2", t2)
    call cfg_get(cfg, "bsm_higgs%t3", t3)
    call cfg_get(cfg, "bsm_higgs%t4", t4)
    call cfg_get(cfg, "bsm_higgs%t5", t5)
    call cfg_get(cfg, "bsm_higgs%t6", t6)
    call cfg_get(cfg, "bsm_higgs%w1", w1)
    call cfg_get(cfg, "bsm_higgs%w2", w2)
    call cfg_get(cfg, "bsm_higgs%w3", w3)
    call cfg_get(cfg, "bsm_higgs%w4", w4)
    call cfg_get(cfg, "bsm_higgs%w5", w5)
    call cfg_get(cfg, "bsm_higgs%c6",    c6_cfg)
    call cfg_get(cfg, "bsm_higgs%ct",    ct_cfg)
    call cfg_get(cfg, "bsm_higgs%cg",    cg_cfg)
    call cfg_get(cfg, "bsm_higgs%c6_sm", c6_sm)
    call cfg_get(cfg, "bsm_higgs%ct_sm", ct_sm)
    call cfg_get(cfg, "bsm_higgs%cg_sm", cg_sm)
    call cfg_get(cfg, "bsm_higgs%cx",    cx)
    call cfg_get(cfg, "bsm_higgs%cx_sm", cx_sm)
    call cfg_get(cfg, "bsm_higgs%mx",    mx)

    bsm_higgs_scenario = trim(adjustl(bsm_higgs_scenario))

    c6_init = c6_cfg(1)
    c6_nval = int( (c6_cfg(2) - c6_cfg(1)) / c6_cfg(3) + 1 )
    c6_step = c6_cfg(3)

    ct_init = ct_cfg(1)
    ct_nval = int( (ct_cfg(2) - ct_cfg(1)) / ct_cfg(3) + 1 )
    ct_step = ct_cfg(3)

    cg_init = cg_cfg(1)
    cg_nval = int( (cg_cfg(2) - cg_cfg(1)) / cg_cfg(3) + 1 )
    cg_step = cg_cfg(3)

    ! [output]
    call cfg_get(cfg, "output%events",   eventfile)
    call cfg_get(cfg, "output%metadata", metadatafile)

    ! [anom_wz]
    call cfg_get(cfg, "anom_wz%enable", anomtgc)
    call cfg_get(cfg, "anom_wz%delg1_z", delg1_z)
    call cfg_get(cfg, "anom_wz%delk_z", delk_z)
    call cfg_get(cfg, "anom_wz%delk_g", delk_g)
    call cfg_get(cfg, "anom_wz%lambda_z", lambda_z)
    call cfg_get(cfg, "anom_wz%lambda_g", lambda_g)
    call cfg_get(cfg, "anom_wz%h1Z", h1Z)
    call cfg_get(cfg, "anom_wz%h1gam", h1gam)
    call cfg_get(cfg, "anom_wz%h2Z", h2Z)
    call cfg_get(cfg, "anom_wz%h2gam", h2gam)
    call cfg_get(cfg, "anom_wz%h3Z", h3Z)
    call cfg_get(cfg, "anom_wz%h3gam", h3gam)
    call cfg_get(cfg, "anom_wz%h4Z", h4Z)
    call cfg_get(cfg, "anom_wz%h4gam", h4gam)
    call cfg_get(cfg, "anom_wz%tevscale", tevscale)

    ! [singletop]

    call cfg_get(cfg, "singletop%c_phiq", c_phiq_in)
    call anomcoup_tbW_set_c1(cmplx(c_phiq_in,0._dp,dp))

    call cfg_get(cfg, "singletop%c_phiphi", c_phiphi_in)
    call anomcoup_tbW_set_c2(cmplx(c_phiphi_in(1), c_phiphi_in(2), dp))

    call cfg_get(cfg, "singletop%c_tw", c_tw_in)
    call anomcoup_tbW_set_c3(cmplx(c_tw_in(1), c_tw_in(2), dp))

    call cfg_get(cfg, "singletop%c_bw", c_bw_in)
    call anomcoup_tbW_set_c4(cmplx(c_bw_in(1), c_bw_in(2), dp))

    call cfg_get(cfg, "singletop%c_tg", c_tg_in)
    call anomcoup_tbW_set_c6(cmplx(c_tg_in(1), c_tg_in(2), dp))

    call cfg_get(cfg, "singletop%c_bg", c_bg_in)
    call anomcoup_tbW_set_c7(cmplx(c_bg_in(1), c_bg_in(2), dp))

    call cfg_get(cfg, "singletop%lambda", lambda_in)
    call anomcoup_tbW_set_lambda(lambda_in)

    call cfg_get(cfg, "singletop%enable_lambda4",  enable_lambda4)
    if (enable_lambda4) then
       enable_lambda2 = .true.
    endif

    call cfg_get(cfg, "singletop%disable_sm", disable_sm)
    call cfg_get(cfg, "singletop%mode_anomcoup", mode_anomcoup)

    ! [extra]
    ! these are flags that are usually set with the runstring
    ! or are optional technical parameters that should usually not be set

    call cfg_get_add(cfg, "extra%debug", debug, .false., "debug")
    call cfg_get(cfg, "extra%verbose", verbose)

    call cfg_get_add(cfg, "extra%toponly", toponly, .false., "")
    call cfg_get(cfg, "extra%new_pspace", new_pspace)
    call cfg_get_add(cfg, "extra%spira", spira, .true., "")

    call cfg_get(cfg, "extra%griddebug", gridDebug)
    call cfg_get(cfg, "extra%dummypdf", dummypdf)

    call cfg_get_add(cfg, "singletop%nnlo_enable_light", singletop2_nnlo_enable_light, .true., "")
    call cfg_get_add(cfg, "singletop%nnlo_enable_heavy_prod", singletop2_nnlo_enable_heavy_prod, .true., "")
    call cfg_get_add(cfg, "singletop%nnlo_enable_heavy_decay", singletop2_nnlo_enable_heavy_decay, .true., "")

    call cfg_get_add(cfg, "singletop%nnlo_enable_interf_lxh", singletop2_nnlo_enable_lxh, .true., "")
    call cfg_get_add(cfg, "singletop%nnlo_enable_interf_lxd", singletop2_nnlo_enable_lxd, .true., "")
    call cfg_get_add(cfg, "singletop%nnlo_enable_interf_hxd", singletop2_nnlo_enable_hxd, .true., "")
    call cfg_get_add(cfg, "singletop%nnlo_fully_inclusive", singletop2_nnlo_fully_inclusive, .false., "")

    call cfg_get_add(cfg, "extra%dosingcheck", dosingcheck, .false., "perform singcheck")

    if (singletop2_nnlo_fully_inclusive) then
        if (rank == 0) then
            write (*,*) "WARNING: Fully inclusive singletop enabled, disabling decay."
            write (*,*) "         Make sure to set ptjetmin = 0, etajetmax = 99.0,"
            write (*,*) "         zerowidth = .true. and removebr = .true. in the input file."
        endif
    endif

    expansionorder = 0
    Wilsonorder = 0

    selectpdfs(:,:) = .true.

    if (cfg_var_configadded(cfg, "extra%pdfchannels")) then
       call cfg_get(cfg, "extra%pdfchannels", pdfchannels)
       pdfcase: select case (pdfchannels)
       case ("noglue")
          selectpdfs(:,0) = .false.
          write (*,*) "WARNING: Gluon PDF flux disabled"
       case ("qq")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 1:5) = .true.
          selectpdfs(2, 1:5) = .true.
          write (*,*) "WARNING: quark-quark channel only!"
       case ("qbqb")
          selectpdfs(:,:) = .false.
          selectpdfs(1, -5:-1) = .true.
          selectpdfs(2, -5:-1) = .true.
          write (*,*) "WARNING: antiquark-antiquark channel only!"
       case ("qg")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 1:5) = .true.
          selectpdfs(2, 0) = .true.
          write (*,*) "WARNING: quark-gluon channel only!"
       case ("qbg")
          selectpdfs(:,:) = .false.
          selectpdfs(1, -5:-1) = .true.
          selectpdfs(2, 0) = .true.
          write (*,*) "WARNING: antiquark-gluon channel only!"
       case ("gq")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 0) = .true.
          selectpdfs(2, 1:5) = .true.
          write (*,*) "WARNING: gluon-quark channel only!"
       case ("gqb")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 0) = .true.
          selectpdfs(2, -5:-1) = .true.
          write (*,*) "WARNING: gluon-antiquark channel only!"
       case ("qqb")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 1:5) = .true.
          selectpdfs(2, -5:-1) = .true.
          write (*,*) "WARNING: quark-antiquark channel only!"
       case ("qbq")
          selectpdfs(:,:) = .false.
          selectpdfs(1, -5:-1) = .true.
          selectpdfs(2, 1:5) = .true.
          write (*,*) "WARNING: antiquark-quark channel only!"
       case ("gg")
          selectpdfs(:,:) = .false.
          selectpdfs(1, 0) = .true.
          selectpdfs(2, 0) = .true.
          write (*,*) "WARNING: gluon-gluon channel only!"
       case default
          read (pdfchannels, *, err=666) pdfchan1, pdfchan2
          selectpdfs(:,:) = .false.
          selectpdfs(1,pdfchan1) = .true.
          selectpdfs(2,pdfchan2) = .true.
          write (*,*) "WARNING: Using only pdf channels", pdfchan1, pdfchan2
          exit pdfcase
666       continue
          write (*,*) "WARNING: Unknown pdf channel selection, doing nothing!"
       end select pdfcase
    endif

    nmin = 1
    nmax = 2

    call cfg_get_add(cfg, "extra%clustering", clustering, .true., "")
    call cfg_get_add(cfg, "extra%realwt", realwt, .false., "")
    call cfg_get_add(cfg, "extra%colourchoice", colourchoice, 0, "")

    call cfg_get_add(cfg, "extra%cutoff", cutoff, 1d-12, "")
    call cfg_get_add(cfg, "extra%cutoff_s", cutoff_s, 1d-5, "")
    call cfg_get_add(cfg, "extra%rtsmin", rtsmin, 1d-8, "")

    call cfg_get_add(cfg, "extra%reweight", enable_reweight_user, .false., &
                   "enable use of reweight_user routine")

    ! [dipoles]
    ! these are also optional
    call cfg_get_add(cfg, "dipoles%aii", aii, 1d0, "aii")
    call cfg_get_add(cfg, "dipoles%aif", aif, 1d0, "aif")
    call cfg_get_add(cfg, "dipoles%afi", afi, 1d0, "afi")
    call cfg_get_add(cfg, "dipoles%aff", aff, 1d0, "aff")
    call cfg_get_add(cfg, "dipoles%bfi", bfi, 1d0, "bfi")
    call cfg_get_add(cfg, "dipoles%bff", bff, 1d0, "bff")

    ! possibility for setting a common value
    if (cfg_var_configadded(cfg, "dipoles%alpha")) then
       call cfg_get_add(cfg, "dipoles%alpha", aii, 1d0, "common alpha parameter")
       aif = aii
       afi = aii
       aff = aii
       bfi = aii
       bff = aii
    endif

    ! if the user sets these to 1.0, we might end up with 1d0 + eps
    if (aii > 1d0) aii = 1d0
    if (aif > 1d0) aif = 1d0
    if (afi > 1d0) afi = 1d0
    if (aff > 1d0) aff = 1d0
    if (bfi > 1d0) bfi = 1d0
    if (bff > 1d0) bff = 1d0

    if (aii /= 1d0) write (*,*) "USING aii = ", aii
    if (afi /= 1d0) write (*,*) "USING afi = ", afi
    if (aif /= 1d0) write (*,*) "USING aif = ", aif
    if (aff /= 1d0) write (*,*) "USING aff = ", aff


!!! resummation settings
    call cfg_get(cfg, "resummation%res_range", res_range)
    call cfg_get(cfg, "resummation%resexp_range", resexp_range)
    call cfg_get(cfg, "resummation%fo_cutoff", res_fo_cutoff)

    call cfg_get_add(cfg, "resummation%ptveto", jetptveto, 1d5, "jet pt veto")
    call cfg_get_add(cfg, "resummation%useBanfid3veto", useBanfid3veto, .true., "Banfi version of d3veto")
    call cfg_get_add(cfg, "resummation%d3vetoR0", d3vetoR0, 1d0, "d3veto R0 (Banfi et al)")
    call cfg_get_add(cfg, "resummation%d3vetokappa", d3vetokappa, 1d0, "d3veto kappa (BNR)")
    call cfg_get_add(cfg, "resummation%gghsinglestep", gghsinglestep, .false., "ggh single step?")

    qtminRes = res_range(1)
    qtmaxRes = res_range(2)

    qtminResexp = resexp_range(1)
    qtmaxResexp = resexp_range(2)

    qtcutoff = res_fo_cutoff

    call cfg_get(cfg, "resummation%transitionswitch", transitionSwitch)
    call cfg_get_add(cfg, "benchmark%enable_fixed_y", enable_fixed_y, .false., "enable fixed y")
    call cfg_get_add(cfg, "benchmark%fixed_y", fixed_y, 0.0d0, "fixed y value")
    call cfg_get_add(cfg, "benchmark%enable_dsigma_dQ", enable_dsigma_dQ, .false., "calculate dsigma/dQ (mV)")
    call cfg_get_add(cfg, "benchmark%generations", generations, 5, "number nf to sum over")
    call cfg_get_add(cfg, "benchmark%fix_alphas_nf5", fix_alphas_nf5, .false., "fix alphas to 5-flavor 3-loop running")

    if (rank == 0) then
        if (resexp_linPC) then
            write (*,*) "******************************************"
            write (*,*) "*** Computing linear power corrections ***"
            write (*,*) "*** for N^kLO coefficient              ***"
            write (*,*) "******************************************"
        endif
    endif

!!! end resummation settings


    ! FINAL SANITY CHECKS AND WARNINGS (AND SETTINGS)
    if ((doscalevar) .and. (kewcorr /= knone)) then
       write(6,*) 'Cannot compute EW corrections and scale variation in a single run'
       stop
    endif

    if (kewcorr == kexact) then
       cutoff = max(cutoff, 1d-5)
    endif

    ! CHOOSER, set-up the variables for the process we wish to consider
    call chooser
    ! additional settings might depend on the chooser output below (n2,n3)

    ! E-M gauge invariance requires that delg1_g=0
    delg1_g=0._dp

    ! for resummation we tuned the non-new_pspace (except for WW,WZ,ZZ), so let's use these
    if (kpart == kresummed) then
       if ( (kcase /= kWWqqbr) .and. (kcase /= kWZbbar) .and. (kcase /= kZZlept)&
       .and.(kcase /= kWHbbar) .and. (kcase /= kZHbbar) ) then
         new_pspace = .false.
       endif
    endif

    !--- check that we have a valid value of 'part'
    if ( (kpart.ne.klord) .and. (kpart.ne.kreal) .and. &
                 (kpart.ne.kvirt) .and. (kpart.ne.ktota) ) then
       if    ( (kpart==ktodk) .and. &
                         ((kcase==kbq_tpq) .or. (kcase==kt_bbar) &
                     .or. (kcase==kW_twdk) .or. (kcase==ktt_bbl) &
                     .or. (kcase==ktt_bbh) .or. (kcase==k4ftwdk) &
                     .or. (kcase==kHWW2lq) .or. (kcase==kqq_ttw) &
                     .or. (kcase==kWWqqbr) .or. (kcase==ktt_bbu) &
                     .or. (kcase==kWHbbar) .or. (kcase==kZHbbar)) ) then
          !--- this is an allowed combination
       elseif ( (kpart==kfrag) .and. &
                       ((kcase==kWgamma) &
                   .or. (kcase==kgamgam) .or. (kcase==kgg2gam) &
                   .or. (kcase==kdirgam) .or. (kcase==kdm_gam) &
                   .or. (kcase==kgmgmjt) .or. (kcase==ktrigam) &
                   .or. (kcase==kfourga) &
                   .or. (kcase==kZ_2gam) .or. (kcase==kZgajet) &
                   .or. (kcase==kW_2gam)) ) then
          !--- this is an allowed combination
       elseif ((kpart==knnlo) .or. (kpart==ksnlo)) then 
          !--- this is an allowed combination
       elseif (kpart == kn3lo) then
           continue
       elseif (kpart==kresummed) then
          continue
       else 
          write(6,*) 'part=',part,' is not a valid option'
          write(6,*) 'for this process number.'
          stop
       endif
    endif


    !--- check that we are not trying to calculate radiation in decay at LO
    if    ( (kpart==klord) .and. &
                   ((kcase==kWWqqdk) .or. (kcase==kHWWdkW) &
               .or. (kcase==ktt_ldk) .or. (kcase==ktt_udk) &
               .or. (kcase==ktt_hdk) .or. (kcase==ktthWdk) &
               .or. (kcase==kttdkay) .or. (kcase==kWtdkay) &
               .or. (kcase==kdk_4ft) .or. (kcase==kttwldk) &
               .or. (kcase==kWHbbdk) .or. (kcase==kZHbbdk)) ) then
       write(6,*) 'This process number cannot be used for'
       write(6,*) 'a LO calculation.'
       stop
    endif

    !--- check that EW corrections are included for this process, if required
    if (kewcorr /= knone) then
       if ( (kcase == kZ_only) .or. (kcase == ktt_tot) &
                 .or. (kcase == ktt_mix) .or. (kcase == ktwojet) &
                 .or. (kcase == ktwo_ew) ) then
          continue
       else
          write(6,*) 'EW corrections not available for this process'
          stop
       endif
    endif

  end subroutine read_config

  subroutine parse_jetalgo()
    implicit none
    include 'clustering.f'

    call cfg_get(cfg, "basicjets%algorithm", algorithm)

    ! Assign choice of jet algorithm to integer variable
    if     (algorithm == 'ktal') then
       jetalgorithm=kt
    elseif (algorithm == 'ankt') then
       jetalgorithm=antikt
    elseif (algorithm == 'cone') then
       jetalgorithm=Rsepcone
    elseif (algorithm == 'hqrk') then
       jetalgorithm=hqrk
    elseif (algorithm == 'none') then
       jetalgorithm=noclustering
    else
       write(6,*) 'Invalid choice of jet algorithm: should be one of'
       write(6,*) 'ktal, ankt, cone, hqrk, none'
       stop
    endif

  end subroutine parse_jetalgo

  subroutine parse_ewcorr()
    implicit none
    include 'ewcorr.f'
    include 'cutoff.f'

    call cfg_get(cfg, "general%ewcorr", ewcorr)

    if     (ewcorr == 'none') then
       kewcorr=knone
    elseif (ewcorr == 'sudakov') then
       kewcorr=ksudakov
    elseif (ewcorr == 'exact') then
       kewcorr=kexact
    else
       write(6,*) 'Unexpected EW correction in input file: ',ewcorr
       stop
    endif

  end subroutine parse_ewcorr

  subroutine parse_seed()
    use cxx11random
    use omp_lib
    use iso_c_binding, only: c_loc
    implicit none
    include 'mpicommon.f'

    integer :: seed, origSeed
    common/seedBlock/seed
!$omp threadprivate(/seedBlock/)

    real(dp), save :: realSeed
!$omp threadprivate(realSeed)

    integer, target, dimension(:), allocatable :: seeds

    integer :: tid

    call cfg_get(cfg, "integration%seed", origSeed)

    allocate(seeds(omp_get_max_threads()))

!$omp parallel do
    do tid=0,omp_get_max_threads()-1
       !random seed if started with special seed value of 0
       if (origSeed == 0) then
          call random_seed()
          call random_number(realSeed)
          seed = -(nint(realSeed * huge(0)) + rank)
          seeds(tid+1) = -seed
       else
          seed = -(origSeed + omp_get_thread_num() + rank*100)
          seeds(tid+1) = -seed
       end if
    enddo
!$omp end parallel do

    call cxx11_init_random(c_loc(seeds))

  end subroutine parse_seed

  subroutine parse_taucut()
    use SCET
    use ptveto, only: usept
    implicit none
    include 'taucut.f'
    include 'kpart.f'
    include 'nproc.f'
    real(dp) xpow

    if (any(kpart==[ksnlo,knnlo,kn3lo]) &
        .or. ((kpart==kresummed) .and. (usept) .and. (kresorder==6))) then
       if (useqt .or. useqt_nnlo) then
         if (any(nproc == [31,32])) then
	     ! useqt is the new implementation up to N3LO from 2207.07056, but only set up for Z
	         useqt = .true.
             useqt_nnlo = .false.
         else
             ! useqt is the implementation from 2202.07738 only up to NNLO
             useqt = .false.
             useqt_nnlo = .true.
         endif

         call setupqt()
         call cfg_get_add(cfg, "nnlo%taucut", taucut, qtcut, "taucut")
         call cfg_get_add(cfg, "nnlo%qtcut", qtcut, qtcut, "qtcut")
         taucut = qtcut
       elseif (usept) then
         qtcut = 0.012_dp
         call cfg_get_add(cfg, "nnlo%taucut", taucut, qtcut, "taucut")
         call cfg_get_add(cfg, "nnlo%qtcut", qtcut, qtcut, "qtcut")
         taucut = qtcut
       else ! jettiness
         call setuptau()
         call cfg_get_add(cfg, "nnlo%qtcut", qtcut, taucut, "qtcut")
         call cfg_get_add(cfg, "nnlo%taucut", taucut, taucut, "taucut")
       endif

       usescet=.true.
       abovecut=.false.
       if (taucut < 0) then
          write(6,*) 'Must specify taucut/qtcut > 0 for SCET calculation'
          stop
       endif
    else
       usescet=.false.
       abovecut=.false.
    endif

    call cfg_get_add(cfg, "nnlo%tauboost", tauboost, .true., "")
    call cfg_get_add(cfg, "nnlo%incpowcorr", incpowcorr, .false., "")
    call cfg_get_add(cfg, "nnlo%onlypowcorr", onlypowcorr, .false., "")


    if (cfg_var_configadded(cfg, "nnlo%tcutarray") .and. usescet) then
       allocate(tcutarray(cfg_var_size(cfg, "nnlo%tcutarray")))
       call cfg_get(cfg, "nnlo%tcutarray", tcutarray)
       if (size(tcutarray) > 0) then
          doMultitaucut = .true.
       endif
    elseif (usescet .and. (kpart /= kresummed)) then
       allocate(tcutarray(5))
       if (useqt .or. useqt_nnlo .or. usept) then
         tcutarray = [taucut*2, taucut*4, taucut*8, taucut*20, taucut*40]
       else
         xpow=sqrt(2._dp)
         tcutarray = [taucut*2**xpow, taucut*4**xpow, taucut*8**xpow, taucut*16**xpow, taucut*32**xpow]
       endif
       doMultitaucut = .true.
    else
       allocate(tcutarray(0))
       doMultitaucut = .false.
    endif

    if (doMultitaucut) then
       if ( (nproc /= 1) .and. (nproc /= 6) .and. &
                         (nproc /= 31) .and. (nproc /= 300) .and. &
                         (nproc /= 111) .and. (nproc /= 112) .and. &
                         (nproc /= 119) .and. (nproc /= 285) .and. &
                         (nproc /= 290) .and. (nproc /= 295) .and. &
                         (nproc < 91) .and. (nproc > 110) ) then
          error stop "this process does not work with multitaucut yet"
       endif
    endif

!$omp parallel
    allocate(scetreweight(size(tcutarray)))
!$omp end parallel

    smallestTaucut = min(minval(tcutarray),taucut)

    ! this must be initialized to .true.
    ! only maketaucut modifies this, and for a non-scet run
    ! realint depends on this being true
!$omp parallel
    includeTaucutgrid(:) = .true.
!$omp end parallel

  end subroutine parse_taucut

  subroutine parse_scales()
    use Scalevar
    use singletop2_scale_m
    use qtResummation_params, only: scalevar_rapidity
    use singletop2_scale_m, only: singletop2_scale_init, use_DDIS
    use ptveto, only: timelikemusq
    implicit none
    include 'nf.f'
    include 'nproc.f'
    include 'initialscales.f'
    include 'facscale.f'
    include 'scale.f'
    include 'stopscales.f'
    include 'dynamicscale.f'
    include 'scalevar.f'
    include 'ewcorr.f'

    call cfg_get(cfg, "scales%renscale", scale)
    call cfg_get(cfg, "scales%facscale", facscale)
    call cfg_get(cfg, "scales%timelikemusq", timelikemusq)

    call cfg_get(cfg, "resummation%scalevar_rapidity", scalevar_rapidity)
    if (scalevar_rapidity) then
        doscalevar = .true.
        extrascalevar = extrascalevar + 2
    endif

    initscale=scale
    initfacscale=facscale

    ! special case for stop+b process
    initrenscale_L=0._dp
    initfacscale_L=0._dp
    initrenscale_H=0._dp
    initfacscale_H=0._dp

    if (((nproc >= 231) .and. (nproc <= 240)) .and. &
                 (scale == 0._dp) .and. (facscale == 0._dp)) then

       call cfg_get(cfg, "scales%renscale_L", initrenscale_L)
       call cfg_get(cfg, "scales%facscale_L", initfacscale_L)
       call cfg_get(cfg, "scales%renscale_H", initrenscale_H)
       call cfg_get(cfg, "scales%facscale_H", initfacscale_H)

       renscale_L=initrenscale_L
       facscale_L=initfacscale_L
       renscale_H=initrenscale_H
       facscale_H=initfacscale_H
       scale=initrenscale_H
       facscale=initfacscale_H        
    endif

    !         dynamic scale

    call cfg_get(cfg, "scales%dynamicscale", dynstring)

    call cfg_get(cfg, "scales%doscalevar", doscalevar)
    call cfg_get(cfg, "scales%maxscalevar", maxscalevar)
    call cfg_get(cfg, "scales%vetoscalevar", vetoscalevar)

    if (doscalevar .and. (.not. any(maxscalevar == [2,6,8]))) then
       error stop "maxscalevar must be 2,6 or 8 when doing scale variation"
    endif

    ! sets use_DDIS flag if dynscale is 'DDIS'
    call singletop2_scale_init()

    if (use_DDIS .and. doscalevar) then
        ! we additionally sample variation around mt as scale 
        maxscalevar = maxscalevar*2 + 1
    endif

    !---  create logical:: variable dynamicscale for use in other routines
    if (  (dynstring == 'no') .or. (dynstring == '.false.') &
             .or. (dynstring == 'none') ) then 
       dynamicscale=.false. 
    else
       dynamicscale=.true. 
    endif




  end subroutine parse_scales

  subroutine parse_part(part)
    implicit none
    include 'kpart.f'
    character(len=*), intent(in) :: part

    coeffonly=.false.
    kpart=0
    if     ((part == 'lo') .or. (part == 'lord')) then
       kpart=klord
    elseif (part == 'virt') then
       coeffonly = .true.
       kpart=kvirt
    elseif (part == 'real') then
       kpart=kreal
    elseif ((part == 'nlo') .or. (part == 'tota') &
               .or. (part == 'nlocoeff') .or. (part == 'totacoeff')) then
       kpart=ktota
    elseif (part == 'frag') then
       kpart=kfrag
    elseif ((part == 'todk') .or. (part == 'nlodk')) then
       kpart=ktodk
    elseif ((part == 'snlo') .or. (part == 'scetnlo') &
               .or. (part == 'snlocoeff') .or. (part == 'scetnlocoeff') &
               .or. part == 'snloR' .or. part == 'snloV') then
       kpart=ksnlo
       if (part == 'snloR') then
          ksnlopart = ksnloR
          coeffonly = .true.
       elseif (part == 'snloV') then
          ksnlopart = ksnloV
          coeffonly = .true.
       else
          ksnlopart = 0
       endif
    elseif ((part == 'nnlo') .or. (part == 'nnlocoeff') &
               .or. (part == 'nnloVV') .or. (part == 'nnloVVcoeff') &
               .or. (part == 'nnloRV') .or. (part == 'nnloRVcoeff') &
               .or. (part == 'nnloRR') .or. (part == 'nnloRRcoeff') ) then
       kpart=knnlo
       if     ((part == 'nnloVV') .or. (part == 'nnloVVcoeff')) then
          knnlopart=knnloVV
       elseif ((part == 'nnloRV') .or. (part == 'nnloRVcoeff')) then
          knnlopart=knnloRV
       elseif ((part == 'nnloRR') .or. (part == 'nnloRRcoeff')) then
          knnlopart=knnloRR
       else
          knnlopart=0
       endif
    elseif (part == 'n3lo' .or. part == 'n3locoeff') then
        kpart=kn3lo
    elseif (part == 'resLO') then
        kpart = kresummed
        krespart = 0
        kresorder = 2
    elseif (part == 'resonlyLO') then
        kpart = kresummed
        krespart = kresonly
        kresorder = 2
    elseif (part == 'resonlyLOp') then
        kpart = kresummed
        krespart = kresonly
        kresorder = 3
    elseif (part == 'resexpNLO') then
       kpart = kresummed
       krespart = kresexp
       kresorder = 4
    elseif (part == 'resonlyNLO') then
       kpart = kresummed
       krespart = kresonly
       kresorder = 4
    elseif (part == 'resaboveNLO') then
       kpart = kresummed
       krespart = kresabove
       kresorder = 4
    elseif (part == 'resmatchcorrNLO') then
       kpart = kresummed
       krespart = kresmatchcorr
       kresorder = 4
    elseif (part == 'resonlyNLOp') then
       kpart = kresummed
       krespart = kresonly
       kresorder = 5
    elseif (part == 'resexpNNLO') then
       kpart = kresummed
       krespart = kresexp
       kresorder = 6
    elseif (part == 'resonlyNNLO') then
       kpart = kresummed
       krespart = kresonly
       kresorder = 6
    elseif (part == 'resaboveNNLO') then
       kpart = kresummed
       krespart = kresabove
       kresorder = 6
    elseif (part == 'resmatchcorrNNLO') then
       kpart = kresummed
       krespart = kresmatchcorr
       kresorder = 6
    elseif (part == 'resLOp') then
       kpart=kresummed
       krespart = 0
       kresorder = 3
    elseif (part == 'resNLO') then
       kpart=kresummed
       krespart = 0
       kresorder = 4
    elseif (part == 'resNLOp') then
       kpart=kresummed
       krespart = 0
       kresorder = 5
    elseif (part == 'resNNLO') then
       kpart=kresummed
       krespart = 0
       kresorder = 6
    elseif (part == 'resNNLOp') then
       kpart=kresummed
       krespart = 0
       kresorder = 7
    elseif (part == 'resonlyNNLOp') then
       kpart = kresummed
       krespart = kresonly
       kresorder = 7
    elseif (part == 'resNNLOp') then
       kpart = kresummed
       krespart = 0
       kresorder = 7
    elseif (part == 'resexpN3LO') then
        kpart = kresummed
        krespart = kresexp
        kresorder = 7
    elseif (part == 'resonlyN3LO') then
        kpart = kresummed
        krespart = kresonly
        kresorder = 8
    endif
    if (index(part,'coeff') > 0) then
       coeffonly=.true.
    endif

    origKpart = kpart

    if (kpart == 0) then
       write(6,*) 'Invalid value of part = ',part
       stop
    endif


  end subroutine parse_part

end module parseinput
