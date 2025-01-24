!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine writereference
      implicit none
      include 'types.f'
      include 'kpart.f'
      include 'nproc.f'
      include 'ewcorr.f'
c      logical:: writerefs
c      common/writerefs/writerefs

      !if (writerefs .eqv. .false.) return

        write(6,*)
        write(6,58) '****************  MCFM references   ****************'
      if (origkpart == kresummed) then
        write(6,58) '*                                                  *'
        write(6,58) '*  Fiducial qT resummation of color-singlet        *'
        write(6,58) '*  processes at N3LL+NNLO (T. Becher, T.Neumann)   *'
        write(6,58) '*    arXiv:2009.11437  (CuTe-MCFM)                 *'
      endif
#ifdef WITH_LIBRARY
        write (6,*) '*                                                  *'
        write (6,*) '*  Accelerating LHC phenomenology with analytic    *'
        write (6,*) '*    one-loop amplitudes: A C++ interface to MCFM  *'
        write (6,*) '*   J.M. Campbell, S. Hoeche, C. Preuss            *'
        write (6,*) '*    arXiv:2107.04472                              *'
#endif
        write(6,58) '*                                                  *'
        write(6,58) '*  Precision phenomenology with MCFM-9.0           *'
        write(6,58) '*   J.M. Campbell, T. Neumann                      *'
        write(6,58) '*    arXiv:1909.09117                              *'
        write(6,58) '*                                                  *'
        write(6,58) '*  An update on vector boson pair production at    *'
        write(6,58) '*    hadron colliders                              *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, PRD60 (1999) 113006 *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Vector boson pair production at the LHC         *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    JHEP 1107 (2011) 018                          *'
        write(6,58) '*                                                  *'
        write(6,58) '****************************************************'

        write(6,58) '**************  Process references   ***************'
        write(6,58) '*                                                  *'

      if (origkpart == kresummed .and. nproc == 285) then
        write(6,58) "* The diphoton spectrum at N3LL'+NNLO              *"
        write(6,58) '*   T. Neumann                                     *'
        write(6,58) '* EPJC 81 (2021) 10, 905, arXiv:2107.12478         *'
        write(6,58) '*                                                  *'
      endif

      if     (nproc == 285) then
        write(6,58) '*  Predictions for diphoton production at the      *'
        write(6,58) '*    LHC through NNLO in QCD                       *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, Ye Li, C. Williams  *'
        write(6,58) '*    JHEP 1607 (2016) 148                          *'
        write(6,58) '*                                                  *'
      endif

      if ((nproc >= 91) .and. (nproc <= 110)) then
        write(6,58) '*  Associated production of a Higgs boson at NNLO  *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    JHEP 1606 (2016) 179                          *'
        write(6,58) '*                                                  *'
      endif

      if (nproc == 41 .and. kpart == knnlo) then
        write (6,58) '* Fiducial Drell-Yan production at the LHC       *'
        write (6,58) '*   improved by qT-resummation at N4LL+N3LO      *'
        write (6,58) '* T. Neumann and J.M. Campbell                   *'
        write (6,58) '*   arXiv: 2207.07056                            *'
        write (6,58) '*                                                *'
        write (6,58) '* Z-boson production in association with a jet   *'
        write (6,58) '*   at NNLO in perturbative QCD                  *'
        write (6,58) '* R. Boughezal, J.M. Campbell, R.K. Ellis,       *'
        write (6,58) '*  C. Focke, W. Giele, X. Liu, F. Petriello      *'
        write (6,58) '* PRL 116, 152001 (2016), arXiv:1512.01291       *'
        write (6,58) '*                                                *'
      endif

      if (any(nproc == [31,32]) .and. kpart == kresummed .or. kpart == kn3lo) then
        write (6,58) '* Fiducial Drell-Yan production at the LHC       *'
        write (6,58) '*   improved by qT-resummation at N4LL+N3LO      *'
        write (6,58) '* T. Neumann and J.M. Campbell                   *'
        write (6,58) '*   arXiv: 2207.07056                            *'
        write (6,58) '*                                                *'
      endif

      if (kpart == knnlo
     &  .and. any(nproc == [1,6,31,32,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,119])) then
        write(6,58) '*  Color singlet production at NNLO in MCFM        *'
        write(6,58) '*   R. Boughezal, J. Campbell, R.K. Ellis,         *'
        write(6,58) '*    C. Focke, W. Giele, X. Liu, F. Petriello,     *'
        write(6,58) '*    C. Williams,  arXiv: 1605.08011               *'
        write(6,58) '*                                                  *'
      endif

      if ((kpart == knnlo) .and. (nproc >= 61) .and. (nproc <= 90)) then
        write(6,58) '*  Non-local slicing approaches for                *'
        write(6,58) '*    NNLO QCD in MCFM                              *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, S. Seth             *'
        write(6,58) '*    arXiv: 2202.07738                             *'
        write(6,58) '*                                                  *'
      endif

        if     (nproc == 287) then
        write(6,58) '*  Triphoton production at hadron colliders        *'
        write(6,58) '*   J.M. Campbell, C. Williams, arXiv:1403.2641    *'
        write(6,58) '*                                                  *'
        elseif (nproc == 289) then
        write(6,58) '*  Four-photon production at the LHC: an           *'
        write(6,58) '*    of 2 -> 4 analytic unitarity                  *'
        write(6,58) '*   T. Dennen, C. Williams, arXiv:1411.3237        *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 20) .or. (nproc == 25)) then
        write(6,58) '*  QCD corrections to the hadronic production of a *'
        write(6,58) '*    heavy quark pair including decay correlations *'
        write(6,58) '*   S. Badger, J.M. Campbell, R.K. Ellis           *'
        write(6,58) '*    arXiv:1011.6647                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 21) .or. (nproc == 26)) then
        write(6,58) '*  Strong radiative correction to Wbb production   *'
        write(6,58) '*    in proton-antiproton collisions               *'
        write(6,58) '*   R.K. Ellis, S. Veseli, hep-ph/9810489          *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 22) .or. (nproc == 27)
     &     .or. (nproc == 44) .or. (nproc == 46)) then
        write(6,58) '*  Next-to-leading order corrections to W+2jet     *'
        write(6,58) '*    and Z+2jet production at hadron colliders     *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, hep-ph/0202176      *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Next-to-leading order QCD predictions for       *'
        write(6,58) '*    W+2jet and Z+2jet production at the CERN LHC  *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, D.L. Rainwater,     *'
        write(6,58) '*    hep-ph/0308195                                *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 50) .and. (nproc <= 56)) then
        write(6,58) '*  Radiative corrections to Zbb production         *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, hep-ph/0006304      *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 123) .and. (nproc <= 127)) then
        write(6,58) '*  Gluon-gluon contributions to W+ W- production   *'
        write(6,58) '*    and Higgs interference effects                *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1107.5569                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 128) .and. (nproc <= 132)) then
        write(6,58) '*  Bounding the Higgs width at the LHC using       *'
        write(6,58) '*    full analytic results for gg -> e-e+mu-mu+    *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1311.3589                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 141) .and. (nproc <= 151)) then
        write(6,58) '*  Top-quark processes at NLO in production        *'
        write(6,58) '*     and decay                                    *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1204.1513     *'
        write(6,58) '*                                                  *'
        elseif (nproc == 1610 .or. nproc == 1650) then
        write(6,58) '*   Single top production in the t-channel         *'
        write(6,58) '*     at NNLO                                      *'
        write(6,58) '*   J.M. Campbell, T. Neumann, Z. Sullivan         *'
        write(6,58) '*    arXiv:2012.01574                              *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 164) .or. (nproc == 169)) then
        write(6,58) '*  Off-shell NLO predictions for                   *'
        write(6,58) '*    single-top-quark production in the            *'
        write(6,58) '*    Standard Model Effective Field Theory         *'
        write(6,58) '*   T. Neumann, Z. Sullivan,  arXiv:1903.11023     *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 161) .and. (nproc <= 177)) then
        write(6,58) '*   Single top production and decay at             *'
        write(6,58) '*     next-to-leading order                        *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, F. Tramontano,      *'
        write(6,58) '*    hep-ph/0408158                                *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 180) .and. (nproc <= 187)) then
        write(6,58) '*  Next-to-leading order corrections to Wt         *'
        write(6,58) '*    production and decay                          *'
        write(6,58) '*   J.M. Campbell, F. Tramontano, hep-ph/0506289   *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 211) .and. (nproc <= 217)) then
        write(6,58) '*  Higgs boson production in weak boson fusion     *'
        write(6,58) '*    at next-to-leading order                      *'
        write(6,58) '*   E.L. Berger, J.M. Campbell, hep-ph/0403194     *'
        write(6,58) '*                                                  *'
        elseif (((nproc >= 220) .and. (nproc <= 229)) .or.
     &          ((nproc >= 2201) .and. (nproc <= 2291))) then
        write(6,58) '*  Higgs constraints from vector boson fusion      *'
        write(6,58) '*    and scattering                                *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1502.02990    *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 270) .and. (nproc <= 279)) then
        write(6,58) '*  Next-to-leading order Higgs + 2 jet production  *'
        write(6,58) '*    via gluon fusion                              *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, G. Zanderighi,      *'
        write(6,58) '*    hep-ph/0608194                                *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Hadronic production of a Higgs boson and        *'
        write(6,58) '*    two jets at next-to-leading order             *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1001.4495                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 290) .and. (nproc <= 299)) then
        write(6,58) '*  The pp -> W (-> l nu) + gamma process at        *'
        write(6,58) '*    next-to-next-to-leading order                 *'
        write(6,58) '*   J.M. Campbell, G. De Laurentis, R.K. Ellis     *'
        write(6,58) '*    and S. Seth, arXiv:2105.00954                 *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 300) .and. (nproc <= 309)) then
        write(6,58) '*  Next-to-leading order predictions for           *'
        write(6,58) '*    Zgam+jet and Zgamgam final states at the LHC  *'
        write(6,58) '*   J.M. Campbell, H. Hartanto, C. Williams        *'
        write(6,58) '*    arXiv:1208.0566                               *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Z+photon production at NNLO                     *'
        write(6,58) '*    including Anomalous Couplings                 *'
        write(6,58) '*   J.M. Campbell, T.Neumann, C. Williams          *'
        write(6,58) '*    arXiv:1802.02981                              *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 500) .and. (nproc <= 516)) then
        write(6,58) '*  ttW production and decay at NLO                 *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1204.5678     *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 540) .and. (nproc <= 569)) then
        write(6,58) '*  Single top production in association with a     *'
        write(6,58) '*    Z boson at the LHC                            *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, R. Rontsch,         *'
        write(6,58) '*    arXiv:1302.3856                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 800) .and. (nproc <= 848)) then
        write(6,58) '*  Next-to-leading order predictions for           *'
        write(6,58) '*    dark matter production at hadron colliders    *'
        write(6,58) '*   P.J. Fox, C. Williams, arXiv:1211.6390         *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 200)) then
        write(6,58) '*  Next-to-leading order predictions for           *'
        write(6,58) '*    Higgs+jet prod. with top-quark mass effects   *'
        write(6,58) '*   T. Neumann, C. Williams, arXiv:1609:00367      *'
        write(6,58) '*   T. Neumann, arXiv:1802.02981                   *'
        write(6,58) '*                                                  *'
        write(6,58) '*   One-loop amplitudes for Higgs + 4 partons      *'
        write(6,58) '*   L. Budge, J.M. Campbell, G. De Laurentis,      *'
        write(6,58) '*   R.K. Ellis, S. Seth, arXiv:2002.04018          *'
        endif

        if (kewcorr /= knone) then
        write(6,58) '*  A study of weak corrections to Drell-Yan,       *'
        write(6,58) '*    top-quark pair and di-jet production at       *'
        write(6,58) '*    high energies with MCFM                       *'
        write(6,58) '*   J.M. Campbell, D. Wackeroth, J. Zhou,          *'
        write(6,58) '*    arXiv:1608.03356                              *'
        write(6,58) '*                                                  *'
        endif

        write(6,58) '****************************************************'

   58 format(a53)

      return
      end

