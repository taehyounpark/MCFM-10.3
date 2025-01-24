
!     The pseudo-rapidity and pT requirements. The pseudo-rapidity
!     requirements are arrays of size 2, corresponding to the forward
!     backward requirement: 1 is the first lepton (hardest), 2 the
!     second lepton, j the jet, b the heavy jet, g the photon, and 1j
!     the combined first lepton and jet.
      real(dp) ::     eta1_min(2), eta1_max(2), pt1_min, pt1_max
      real(dp) ::     eta2_min(2), eta2_max(2), pt2_min, pt2_max
      real(dp) ::     etaj_min(2), etaj_max(2), ptj_min, ptj_max
      real(dp) ::     etab_min(2), etab_max(2), ptb_min, ptb_max
      real(dp) ::     etag_min(2), etag_max(2), ptg_min, ptg_max
      real(dp) ::     pt1j_min, pt1j_max
!     The variables stored for histogramming. All variables are arrays
!     of size 2 corresponding to forward and backward.
      real(dp) ::     eta1(2), eta2(2), etaj(2), etab(2), etag(2)
      real(dp) ::     pt1(2), pt2(2), ptj(2), ptb(2), ptg(2), pt1j(2), m1j(2)
!     The minimum number of letpons, jets, and heavy jets.
      integer ::     nl_min, nj_min, nb_min
!     The cut mode and direction mode. The cut modes are as follows.
!     0: MCFM cuts.
!     1: LHCb-ANA-2014-076, cut on hardest lepton and jet in acceptance.
!        If the event contains heavy flavor and hard_bjet is true the hardest
!        jet must be heavy.
!     2: Same as (1), but ignore the lepton from the W decay in nprocs 71,
!        76, 181, and 186.
!     3: Same as (1), but ignore the leptons from the Z decay in nprocs 71
!        and 76.
!     The direction mode specifies how to handle asymmetric
!     pseudo-rapidity.
!     0: Use the event and flipped (in pseudo-rapidity) event.
!     1: Use only the event.
!     2: Use only the flipped event.
      integer :: cut_mode, dir_mode
!     Logical if the event passes the forward and backward requirements,
!     and logical if the hardest jet must be heavy.
      logical ::  lhcb_pass(2), hard_bjet
      common/lhcb1/ eta1_min, eta1_max, pt1_min, pt1_max, eta2_min, eta2_max, pt2_min, pt2_max
!$omp threadprivate(/lhcb1/)
      common/lhcb2/ etaj_min, etaj_max, ptj_min, ptj_max, etab_min, etab_max, ptb_min, ptb_max
!$omp threadprivate(/lhcb2/)
      common/lhcb3/ etag_min, etag_max, ptg_min, ptg_max, pt1j_min, pt1j_max
!$omp threadprivate(/lhcb3/)
      common/lhcb4/ eta1, eta2, etaj, etab, etag, pt1, pt2, ptj, ptb, ptg, pt1j, m1j
!$omp threadprivate(/lhcb4/)
      common/lhcb5/ nl_min, nj_min, nb_min, cut_mode, dir_mode, lhcb_pass, hard_bjet
!$omp threadprivate(/lhcb5/)
