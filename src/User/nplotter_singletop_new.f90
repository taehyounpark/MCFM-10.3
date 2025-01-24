!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
module nplotter_singletop
      use MCFMPlotting
      implicit none

      public :: setup, book
      public :: tag_bjets

      ! save b-quantum number for jets
      integer, save, public, protected :: jet_bness(4)
!$omp threadprivate(jet_bness)
      ! saves whether jet contains main b-quark from top-decay (particle 5)
      logical, save, public, protected :: jet_has_mainb(4)
!$omp threadprivate(jet_has_mainb)

      private

      integer, save, allocatable :: histos(:)

      contains

      subroutine setup()
          use types
          implicit none

          allocate(histos(33))

          histos(1) = plot_setup_uniform(0._dp,300._dp,10._dp,'pttop')
          histos(2) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptlj1')
          histos(3) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptlj2')
          histos(4) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptlj3')

          histos(5) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptj1')
          histos(6) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptj2')

          histos(7) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptb1')
          histos(8) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptb2')
          histos(9) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptb3')

          histos(10) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptl')
          histos(11) = plot_setup_uniform(0._dp,300._dp,10._dp,'ptbl')

          histos(12) = plot_setup_uniform(0._dp,5._dp,1._dp,'jetbins')
          histos(13) = plot_setup_uniform(0._dp,5._dp,1._dp,'bjetbins')

          histos(14) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'eta_l')
          histos(15) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'ytop')

          !histos(16) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etalj1')
          histos(16) = plot_setup_uniform(-5.0_dp,5.0_dp,0.5_dp,'etalj1')

          histos(17) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etalj2')
          histos(18) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etalj3')

          !histos(19) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etaj1')
          histos(19) = plot_setup_uniform(-5.0_dp,5.0_dp,0.5_dp,'etaj1')
          histos(20) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etaj2')

          histos(21) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etab1')
          histos(22) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etab2')
          histos(23) = plot_setup_uniform(-4.6_dp,4.6_dp,0.4_dp,'etab3')

          histos(24) = plot_setup_uniform(0._dp,200.0_dp,10._dp,'mbl')

          histos(25) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'cosxprod')
          histos(26) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'cosyprod')
          histos(27) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'coszprod') ! costheta

          histos(28) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'coslstar')
          histos(29) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'coslN')
          histos(30) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'coslT')

          histos(31) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'cosphiNstar')
          histos(32) = plot_setup_uniform(-1._dp,1._dp,0.1_dp,'cosphiTstar')

          histos(33) = plot_setup_uniform(0._dp,5._dp,1._dp,'ljetbins')



      end subroutine

      subroutine tag_bjets(p)
          use types
          use singletop2_nnlo_vars, only: max_bcontrib, b_contrib, &
              currentContrib
          use jettagging
          use SCET, only: currentNd
          use singletop2_nnlo_vars, only: singletop2_nnlo_fully_inclusive
          implicit none
          include 'mxpart.f'
          include 'nwz.f'
          include 'ptilde.f'
          include 'ipsgen.f'
          include 'kpart.f'

          real(dp), intent(in) :: p(mxpart,4)


          ! b-ness of individual partons, used to update b-ness of jets
          ! first: parton number (5,6,7,8), then currentContrib piece, then b_contrib
          integer :: parton_bness(8, 6, max_bcontrib)
          integer :: j,k
          integer :: njets

          logical, parameter :: infraredsafe = .true.

          real(dp) :: massvec 
          integer :: mainb
          integer :: nextjet
          integer, allocatable :: jetindices(:)

          jet_has_mainb(:) = .false.
          parton_bness(:,:,:) = 0
          jet_bness(:) = 0


          if (singletop2_nnlo_fully_inclusive) then
              allocate(jetindices(3), source=[6,7,8])
              njets = count(p(jetindices,4) /= 0._dp)
              nextjet = 6
          else
              allocate(jetindices(4), source=[5,6,7,8])
              njets = count(p(jetindices,4) /= 0._dp)
              nextjet = 5
          endif

          if (nwz == -1) then
              stop "todo: anti-top b-tagging"
              ! likely just flipping sign of parton_bness
          endif

          ! these are the g -> b b~ splitting pieces
          ! that should be ignored in the naive scheme
#define IPSDECAY 3
#define IPSHEAVY 2
#define IPSLIGHT 1
#define IPSLXH 4
#define IPSHXD 6

#define ONLY_MAIN_B 1
#define WITH_B 2
#define WITH_BBAR 3
#define WITH_B_BBAR 4
#define WITH_BBAR_BBAR 5

          ! main-b, always present
          ! Comment then for fully inclusive stable-top
          parton_bness(5,:,:) = 1

          ! is eliminated below
          parton_bness(7, IPSDECAY, WITH_B_BBAR) = 1
          parton_bness(8, IPSDECAY, WITH_B_BBAR) = -1

          parton_bness(7, IPSLXH, WITH_BBAR) = -1

          if (currentContrib == IPSHXD) then
              if (kpart == kreal .and. ipsgen == 9) then ! RV
                  parton_bness(7, IPSHXD, WITH_BBAR) = -1
              elseif (kpart == kreal .and. ipsgen == 8) then ! RR
                  parton_bness(8, IPSHXD, WITH_BBAR) = -1
              endif
          endif

          parton_bness(7, IPSHEAVY, WITH_BBAR) = -1

          ! is eliminated below
          parton_bness(7, IPSHEAVY, WITH_B_BBAR) = -1
          parton_bness(8, IPSHEAVY, WITH_B_BBAR) = 1

          parton_bness(7, IPSHEAVY, WITH_BBAR_BBAR) = -1
          parton_bness(8, IPSHEAVY, WITH_BBAR_BBAR) = -1


          parton_bness(7, IPSLIGHT, WITH_B) = 1
          parton_bness(7, IPSLIGHT, WITH_BBAR) = -1

          ! is eliminated below
          parton_bness(7, IPSLIGHT, WITH_B_BBAR) = 1
          parton_bness(8, IPSLIGHT, WITH_B_BBAR) = -1

          mainb = 5
          ! remove flavor assignment of g -> b b~ splitting
          if (b_contrib == WITH_B_BBAR) then
              if (currentContrib == IPSDECAY) then
                  ! in case of decay we have to include the main b into the
                  ! decision, otherwise 7 and 8 will always be the b b~ pair
                  if (massvec(ptilde(currentNd,5,:)+ptilde(currentNd,8,:)) > &
                      massvec(ptilde(currentNd,7,:)+ptilde(currentNd,8,:))) then

                      parton_bness(7, currentContrib, b_contrib) = 0
                      parton_bness(8, currentContrib, b_contrib) = 0
                      mainb = 5
                  else
                      parton_bness(5, currentContrib, b_contrib) = 0
                      parton_bness(8, currentContrib, b_contrib) = 0
                      mainb = 7
                  endif
              else
                  parton_bness(7, currentContrib, b_contrib) = 0
                  parton_bness(8, currentContrib, b_contrib) = 0
              endif
          endif

          do j=1,njets
            do k=nextjet,8
                if (iand(bitflag(k), jetcontent(jetindices(j),currentNd)) > 0) then
                    if (k == mainb) jet_has_mainb(j) = .true.
                    jet_bness(j) = jet_bness(j) + parton_bness(k,currentContrib,b_contrib)
                endif
            enddo
          enddo

      end subroutine

      subroutine book(p,wt,ids,vals,wts)
          use types
          use singletop2_nnlo_vars, only: max_bcontrib, b_contrib, corr_on_beam
          use jettagging
          use SCET, only: currentNd
          use qsort_m
          use singletop2_nnlo_vars
          use mathfun, only: crossprod, euclideanNorm
          use ieee_arithmetic
          implicit none
          include 'mxpart.f'
          include 'ipsgen.f'
          include 'nwz.f'
          include 'ptilde.f'
          include 'masses.f'

          real(dp), intent(in) :: p(mxpart,4), wt
          integer, allocatable, intent(out) :: ids(:)
          real(dp), allocatable, intent(out) :: vals(:)
          real(dp), allocatable, intent(out) :: wts(:)

          integer :: numb, numnonb
          integer :: njets

          real(dp) :: ptpure, etarappure, yrappure
          real(dp) :: ptop(4)

          real(dp) :: pttop, ptlj1, ptlj2, ptlj3, ptb1, ptb2, ptb3, ptl, ptbl
          real(dp) :: ptj1, ptj2, etaj1, etaj2
          real(dp) :: jetbin, bjetbin, ljetbin, puremass
          real(dp) :: etal,ytop,etalj1,etalj2,etalj3,etab1,etab2,etab3
          real(dp) :: mbl
          real(dp) :: cosxprod, cosyprod, coszprod
          real(dp) :: coslstar, coslN, coslT, cosphiNstar, cosphiTstar


          real(dp) :: pwrest_w(4), pspecrest(4)
          real(dp) :: plepwrest(4), beam(4), beamrest(4)

          real(dp) :: st(3), x(3), y(3), z(3)
          real(dp) :: xprod(3), yprod(3), zprod(3), plepWrestproj(3)
          integer :: j,k

          integer :: bjetindices(4)
          integer :: nonbjetindices(4)

          real(dp) :: plep(4), plight(4), pleprest(4), plightrest(4), &
                ptoprest(4), pw(4), pwrest(4)

          integer :: lastnonjet

          call tag_bjets(p)

          if (singletop2_nnlo_fully_inclusive) then
              njets = count(p([6,7,8],4) /= 0._dp)
          else
              njets = count(p([5,6,7,8],4) /= 0._dp)
          endif

          numb = count(jet_bness(1:njets) /= 0)
          numnonb = count(jet_bness(1:njets) == 0)

          ! find b-jets
          if (numb > 0) then
              bjetindices(:) = 0._dp
              k = 1
              do j=1,njets
                  if (jet_bness(j) /= 0) then
                      bjetindices(k) = j
                      k = k + 1
                  endif
              enddo
          endif

          ! find non-b-jets
          nonbjetindices(:) = 0._dp
          k = 1
          do j=1,njets
              if (jet_bness(j) == 0) then
                  nonbjetindices(k) = j
                  k = k + 1
              endif
          enddo

          if (singletop2_nnlo_fully_inclusive) then
              lastnonjet = 5
          else
              lastnonjet = 4
          endif


          if (singletop2_nnlo_fully_inclusive) then
              ptop(:) = ptilde(currentNd,3,:) + ptilde(currentNd,4,:) + ptilde(currentNd,5,:)
          else
          ! top reconstruction with leading b bjet
              ptop(:) = ptilde(currentNd,3,:) + ptilde(currentNd,4,:) + p(bjetindices(1)+lastnonjet,:)
          endif

          pttop = ptpure(ptop)


          if (numnonb > 0) then
              ptlj1 = ptpure(p(nonbjetindices(1)+lastnonjet,:))
          else
              ptlj1 = -1._dp
          endif
          if (numnonb > 1) then
              ptlj2 = ptpure(p(nonbjetindices(2)+lastnonjet,:))
          else
              ptlj2 = -1._dp
          endif
          if (numnonb > 2) then
              ptlj3 = ptpure(p(nonbjetindices(3)+lastnonjet,:))
          else
              ptlj3 = -1._dp
          endif

          ptj1 = ptpure(p(lastnonjet+1,:))
          ptj2 = ptpure(p(lastnonjet+2,:))

          if (numb > 0) then
              ptb1 = ptpure(p(bjetindices(1)+lastnonjet,:))
          else
              ptb1 = -1._dp
          endif
          if (numb > 1) then
              ptb2 = ptpure(p(bjetindices(2)+lastnonjet,:))
          else
              ptb2 = -1._dp
          endif
          if (numb > 2) then
              ptb3 = ptpure(p(bjetindices(3)+lastnonjet,:))
          else
              ptb3 = -1._dp
          endif

          ptl = ptpure(ptilde(currentNd,4,:))

          if (numb > 0) then
              ptbl = ptpure(p(bjetindices(1)+lastnonjet,:) + ptilde(currentNd,4,:))
          else
              ptbl = -1._dp
          endif

          jetbin = real(njets,dp) + 0.5_dp
          bjetbin = real(numb,dp) + 0.5_dp
          ljetbin = real(numnonb,dp) + 0.5_dp

          etal = etarappure(ptilde(currentNd,4,:))
          ytop = yrappure(ptop)
          if (numnonb > 0) then
              etalj1 = etarappure(p(nonbjetindices(1)+lastnonjet,:))
          else
              etalj1 = -100._dp
          endif
          if (numnonb > 1) then
              etalj2 = etarappure(p(nonbjetindices(2)+lastnonjet,:))
          else
              etalj2 = -100._dp
          endif
          if (numnonb > 2) then
              etalj3 = etarappure(p(nonbjetindices(3)+lastnonjet,:))
          else
              etalj3 = -100._dp
          endif

          etaj1 = etarappure(p(lastnonjet+1,:))
          etaj2 = etarappure(p(lastnonjet+2,:))

          if (numb > 0) then
              etab1 = etarappure(p(bjetindices(1)+lastnonjet,:))
          else
              etab1 = -100._dp
          endif
          if (numb > 1) then
              etab2 = etarappure(p(bjetindices(2)+lastnonjet,:))
          else
              etab2 = -100._dp
          endif
          if (numb > 2) then
              etab3 = etarappure(p(bjetindices(3)+lastnonjet,:))
          else
              etab3 = -100._dp
          endif

          if (numb > 0) then
              mbl = puremass(p(bjetindices(1)+lastnonjet,:) + ptilde(currentNd,4,:))
          else
              mbl = -1._dp
          endif


!         costheta = 0._dp

          if (numb > 0 .and. numnonb > 0) then

          ptoprest(:) = 0._dp
          ptoprest(4) = mt

          plep(:) = ptilde(currentNd,4,:)
          plight(:) = p(nonbjetindices(1) + lastnonjet,:)

          call boostx(plep,ptop,ptoprest,pleprest)
          call boostx(plight,ptop,ptoprest,plightrest)

!         costheta = sum(pleprest(1:3)*plightrest(1:3)) &
!               /sqrt(sum(pleprest(1:3)**2))/sqrt(sum(plightrest(1:3)**2))


          !!!
          !!! FANCY ANGLES
          !!!

          ! construction of coordinate system in top rest frame
          ! the following observables probe the decay part (1005.5382)

          pw(:) = p(3,:) + p(4,:)
          pwrest_W(:) = 0._dp
          pwrest_W(4) = puremass(pw)

          ! pwrest is W vector in top rest frame
          ! pwrest_W is W vector in W rest fram

          ! pwrest defines z axis
          call boostx(pw,ptop,ptoprest,pwrest)

          ! normalize axis
          z(1:3) = pwrest(1:3)
          z(:) = z(:) / euclideanNorm(z)

          ! pspecrest defines \hat{s}_t
          call boostx(p(nonbjetindices(1)+lastnonjet,:), ptop, ptoprest, pspecrest)

          st(1:3) = pspecrest(1:3)
          st(:) = st(:) / euclideanNorm(st)

          y(:) = crossprod(z, st)
          y(:) = y(:) / euclideanNorm(y)

          x(:) = crossprod(y,z)
          x(:) = x(:) / euclideanNorm(x)

          ! lepton in W rest frame
          call boostx(plep, pw, pwrest_W, plepWrest)

          coslstar = sum(z(:)*plepWrest(1:3))/ euclideanNorm(plepWrest(1:3))

          coslN = sum(-y(:)*plepWrest(1:3)) / euclideanNorm(plepWrest(1:3))

          coslT = sum(x(:)*plepWrest(1:3)) / euclideanNorm(plepWrest(1:3))

          plepWrestproj(:) = plepWrest(1:3) - sum(plepWrest(1:3)*z(1:3)) * z(1:3)

          ! angle between N (-y) and the projection
          cosphiNstar = sum(-y(:)*plepWrestproj(:)) / euclideanNorm(plepWrestproj)

          ! angle between T (x) and the projection
          cosphiTstar = sum(x(:)*plepWrestproj(:)) / euclideanNorm(plepWrestproj)


          ! construction of the coordinate system that probes the production part
          ! 1404.1585

          zprod(:) = st(:)

          ! opposite sign for p1, since MCFM uses negative incoming momenta
          if ( (-p(1,3)) * p(nonbjetindices(1)+lastnonjet,3) > 0d0 ) then
              beam(:) = p(1,:)
          else
              beam(:) = p(2,:)
          endif

          call boostx(beam, ptop, ptoprest, beamrest)

          yprod(:) = crossprod(zprod, beamrest(1:3))
          yprod(:) = yprod(:) / euclideanNorm(yprod)

          xprod(:) = crossprod(yprod, zprod)

          ! angles of lepton w.r.t xprod,yprod,zprod

          cosxprod = sum(xprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))
          cosyprod = sum(yprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))
          coszprod = sum(zprod(:)*pleprest(1:3)) / euclideanNorm(pleprest(1:3))

          !!! FINALIZATION

          endif

          ids = histos
          vals = [pttop,ptlj1,ptlj2,ptlj3,ptj1,ptj2,ptb1,ptb2,ptb3,ptl,ptbl, &
                    jetbin,bjetbin,etal,ytop,etalj1,etalj2,etalj3, &
                etaj1,etaj2,etab1,etab2,etab3,mbl, &
                cosxprod, cosyprod, coszprod, &
                coslstar, coslN, coslT, &
                cosphiNstar, cosphiTstar, ljetbin ]

          wts =  [wt,wt,wt,wt,wt,wt,wt,wt,wt,wt, &
                  wt,wt,wt,wt,wt,wt,wt,wt,wt,wt, &
                  wt,wt,wt,wt,wt,wt,wt,wt,wt,wt, &
                  wt,wt,wt]

          do j=1,33
              if (ieee_is_nan(vals(j)) .or. (.not. ieee_is_finite(vals(j)))) then
                  vals(j) = -100d0
              endif
          enddo

      end subroutine

end module
