!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c This is just a normal MCFM-style interface to the virtual matrix elements
      subroutine qqb_zaj_v_new_interface(p_mcfm,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp), intent(in) :: p_mcfm(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)
      real(dp) :: msq0(-nf:nf, -nf:nf)

      call qqb_zaj_v_new(p_mcfm,msq,msq0)

      end subroutine

      subroutine qqb_zaj_v_new(p_mcfm,msq,msq0)
          use helicities
          use zaj_virtamps_l_m, only : zaj_virtamp_l1, zaj_virtamp_l2
          use zaj_virtamps_m, only : zaj_virtamp_q1, zaj_virtamp_q2, zaj_virtamp_q3, zaj_virtamp_q4
          use zaj_treeamps_m, only : zaj_tree_isr_pp, zaj_tree_isr_pm,
     &               zaj_tree_isr_mp, zaj_tree_fsr_pp, zaj_tree_fsr_pm,
     &               zaj_tree_anomZZ_pp, zaj_tree_anomZZ_pm,
     &               zaj_tree_anomZZ_mp, zaj_tree_anomZa_pp,
     &               zaj_tree_anomZa_pm, zaj_tree_anomZa_mp,
     &               zaj_tree_anomaZ_pp, zaj_tree_anomaZ_pm,
     &               zaj_tree_anomaZ_mp
          use omega_ee3jet_m, only : amp_zaj_fsr_pp, amp_zaj_fsr_pm,
     &            amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm,
     &            amp_zaj_anomZa_pp, amp_zaj_anomZa_pm,
     &            amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm,
     &            zaj_crossings_l_new,
     &            zaj_crossings_anom_new
          use VVconfig_m

           implicit none
          include 'types.f'
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'qcdcouple.f'
          ! we use q1 here to determine if we need to compute the
          ! final state radiation pieces, they are cheap; but still...
          include 'zcouple.f'! only q1
          include 'epinv.f'! for subuv piece
          include 'nflav.f'! for subuv piec
          include 'scheme.f'
          include 'anomcoup.f'! logical :: anomtgc
          include 'blha.f'


          real(dp), intent(in) :: p_mcfm(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf, -nf:nf)
          real(dp), intent(out) :: msq0(-nf:nf, -nf:nf)

          complex(dp) :: amps_q1(2,2,2,2), amps_q2(2,2,2,2), amps_q3(2,2,2,2), amps_q4(2,2,2,2)
          complex(dp) :: amps_l2(2,2,2,2)
          complex(dp) :: ampsdress(2,2,2,2)

          complex(dp) :: ampsTree_q(2,2,2,2)
          complex(dp) :: ampsTreedress(2,2,2,2)

          complex(dp) :: amps_l1_new(0:1, 2,2,2,2)
          complex(dp) :: amps_anomZZ_new(0:1, 2,2,2,2)
          complex(dp) :: amps_anomZa_new(0:1, 2,2,2,2)
          complex(dp) :: amps_anomaZ_new(0:1, 2,2,2,2)

          integer (kind(decayElAntiEl)) :: vDecay

          complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          real(dp) :: p(mxpart,4)
          real(dp) :: subuv
          integer :: qi, i, si

          p(1,:) = p_mcfm(2,:) ! quark
          p(2,:) = p_mcfm(1,:) ! antiquark
          p(3,:) = p_mcfm(6,:) ! gluon
          p(4,:) = p_mcfm(5,:) ! photon
          p(5,:) = p_mcfm(4,:) ! antilepton
          p(6,:) = p_mcfm(3,:) ! lepton
          call spinoru(6,p,za,zb)

          vDecay = decayChannel()

          amps_q1 = 0._dp
          amps_q2 = 0._dp
          amps_q3 = 0._dp
          amps_q4 = 0._dp
          amps_l1_new = 0._dp
          amps_l2 = 0._dp
          ampsTree_q = 0._dp
          amps_anomZZ_new = 0._dp
          amps_anomZa_new = 0._dp
          amps_anomaZ_new = 0._dp
          msq = 0._dp
          msq0 = 0._dp

          if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype == 1))) then

          call zaj_virtamp_q1(1,2,3,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(1,2,3,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(1,2,3,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(1,2,3,4,5,6,za,zb,amps_q4)
          call zaj_crossings(1,2,3,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 1,2,3,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(1,2,3,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 1,2,3,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(qi,-qi) = aveqq*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(qi,-qi) = aveqq*sum(abs(ampsTreedress)**2)
          enddo

          endif

          if (useblha == 0) then

          call zaj_virtamp_q1(2,1,3,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(2,1,3,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(2,1,3,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(2,1,3,4,5,6,za,zb,amps_q4)
          call zaj_crossings(2,1,3,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 2,1,3,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(2,1,3,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 2,1,3,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(-qi,qi) = aveqq*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(-qi,qi) = aveqq*sum(abs(ampsTreedress)**2)
          enddo

          endif


          if (useblha == 0) then

          call zaj_virtamp_q1(3,2,1,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(3,2,1,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(3,2,1,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(3,2,1,4,5,6,za,zb,amps_q4)
          call zaj_crossings(3,2,1,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 3,2,1,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(3,2,1,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 3,2,1,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(qi,0) = aveqg*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(qi,0) = aveqg*sum(abs(ampsTreedress)**2)
          enddo

          endif


          if (useblha == 0) then

          call zaj_virtamp_q1(2,3,1,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(2,3,1,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(2,3,1,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(2,3,1,4,5,6,za,zb,amps_q4)
          call zaj_crossings(2,3,1,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 2,3,1,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(2,3,1,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 2,3,1,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(-qi,0) = aveqg*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(-qi,0) = aveqg*sum(abs(ampsTreedress)**2)
          enddo

          endif


          if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype == 2))) then

          call zaj_virtamp_q1(3,1,2,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(3,1,2,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(3,1,2,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(3,1,2,4,5,6,za,zb,amps_q4)
          call zaj_crossings(3,1,2,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 3,1,2,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(3,1,2,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 3,1,2,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(0,qi) = aveqg*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(0,qi) = aveqg*sum(abs(ampsTreedress)**2)
          enddo

          endif


          if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype == 3))) then

          call zaj_virtamp_q1(1,3,2,4,5,6,za,zb,amps_q1)
          call zaj_virtamp_q2(1,3,2,4,5,6,za,zb,amps_q2)
          call zaj_virtamp_q3(1,3,2,4,5,6,za,zb,amps_q3)
          call zaj_virtamp_q4(1,3,2,4,5,6,za,zb,amps_q4)
          call zaj_crossings(1,3,2,4,5,6,za,zb,
     &       zaj_tree_isr_pp, zaj_tree_isr_pm, zaj_tree_isr_mp, ampsTree_q)
          if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
            call zaj_crossings_l_new(p, 1,3,2,4,5,6,za,zb,amps_l1_new)
            call zaj_virtamp_l2(1,3,2,4,5,6,za,zb,amps_l2)
          endif
          if (anomtgc) then
            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomZZ_pp, amp_zaj_anomZZ_pm, amps_anomZZ_new)
            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomZa_pp, amp_zaj_anomZa_pm, amps_anomZa_new)
            call zaj_crossings_anom_new(p, 1,3,2,4,5,6,za,zb,
     &          amp_zaj_anomaZ_pp, amp_zaj_anomaZ_pm, amps_anomaZ_new)
          endif

          do qi=1,2
            call zaj_treeamps_dress(qi,za,zb,ampsTree_q,
     &              amps_l1_new(0,:,:,:,:),
     &              amps_anomZZ_new(0,:,:,:,:), amps_anomZa_new(0,:,:,:,:),
     &              amps_anomaZ_new(0,:,:,:,:), ampsTreedress)
            call zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &              amps_l1_new(1,:,:,:,:),amps_l2,
     &              amps_anomZZ_new(1,:,:,:,:), amps_anomZa_new(1,:,:,:,:),
     &              amps_anomaZ_new(1,:,:,:,:), ampsdress)
            msq(0,-qi) = aveqg*2*sum(real(conjg(ampsTreedress)*ampsdress))
            msq0(0,-qi) = aveqg*sum(abs(ampsTreedress)**2)
          enddo

          endif


          do i=3,Nf
            si = 2 - modulo(i,2) ! this is the sequence 1,2,1,2,1,..

            msq(-i,i) = msq(-si,si)
            msq(i,-i) = msq(si,-si)
            msq(-i,0) = msq(-si,0)
            msq(i,0) = msq(si,0)
            msq(0,-i) = msq(0,-si)
            msq(0,i) = msq(0,si)

            msq0(-i,i) = msq0(-si,si)
            msq0(i,-i) = msq0(si,-si)
            msq0(-i,0) = msq0(-si,0)
            msq0(i,0) = msq0(si,0)
            msq0(0,-i) = msq0(0,-si)
            msq0(0,i) = msq0(0,si)
          enddo

          msq0(:,:) = msq0(:,:) * 8._dp
          msq(:,:) = msq(:,:) * 8._dp * ason2pi / 2._dp

          scheme = 'dred'
          subuv=ason2pi*xn*(epinv*(11-2*nflav/xn)-1)/6
          msq(:,:) = msq(:,:) - subuv*msq0(:,:)

      end subroutine

c     dresses raw helicity amps obtained from zaj_treeamps with
c       coupling factors, Z propagator, given quark flavor qi
      subroutine zaj_virtamps_dress(qi,za,zb,amps_q1,amps_q2,amps_q3,amps_q4,
     &             amps_l1,amps_l2,amps_anomZZ,amps_anomZa,amps_anomaZ,amps)
          use helicities

          implicit none
          include 'types.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'nf.f'
          include 'zcouple.f'
          include 'zcouple_cms.f'
          include 'ewcouple.f'
          include 'ewcharge.f'
          include 'qcdcouple.f'
          include 'masses.f'

          integer, intent(in) :: qi
          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          complex(dp), intent(in) :: amps_q1(2,2,2,2), amps_q2(2,2,2,2),
     &                               amps_q3(2,2,2,2), amps_q4(2,2,2,2),
     &                               amps_l1(2,2,2,2), amps_l2(2,2,2,2)
          complex(dp), intent(in) :: amps_anomZZ(2,2,2,2), amps_anomZa(2,2,2,2)
          complex(dp), intent(in) :: amps_anomaZ(2,2,2,2)
          complex(dp), intent(out) :: amps(2,2,2,2)

          integer :: h1,h2,h3,h4
          complex(dp) :: vq, vl
          complex(dp) :: propzQQ, propzQL
          complex(dp) :: propZ_56, propZ_456, propA_56, propA_456
          real(dp) :: sumQiNf
          complex(dp) :: sumViNf, sinw, cosw, sumQiTwoTF

          integer :: i1,i2,i3
          real(dp) :: s,t
          real(dp), parameter :: TwoTF(1:nf)=(/ -1._dp, +1._dp, -1._dp, +1._dp, -1._dp /)
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
          t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

          propzQQ = s(5,6)/(s(5,6) - zmass**2 + im*zmass*zwidth)
          propzQL =t(4,5,6)/(t(4,5,6) - zmass**2 + im*zmass*zwidth)

          propZ_56 = 1/(s(5,6) - zmass**2 + im*zmass*zwidth)
          propZ_456 = 1/(t(4,5,6) - zmass**2 + im*zmass*zwidth)
          propA_56 = 1/s(5,6)
          propA_456 = 1/t(4,5,6)

          sumQiNf = sum(Q(1:Nf)**2)
          sumViNf = sum((zl(1:Nf)+zr(1:Nf))*Q(1:Nf))
          sumQiTwoTF = sum(Q(1:Nf)*TwoTF(1:Nf))
          sinw = sqrt(zxw)
          cosw = sqrt(1._dp - zxw)


          amps(:,:,:,:) = 0._dp

          ! (h3,h3,h2,h1) = helicities of (q,l,g,gam)
          do h1=1,2
            do h2=1,2
              do h3=1,2
                if (h3 == helLeft) then
                  vl = zl1
                elseif (h3 == helRight) then
                  vl = zr1
                endif

                do h4=1,2
                  if (h4 == helLeft) then
                    vq = zl(qi)
                  elseif (h4 == helRight) then
                    vq = zr(qi)
                  endif

c   Axial quark loop contributions are included in the limit mt -> infinity
                  amps(h4,h3,h2,h1) =
     &              twort2*sqrt(zesq)**3*sqrt(gsq)*(
     &              amps_q1(h4,h3,h2,h1) * Q(qi)*(q1*Q(qi) + vl*vq*propzQQ)
     &            + amps_q2(h4,h3,h2,h1) * (q1*sumQiNf + vl*sumViNf*propzQQ/2)
     &            + amps_q3(h4,h3,h2,h1) * (vl/2/sinw/cosw*propzQQ)*sumQiTwoTF
     &            + amps_q4(h4,h3,h2,h1) * Q(qi)*(vl/2/sinw/cosw*propzQQ)
     &            + amps_l1(h4,h3,h2,h1) * q1*(q1*Q(qi)*propA_456 + vl*vq*propZ_456)
     &            + amps_l2(h4,h3,h2,h1) * q1*(vl/2/sinw/cosw*propzQL)
     &            + amps_anomZZ(h4,h3,h2,h1) * vl*vq * propZ_456 * propZ_56
     &            + amps_anomZa(h4,h3,h2,h1) * q1*vq * propZ_456 * propA_56
     &            + amps_anomaZ(h4,h3,h2,h1) * Q(qi)*vl * propA_456 * propZ_56
     &              )
                enddo
              enddo
            enddo
          enddo


      end subroutine

