!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module singletop_jet3
      use ieee_arithmetic
      use types
      use singletop2_scale_m
      use singletop_jet2, only: Wtoponshell_gen
      use singletop_jetdeps

      public :: singletop_jet_decay
      public :: singletop_jet_decay_all
      public :: singletop_jet_decay_virt_all
      public :: singletop_jet_decay_real_all
      public :: singletop_jet_decay_gvec
      public :: singletop_jet_decay_gs

      private :: tdecay_gen,adecay_gen

      private :: singletop_jet_decay_intdip

      private

      contains

#define WITH_B_BBAR 4

      subroutine singletop_jet_decay_gvec(p,n,ip,msq)
          use types
          use singletop2_nnlo_vars
          use singletop_jet2_deps, only: wtgvecn
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'nwz.f'
          include 'masses.f'
          include 'constants.f'
          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(in) :: n(4)
          integer, intent(in) :: ip
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

          integer:: i3,i4,i1,i6,iq,ibeam1,ibeam2
          real(dp):: msq_q_b,msq_qb_b,
     &     s16,dot,fac

          integer :: csave

          if (ip /=7) then
            write(6,*) 'Error in singletop_jet_decay_gvec, in is not 7 :',ip
            error stop
          endif

          msq(:,:)=0._dp

          csave = corr_on_beam

c loop over corrections for initial-state b-quark in both beams
          do corr_on_beam=1,max_corr_on_beam

          if (corr_on_beam == 1) then
            ibeam1=2
            ibeam2=1
          elseif (corr_on_beam == 2) then
            ibeam1=1
            ibeam2=2
          endif

c--- set up lepton variables depending on whether it's t or tbar
          if     (nwz == -1) then
            i3=3
            i4=4
            i1=6
            i6=ibeam1
            iq=-1
            write(6,*) 'singletop_jet_decay_gvec not tested for t-bar yet!'
            stop
          elseif (nwz == +1) then
            i3=4
            i4=3
            i1=ibeam1
            i6=6
            iq=+1
          else
            write(6,*) 'Error in singletop_jet_decay_gvec, nwz is not +1 or -1 :',nwz
            stop
          endif

c--- fac adds additional color factor and corrects for
c--- (unnecessary) W width in propagator
          s16=two*dot(p,i1,i6)
          fac=xn*((s16-wmass**2)**2+(wmass*wwidth)**2)/(s16-wmass**2)**2

          msq_q_b=aveqq*fac*wtgvecn(mt,twidth,7,5,i3,i4,i1,i6,ibeam2,p,n)
          msq_qb_b=aveqq*fac*wtgvecn(mt,twidth,7,5,i3,i4,i6,i1,ibeam2,p,n)

          if (corr_on_beam == 2) then
            msq([2,4],5*iq)=msq_q_b
            msq([-1,-3],5*iq)=msq_qb_b
          else
            msq(5*iq,[2,4])=msq_q_b
            msq(5*iq,[-1,-3])=msq_qb_b
          endif

          enddo

          corr_on_beam = csave

      end subroutine

      subroutine singletop_jet_decay_intdip(p,msq)
          use singletop2_nnlo_vars
          use singletop2_scale_m
          use types
          use singletop2_scet_heavy_decay, only: intdip_fi_mt, intdip_fi_mt_gg
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'scheme.f'
          include 'colstruc.f'
          include 'nflav.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          real(dp) :: ff_qq, ff_gg, ff_gq, dotvec, xl57

          real(dp) :: msqborn(-nf:nf,-nf:nf, max_bcontrib,max_corr_on_beam)
          real(dp) :: beam2,beam1


          integer, parameter :: noglue(4) = [-3,-1,2,4]

          scheme = 'tH-V'
          caonly = .true.

          call singletop_jet_decay_all(p,msqborn)

          msq = 0._dp
          beam1 = 0._dp
          beam2 = 0._dp
          corr_on_beam = 2

          !! gg contribution

          ! for 87top with 5 recoil and 78top with 5 recoil, 1/2 identical particle
          ! fixed by log(1 - p345^2/mt^2)
          beam2 = beam2 + 2*ca*(1._dp/2._dp) * (intdip_fi_mt_gg(p, 5) * 2._dp)

          ! for 58top with 7 recoil and 57top with 8 recoil, 1/2 identical particle
          ! fixed by log(1 - p347^2/mt^2)
          beam2 = beam2 - cf*((ca-2*cf)/cf/2._dp)  * (intdip_fi_mt(p, 7) * 2._dp)

          ! All factors below this point understood

          ! following two fixed by log(p57^2/mu^2)
          xl57=log(2._dp*dotvec(p(5,:),p(7,:))/renscale_beam2_isheavy_onheavy**2)

          beam2 = beam2 + ca*(1._dp/2._dp) * ff_gg(0._dp,xl57,1)

          ! 587 + 578
          beam2 = beam2 + 2._dp*cf*(ca/cf/2._dp) * ff_qq(0._dp,xl57,1)

c         !! bb contribution
          beam2 = beam2 + 2._dp*tr*(1._dp/2._dp) * 2._dp
     &              * ff_gq(0._dp,xl57,1) * (ca / nflav / (2._dp*tr))

c         !! qq contrib
          beam2 = beam2 + 2._dp*tr*(nf-1)*(1._dp/2._dp) * 2._dp
     &              * ff_gq(0._dp,xl57,1) * (ca / nflav / (2._dp*tr))

          msq(noglue,5, 1,2) = msqborn(noglue,5, 1,2) * beam2 * as_heavy_beam2/4._dp/pi


          corr_on_beam = 1

          ! for 87top with 5 recoil and 78top with 5 recoil, 1/2 identical particle
          ! fixed by log(1 - p345^2/mt^2)
          beam1 = beam1 + 2*ca*(1._dp/2._dp) * (intdip_fi_mt_gg(p, 5) * 2._dp)

          ! for 58top with 7 recoil and 57top with 8 recoil, 1/2 identical particle
          ! fixed by log(1 - p347^2/mt^2)
          beam1 = beam1 - cf*((ca-2*cf)/cf/2._dp)  * (intdip_fi_mt(p, 7) * 2._dp)

          ! All factors below this point understood

          ! following two fixed by log(p57^2/mu^2)
          xl57=log(2._dp*dotvec(p(5,:),p(7,:))/renscale_beam1_isheavy_onheavy**2)

          beam1 = beam1 + ca*(1._dp/2._dp) * ff_gg(0._dp,xl57,1)

          ! 587 + 578
          beam1 = beam1 + 2._dp*cf*(ca/cf/2._dp) * ff_qq(0._dp,xl57,1)

c         !! bb contribution
          beam1 = beam1 + 2._dp*tr*(1._dp/2._dp) * 2._dp
     &              * ff_gq(0._dp,xl57,1) * (ca / nflav / (2._dp*tr))

c         !! qq contrib
          beam1 = beam1 + 2._dp*tr*(nf-1)*(1._dp/2._dp) * 2._dp
     &              * ff_gq(0._dp,xl57,1) * (ca / nflav / (2._dp*tr))


          msq(5,noglue, 1,1) = msqborn(5,noglue, 1,1) * beam1 * as_heavy_beam1/4._dp/pi

      end subroutine singletop_jet_decay_intdip

      subroutine singletop_jet_decay_gs(p,ndmx,msqc)
          use singletop2_nnlo_vars
          use singletop2_scale_m
          use types
          use singletop2_scet_heavy_decay, only: dips_fi_mt, dips_fi_mt_gg, dips_ff_ident
          use singletop2_scet_light, only: singletop2_scet_tree_ub,
     &        singletop2_scet_tree_bu
          implicit none
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'qqgg.f'

          real(dp), intent(in) :: p(mxpart,4)
          integer, intent(in) :: ndmx
          real(dp), intent(out) :: msqc(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          real(dp) :: msqv(-nf:nf,-nf:nf), subv
          real(dp) :: msqborn(-nf:nf,-nf:nf), sub(4)

          integer, parameter :: noglue(4) = [-3,-1,2,4]

          external donothing_gvec

          msqc = 0._dp

          corr_islight = .false. ! for singletop2_set_dipscale
          corr_beam1 = .false. ! for singletop2_set_dipscale
          corr_on_beam = 2

          ! tested cases for g g piece
          ! 7 soft, 8 soft, 5||7, 5||8, 7||8

          ! All factors except those related to dips_fi_mt* understood

          call dips_fi_mt_gg(1, p, singletop_jet_decay, singletop_jet_decay_gvec, msqborn, 8,7,5)
          msqc(1,noglue,5, 1,2) = -4*ca*msqborn(noglue,5) * (1._dp/2._dp) / 2._dp

          call dips_fi_mt_gg(2, p, singletop_jet_decay, singletop_jet_decay_gvec, msqborn, 7,8,5)
          msqc(2,noglue,5, 1,2) = -4*ca*msqborn(noglue,5) * (1._dp/2._dp) / 2._dp

          call dips_fi_mt(3, p, singletop_jet_decay, msqborn, 5,8,7)
          msqc(3,noglue,5, 1,2) = -cf*msqborn(noglue,5) * (ca-2*cf)/cf/2._dp / 2._dp

          call dips(4, p, 5,8,7, sub, subv, msqborn, msqv, singletop_jet_decay, donothing_gvec)
          ! identical particle factor of 1/2
          msqc(4,noglue,5, 1,2) = 2._dp*cf*msqborn(noglue,5)*sub(qq) * (ca/cf/2._dp) / 2._dp

          call dips_fi_mt(5, p, singletop_jet_decay, msqborn, 5,7,8)
          msqc(5,noglue,5, 1,2) = -cf*msqborn(noglue,5) * (ca-2*cf)/cf/2._dp / 2._dp

          call dips(6, p, 5,7,8, sub, subv, msqborn, msqv, singletop_jet_decay, donothing_gvec)
          ! identical particle factor of 1/2
          msqc(6,noglue,5, 1,2) = 2._dp*cf*msqborn(noglue,5)*sub(qq) * (ca/cf/2._dp) / 2._dp

          ! identical particle factor of 1/2
          call dips(7, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          msqc(7,noglue,5, 1,2) = 2._dp*ca*(msqborn(noglue,5)*sub(gg) + msqv(noglue,5)*subv) * (1._dp/2._dp) * (1._dp/2._dp)

          ! g -> b b~
          call dips_ff_ident(8, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! identical particle factor of 1/2 cancelled by factor of two because, for g->qqb splitting,
          ! we include only 78_5, not 78_5 + 78_top
          msqc(8,noglue,5, 1,2) = 2._dp*tr*(msqborn(noglue,5)*sub(gq) - msqv(noglue,5)*subv) * (1._dp/2._dp)

          call dips_ff_ident(10, p, 5,8,7, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! identical particle factor of 1/2 cancelled by factor of two because, for g->qqb splitting,
          ! we include only 78_5, not 78_5 + 78_top
          msqc(10,noglue,5, 1,2) = 2._dp*tr*(msqborn(noglue,5)*sub(gq) - msqv(noglue,5)*subv) * (1._dp/2._dp)

          ! g -> q q~
          call dips(9, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! factor of two because, for g->qqb splitting, we include only 78_5, not 78_5 + 78_top
          msqc(9,noglue,5, 1,2) = 2._dp*tr*(nf-1)*(msqborn(noglue,5)*sub(gq) - msqv(noglue,5)*subv) * (1._dp/2._dp) * 2._dp


          corr_islight = .false. ! for singletop2_set_dipscale
          corr_beam1 = .true. ! for singletop2_set_dipscale
          corr_on_beam = 1

          call dips_fi_mt_gg(11, p, singletop_jet_decay, singletop_jet_decay_gvec, msqborn, 8,7,5)
          msqc(11,5,noglue, 1,1) = -4*ca*msqborn(5,noglue) * (1._dp/2._dp) / 2._dp

          call dips_fi_mt_gg(12, p, singletop_jet_decay, singletop_jet_decay_gvec, msqborn, 7,8,5)
          msqc(12,5,noglue, 1,1) = -4*ca*msqborn(5,noglue) * (1._dp/2._dp) / 2._dp

          call dips_fi_mt(13, p, singletop_jet_decay, msqborn, 5,8,7)
          msqc(13,5,noglue, 1,1) = -cf*msqborn(5,noglue) * (ca-2*cf)/cf/2._dp / 2._dp

          call dips(14, p, 5,8,7, sub, subv, msqborn, msqv, singletop_jet_decay, donothing_gvec)
          ! identical particle factor of 1/2
          msqc(14,5,noglue, 1,1) = 2._dp*cf*msqborn(5,noglue)*sub(qq) * (ca/cf/2._dp) / 2._dp

          call dips_fi_mt(15, p, singletop_jet_decay, msqborn, 5,7,8)
          msqc(15,5,noglue, 1,1) = -cf*msqborn(5,noglue) * (ca-2*cf)/cf/2._dp / 2._dp

          call dips(16, p, 5,7,8, sub, subv, msqborn, msqv, singletop_jet_decay, donothing_gvec)
          ! identical particle factor of 1/2
          msqc(16,5,noglue, 1,1) = 2._dp*cf*msqborn(5,noglue)*sub(qq) * (ca/cf/2._dp) / 2._dp

          ! identical particle factor of 1/2
          call dips(17, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          msqc(17,5,noglue, 1,1) = 2._dp*ca*(msqborn(5,noglue)*sub(gg) + msqv(5,noglue)*subv) * (1._dp/2._dp) * (1._dp/2._dp)

          ! g -> b b~
          call dips_ff_ident(18, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! identical particle factor of 1/2 cancelled by factor of two because, for g->qqb splitting,
          ! we include only 78_5, not 78_5 + 78_top
          msqc(18,5,noglue, 1,1) = 2._dp*tr*(msqborn(5,noglue)*sub(gq) - msqv(5,noglue)*subv) * (1._dp/2._dp)

          call dips_ff_ident(20, p, 5,8,7, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! identical particle factor of 1/2 cancelled by factor of two because, for g->qqb splitting,
          ! we include only 78_5, not 78_5 + 78_top
          msqc(20,5,noglue, 1,1) = 2._dp*tr*(msqborn(5,noglue)*sub(gq) - msqv(5,noglue)*subv) * (1._dp/2._dp)

          ! g -> q q~
          call dips(19, p, 7,8,5, sub, subv, msqborn, msqv, singletop_jet_decay, singletop_jet_decay_gvec)
          ! factor of two because, for g->qqb splitting, we include only 78_5, not 78_5 + 78_top
          msqc(19,5,noglue, 1,1) = 2._dp*tr*(nf-1)*(msqborn(5,noglue)*sub(gq) - msqv(5,noglue)*subv) * (1._dp/2._dp) * 2._dp

      end subroutine singletop_jet_decay_gs

      ! copied from Wt/tree.f
      subroutine Wttree(mq,ig,is,ie,in,it,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: is,ig,ie,in,it,i,j
      complex(dp):: amp(2,2)
      real(dp):: mq,prop

      prop=sqrt((real(za(ie,in)*zb(in,ie))-wmass**2)**2
     &          +(wmass*wwidth)**2)

c---- helitities: amp(ht,hg)
c--- labels on amplitudes represent helicities for the heavy quark
c--- and the gluon respectively
c---   1 = negative helitity, 2 = positive helitity
c--- heavy quark momentum is made massless (it) with the gluon momentum ig
      amp(1,2)=za(ie,it)
     & /za(ig,is)/za(ig,it)*(za(is,it)*zb(is,in)+za(ig,it)*zb(ig,in))
      amp(2,2)=-mq/za(ig,it)*za(ig,ie)
     & /za(ig,is)/za(ig,it)*(za(is,it)*zb(is,in)+za(ig,it)*zb(ig,in))
      amp(1,1)=-(za(ig,ie)*zb(ig,is)+za(ie,it)*zb(is,it))*zb(is,in)
     & /zb(ig,is)/zb(ig,it)
      amp(2,1)=-mq/za(ig,it)/zb(it,ig)*za(ig,ie)*zb(is,in)*zb(is,it)
     & /zb(ig,is)

      do i=1,2
      do j=1,2
      amp(i,j)=amp(i,j)/prop
      enddo
      enddo

      return
      end

      subroutine singletop_jet_decay_all(p,msq)
          use types
          use singletop2_nnlo_vars
          !use NNTopDec, only: top23treesqamp
      implicit none
c Crossed from routines representing radiation in production
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
      include 'zprods_com.f'
      integer:: ht,hb,h2,i1,i2

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp):: fac,msq_q_b,msq_qb_b
      complex(dp):: prop
      complex(dp):: mrad(2,2,2),
     & mtotq_b(2,2),mtotqb_b(2,2),ampq_b(2),ampqb_b(2)

      msq = 0._dp

      do corr_on_beam=1,max_corr_on_beam

          mtotq_b(:,:)=czip
          mtotqb_b(:,:)=czip

          if (corr_on_beam == 1) then
            i1=2
            i2=1
          elseif (corr_on_beam == 2) then
            i1=1
            i2=2
          endif

          if (nwz == +1) then
             call Watoponshell_gen(5,7,3,4,p,0,mrad)
c 3 is the vector used to make the top quark massless in Watoponshell_gen
             call adecay_gen(p,i2,6,i1,3,ampq_b)
             call adecay_gen(p,i2,i1,6,3,ampqb_b)

             do hb=1,1; do h2=1,2; do ht=1,2
                 mtotq_b(hb,h2)=mtotq_b(hb,h2)+ampq_b(ht)*mrad(hb,h2,ht)
                 mtotqb_b(hb,h2)=mtotqb_b(hb,h2)+ampqb_b(ht)*mrad(hb,h2,ht)
             enddo; enddo; enddo
          else
              error stop __FILE__//": nwz=-1 not implemented"
          endif

          prop=cplx2(zip,mt*twidth)

          ! this does not really matter since we're in the decay piece
          ! which has a completely separate scale (see singletop2_scale_setup)
          if (corr_on_beam == 2) then
              fac=V*xn*gwsq**4*(4*pi*as_heavy_beam2)/abs(prop)**2
          else
              fac=V*xn*gwsq**4*(4*pi*as_heavy_beam1)/abs(prop)**2
          endif

          msq_q_b = sum(abs(mtotq_b)**2)
          msq_qb_b = sum(abs(mtotqb_b)**2)

          if (corr_on_beam == 2) then
            msq([2,4],5, 1,2)=aveqq*fac*msq_q_b
            msq([-1,-3],5, 1,2)=aveqq*fac*msq_qb_b
          else
            msq(5,[2,4], 1,1)=aveqq*fac*msq_q_b
            msq(5,[-1,-3], 1,1)=aveqq*fac*msq_qb_b
          endif


      enddo

      !write (*,*) "ub = ", msq(2,5, 1,2)
      !call top23treesqamp(p, as_heavy_beam2, newub)
      !write (*,*) "ubnew = ", newub / msq(2,5, 1,2)

      !pause

      end subroutine singletop_jet_decay_all


      ! for use with dipole subtractions, gives both corr_on_beam pieces
      subroutine singletop_jet_decay(p,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'types.f'
c Crossed from routines representing radiation in production
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
      include 'zprods_com.f'

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

      integer:: ht,hb,h2,i1,i2
      real(dp):: fac,msq_q_b,msq_qb_b
      complex(dp):: prop
      complex(dp):: mrad(2,2,2),
     & mtotq_b(2,2),mtotqb_b(2,2),ampq_b(2),ampqb_b(2)
      integer :: csave

      msq = 0._dp

      csave = corr_on_beam

      do corr_on_beam=1,max_corr_on_beam
c---initialize

      if (corr_on_beam == 1) then
        i1=2
        i2=1
      elseif (corr_on_beam == 2) then
        i1=1
        i2=2
      endif

      if (nwz == +1) then
          mtotq_b(:,:)=czip
          mtotqb_b(:,:)=czip
          mrad = 0._dp

         call Watoponshell_gen(5,7,3,4,p,0,mrad)
c 3 is the vector used to make the top quark massless in Watoponshell_gen
         call adecay_gen(p,i2,6,i1,3,ampq_b)
         call adecay_gen(p,i2,i1,6,3,ampqb_b)

         do hb=1,1; do h2=1,2; do ht=1,2
             mtotq_b(hb,h2)=mtotq_b(hb,h2)+ampq_b(ht)*mrad(hb,h2,ht)
             mtotqb_b(hb,h2)=mtotqb_b(hb,h2)+ampqb_b(ht)*mrad(hb,h2,ht)
         enddo; enddo; enddo
      else
          error stop __FILE__//": nwz=-1 not implemented"
      endif

      prop=cplx2(zip,mt*twidth)
      if (corr_on_beam == 2) then
          fac=V*xn*gwsq**4*(4*pi*as_heavy_beam2)/abs(prop)**2
      else
          fac=V*xn*gwsq**4*(4*pi*as_heavy_beam1)/abs(prop)**2
      endif

      msq_q_b=zero
      msq_qb_b=zero
      do hb=1,2
      do h2=1,2
      msq_q_b=msq_q_b+abs(mtotq_b(hb,h2))**2
      msq_qb_b=msq_qb_b+abs(mtotqb_b(hb,h2))**2
      enddo
      enddo

      if (corr_on_beam == 2) then

        msq(2,5)=aveqq*fac*msq_q_b
        msq(-1,5)=aveqq*fac*msq_qb_b

        msq(4,5)=msq(2,5)
        msq(-3,5)=msq(-1,5)

      else

        msq(5,2)=aveqq*fac*msq_q_b
        msq(5,-1)=aveqq*fac*msq_qb_b

        msq(5,4)=msq(5,2)
        msq(5,-3)=msq(5,-1)

      endif

      enddo

      corr_on_beam = csave

      return
      end

      subroutine singletop_jet_decay_virt_all(p,msq)
          use types
          use singletop2_nnlo_vars
          use singletop_jet2_deps, only: virt_pp, virt_mm, virt_pm, virt_mp
      implicit none
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

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp) :: msqintdip(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      integer:: ibeam1,ibeam2,i6,i1,i3,i4,iq,ig,ib,k1,k6,hg
      real(dp) :: q(mxpart,4),fac,r(mxpart,4),wprop,msqbit_q_b,msqbit_qb_b
      complex(dp):: ampl0_q_b(2),amplv_q_b(2),ampl0_qb_b(2),amplv_qb_b(2),
     & amplp(2,2),ampld_q_b(2),ampld_qb_b(2)
      real(dp):: twotDg,dot
      complex(dp):: s_ft(2,2)

      real(dp) :: musq

      scheme = 'tH-V'

      msq = 0._dp

      do corr_on_beam=1,max_corr_on_beam

          if (corr_on_beam == 1) then
            ibeam1=2
            ibeam2=1
          elseif (corr_on_beam == 2) then
            ibeam1=1
            ibeam2=2
          endif

c--- set up lepton variables depending on whether it's t or tbar
          if     (nwz == -1) then
c anti top quark, so need crossed amplitudes for top quark
            k6=3
            k1=4
            i3=6
            i4=ibeam1
            iq=-1 ! antitop quark
          elseif (nwz == +1) then
c top quark, so need crossed amplitudes for antitop quark
            k6=4
            k1=3
            i3=ibeam1
            i4=6
            iq=1 ! top quark
          else
            write(6,*) 'Error in singletop_jet_decay_virt, nwz is not +1 or -1 :   ',nwz
            stop
          endif

c--- overall factor contained in diag.prc
          if (corr_on_beam == 2) then
              fac=aveqq*(as_heavy_beam2/2._dp/pi)*(V/two)*xn*two*(4*pi*as_heavy_beam2)*gwsq**4
          else
              fac=aveqq*(as_heavy_beam2/2._dp/pi)*(V/two)*xn*two*(4*pi*as_heavy_beam1)*gwsq**4
          endif

          ig=7
          ib=5
          i1=k1
          i6=k6

c-- calculate auxiliary momentum array - gq case
          twotDg=two*(+dot(p,1,ig)+dot(p,2,ig)+dot(p,6,ig))
          q(1:7,:)=p(1:7,:)
          q(8,:)=+p(1,:)+p(2,:)+p(6,:)-p(ig,:)*mt**2/twotDg

c---fill matrices of spinor products
          call spinoru(8,q,za,zb)

c--- Note: call to tree now passes top mass as a parameter
          call Wttree(mt,ig,ib,i6,i1,8,amplp)

c ig is the vector used to make the top quark massless in virt_pp etc.
          call tdecay_gen(p,ibeam2,i3,i4,ig,ampld_q_b)
          call tdecay_gen(p,ibeam2,i4,i3,ig,ampld_qb_b)
          ampld_q_b(:)=ampld_q_b(:)/(mt*twidth)
          ampld_qb_b(:)=ampld_qb_b(:)/(mt*twidth)

          r(1,:)=p(ig,:)
          r(2,:)=p(ib,:)
          r(3,:)=p(i6,:)
          r(4,:)=p(i1,:)
          r(5,:)=-r(1,:)-r(2,:)-r(3,:)-r(4,:)

          if (corr_on_beam == 2) then
              musq = renscale_beam2_isheavy_onheavy**2
          else
              musq = renscale_beam1_isheavy_onheavy**2
          endif

          s_ft(2,2)=virt_pp(mt,1,2,3,4,5,r,musq)
          s_ft(2,1)=virt_pm(mt,1,2,3,4,5,r,musq)
          s_ft(1,1)=virt_mm(mt,1,2,3,4,5,r,musq)
          s_ft(1,2)=virt_mp(mt,1,2,3,4,5,r,musq)

          wprop=sqrt((s(i1,i6)-wmass**2)**2+(wmass*wwidth)**2)

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark, for Born and virtual

          do hg=1,2
c--- lowest order amplitudes
              ampl0_q_b(hg)=amplp(1,hg)*ampld_q_b(1)+amplp(2,hg)*ampld_q_b(2)
              ampl0_qb_b(hg)=amplp(1,hg)*ampld_qb_b(1)+amplp(2,hg)*ampld_qb_b(2)
c--- virtual amplitudes
              amplv_q_b(hg)=(s_ft(1,hg)/wprop)*ampld_q_b(1)+(s_ft(2,hg)/wprop)*ampld_q_b(2)
              amplv_qb_b(hg)=(s_ft(1,hg)/wprop)*ampld_qb_b(1)+(s_ft(2,hg)/wprop)*ampld_qb_b(2)
          enddo

          msqbit_q_b=zip
          msqbit_qb_b=zip
          do hg=1,2
              msqbit_q_b=msqbit_q_b+real(amplv_q_b(hg)*conjg(ampl0_q_b(hg)))
              msqbit_qb_b=msqbit_qb_b+real(amplv_qb_b(hg)*conjg(ampl0_qb_b(hg)))
c            msqbit_q_b=msqbit_q_b+real(ampl0_q_b(hg)*conjg(ampl0_q_b(hg))) ! To check LO
c            msqbit_qb_b=msqbit_qb_b+real(ampl0_qb_b(hg)*conjg(ampl0_qb_b(hg))) ! To check LO
          enddo
          msqbit_q_b=msqbit_q_b*fac
          msqbit_qb_b=msqbit_qb_b*fac

          if (corr_on_beam == 2) then
              msq([2,4]*iq,5, 1,corr_on_beam) = msqbit_q_b
              msq([-1,-3]*iq,5, 1,corr_on_beam) = msqbit_qb_b
          else
              msq(5, [2,4]*iq, 1,corr_on_beam) = msqbit_q_b
              msq(5, [-1,-3]*iq, 1,corr_on_beam) = msqbit_qb_b
          endif
      enddo

      call singletop_jet_decay_intdip(p,msqintdip)
      msq = msq + msqintdip


      end subroutine singletop_jet_decay_virt_all


      subroutine singletop_jet_decay_virt(p,msq)
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
      include 'qcdcouple.f'
      include 'nwz.f'
      include 'scheme.f'
      integer:: ibeam1,ibeam2,i6,i1,i3,i4,iq,ig,ib,k1,k6,hg
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & q(mxpart,4),fac,r(mxpart,4),wprop,msqbit_q_b,msqbit_qb_b
      complex(dp):: ampl0_q_b(2),amplv_q_b(2),ampl0_qb_b(2),amplv_qb_b(2),
     & amplp(2,2),ampld_q_b(2),ampld_qb_b(2)
      real(dp):: twotDg,dot
      complex(dp):: s_ft(2,2),virt_pp,virt_pm,virt_mp,virt_mm

      scheme = 'tH-V'

      error stop "singletop_jet_decay_virt should not be used"

      msq(:,:)=zero

      if (corr_on_beam == 1) then
        ibeam1=2
        ibeam2=1
      elseif (corr_on_beam == 2) then
        ibeam1=1
        ibeam2=2
      endif

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
c antitop quark, so need crossed amplitudes for top quark
        k6=3
        k1=4
        i3=6
        i4=ibeam1
        iq=-1 ! antitop quark
      elseif (nwz == +1) then
c top quark, so need crossed amplitudes for antitop quark
        k6=4
        k1=3
        i3=ibeam1
        i4=6
        iq=1 ! top quark
      else
        write(6,*) 'Error in singletop_jet_decay_virt, nwz is not +1 or -1 :   ',nwz
        stop
      endif

c--- overall factor contained in diag.prc
      fac=aveqq*ason2pi*(V/two)*xn*two*gsq*gwsq**4

      ig=7
      ib=5
      i1=k1
      i6=k6

c-- calculate auxiliary momentum array - gq case
      twotDg=two*(+dot(p,1,ig)+dot(p,2,ig)+dot(p,6,ig))
      q(1:7,:)=p(1:7,:)
      q(8,:)=+p(1,:)+p(2,:)+p(6,:)-p(ig,:)*mt**2/twotDg

c---fill matrices of spinor products
      call spinoru(8,q,za,zb)

c--- Note: call to tree now passes top mass as a parameter
      call Wttree(mt,ig,ib,i6,i1,8,amplp)

c ig is the vector used to make the top quark massless in virt_pp etc.
      call tdecay_gen(p,ibeam2,i3,i4,ig,ampld_q_b)
      call tdecay_gen(p,ibeam2,i4,i3,ig,ampld_qb_b)
      ampld_q_b(:)=ampld_q_b(:)/(mt*twidth)
      ampld_qb_b(:)=ampld_qb_b(:)/(mt*twidth)

      r(1,:)=p(ig,:)
      r(2,:)=p(ib,:)
      r(3,:)=p(i6,:)
      r(4,:)=p(i1,:)
      r(5,:)=-r(1,:)-r(2,:)-r(3,:)-r(4,:)

      s_ft(2,2)=virt_pp(mt,1,2,3,4,5,r)
      s_ft(2,1)=virt_pm(mt,1,2,3,4,5,r)
      s_ft(1,1)=virt_mm(mt,1,2,3,4,5,r)
      s_ft(1,2)=virt_mp(mt,1,2,3,4,5,r)

      wprop=sqrt((s(i1,i6)-wmass**2)**2+(wmass*wwidth)**2)

c--- Construct factored form of amplitudes by adding the helicities of
c--- the heavy quark, for Born and virtual

      do hg=1,2
c--- lowest order amplitudes
      ampl0_q_b(hg)=amplp(1,hg)*ampld_q_b(1)+amplp(2,hg)*ampld_q_b(2)
      ampl0_qb_b(hg)=amplp(1,hg)*ampld_qb_b(1)+amplp(2,hg)*ampld_qb_b(2)
c--- virtual amplitudes
      amplv_q_b(hg)=(s_ft(1,hg)/wprop)*ampld_q_b(1)+(s_ft(2,hg)/wprop)*ampld_q_b(2)
      amplv_qb_b(hg)=(s_ft(1,hg)/wprop)*ampld_qb_b(1)+(s_ft(2,hg)/wprop)*ampld_qb_b(2)
      enddo

      msqbit_q_b=zip
      msqbit_qb_b=zip
      do hg=1,2
      msqbit_q_b=msqbit_q_b+real(amplv_q_b(hg)*conjg(ampl0_q_b(hg)))
      msqbit_qb_b=msqbit_qb_b+real(amplv_qb_b(hg)*conjg(ampl0_qb_b(hg)))
c      msqbit_q_b=msqbit_q_b+real(ampl0_q_b(hg)*conjg(ampl0_q_b(hg)))/ason2pi ! To check LO
c      msqbit_qb_b=msqbit_qb_b+real(ampl0_qb_b(hg)*conjg(ampl0_qb_b(hg)))/ason2pi ! To check LO
      enddo
      msqbit_q_b=msqbit_q_b*fac
      msqbit_qb_b=msqbit_qb_b*fac

      msq(2*iq,5) = msqbit_q_b
      msq(-1*iq,5) = msqbit_qb_b

      msq(4*iq,5)=msq(2*iq,5)
      msq(-3*iq,5)=msq(-1*iq,5)

      if (corr_on_beam == 1) then
        msq = transpose(msq)
      endif

      return
      end

      subroutine singletop_jet_decay_real_all(p,msqall)
          use types
          use singletop2_nnlo_vars
          !use NNTopDec, only: top24treebb, top24treegg, top24treeqq
          implicit none
c Heavily based on routine for t-channel single top + b production (qg_tbqdk_g)
          include 'constants.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'cplx.h'
          include 'masses.f'
          include 'ckm.f'
          include 'noglue.f'
          include 'nwz.f'
          include 'stopscales.f'
          include 'ewcouple.f'

          real(dp), intent(in) :: p(mxpart,4)
          real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

          real(dp):: msq_qb,msq_qbarb,msq_qb_qq,msq_qbarb_qq,msq_qb_bb,msq_qbarb_bb
          integer:: i1,ibeam1,ibeam2,i3,i4,i6,i7,i8,iq,m

          msqall = 0._dp

          do m=1,maxbeams
              corr_on_beam = beams_enabled(m)

              if (corr_on_beam == 1) then
                ibeam1=2
                ibeam2=1
              else
                ibeam1=1
                ibeam2=2
              endif

              if (nwz == -1) then
c antitop quark, so need crossed amplitudes for top quark
                  i3=3
                  i4=4
                  i1=ibeam1
                  i6=6
                  i7=7
                  i8=8
                  iq=-1 ! antitop quark
              elseif (nwz == +1) then
c top quark, so need crossed amplitudes for antitop quark
                  i3=4
                  i4=3
                  i1=6
                  i6=ibeam1
                  i7=8
                  i8=7
                  iq=1 ! top quark
              else
                  write(6,*) 'Error in singletop_jet_decay_real, nwz is not +1 or -1 :   ',nwz
                  stop
              endif

              msq_qb = 0._dp
              msq_qbarb = 0._dp
              msq_qb_qq = 0._dp
              msq_qbarb_qq = 0._dp
              msq_qb_bb = 0._dp
              msq_qbarb_bb = 0._dp

              call interdk_gen(p,i4,8,i6,i1,ibeam2,5,i3,7,mt,msq_qb)              ! d  b~ -> t~(e- nu~ b~ g g) u
              call interdk_gen(p,i4,8,i1,i6,ibeam2,5,i3,7,mt,msq_qbarb)           ! u~ b~ -> t~(e- nu~ b~ g g) d~

              call interdk_qq_gen(p,i4,i8,i6,i1,ibeam2,5,i3,i7,mt,msq_qb_qq)        ! d  b~ -> t~(e- nu~ b~ q q~) u
              call interdk_qq_gen(p,i4,i8,i1,i6,ibeam2,5,i3,i7,mt,msq_qbarb_qq)     ! u~ b~ -> t~(e- nu~ b~ q q~) d~

              call interdk_qqid_gen(p,i4,i8,i6,i1,ibeam2,5,i3,i7,mt,msq_qb_bb)      ! d  b~ -> t~(e- nu~ b~ b b~) u
              call interdk_qqid_gen(p,i4,i8,i1,i6,ibeam2,5,i3,i7,mt,msq_qbarb_bb)   ! u~ b~ -> t~(e- nu~ b~ b b~) d~

              if (corr_on_beam == 2) then
                  !call top24treebb(p,as_heavy_beam2,msq_qb_bb)
                  !msq_qb_bb = msq_qb_bb / half/aveqq

                  !call top24treeqq(p,as_heavy_beam2,msq_qb_qq)
                  !msq_qb_qq = msq_qb_qq / aveqq / (nf-1._dp)

                  msqall([2,4]*iq,5, 1,2) = aveqq*(half*msq_qb+(nf-1)*msq_qb_qq)
                  msqall([2,4]*iq,5, 4,2) = aveqq*half*msq_qb_bb
                  msqall([-1,-3]*iq,5, 1,2) = aveqq*(half*msq_qbarb+(nf-1)*msq_qbarb_qq)
                  msqall([-1,-3]*iq,5, 4,2) = aveqq*half*msq_qbarb_bb

                  !write (*,*) "msq qq~ = ", aveqq*(nf-1)*msq_qb_qq
                  !call top24treeqq(p,as_heavy_beam2,ubnew)
                  !write (*,*) "ratio qq~ ", ubnew / (aveqq*(nf-1)*msq_qb_qq)
                  !write (*,*) ""
                  !write (*,*) "msq bb~ = ", aveqq*half*msq_qb_bb
                  !call top24treebb(p,as_heavy_beam2,ubnew)
                  !write (*,*) "ratio bb~ ", ubnew / (aveqq*half*msq_qb_bb)
                  !write (*,*) ""
                  !write (*,*) "msq gg = ", aveqq*half*msq_qb
                  !call top24treegg(p,as_heavy_beam2,ubnew)
                  !write (*,*) "ratio gg ", ubnew / (aveqq*half*msq_qb)
                  !write (*,*) ""

                  !pause
              else
                  msqall(5, [2,4]*iq, 1,1) = aveqq*(half*msq_qb+(nf-1)*msq_qb_qq)
                  msqall(5, [2,4]*iq, 4,1) = aveqq*half*msq_qb_bb
                  msqall(5, [-1,-3]*iq, 1,1) = aveqq*(half*msq_qbarb+(nf-1)*msq_qbarb_qq)
                  msqall(5, [-1,-3]*iq, 4,1) = aveqq*half*msq_qbarb_bb
              endif
          enddo


      end subroutine singletop_jet_decay_real_all


      subroutine singletop_jet_decay_real(p,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'types.f'
c Heavily based on routine for t-channel single top + b production (qg_tbqdk_g)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'noglue.f'
      include 'nwz.f'
      include 'stopscales.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4),msq_qb,msq_qbarb,msq_qb_qq,msq_qbarb_qq,msq_qb_bb,msq_qbarb_bb
      real(dp):: msq(-nf:nf,-nf:nf)
      integer:: i1,ibeam1,ibeam2,i3,i4,i6,i7,i8,iq

      error stop "singletop_jet_decay_real should not be used"

c--- initialize
      msq(:,:)=0._dp

c--- initialize extra variables
      as_L=as; as_H=as

      if (corr_on_beam == 1) then
        ibeam1=2
        ibeam2=1
      else
        ibeam1=1
        ibeam2=2
      endif

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
c antitop quark, so need crossed amplitudes for top quark
        i3=3
        i4=4
        i1=ibeam1
        i6=6
        i7=7
        i8=8
        iq=-1 ! antitop quark
      elseif (nwz == +1) then
c top quark, so need crossed amplitudes for antitop quark
        i3=4
        i4=3
        i1=6
        i6=ibeam1
        i7=8
        i8=7
        iq=1 ! top quark
      else
        write(6,*) 'Error in singletop_jet_decay_real, nwz is not +1 or -1 :   ',nwz
        stop
      endif

      call interdk_gen(p,i4,8,i6,i1,ibeam2,5,i3,7,mt,msq_qb)              ! d  b~ -> t~(e- nu~ b~ g g) u
      call interdk_gen(p,i4,8,i1,i6,ibeam2,5,i3,7,mt,msq_qbarb)           ! u~ b~ -> t~(e- nu~ b~ g g) d~

      call interdk_qq_gen(p,i4,i8,i6,i1,ibeam2,5,i3,i7,mt,msq_qb_qq)        ! d  b~ -> t~(e- nu~ b~ q q~) u
      call interdk_qq_gen(p,i4,i8,i1,i6,ibeam2,5,i3,i7,mt,msq_qbarb_qq)     ! u~ b~ -> t~(e- nu~ b~ q q~) d~

      call interdk_qqid_gen(p,i4,i8,i6,i1,ibeam2,5,i3,i7,mt,msq_qb_bb)      ! d  b~ -> t~(e- nu~ b~ b b~) u
      call interdk_qqid_gen(p,i4,i8,i1,i6,ibeam2,5,i3,i7,mt,msq_qbarb_bb)   ! u~ b~ -> t~(e- nu~ b~ b b~) d~

c      write(6,*) 'msq_qb   ',msq_qb
c      write(6,*) 'msq_qb_qq',msq_qb_qq
c      write(6,*) 'msq_qb_bb',msq_qb_bb

      msq(2*iq,5) = aveqq*(half*msq_qb+(nf-1)*msq_qb_qq+msq_qb_bb)
      msq(-1*iq,5) = aveqq*(half*msq_qbarb+(nf-1)*msq_qbarb_qq+msq_qbarb_bb)

      msq(4*iq,:) = msq(2*iq,:)
      msq(-3*iq,:) = msq(-1*iq,:)

      if (corr_on_beam == 1) then
        msq = transpose(msq)
      endif

      return
      end




      subroutine tdecay_gen(p,jb,jn,je,jeta,amp)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: J.M.Campbell, January 2020                               *
c                                                                      *
c     Amplitudes for the process t -> b(jb) + nu (jn) + e+ (je)        *
c                                                                      *
c     where the bottom quark is massless and the top quark is made     *
c     massless with respect to vector jeta                             *
c                                                                      *
c     Taken from arXiv: 1204.1513, Section 4                           *
c     with W propagator added, only including a width when appropriate *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer jb,jn,je,jeta,kb,kn,ke,kt,keta
      real(dp):: p(mxpart,4),q(mxpart,4),dot,alpha,sw
      complex(dp):: amp(2),cprop
      parameter(kb=1,kn=2,ke=3,kt=4,keta=5)

      q(kb,:)=p(jb,:)
      q(kn,:)=p(jn,:)
      q(ke,:)=p(je,:)
      q(keta,:)=p(jeta,:)
      q(kt,:)=q(kb,:)+q(kn,:)+q(ke,:)
      alpha=mt**2/(2._dp*dot(q,kt,keta))
      q(kt,:)=q(kt,:)-alpha*q(keta,:)

      call spinoru(5,q,za,zb)

      sw=s(kn,ke)

      if (sw < zip) then
        cprop=cplx2(sw-wmass**2,zip)
      else
        cprop=cplx2(sw-wmass**2,wmass*wwidth)
      endif

      amp(1) = za(kb,kn)*zb(ke,kt)

      amp(2) = mt*za(kb,kn)*zb(ke,keta)/zb(kt,keta)

      amp(:)=amp(:)/cprop

      return
      end




      subroutine adecay_gen(p,jb,je,jn,jeta,amp)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: J.M.Campbell, January 2020                               *
c                                                                      *
c     Amplitudes for the process t~ -> b~(jb) + e- (je) + nu~ (jn)     *
c     with W propagator added, only including a width when appropriate *
c                                                                      *
c     where the bottom quark is massless and the top quark is made     *
c     massless with respect to vector jeta                             *
c                                                                      *
c     Taken from arXiv: 1204.1513, Section 4                           *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer jb,jn,je,jeta,kb,kn,ke,kt,keta
      real(dp):: p(mxpart,4),q(mxpart,4),dot,alpha,sw
      complex(dp):: amp(2),cprop
      parameter(kb=1,ke=2,kn=3,kt=4,keta=5)

      q(kb,:)=p(jb,:)
      q(kn,:)=p(jn,:)
      q(ke,:)=p(je,:)
      q(keta,:)=p(jeta,:)
      q(kt,:)=q(kb,:)+q(kn,:)+q(ke,:)
      alpha=mt**2/(2._dp*dot(q,kt,keta))
      q(kt,:)=q(kt,:)-alpha*q(keta,:)

      call spinoru(5,q,za,zb)

      sw=s(kn,ke)

      if (sw < zip) then
        cprop=cplx2(sw-wmass**2,zip)
      else
        cprop=cplx2(sw-wmass**2,wmass*wwidth)
      endif

      amp(1) = -mt*zb(kb,kn)*za(ke,keta)/za(kt,keta)

      amp(2) = -zb(kb,kn)*za(ke,kt)

      amp(:)=amp(:)/cprop

      return
      end
      end module


