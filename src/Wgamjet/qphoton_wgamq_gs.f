!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qphoton_wgamq_gs(p,msq)
      implicit none
      include 'types.f'
c---QED matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+photon(-p2) -->  l(p3)+a(p4)+gamma(p5)+q(p6)
c   positively charged W only

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'qcdcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'nwz.f'
      integer:: j,k,i4
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: dummyv(-nf:nf,-nf:nf),dsubv,
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),sub16_2(4),sub26_1(4)
      external qqb_wgam,donothing_gvec

      msq(:,:,:)=zip

      ndmax=2

c position of charged lepton in W decay
c not referenced ?
      if (nwz == +1) then
        i4=4
      else
        i4=3
      endif

c Note that for this process we always demand that particle 5 be a hard, isolated photon.
c As a result, no singularities associated with 5 being soft or collinear are included
c but we include a factor of two to count this possibility

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub16_2(:)=sub16_2(:)/gsq*abs(zesq)
      call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub26_1(:)=sub26_1(:)/gsq*abs(zesq)

      do j=-nf,nf
      do k=-nf,nf

      if ((j==0) .and. (k==0)) cycle

      if     ((j == 0) .and. (k  /=  0)) then
       msq(1,0,k)=two*xn*sub16_2(qg)*(
     &  Q(1)**2*(msq16_2(+1,k)+msq16_2(+3,k)+msq16_2(+5,k))
     & +Q(2)**2*(msq16_2(+2,k)+msq16_2(+4,k))
     & +Q(1)**2*(msq16_2(-1,k)+msq16_2(-3,k)+msq16_2(-5,k))
     & +Q(2)**2*(msq16_2(-2,k)+msq16_2(-4,k)))

      elseif ((j  /=  0) .and. (k == 0)) then
       msq(2,j,0)=two*xn*sub26_1(qg)*(
     &  Q(1)**2*(msq26_1(j,+1)+msq26_1(j,+3)+msq26_1(j,+5))
     & +Q(2)**2*(msq26_1(j,+2)+msq26_1(j,+4))
     & +Q(1)**2*(msq26_1(j,-1)+msq26_1(j,-3)+msq26_1(j,-5))
     & +Q(2)**2*(msq26_1(j,-2)+msq26_1(j,-4)))
      endif

      enddo
      enddo

      return
      end


