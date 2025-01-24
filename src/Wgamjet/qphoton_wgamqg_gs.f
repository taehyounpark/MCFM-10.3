!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qphoton_wgamqg_gs(p,msq)
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
      integer:: i4
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: dummyv(-nf:nf,-nf:nf),dsubv,
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),sub16_2(4),sub26_1(4),
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),sub17_2(4),sub27_1(4),
     & msq16_2qcd(-nf:nf,-nf:nf),msq26_1qcd(-nf:nf,-nf:nf),sub16_2qcd(4),sub26_1qcd(4),
     & msq17_2qcd(-nf:nf,-nf:nf),msq27_1qcd(-nf:nf,-nf:nf),sub17_2qcd(4),sub27_1qcd(4)
       external qqb_wgamg,qphoton_wgamq,donothing_gvec

      msq(:,:,:)=zip

      ndmax=4

c position of charged lepton in W decay
c not referenced?
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
      call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub16_2(:)=sub16_2(:)/gsq*abs(zesq)
      call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub26_1(:)=sub26_1(:)/gsq*abs(zesq)
      call dips(3,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub17_2(:)=sub17_2(:)/gsq*abs(zesq)
      call dips(4,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub27_1(:)=sub27_1(:)/gsq*abs(zesq)

      call dips(1,p,1,6,2,sub16_2qcd,dsubv,msq16_2qcd,dummyv,qphoton_wgamq,donothing_gvec)
      call dips(2,p,2,6,1,sub26_1qcd,dsubv,msq26_1qcd,dummyv,qphoton_wgamq,donothing_gvec)
      call dips(3,p,1,7,2,sub17_2qcd,dsubv,msq17_2qcd,dummyv,qphoton_wgamq,donothing_gvec)
      call dips(4,p,2,7,1,sub27_1qcd,dsubv,msq27_1qcd,dummyv,qphoton_wgamq,donothing_gvec)

c Photon-gluon subtractions
c Initial-state photon singularities
      msq(1,1,0) = two*xn*sub16_2(qg)*(
     &  Q(1)**2*(msq16_2(-1,0)+msq16_2(-3,0)+msq16_2(-5,0))
     & +Q(2)**2*(msq16_2(-2,0)+msq16_2(-4,0)))

      msq(3,1,0) = two*xn*sub17_2(qg)*(
     &  Q(1)**2*(msq17_2(+1,0)+msq17_2(+3,0)+msq17_2(+5,0))
     & +Q(2)**2*(msq17_2(+2,0)+msq17_2(+4,0)))

c Initial-state gluon singularities
      msq(2,1,0) = two*tr*sub26_1qcd(qg)*(
     &  msq26_1qcd(0,-1)+msq26_1qcd(0,-3)+msq26_1qcd(0,-5)
     & +msq26_1qcd(0,-2)+msq26_1qcd(0,-4))

      msq(4,1,0) = two*tr*sub27_1qcd(qg)*(
     &  msq27_1qcd(0,+1)+msq27_1qcd(0,+3)+msq27_1qcd(0,+5)
     & +msq27_1qcd(0,+2)+msq27_1qcd(0,+4))

c Gluon-photon subtractions
c Initial-state photon singularities
      msq(2,0,1) = two*xn*sub26_1(qg)*(
     &  Q(1)**2*(msq26_1(0,-1)+msq26_1(0,-3)+msq26_1(0,-5))
     & +Q(2)**2*(msq26_1(0,-2)+msq26_1(0,-4)))

      msq(4,0,1) = two*xn*sub27_1(qg)*(
     &  Q(1)**2*(msq27_1(0,+1)+msq27_1(0,+3)+msq27_1(0,+5))
     & +Q(2)**2*(msq27_1(0,+2)+msq27_1(0,+4)))

c Initial-state gluon singularities
      msq(1,0,1) = two*tr*sub16_2qcd(qg)*(
     &  msq16_2qcd(-1,0)+msq16_2qcd(-3,0)+msq16_2qcd(-5,0)
     & +msq16_2qcd(-2,0)+msq16_2qcd(-4,0))

      msq(3,0,1) = two*tr*sub17_2qcd(qg)*(
     &  msq17_2qcd(+1,0)+msq17_2qcd(+3,0)+msq17_2qcd(+5,0)
     & +msq17_2qcd(+2,0)+msq17_2qcd(+4,0))

c      do j=-nf,nf
c      do k=-nf,nf

c      if ((j==0) .and. (k==0)) cycle

c      if     ((j == 0) .and. (k  /=  0)) then
c!       msq(1,0,k)=xn*sub16_2(qg)*(
c!     &  Q(1)**2*(msq16_2(+1,k)+msq16_2(+3,k)+msq16_2(+5,k))
c!     & +Q(2)**2*(msq16_2(+2,k)+msq16_2(+4,k))
c!     & +Q(1)**2*(msq16_2(-1,k)+msq16_2(-3,k)+msq16_2(-5,k))
c!     & +Q(2)**2*(msq16_2(-2,k)+msq16_2(-4,k)))

c!       msq(3,0,k)=xn*sub17_2(qg)*(
c!     &  Q(1)**2*(msq17_2(+1,k)+msq17_2(+3,k)+msq17_2(+5,k))
c!     & +Q(2)**2*(msq17_2(+2,k)+msq17_2(+4,k))
c!     & +Q(1)**2*(msq17_2(-1,k)+msq17_2(-3,k)+msq17_2(-5,k))
c!     & +Q(2)**2*(msq17_2(-2,k)+msq17_2(-4,k)))

c      elseif ((j  /=  0) .and. (k == 0)) then
c!       msq(2,j,0)=xn*sub26_1(qg)*(
c!     &  Q(1)**2*(msq26_1(j,+1)+msq26_1(j,+3)+msq26_1(j,+5))
c!     & +Q(2)**2*(msq26_1(j,+2)+msq26_1(j,+4))
c!     & +Q(1)**2*(msq26_1(j,-1)+msq26_1(j,-3)+msq26_1(j,-5))
c!     & +Q(2)**2*(msq26_1(j,-2)+msq26_1(j,-4)))

c!       msq(4,j,0)=xn*sub27_1(qg)*(
c!     &  Q(1)**2*(msq27_1(j,+1)+msq27_1(j,+3)+msq27_1(j,+5))
c!     & +Q(2)**2*(msq27_1(j,+2)+msq27_1(j,+4))
c!     & +Q(1)**2*(msq27_1(j,-1)+msq27_1(j,-3)+msq27_1(j,-5))
c!     & +Q(2)**2*(msq27_1(j,-2)+msq27_1(j,-4)))
c      endif

c      enddo
c      enddo

      return
      end


