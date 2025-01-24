!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qphoton_wq_gs(p,msq)
      implicit none
      include 'types.f'
c---QED matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+photon(-p2) -->  l(p3)+a(p4)+gamma(p5)+q(p5)
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
      real(dp):: dummyv(-nf:nf,-nf:nf),dsubv,sub15_2(4),sub25_1(4),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf)
      external qqb_w,donothing_gvec

      msq(:,:,:)=zip

      ndmax=2

c position of charged lepton in W decay
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
      call dips(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_2(:)=sub15_2(:)/gsq*abs(zesq)
      call dips(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_1(:)=sub25_1(:)/gsq*abs(zesq)

      do j=-nf,nf
      do k=-nf,nf

      if ((j==0) .and. (k==0)) cycle

      if     ((j == 0) .and. (k  /=  0)) then
       msq(1,0,k)=two*xn*sub15_2(qg)*(
     &  Q(1)**2*(msq15_2(+1,k)+msq15_2(+3,k)+msq15_2(+5,k))
     & +Q(2)**2*(msq15_2(+2,k)+msq15_2(+4,k))
     & +Q(1)**2*(msq15_2(-1,k)+msq15_2(-3,k)+msq15_2(-5,k))
     & +Q(2)**2*(msq15_2(-2,k)+msq15_2(-4,k)))

      elseif ((j  /=  0) .and. (k == 0)) then
       msq(2,j,0)=two*xn*sub25_1(qg)*(
     &  Q(1)**2*(msq25_1(j,+1)+msq25_1(j,+3)+msq25_1(j,+5))
     & +Q(2)**2*(msq25_1(j,+2)+msq25_1(j,+4))
     & +Q(1)**2*(msq25_1(j,-1)+msq25_1(j,-3)+msq25_1(j,-5))
     & +Q(2)**2*(msq25_1(j,-2)+msq25_1(j,-4)))
      endif

      enddo
      enddo

      return
      end


