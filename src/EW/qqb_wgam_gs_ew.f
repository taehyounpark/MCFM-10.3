!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wgam_gs_ew(p,msq)
      implicit none
      include 'types.f'
c---QED matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+photon(-p2) -->  l(p3)+a(p4)+gamma(p5)+gamma(p6)
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
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),sub16_2(4),sub26_1(4),
     & msq16_4(-nf:nf,-nf:nf),msq26_4(-nf:nf,-nf:nf),sub16_4(4),sub26_4(4),
     & msq46_2(-nf:nf,-nf:nf),msq46_1(-nf:nf,-nf:nf),sub46_2(4),sub46_1(4),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),sub15_2(4),sub25_1(4),
     & msq15_4(-nf:nf,-nf:nf),msq25_4(-nf:nf,-nf:nf),sub15_4(4),sub25_4(4),
     & msq45_2(-nf:nf,-nf:nf),msq45_1(-nf:nf,-nf:nf),sub45_2(4),sub45_1(4)
      external qqb_wgam,donothing_gvec

      msq(:,:,:)=zip

      ndmax=8

c position of charged lepton in W decay
      if (nwz == +1) then
        i4=4
      else
        i4=3
      endif

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub16_2(:)=sub16_2(:)/gsq*abs(zesq)
      call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub26_1(:)=sub26_1(:)/gsq*abs(zesq)

      call dips(3,p,i4,6,1,sub46_1,dsubv,msq46_1,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub46_1(:)=sub46_1(:)/gsq*abs(zesq)
      call dips(3,p,1,6,i4,sub16_4,dsubv,msq16_4,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub16_4(:)=sub16_4(:)/gsq*abs(zesq)

      call dips(4,p,i4,6,2,sub46_2,dsubv,msq46_2,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub46_2(:)=sub46_2(:)/gsq*abs(zesq)
      call dips(4,p,2,6,i4,sub26_4,dsubv,msq26_4,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub26_4(:)=sub26_4(:)/gsq*abs(zesq)

c--- now repeat with 6 <-> 5
      call dips(5,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_2(:)=sub15_2(:)/gsq*abs(zesq)
      call dips(6,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_1(:)=sub25_1(:)/gsq*abs(zesq)

      call dips(7,p,i4,5,1,sub45_1,dsubv,msq45_1,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_1(:)=sub45_1(:)/gsq*abs(zesq)
      call dips(7,p,1,5,i4,sub15_4,dsubv,msq15_4,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_4(:)=sub15_4(:)/gsq*abs(zesq)

      call dips(8,p,i4,5,2,sub45_2,dsubv,msq45_2,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_2(:)=sub45_2(:)/gsq*abs(zesq)
      call dips(8,p,2,5,i4,sub25_4,dsubv,msq25_4,dummyv,qqb_wgam,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_4(:)=sub25_4(:)/gsq*abs(zesq)

      do j=-nf,nf
      do k=-nf,nf

      if ((j==0) .and. (k==0)) cycle

      if    ((j > 0) .and. (k < 0)
     &   .or.(j < 0) .and. (k > 0)) then
       msq(1,j,k)=Q(j)*Q(-k)*sub16_2(qq)*msq16_2(j,k)
       msq(2,j,k)=Q(j)*Q(-k)*sub26_1(qq)*msq26_1(j,k)
       msq(3,j,k)=-Q(+j)*Qe*nwz*(sub16_4(qq)+sub46_1(qq))*msq16_4(j,k)
       msq(4,j,k)=+Q(-k)*Qe*nwz*(sub26_4(qq)+sub46_2(qq))*msq26_4(j,k)

       msq(5,j,k)=Q(j)*Q(-k)*sub15_2(qq)*msq15_2(j,k)
       msq(6,j,k)=Q(j)*Q(-k)*sub25_1(qq)*msq25_1(j,k)
       msq(7,j,k)=-Q(+j)*Qe*nwz*(sub15_4(qq)+sub45_1(qq))*msq15_4(j,k)
       msq(8,j,k)=+Q(-k)*Qe*nwz*(sub25_4(qq)+sub45_2(qq))*msq25_4(j,k)

c These are included as part of the Wgaj_a process

c      elseif ((j == 0) .and. (k  /=  0)) then
c       msq(1,0,k)=two*xn*sub16_2(qg)*(
c     &  Q(1)**2*(msq16_2(+1,k)+msq16_2(+3,k)+msq16_2(+5,k))
c     & +Q(2)**2*(msq16_2(+2,k)+msq16_2(+4,k))
c     & +Q(1)**2*(msq16_2(-1,k)+msq16_2(-3,k)+msq16_2(-5,k))
c     & +Q(2)**2*(msq16_2(-2,k)+msq16_2(-4,k)))

c      elseif ((j  /=  0) .and. (k == 0)) then
c       msq(2,j,0)=two*xn*sub26_1(qg)*(
c     &  Q(1)**2*(msq26_1(j,+1)+msq26_1(j,+3)+msq26_1(j,+5))
c     & +Q(2)**2*(msq26_1(j,+2)+msq26_1(j,+4))
c     & +Q(1)**2*(msq26_1(j,-1)+msq26_1(j,-3)+msq26_1(j,-5))
c     & +Q(2)**2*(msq26_1(j,-2)+msq26_1(j,-4)))
      endif

      enddo
      enddo

      return
      end


