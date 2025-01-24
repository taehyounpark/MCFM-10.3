!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_w_ew_gs(p,msq)
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
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),sub15_2(4),sub25_1(4),
     & msq15_4(-nf:nf,-nf:nf),msq25_4(-nf:nf,-nf:nf),sub15_4(4),sub25_4(4),
     & msq45_2(-nf:nf,-nf:nf),msq45_1(-nf:nf,-nf:nf),sub45_2(4),sub45_1(4)
      external qqb_w,donothing_gvec

      msq(:,:,:)=zip

      ndmax=4

c position of charged lepton in W decay
      if (nwz == +1) then
        i4=4
      else
        i4=3
      endif

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies

      call dips(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_2(:)=sub15_2(:)/gsq*abs(zesq)
      call dips(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_1(:)=sub25_1(:)/gsq*abs(zesq)

      call dips(3,p,i4,5,1,sub45_1,dsubv,msq45_1,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_1(:)=sub45_1(:)/gsq*abs(zesq)
      call dips(3,p,1,5,i4,sub15_4,dsubv,msq15_4,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_4(:)=sub15_4(:)/gsq*abs(zesq)

      call dips(4,p,i4,5,2,sub45_2,dsubv,msq45_2,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_2(:)=sub45_2(:)/gsq*abs(zesq)
      call dips(4,p,2,5,i4,sub25_4,dsubv,msq25_4,dummyv,qqb_w,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_4(:)=sub25_4(:)/gsq*abs(zesq)

      do j=-nf,nf
      do k=-nf,nf

      if ((j==0) .and. (k==0)) cycle

      if    ((j > 0) .and. (k < 0)
     &   .or.(j < 0) .and. (k > 0)) then
       msq(1,j,k)=Q(j)*Q(-k)*two*sub15_2(qq)*msq15_2(j,k)
       msq(2,j,k)=Q(j)*Q(-k)*two*sub25_1(qq)*msq25_1(j,k)
       msq(3,j,k)=-Q(+j)*Qe*nwz*two*(sub15_4(qq)+sub45_1(qq))*msq15_4(j,k)
       msq(4,j,k)=+Q(-k)*Qe*nwz*two*(sub25_4(qq)+sub45_2(qq))*msq25_4(j,k)

      endif

      enddo
      enddo

      return
      end


