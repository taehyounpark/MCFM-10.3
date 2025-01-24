!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wgamg_gs_ew(p,msq)
      implicit none
      include 'types.f'
c---QED matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+photon(-p2) -->  l(p3)+a(p4)+gamma(p5)+g(p6)+gamma(p7)
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
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),Qopp
      real(dp):: dummyv(-nf:nf,-nf:nf),dsubv,
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),sub17_2(4),sub27_1(4),
     & msq17_4(-nf:nf,-nf:nf),msq27_4(-nf:nf,-nf:nf),sub17_4(4),sub27_4(4),
     & msq47_2(-nf:nf,-nf:nf),msq47_1(-nf:nf,-nf:nf),sub47_2(4),sub47_1(4),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),sub15_2(4),sub25_1(4),
     & msq15_4(-nf:nf,-nf:nf),msq25_4(-nf:nf,-nf:nf),sub15_4(4),sub25_4(4),
     & msq45_2(-nf:nf,-nf:nf),msq45_1(-nf:nf,-nf:nf),sub45_2(4),sub45_1(4),
     & msq17_6(-nf:nf,-nf:nf),msq67_1(-nf:nf,-nf:nf),sub17_6(4),sub67_1(4),
     & msq27_6(-nf:nf,-nf:nf),msq67_2(-nf:nf,-nf:nf),sub27_6(4),sub67_2(4),
     & msq67_4(-nf:nf,-nf:nf),msq47_6(-nf:nf,-nf:nf),sub67_4(4),sub47_6(4),
     & msq15_6(-nf:nf,-nf:nf),msq65_1(-nf:nf,-nf:nf),sub15_6(4),sub65_1(4),
     & msq25_6(-nf:nf,-nf:nf),msq65_2(-nf:nf,-nf:nf),sub25_6(4),sub65_2(4),
     & msq65_4(-nf:nf,-nf:nf),msq45_6(-nf:nf,-nf:nf),sub65_4(4),sub45_6(4)
      external qqb_wgamg,qqb_wgamg_swap,donothing_gvec

      msq(:,:,:)=zip

      ndmax=16

c position of charged lepton in W decay
      if (nwz == +1) then
        i4=4
      else
        i4=3
      endif

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub17_2(:)=sub17_2(:)/gsq*abs(zesq)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub27_1(:)=sub27_1(:)/gsq*abs(zesq)

      call dips(3,p,i4,7,1,sub47_1,dsubv,msq47_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub47_1(:)=sub47_1(:)/gsq*abs(zesq)
      call dips(3,p,1,7,i4,sub17_4,dsubv,msq17_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub17_4(:)=sub17_4(:)/gsq*abs(zesq)

      call dips(4,p,i4,7,2,sub47_2,dsubv,msq47_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub47_2(:)=sub47_2(:)/gsq*abs(zesq)
      call dips(4,p,2,7,i4,sub27_4,dsubv,msq27_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub27_4(:)=sub27_4(:)/gsq*abs(zesq)

c--- now repeat with 7 <-> 5
      call dips(5,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_2(:)=sub15_2(:)/gsq*abs(zesq)
      call dips(6,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_1(:)=sub25_1(:)/gsq*abs(zesq)

      call dips(7,p,i4,5,1,sub45_1,dsubv,msq45_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_1(:)=sub45_1(:)/gsq*abs(zesq)
      call dips(7,p,1,5,i4,sub15_4,dsubv,msq15_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_4(:)=sub15_4(:)/gsq*abs(zesq)

      call dips(8,p,i4,5,2,sub45_2,dsubv,msq45_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_2(:)=sub45_2(:)/gsq*abs(zesq)
      call dips(8,p,2,5,i4,sub25_4,dsubv,msq25_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_4(:)=sub25_4(:)/gsq*abs(zesq)

c Extra for gq configurations
      call dips(9,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub17_6(:)=sub17_6(:)/gsq*abs(zesq)
      call dips(9,p,6,7,1,sub67_1,dsubv,msq67_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub67_1(:)=sub67_1(:)/gsq*abs(zesq)
      call dips(10,p,6,7,i4,sub67_4,dsubv,msq67_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub67_4(:)=sub67_4(:)/gsq*abs(zesq)
      call dips(11,p,i4,7,6,sub47_6,dsubv,msq47_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub47_6(:)=sub47_6(:)/gsq*abs(zesq)
      call dips(12,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub27_6(:)=sub27_6(:)/gsq*abs(zesq)
      call dips(12,p,6,7,2,sub67_2,dsubv,msq67_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub67_2(:)=sub67_2(:)/gsq*abs(zesq)

c now repeat extra ones with 7 <-> 5
      call dips(13,p,1,5,6,sub15_6,dsubv,msq15_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub15_6(:)=sub15_6(:)/gsq*abs(zesq)
      call dips(13,p,6,5,1,sub65_1,dsubv,msq65_1,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub65_1(:)=sub65_1(:)/gsq*abs(zesq)
      call dips(14,p,6,5,i4,sub65_4,dsubv,msq65_4,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub65_4(:)=sub65_4(:)/gsq*abs(zesq)
      call dips(15,p,i4,5,6,sub45_6,dsubv,msq45_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub45_6(:)=sub45_6(:)/gsq*abs(zesq)
      call dips(16,p,2,5,6,sub25_6,dsubv,msq25_6,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub25_6(:)=sub25_6(:)/gsq*abs(zesq)
      call dips(16,p,6,5,2,sub65_2,dsubv,msq65_2,dummyv,qqb_wgamg,donothing_gvec)
c Change from QCD splitting to QED one (use gsq now, as set in dipole routine)
      sub65_2(:)=sub65_2(:)/gsq*abs(zesq)


      do j=-nf,nf
      do k=-nf,nf

      if ((j==0) .and. (k==0)) cycle

      if    ((j > 0) .and. (k < 0)
     &   .or.(j < 0) .and. (k > 0)) then
       msq(1,j,k)=Q(j)*Q(-k)*sub17_2(qq)*msq17_2(j,k)
       msq(2,j,k)=Q(j)*Q(-k)*sub27_1(qq)*msq27_1(j,k)
       msq(3,j,k)=-Q(+j)*Qe*nwz*(sub17_4(qq)+sub47_1(qq))*msq17_4(j,k)
       msq(4,j,k)=+Q(-k)*Qe*nwz*(sub27_4(qq)+sub47_2(qq))*msq27_4(j,k)

       msq(5,j,k)=Q(j)*Q(-k)*sub15_2(qq)*msq15_2(j,k)
       msq(6,j,k)=Q(j)*Q(-k)*sub25_1(qq)*msq25_1(j,k)
       msq(7,j,k)=-Q(+j)*Qe*nwz*(sub15_4(qq)+sub45_1(qq))*msq15_4(j,k)
       msq(8,j,k)=+Q(-k)*Qe*nwz*(sub25_4(qq)+sub45_2(qq))*msq25_4(j,k)

c Need to add 7<->5 exchanges below
       elseif ((j > 0) .and. (k == 0)) then
       Qopp=one/three-Q(j)
       msq(9,j,k)=Q(j)*Qopp*(sub17_6(qq)+sub67_1(qq))*msq17_6(j,k)
       msq(3,j,k)=-Q(j)*Qe*nwz*(sub17_4(qq)+sub47_1(qq))*msq17_4(j,k)
       msq(10,j,k)=+Qopp*Qe*nwz*sub67_4(qq)*msq67_4(j,k)
       msq(11,j,k)=+Qopp*Qe*nwz*sub47_6(qq)*msq47_6(j,k)

       msq(13,j,k)=Q(j)*Qopp*(sub15_6(qq)+sub65_1(qq))*msq15_6(j,k)
       msq(7,j,k)=-Q(j)*Qe*nwz*(sub15_4(qq)+sub45_1(qq))*msq15_4(j,k)
       msq(14,j,k)=+Qopp*Qe*nwz*sub65_4(qq)*msq65_4(j,k)
       msq(15,j,k)=+Qopp*Qe*nwz*sub45_6(qq)*msq45_6(j,k)

       elseif ((j < 0) .and. (k == 0)) then
       Qopp=one/three-Q(-j)
       msq(9,j,k)=Qopp*Q(-j)*(sub17_6(qq)+sub67_1(qq))*msq17_6(j,k)
       msq(3,j,k)=+Q(-j)*Qe*nwz*(sub17_4(qq)+sub47_1(qq))*msq17_4(j,k)
       msq(10,j,k)=-Qopp*Qe*nwz*sub67_4(qq)*msq67_4(j,k)
       msq(11,j,k)=-Qopp*Qe*nwz*sub47_6(qq)*msq47_6(j,k)

       msq(13,j,k)=Qopp*Q(-j)*(sub15_6(qq)+sub65_1(qq))*msq15_6(j,k)
       msq(7,j,k)=+Q(-j)*Qe*nwz*(sub15_4(qq)+sub45_1(qq))*msq15_4(j,k)
       msq(14,j,k)=-Qopp*Qe*nwz*sub65_4(qq)*msq65_4(j,k)
       msq(15,j,k)=-Qopp*Qe*nwz*sub45_6(qq)*msq45_6(j,k)

       elseif ((j == 0) .and. (k > 0)) then
       Qopp=one/three-Q(k)
       msq(12,j,k)=Q(k)*Qopp*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
       msq(4,j,k)=-Q(k)*Qe*nwz*(sub27_4(qq)+sub47_2(qq))*msq27_4(j,k)
       msq(10,j,k)=+Qopp*Qe*nwz*sub67_4(qq)*msq67_4(j,k)
       msq(11,j,k)=+Qopp*Qe*nwz*sub47_6(qq)*msq47_6(j,k)

       msq(16,j,k)=Q(k)*Qopp*(sub25_6(qq)+sub65_2(qq))*msq25_6(j,k)
       msq(8,j,k)=-Q(k)*Qe*nwz*(sub25_4(qq)+sub45_2(qq))*msq25_4(j,k)
       msq(14,j,k)=+Qopp*Qe*nwz*sub65_4(qq)*msq65_4(j,k)
       msq(15,j,k)=+Qopp*Qe*nwz*sub45_6(qq)*msq45_6(j,k)

       elseif ((j == 0) .and. (k < 0)) then
       Qopp=one/three-Q(-k)
       msq(12,j,k)=Qopp*Q(-k)*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
       msq(4,j,k)=+Q(-k)*Qe*nwz*(sub27_4(qq)+sub47_2(qq))*msq27_4(j,k)
       msq(10,j,k)=-Qopp*Qe*nwz*sub67_4(qq)*msq67_4(j,k)
       msq(11,j,k)=-Qopp*Qe*nwz*sub47_6(qq)*msq47_6(j,k)

       msq(16,j,k)=Qopp*Q(-k)*(sub25_6(qq)+sub65_2(qq))*msq25_6(j,k)
       msq(8,j,k)=+Q(-k)*Qe*nwz*(sub25_4(qq)+sub45_2(qq))*msq25_4(j,k)
       msq(14,j,k)=-Qopp*Qe*nwz*sub65_4(qq)*msq65_4(j,k)
       msq(15,j,k)=-Qopp*Qe*nwz*sub45_6(qq)*msq45_6(j,k)

      endif

      enddo
      enddo

      return
      end
