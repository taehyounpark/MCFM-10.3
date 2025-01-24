!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgam_gvec(p,n,in,msq)
      implicit none
c----Matrix element for Wgamma production
c----averaged over initial colours and spins
c    contracted with the vector v(mu)
c For nwz=-1
c  d(-p1)+ubar(-p2)--> W^-(e-(p3)+ve~(p4))+ Z(mu^-(p5)+mu^+(p6))+g(p7)
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'zcouple_cms.f'
      include 'pchoice.f'
      include 'nwz.f'
      integer in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),fac,fac1
      complex(dp):: d_ub(2),ub_d(2),d_g(2),g_d(2),ub_g(2),g_ub(2),amp(2)

      msq(:,:)=zip
      fac1=V/4._dp*gsq*abs((zesq/zxw)**2*zesq)

      call checkndotp(p,n,in)

      call spinoru(6,p,za,zb)

      if (nwz == -1) then
      if     (in  ==  6) then
        call a6Wgamn(1,2,3,4,5,6,p,n,za,zb,d_ub)
        call a6Wgamn(2,1,3,4,5,6,p,n,za,zb,ub_d)
      elseif (in  ==  1) then
        call a6Wgamn(6,2,3,4,5,1,p,n,za,zb,g_ub)
        call a6Wgamn(2,6,3,4,5,1,p,n,za,zb,g_d)
      elseif (in  ==  2) then
        call a6Wgamn(1,6,3,4,5,2,p,n,za,zb,d_g)
        call a6Wgamn(6,1,3,4,5,2,p,n,za,zb,ub_g)
      else
        write(6,*) 'Invalid value of in=',in,' in qqb_wzg_gvec.f!'
        stop
      endif
      elseif (nwz == +1) then
      if     (in  ==  6) then
        call a6Wgamn(2,1,4,3,5,6,p,n,zb,za,d_ub)
        call a6Wgamn(1,2,4,3,5,6,p,n,zb,za,ub_d)
      elseif (in  ==  1) then
        call a6Wgamn(2,6,4,3,5,1,p,n,zb,za,g_ub)
        call a6Wgamn(6,2,4,3,5,1,p,n,zb,za,g_d)
      elseif (in  ==  2) then
        call a6Wgamn(6,1,4,3,5,2,p,n,zb,za,d_g)
        call a6Wgamn(1,6,4,3,5,2,p,n,zb,za,ub_g)
      endif
      endif

      do j=-2,2
      if (((j==+1).or.(j==-2)) .and. (nwz == +1)) cycle
      if (((j==-1).or.(j==+2)) .and. (nwz == -1)) cycle
      do k=-2,2
      if (((k==+1).or.(k==-2)) .and. (nwz == +1)) cycle
      if (((k==-1).or.(k==+2)) .and. (nwz == -1)) cycle
      if (j == k)  cycle

      if     ((j > 0) .and. (k < 0)) then
        fac=aveqq*fac1
        amp(:)=d_ub(:)
      elseif ((j < 0) .and. (k > 0)) then
        fac=aveqq*fac1
        amp(:)=ub_d(:)
      elseif ((j > 0) .and. (k == 0)) then
        fac=aveqg*fac1
        amp(:)=d_g(:)
      elseif ((j < 0) .and. (k == 0)) then
        fac=aveqg*fac1
        amp(:)=ub_g(:)
      elseif ((j == 0) .and. (k > 0)) then
        fac=aveqg*fac1
        amp(:)=g_d(:)
      elseif ((j == 0) .and. (k < 0)) then
        fac=aveqg*fac1
        amp(:)=g_ub(:)
      endif

      msq(j,k)=fac*(cdabs(amp(1))**2+cdabs(amp(2))**2)

      enddo
      enddo

c Extend to other flavours
      do j=3,4
      msq(j,0)=msq((j-2),0)
      msq(0,j)=msq(0,(j-2))
      msq(-j,0)=msq(-(j-2),0)
      msq(0,-j)=msq(0,-(j-2))
      enddo
      msq(+3,-4)=msq(+1,-2)
      msq(-4,+3)=msq(-2,+1)
      msq(+4,-3)=msq(+2,-1)
      msq(-3,+4)=msq(-1,+2)

      return
      end

