!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgamg_gs(p,msq)
c***********************************************************************
c     Author: J.M.Campbell                                             *
c     May, 2020.                                                       *
c***********************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W gamma(p5) + parton(p6) + parton(p7)
c                           |
c                            -->l(p3)+a(p4)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'nwz.f'
      include 'nflav.f'
      integer j,k,nd,i1,i2,i3,i4
c --- remember: nd will count the dipoles

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq16_7(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),msq26_7(-nf:nf,-nf:nf),
     & msq67_1v(-nf:nf,-nf:nf),msq67_2v(-nf:nf,-nf:nf),
     & msq27_6v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & msq16_7v(-nf:nf,-nf:nf),msq17_2v(-nf:nf,-nf:nf),
     & msq17_6v(-nf:nf,-nf:nf),msq26_7v(-nf:nf,-nf:nf),
     & msq16_2v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),sub17_2(4),sub27_1(4),
     & sub16_7(4),sub17_6(4),sub26_7(4),sub27_6(4),
     & sub67_1(4),sub67_2(4),sub67_1v,sub67_2v,
     & sub27_6v,sub27_1v,sub17_6v,sub17_2v,sub16_2v,sub16_7v,sub26_7v,
     & sub26_1v
      external qqb_wgamg,qqb_wgam_gvec

      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(2,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(3,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(4,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & qqb_wgamg,qqb_wgam_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,6,7,sub16_7,sub16_7v,msq16_7,msq16_7v,
     & qqb_wgamg,qqb_wgam_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,6,7,1,sub67_1,sub67_1v,dummy,msq67_1v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(5,p,1,7,6,sub17_6,sub17_6v,msq17_6,msq17_6v,
     & qqb_wgamg,qqb_wgam_gvec)

      call dips(6,p,2,7,6,sub27_6,sub27_6v,msq27_6,msq27_6v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(6,p,6,7,2,sub67_2,sub67_2v,dummy,msq67_2v,
     & qqb_wgamg,qqb_wgam_gvec)
      call dips(6,p,2,6,7,sub26_7,sub26_7v,msq26_7,msq26_7v,
     & qqb_wgamg,qqb_wgam_gvec)

      msq(:,:,:)=0._dp

c--- DEBUG: 4-quark only for now
c      goto 33

c--- subtractions for 2-quark pieces
      do j=-2,2
      if (((j==+1).or.(j==-2)) .and. (nwz == +1)) cycle
      if (((j==-1).or.(j==+2)) .and. (nwz == -1)) cycle
      do k=-2,2
      if (((k==+1).or.(k==-2)) .and. (nwz == +1)) cycle
      if (((k==-1).or.(k==+2)) .and. (nwz == -1)) cycle
      if ((j == k) .and. (abs(j) > 0)) cycle

c--- do only q-qb and qb-q cases
      if (  ((j  >  0).and.(k  <  0))
     & .or. ((j  <  0).and.(k  >  0))) then
      msq(1,j,k)=-msq16_2(j,k)*sub16_2(qq)/xn
      msq(2,j,k)=-msq26_1(j,k)*sub26_1(qq)/xn
      msq(3,j,k)=-msq17_2(j,k)*sub17_2(qq)/xn
      msq(4,j,k)=-msq27_1(j,k)*sub27_1(qq)/xn
      msq(5,j,k)=xn*(
     &  +msq16_7(j,k)*(sub16_7(qq)+0.5_dp*sub67_1(gg))
     &  +0.5_dp*msq67_1v(j,k)*sub67_1v
     &  +msq17_6(j,k)*(sub17_6(qq)+0.5_dp*sub67_1(gg))
     &  +0.5_dp*msq67_1v(j,k)*sub67_1v)
      msq(6,j,k)=xn*(
     &  (msq27_6(j,k)*(sub27_6(qq)+0.5_dp*sub67_2(gg))
     &   +0.5_dp*msq67_2v(j,k)*sub67_2v)
     & +(msq26_7(j,k)*(sub26_7(qq)+0.5_dp*sub67_2(gg))
     &   +0.5_dp*msq67_2v(j,k)*sub67_2v))

c--- note statistical factor of one half for two gluons in the final state
      do nd=1,ndmax
        msq(nd,j,k)=half*msq(nd,j,k)
      enddo

      elseif ((k  ==  0).and. (j  /=  0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2d0*tr*(msq26_1(j,-4)+msq26_1(j,-3)
     &                  +msq26_1(j,-2)+msq26_1(j,-1)+msq26_1(j,+1)
     &                  +msq26_1(j,+2)+msq26_1(j,+3)+msq26_1(j,+4)
     &                  )*sub26_1(qg)
      msq(3,j,k)=xn*msq17_2(j,k)*sub17_2(qq)
      msq(4,j,k)=xn*(msq27_1(j,k)*sub27_1(gg)+msq27_1v(j,k)*sub27_1v)
      msq(5,j,k)=-msq17_6(j,k)*(sub17_6(qq)+sub67_1(qq))/xn
      msq(6,j,k)=xn*(msq27_6(j,k)*sub27_6(gg)+msq27_6v(j,k)*sub27_6v
     &              +msq27_6(j,k)*sub67_2(qq))

      elseif ((j  ==  0).and.(k /= 0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2d0*tr*(msq16_2(-4,k)+msq16_2(-3,k)
     &                  +msq16_2(-2,k)+msq16_2(-1,k)+msq16_2(+1,k)
     &                  +msq16_2(+2,k)+msq16_2(+3,k)+msq16_2(+4,k)
     &                  )*sub16_2(qg)
      msq(3,j,k)=xn*(msq17_2(j,k)*sub17_2(gg)+msq17_2v(j,k)*sub17_2v)
      msq(4,j,k)=xn*msq27_1(j,k)*sub27_1(qq)
      msq(5,j,k)=xn*(msq17_6(j,k)*sub17_6(gg)+msq17_6v(j,k)*sub17_6v
     &              +msq17_6(j,k)*sub67_1(qq))
      msq(6,j,k)=-msq27_6(j,k)*(sub27_6(qq)+sub67_2(qq))/xn

      elseif ((j  ==  0).and.(k  ==  0)) then
c--- g-g case
c--- note g,g = 1,2 and q=6, qb=7 so (16),(26)-->qb and (17),(27)-->q
      msq(1,j,k)=(msq16_2(-1,k)+msq16_2(-2,k)+msq16_2(-3,k)
     &           +msq16_2(-4,k))*sub16_2(qg)*2d0*tr
      msq(2,j,k)=(msq26_1(k,-1)+msq26_1(k,-2)+msq26_1(k,-3)
     &           +msq26_1(k,-4))*sub26_1(qg)*2d0*tr
      msq(3,j,k)=(msq17_2(+4,k)+msq17_2(+3,k)
     &           +msq17_2(+2,k)+msq17_2(+1,k))*sub17_2(qg)*2d0*tr
      msq(4,j,k)=(msq27_1(k,+4)+msq27_1(k,+3)
     &           +msq27_1(k,+2)+msq27_1(k,+1))*sub27_1(qg)*2d0*tr

      endif

      enddo
      enddo

      msq(:,+3,-4)=msq(:,+1,-2)
      msq(:,-4,+3)=msq(:,-2,+1)
      msq(:,-3,+4)=msq(:,-1,+2)
      msq(:,+4,-3)=msq(:,+2,-1)
      do j=3,4
      msq(:,j,0)=msq(:,(j-2),0)
      msq(:,0,j)=msq(:,0,(j-2))
      msq(:,-j,0)=msq(:,-(j-2),0)
      msq(:,0,-j)=msq(:,0,-(j-2))
      enddo

c--- DEBUG: 2-quark contribution only for now
c      return

c   33 continue

      if (nwz == -1) then
        i1=1
        i2=2
        i3=3
        i4=4
      else
        i1=2
        i2=1
        i3=4
        i4=3
      endif

c (16) collinear singularities
      msq(1,i2,i1)=msq(1,i2,i1)+half*(xn-1._dp/xn)
     &    *(msq16_2(0,i1)*sub16_2(gq)+msq16_2v(0,i1)*sub16_2v)
      msq(1,i1,i3)=msq(1,i1,i3)+(xn-1._dp/xn)
     &    *(msq16_2(0,i3)*sub16_2(gq)+msq16_2v(0,i3)*sub16_2v)
      msq(1,i2,i3)=msq(1,i2,i3)+(xn-1._dp/xn)
     &    *(msq16_2(0,i3)*sub16_2(gq)+msq16_2v(0,i3)*sub16_2v)

      msq(1,-i1,-i2)=msq(1,-i1,-i2)+half*(xn-1._dp/xn)
     &    *(msq16_2(0,-i2)*sub16_2(gq)+msq16_2v(0,-i2)*sub16_2v)
      msq(1,-i2,-i4)=msq(1,-i2,-i4)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i4)*sub16_2(gq)+msq16_2v(0,-i4)*sub16_2v)
      msq(1,-i1,-i4)=msq(1,-i1,-i4)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i4)*sub16_2(gq)+msq16_2v(0,-i4)*sub16_2v)

      msq(1,i1,-i2)=msq(1,i1,-i2)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i2)*sub16_2(gq)+msq16_2v(0,-i2)*sub16_2v)
      msq(1,i1,-i4)=msq(1,i1,-i4)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i4)*sub16_2(gq)+msq16_2v(0,-i4)*sub16_2v)
      msq(1,i2,-i2)=msq(1,i2,-i2)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i2)*sub16_2(gq)+msq16_2v(0,-i2)*sub16_2v)
      msq(1,i2,-i4)=msq(1,i2,-i4)+(xn-1._dp/xn)
     &    *(msq16_2(0,-i4)*sub16_2(gq)+msq16_2v(0,-i4)*sub16_2v)

      msq(1,-i2,i1)=msq(1,-i2,i1)+(xn-1._dp/xn)
     &    *(msq16_2(0,i1)*sub16_2(gq)+msq16_2v(0,i1)*sub16_2v)
      msq(1,-i2,i3)=msq(1,-i2,i3)+(xn-1._dp/xn)
     &    *(msq16_2(0,i3)*sub16_2(gq)+msq16_2v(0,i3)*sub16_2v)
      msq(1,-i1,i1)=msq(1,-i1,i1)+(xn-1._dp/xn)
     &    *(msq16_2(0,i1)*sub16_2(gq)+msq16_2v(0,i1)*sub16_2v)
      msq(1,-i1,i3)=msq(1,-i1,i3)+(xn-1._dp/xn)
     &    *(msq16_2(0,i3)*sub16_2(gq)+msq16_2v(0,i3)*sub16_2v)

c (27) collinear singularities
      msq(4,i1,i1)=msq(4,i1,i1)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,i2)=msq(4,i1,i2)+half*(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,i3)=msq(4,i1,i3)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,i4)=msq(4,i1,i4)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)

      msq(4,-i2,-i2)=msq(4,-i2,-i2)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,-i1)=msq(4,-i2,-i1)+half*(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,-i4)=msq(4,-i2,-i4)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,-i3)=msq(4,-i2,-i3)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)

      msq(4,i1,-i4)=msq(4,i1,-i4)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,-i3)=msq(4,i1,-i3)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,-i2)=msq(4,i1,-i2)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i1,-i1)=msq(4,i1,-i1)+(xn-1._dp/xn)
     &    *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)

      msq(4,-i2,i4)=msq(4,-i2,i4)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,i3)=msq(4,-i2,i3)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,i2)=msq(4,-i2,i2)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i2,i1)=msq(4,-i2,i1)+(xn-1._dp/xn)
     &    *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)

c (17) collinear singularities
      msq(3,i1,i1)=msq(3,i1,i1)+(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(3,i2,i1)=msq(3,i2,i1)+half*(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)

      msq(3,-i2,-i2)=msq(3,-i2,-i2)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(3,-i1,-i2)=msq(3,-i1,-i2)+half*(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)

c (26) collinear singularities
      msq(2,i1,i2)=msq(2,i1,i2)+half*(xn-1._dp/xn)
     &     *(msq26_1(i1,0)*sub26_1(gq)+msq26_1v(i1,0)*sub26_1v)

      msq(2,-i2,-i1)=msq(2,-i2,-i1)+half*(xn-1._dp/xn)
     &     *(msq26_1(-i2,0)*sub26_1(gq)+msq26_1v(-i2,0)*sub26_1v)

c (67) collinear singularities
      msq(5,i1,-i2)=msq(5,i1,-i2)+tr*nflav
     &     *(msq17_6(i1,-i2)*sub67_1(gq)-msq67_1v(i1,-i2)*sub67_1v)
      msq(6,i1,-i2)=msq(6,i1,-i2)+tr*nflav
     &     *(msq27_6(i1,-i2)*sub67_2(gq)-msq67_2v(i1,-i2)*sub67_2v)

      msq(5,-i2,i1)=msq(5,-i2,i1)+tr*nflav
     &     *(msq17_6(-i2,i1)*sub67_1(gq)-msq67_1v(-i2,i1)*sub67_1v)
      msq(6,-i2,i1)=msq(6,-i2,i1)+tr*nflav
     &     *(msq27_6(-i2,i1)*sub67_2(gq)-msq67_2v(-i2,i1)*sub67_2v)

c Relations to get other flavor combinations
      msq(:,i3,i1)=msq(:,i1,i3)
      msq(:,i3,i2)=msq(:,i1,i4)
      msq(:,i3,i3)=msq(:,i1,i1)
      msq(:,i3,i4)=msq(:,i1,i2)
      msq(:,i4,i1)=msq(:,i2,i3)
      msq(:,i4,i3)=msq(:,i2,i1)

      msq(:,i3,-i4)=msq(:,i1,-i2)
      msq(:,i3,-i3)=msq(:,i1,-i1)
      msq(:,i3,-i2)=msq(:,i1,-i4)
      msq(:,i3,-i1)=msq(:,i1,-i3)
      msq(:,i4,-i2)=msq(:,i2,-i4)
      msq(:,i4,-i4)=msq(:,i2,-i2)

      msq(:,-i4,i1)=msq(:,-i2,i3)
      msq(:,-i4,i2)=msq(:,-i2,i4)
      msq(:,-i4,i3)=msq(:,-i2,i1)
      msq(:,-i4,i4)=msq(:,-i2,i2)
      msq(:,-i3,i1)=msq(:,-i1,i3)
      msq(:,-i3,i3)=msq(:,-i1,i1)

      msq(:,-i4,-i4)=msq(:,-i2,-i2)
      msq(:,-i4,-i3)=msq(:,-i2,-i1)
      msq(:,-i4,-i2)=msq(:,-i2,-i4)
      msq(:,-i4,-i1)=msq(:,-i2,-i3)
      msq(:,-i3,-i4)=msq(:,-i1,-i2)
      msq(:,-i3,-i2)=msq(:,-i1,-i4)

      if (nflav == 5) then
      msq(4,i1,[5,-5])=(xn-1._dp/xn)*(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)
      msq(4,i3,[5,-5])=(xn-1._dp/xn)*(msq27_1(i3,0)*sub27_1(gq)+msq27_1v(i3,0)*sub27_1v)
      msq(4,-i2,[5,-5])=(xn-1._dp/xn)*(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)
      msq(4,-i4,[5,-5])=(xn-1._dp/xn)*(msq27_1(-i4,0)*sub27_1(gq)+msq27_1v(-i4,0)*sub27_1v)

      msq(1,[5,-5],i1)=(xn-1._dp/xn)*(msq16_2(0,i1)*sub16_2(gq)+msq16_2v(0,i1)*sub16_2v)
      msq(1,[5,-5],i3)=(xn-1._dp/xn)*(msq16_2(0,i3)*sub16_2(gq)+msq16_2v(0,i3)*sub16_2v)
      msq(1,[5,-5],-i2)=(xn-1._dp/xn)*(msq16_2(0,-i2)*sub16_2(gq)+msq16_2v(0,-i2)*sub16_2v)
      msq(1,[5,-5],-i4)=(xn-1._dp/xn)*(msq16_2(0,-i4)*sub16_2(gq)+msq16_2v(0,-i4)*sub16_2v)
      endif

      return
      end

