!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wzg_gs(p,msq)
c***********************************************************************
c     Author: J.M.Campbell                                             *
c     April, 2020.                                                     *
c***********************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W Z parton(p7) + parton(p8) + g(p7)
c                           | |
c                           | ---> l(p5)+a(p6)
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
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq18_2(-nf:nf,-nf:nf),msq28_1(-nf:nf,-nf:nf),
     & msq17_8(-nf:nf,-nf:nf),msq28_7(-nf:nf,-nf:nf),
     & msq18_7(-nf:nf,-nf:nf),msq27_8(-nf:nf,-nf:nf),
     & msq78_1v(-nf:nf,-nf:nf),msq78_2v(-nf:nf,-nf:nf),
     & msq28_7v(-nf:nf,-nf:nf),msq28_1v(-nf:nf,-nf:nf),
     & msq17_8v(-nf:nf,-nf:nf),msq18_2v(-nf:nf,-nf:nf),
     & msq18_7v(-nf:nf,-nf:nf),msq27_8v(-nf:nf,-nf:nf),
     & msq17_2v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),sub18_2(4),sub28_1(4),
     & sub17_8(4),sub18_7(4),sub27_8(4),sub28_7(4),
     & sub78_1(4),sub78_2(4),sub78_1v,sub78_2v,
     & sub28_7v,sub28_1v,sub18_7v,sub18_2v,sub17_2v,sub17_8v,sub27_8v,
     & sub27_1v
      external qqb_wzg,qqb_wzg_gvec

      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(2,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(3,p,1,8,2,sub18_2,sub18_2v,msq18_2,msq18_2v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(4,p,2,8,1,sub28_1,sub28_1v,msq28_1,msq28_1v,
     & qqb_wzg,qqb_wzg_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,7,8,sub17_8,sub17_8v,msq17_8,msq17_8v,
     & qqb_wzg,qqb_wzg_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,7,8,1,sub78_1,sub78_1v,dummy,msq78_1v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(5,p,1,8,7,sub18_7,sub18_7v,msq18_7,msq18_7v,
     & qqb_wzg,qqb_wzg_gvec)

      call dips(6,p,2,8,7,sub28_7,sub28_7v,msq28_7,msq28_7v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(6,p,7,8,2,sub78_2,sub78_2v,dummy,msq78_2v,
     & qqb_wzg,qqb_wzg_gvec)
      call dips(6,p,2,7,8,sub27_8,sub27_8v,msq27_8,msq27_8v,
     & qqb_wzg,qqb_wzg_gvec)

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
      msq(1,j,k)=-msq17_2(j,k)*sub17_2(qq)/xn
      msq(2,j,k)=-msq27_1(j,k)*sub27_1(qq)/xn
      msq(3,j,k)=-msq18_2(j,k)*sub18_2(qq)/xn
      msq(4,j,k)=-msq28_1(j,k)*sub28_1(qq)/xn
      msq(5,j,k)=xn*(
     &  +msq17_8(j,k)*(sub17_8(qq)+0.5_dp*sub78_1(gg))
     &  +0.5_dp*msq78_1v(j,k)*sub78_1v
     &  +msq18_7(j,k)*(sub18_7(qq)+0.5_dp*sub78_1(gg))
     &  +0.5_dp*msq78_1v(j,k)*sub78_1v)
      msq(6,j,k)=xn*(
     &  (msq28_7(j,k)*(sub28_7(qq)+0.5_dp*sub78_2(gg))
     &   +0.5_dp*msq78_2v(j,k)*sub78_2v)
     & +(msq27_8(j,k)*(sub27_8(qq)+0.5_dp*sub78_2(gg))
     &   +0.5_dp*msq78_2v(j,k)*sub78_2v))

c--- note statistical factor of one half for two gluons in the final state
      do nd=1,ndmax
        msq(nd,j,k)=half*msq(nd,j,k)
      enddo

      elseif ((k  ==  0).and. (j  /=  0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2d0*tr*(msq27_1(j,-4)+msq27_1(j,-3)
     &                  +msq27_1(j,-2)+msq27_1(j,-1)+msq27_1(j,+1)
     &                  +msq27_1(j,+2)+msq27_1(j,+3)+msq27_1(j,+4)
     &                  )*sub27_1(qg)
      msq(3,j,k)=xn*msq18_2(j,k)*sub18_2(qq)
      msq(4,j,k)=xn*(msq28_1(j,k)*sub28_1(gg)+msq28_1v(j,k)*sub28_1v)
      msq(5,j,k)=-msq18_7(j,k)*(sub18_7(qq)+sub78_1(qq))/xn
      msq(6,j,k)=xn*(msq28_7(j,k)*sub28_7(gg)+msq28_7v(j,k)*sub28_7v
     &              +msq28_7(j,k)*sub78_2(qq))

      elseif ((j  ==  0).and.(k /= 0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2d0*tr*(msq17_2(-4,k)+msq17_2(-3,k)
     &                  +msq17_2(-2,k)+msq17_2(-1,k)+msq17_2(+1,k)
     &                  +msq17_2(+2,k)+msq17_2(+3,k)+msq17_2(+4,k)
     &                  )*sub17_2(qg)
      msq(3,j,k)=xn*(msq18_2(j,k)*sub18_2(gg)+msq18_2v(j,k)*sub18_2v)
      msq(4,j,k)=xn*msq28_1(j,k)*sub28_1(qq)
      msq(5,j,k)=xn*(msq18_7(j,k)*sub18_7(gg)+msq18_7v(j,k)*sub18_7v
     &              +msq18_7(j,k)*sub78_1(qq))
      msq(6,j,k)=-msq28_7(j,k)*(sub28_7(qq)+sub78_2(qq))/xn

      elseif ((j  ==  0).and.(k  ==  0)) then
c--- g-g case
c--- note g,g = 1,2 and q=7, qb=8 so (17),(27)-->qb and (18),(28)-->q
      msq(1,j,k)=(msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)
     &           +msq17_2(-4,k))*sub17_2(qg)*2d0*tr
      msq(2,j,k)=(msq27_1(k,-1)+msq27_1(k,-2)+msq27_1(k,-3)
     &           +msq27_1(k,-4))*sub27_1(qg)*2d0*tr
      msq(3,j,k)=(msq18_2(+4,k)+msq18_2(+3,k)
     &           +msq18_2(+2,k)+msq18_2(+1,k))*sub18_2(qg)*2d0*tr
      msq(4,j,k)=(msq28_1(k,+4)+msq28_1(k,+3)
     &           +msq28_1(k,+2)+msq28_1(k,+1))*sub28_1(qg)*2d0*tr

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

c (17) collinear singularities
      msq(1,i2,i1)=msq(1,i2,i1)+half*(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(1,i1,i3)=msq(1,i1,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)
      msq(1,i2,i3)=msq(1,i2,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)

      msq(1,-i1,-i2)=msq(1,-i1,-i2)+half*(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(1,-i2,-i4)=msq(1,-i2,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)
      msq(1,-i1,-i4)=msq(1,-i1,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)

      msq(1,i1,-i2)=msq(1,i1,-i2)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(1,i1,-i4)=msq(1,i1,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)
      msq(1,i2,-i2)=msq(1,i2,-i2)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(1,i2,-i4)=msq(1,i2,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)

      msq(1,-i2,i1)=msq(1,-i2,i1)+(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(1,-i2,i3)=msq(1,-i2,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)
      msq(1,-i1,i1)=msq(1,-i1,i1)+(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(1,-i1,i3)=msq(1,-i1,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)

      if (nflav == 5) then
      msq(1,5,i1)=msq(1,5,i1)+(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(1,5,i3)=msq(1,5,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)
      msq(1,5,-i2)=msq(1,5,-i2)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(1,5,-i4)=msq(1,5,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)
      msq(1,-5,i1)=msq(1,-5,i1)+(xn-1._dp/xn)
     &    *(msq17_2(0,i1)*sub17_2(gq)+msq17_2v(0,i1)*sub17_2v)
      msq(1,-5,i3)=msq(1,-5,i3)+(xn-1._dp/xn)
     &    *(msq17_2(0,i3)*sub17_2(gq)+msq17_2v(0,i3)*sub17_2v)
      msq(1,-5,-i2)=msq(1,-5,-i2)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i2)*sub17_2(gq)+msq17_2v(0,-i2)*sub17_2v)
      msq(1,-5,-i4)=msq(1,-5,-i4)+(xn-1._dp/xn)
     &    *(msq17_2(0,-i4)*sub17_2(gq)+msq17_2v(0,-i4)*sub17_2v)
      endif

c (28) collinear singularities
      msq(4,i1,i1)=msq(4,i1,i1)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,i2)=msq(4,i1,i2)+half*(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,i3)=msq(4,i1,i3)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,i4)=msq(4,i1,i4)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)

      msq(4,-i2,-i2)=msq(4,-i2,-i2)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,-i1)=msq(4,-i2,-i1)+half*(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,-i4)=msq(4,-i2,-i4)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,-i3)=msq(4,-i2,-i3)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)

      msq(4,i1,-i4)=msq(4,i1,-i4)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,-i3)=msq(4,i1,-i3)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,-i2)=msq(4,i1,-i2)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i1,-i1)=msq(4,i1,-i1)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)

      msq(4,-i2,i4)=msq(4,-i2,i4)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,i3)=msq(4,-i2,i3)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,i2)=msq(4,-i2,i2)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i2,i1)=msq(4,-i2,i1)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)

      if (nflav == 5) then
      msq(4,i1,5)=msq(4,i1,5)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i3,5)=msq(4,i3,5)+(xn-1._dp/xn)
     &    *(msq28_1(i3,0)*sub28_1(gq)+msq28_1v(i3,0)*sub28_1v)
      msq(4,-i2,5)=msq(4,-i2,5)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i4,5)=msq(4,-i4,5)+(xn-1._dp/xn)
     &    *(msq28_1(-i4,0)*sub28_1(gq)+msq28_1v(-i4,0)*sub28_1v)
      msq(4,i1,-5)=msq(4,i1,-5)+(xn-1._dp/xn)
     &    *(msq28_1(i1,0)*sub28_1(gq)+msq28_1v(i1,0)*sub28_1v)
      msq(4,i3,-5)=msq(4,i3,-5)+(xn-1._dp/xn)
     &    *(msq28_1(i3,0)*sub28_1(gq)+msq28_1v(i3,0)*sub28_1v)
      msq(4,-i2,-5)=msq(4,-i2,-5)+(xn-1._dp/xn)
     &    *(msq28_1(-i2,0)*sub28_1(gq)+msq28_1v(-i2,0)*sub28_1v)
      msq(4,-i4,-5)=msq(4,-i4,-5)+(xn-1._dp/xn)
     &    *(msq28_1(-i4,0)*sub28_1(gq)+msq28_1v(-i4,0)*sub28_1v)
      endif

c (18) collinear singularities
      msq(3,i1,i1)=msq(3,i1,i1)+(xn-1._dp/xn)
     &    *(msq18_2(0,i1)*sub18_2(gq)+msq18_2v(0,i1)*sub18_2v)
      msq(3,i2,i1)=msq(3,i2,i1)+half*(xn-1._dp/xn)
     &    *(msq18_2(0,i1)*sub18_2(gq)+msq18_2v(0,i1)*sub18_2v)

      msq(3,-i2,-i2)=msq(3,-i2,-i2)+(xn-1._dp/xn)
     &    *(msq18_2(0,-i2)*sub18_2(gq)+msq18_2v(0,-i2)*sub18_2v)
      msq(3,-i1,-i2)=msq(3,-i1,-i2)+half*(xn-1._dp/xn)
     &    *(msq18_2(0,-i2)*sub18_2(gq)+msq18_2v(0,-i2)*sub18_2v)

c (27) collinear singularities
      msq(2,i1,i2)=msq(2,i1,i2)+half*(xn-1._dp/xn)
     &     *(msq27_1(i1,0)*sub27_1(gq)+msq27_1v(i1,0)*sub27_1v)

      msq(2,-i2,-i1)=msq(2,-i2,-i1)+half*(xn-1._dp/xn)
     &     *(msq27_1(-i2,0)*sub27_1(gq)+msq27_1v(-i2,0)*sub27_1v)

c (78) collinear singularities
      msq(5,i1,-i2)=msq(5,i1,-i2)+tr*nflav
     &     *(msq18_7(i1,-i2)*sub78_1(gq)-msq78_1v(i1,-i2)*sub78_1v)
      msq(6,i1,-i2)=msq(6,i1,-i2)+tr*nflav
     &     *(msq28_7(i1,-i2)*sub78_2(gq)-msq78_2v(i1,-i2)*sub78_2v)

      msq(5,-i2,i1)=msq(5,-i2,i1)+tr*nflav
     &     *(msq18_7(-i2,i1)*sub78_1(gq)-msq78_1v(-i2,i1)*sub78_1v)
      msq(6,-i2,i1)=msq(6,-i2,i1)+tr*nflav
     &     *(msq28_7(-i2,i1)*sub78_2(gq)-msq78_2v(-i2,i1)*sub78_2v)

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

      return
      end

