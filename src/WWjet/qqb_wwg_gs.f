!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wwg_gs(p,msq)
c***********************************************************************
c     Author: R.K.Ellis                                                *
c     August, 2008.                                                    *
c***********************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W^+ W^-+ parton(p7) + parton(p8) + g(p7)
c                           |   |
c                           |   ---> l(p5)+a(p6)
c                           |
c                            -->l(p3)+a(p4)
c   positively charged W only

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd
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
      external qqb_wwg,qqb_wwg_gvec

      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(2,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(3,p,1,8,2,sub18_2,sub18_2v,msq18_2,msq18_2v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(4,p,2,8,1,sub28_1,sub28_1v,msq28_1,msq28_1v,
     & qqb_wwg,qqb_wwg_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,7,8,sub17_8,sub17_8v,msq17_8,msq17_8v,
     & qqb_wwg,qqb_wwg_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,7,8,1,sub78_1,sub78_1v,dummy,msq78_1v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(5,p,1,8,7,sub18_7,sub18_7v,msq18_7,msq18_7v,
     & qqb_wwg,qqb_wwg_gvec)

      call dips(6,p,2,8,7,sub28_7,sub28_7v,msq28_7,msq28_7v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(6,p,7,8,2,sub78_2,sub78_2v,dummy,msq78_2v,
     & qqb_wwg,qqb_wwg_gvec)
      call dips(6,p,2,7,8,sub27_8,sub27_8v,msq27_8,msq27_8v,
     & qqb_wwg,qqb_wwg_gvec)

      msq(:,:,:)=0d0

c--- DEBUG: 4-quark only for now
c      goto 33

c--- subtractions for 2-quark pieces
      do j=-4,4
      do k=-4,4
c--- do only q-qb and qb-q cases
      if (  ((j  >  0).and.(k  <  0))
     & .or. ((j  <  0).and.(k  >  0))) then
      msq(1,j,k)=-msq17_2(j,k)*sub17_2(qq)/xn
      msq(2,j,k)=-msq27_1(j,k)*sub27_1(qq)/xn
      msq(3,j,k)=-msq18_2(j,k)*sub18_2(qq)/xn
      msq(4,j,k)=-msq28_1(j,k)*sub28_1(qq)/xn
      msq(5,j,k)=xn*(
     &  +msq17_8(j,k)*(sub17_8(qq)+0.5d0*sub78_1(gg))
     &  +0.5d0*msq78_1v(j,k)*sub78_1v
     &  +msq18_7(j,k)*(sub18_7(qq)+0.5d0*sub78_1(gg))
     &  +0.5d0*msq78_1v(j,k)*sub78_1v)
      msq(6,j,k)=xn*(
     &  (msq28_7(j,k)*(sub28_7(qq)+0.5d0*sub78_2(gg))
     &   +0.5d0*msq78_2v(j,k)*sub78_2v)
     & +(msq27_8(j,k)*(sub27_8(qq)+0.5d0*sub78_2(gg))
     &   +0.5d0*msq78_2v(j,k)*sub78_2v))

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

c--- DEBUG: 2-quark contribution only for now
c      return

c   33 continue

c--- subtractions for 4-quark pieces
      do j=-4,4
      do k=-4,4

      if ((j  >  0) .and. (k  >  0)) then
        if (j  /=  k) then
c--- Q Q - different flavours
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &  *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
        else
c--- Q Q - same flavours
        msq(1,j,k)=msq(1,j,k)+half*(xn-1d0/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(2,j,k)=msq(2,j,k)+half*(xn-1d0/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
        msq(3,j,k)=msq(3,j,k)+half*(xn-1d0/xn)
     &  *(msq18_2(0,k)*sub18_2(gq)+msq18_2v(0,k)*sub18_2v)
        msq(4,j,k)=msq(4,j,k)+half*(xn-1d0/xn)
     &  *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
        endif

       elseif ((j  <  0).and.(k  <  0)) then
        if (j  /=  k) then
c--- QBAR QBAR - different flavours
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &  *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
        else
c--- QBAR QBAR - same flavours
        msq(1,j,k)=msq(1,j,k)+half*(xn-1d0/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(2,j,k)=msq(2,j,k)+half*(xn-1d0/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
        msq(3,j,k)=msq(3,j,k)+half*(xn-1d0/xn)
     &  *(msq18_2(0,k)*sub18_2(gq)+msq18_2v(0,k)*sub18_2v)
        msq(4,j,k)=msq(4,j,k)+half*(xn-1d0/xn)
     &  *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
        endif

      elseif ((j  >  0).and.(k  <  0)) then
c--- Q QBAR
        if (j  ==  -k) then
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &    *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &    *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
c nf=4d0
        msq(5,j,k)=msq(5,j,k)+tr*(4d0)
     &    *(msq18_7(j,k)*sub78_1(gq)-msq78_1v(j,k)*sub78_1v)
c nf=4d0
        msq(6,j,k)=msq(6,j,k)+tr*(4d0)
     &    *(msq28_7(j,k)*sub78_2(gq)-msq78_2v(j,k)*sub78_2v)
        else
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &    *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &    *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
        endif
c--QBAR Q
      elseif ((j  <  0).and.(k  >  0)) then
        if (j  ==  -k) then
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &    *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &    *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
c nf=4d0
        msq(5,j,k)=msq(5,j,k)+tr*(4d0)
     &    *(msq18_7(j,k)*sub78_1(gq)-msq78_1v(j,k)*sub78_1v)
c nf=4d0
        msq(6,j,k)=msq(6,j,k)+tr*(4d0)
     &    *(msq28_7(j,k)*sub78_2(gq)-msq78_2v(j,k)*sub78_2v)
        else
        msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     &    *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
        msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     &    *(msq28_1(j,0)*sub28_1(gq)+msq28_1v(j,0)*sub28_1v)
      endif

      endif

      enddo
      enddo

      return
      end

