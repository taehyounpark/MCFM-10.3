!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_w_twdk_gs(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell                                            *
c     January, 2005.                                                   *
c***********************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W + t [or t~](p5) + g(p6)
c                           |
c                            -->l(p3)+a(p4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'nores.f'
      include 'nwz.f'
      integer:: j,k,nd,iq
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq26_5(-nf:nf,-nf:nf),msq16_5(-nf:nf,-nf:nf),
     & msq56_1(-nf:nf,-nf:nf),msq56_2(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq16_2v(-nf:nf,-nf:nf),msq16_5v(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),sub16_5(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub26_1v,sub16_5v,sub16_2v
      external qqb_w_twdk,qqb_w_twdk_gvec,donothing_gvec

c--- Note that the subtractions here are very similar to the ones
c--- for W+1 jet, except that the top mass means that we must
c--- separate the initial-final and final-initial dipoles

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        iq=+1 ! quark initial state
      elseif (nwz == +1) then
        iq=-1 ! antiquark initial state
      else
        write(6,*) 'Error in qqb_w_twdk_gs, nwz is not +1 or -1 :  ',nwz
      stop
      endif

      ndmax=6

      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

c--- calculate the initial-initial dipoles
      call dips_mass(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_w_twdk,qqb_w_twdk_gvec)
      call dips_mass(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_w_twdk,qqb_w_twdk_gvec)

c--- now the initial-final ones
      call dips_mass(5,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & qqb_w_twdk,qqb_w_twdk_gvec)
      call dips_mass(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qqb_w_twdk,qqb_w_twdk_gvec)

c--- now the final-initial ones
      call dips_mass(1,p,5,6,1,sub56_1,sub56_1v,msq56_1,dummyv,
     & qqb_w_twdk,donothing_gvec)
      call dips_mass(2,p,5,6,2,sub56_2,sub56_2v,msq56_2,dummyv,
     & qqb_w_twdk,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=zero
      enddo
      enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf

      if     ((k == 0).and. (j  /=  0)) then
c--- q-g and qb-g cases
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-msq16_5(j,k)*sub16_5(qq)/xn
      msq(1,j,k)=-msq56_1(j,k)*sub56_1(qq)/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v)
      msq(2,j,k)=xn*msq56_2(j,k)*sub56_2(qq)

      elseif ((j == 0).and.(k /= 0)) then
c--- g-q and g-qb cases
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v)
      msq(1,j,k)=xn*msq56_1(j,k)*sub56_1(qq)
      msq(6,j,k)=-msq26_5(j,k)*sub26_5(qq)/xn
      msq(2,j,k)=-msq56_2(j,k)*sub56_2(qq)/xn

      elseif ((j == 0).and.(k == 0)) then
      if (nores .eqv. .false.) then
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(+5,k))*sub16_2(qg)*two*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,+5))*sub26_1(qg)*two*tr
      endif

      endif

      enddo
      enddo

      do j=1,4
c--- subtraction terms for 4-quark matrix elements
      msq(3,j,+5*iq)=two*cf*(
     &          msq16_2(0,+5*iq)*sub16_2(gq)+msq16_2v(0,+5*iq)*sub16_2v)
      msq(4,+5*iq,j)=two*cf*(
     &          msq26_1(+5*iq,0)*sub26_1(gq)+msq26_1v(+5*iq,0)*sub26_1v)

      msq(3,-j,+5*iq)=msq(3,+j,+5*iq)
      msq(4,+5*iq,-j)=msq(4,+5*iq,+j)
      enddo

      return
      end
