!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_tbbdk_gs(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell                                            *
c     February, 2012.                                                  *
c     calculate the subtraction terms for the process                  *
c                                                                      *
c     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)                *
c     or                                                               *
c     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)                *
c                                                                      *
c     Top (antitop) is kept strictly on-shell                          *
c     although all spin correlations are retained.                     *
c                                                                      *
c     NOTE: this routine is a replacement for qqb_tbb_gs.f that        *
c           includes the effect of the b-quark mass.                   *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'masses.f'
      include 'breit.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq35_4(-nf:nf,-nf:nf),msq45_3(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),
     & sub35_4(4),sub45_3(4),
     & dummyv(-nf:nf,-nf:nf),dsubv
      real(dp):: oldmass2
      external qqb_tbbdk,donothing_gvec

      ndmax=4

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips_mass(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & qqb_tbbdk,donothing_gvec)
      call dips_mass(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     & qqb_tbbdk,donothing_gvec)

c--- dipoles for the final state
      qqproc=.true.
      oldmass2=mass2
      mass2=mt
      call dips_mass(3,p,3,5,4,sub35_4,dsubv,msq35_4,dummyv,
     & qqb_tbbdk,donothing_gvec)
      mass2=mb
      call dips_mass(4,p,4,5,3,sub45_3,dsubv,msq45_3,dummyv,
     & qqb_tbbdk,donothing_gvec)
      mass2=oldmass2

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if      ((j > 0) .and. (k < 0)
     &    .or. (j < 0) .and. (k > 0)) then
        msq(1,j,k)=2._dp*cf*sub15_2(qq)*msq15_2(j,k)
        msq(2,j,k)=2._dp*cf*sub25_1(qq)*msq25_1(j,k)
        msq(3,j,k)=2._dp*cf*sub35_4(qq)*msq35_4(j,k)
        msq(4,j,k)=2._dp*cf*sub45_3(qq)*msq45_3(j,k)
      elseif ((j  /=  0) .and. (k == 0)) then
        msq(2,j,k)=2._dp*tr*sub25_1(qg)*(
     &   msq25_1(j,+1)+msq25_1(j,+2)+msq25_1(j,+3)+msq25_1(j,+4)
     &  +msq25_1(j,+5)+msq25_1(j,-1)+msq25_1(j,-2)+msq25_1(j,-3)
     &  +msq25_1(j,-4)+msq25_1(j,-5))
      elseif ((j == 0) .and. (k  /=  0)) then
        msq(1,j,k)=2._dp*tr*sub15_2(qg)*(
     &   msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)+msq15_2(+4,k)
     &  +msq15_2(+5,k)+msq15_2(-1,k)+msq15_2(-2,k)+msq15_2(-3,k)
     &  +msq15_2(-4,k)+msq15_2(-5,k))
      endif

      enddo
      enddo

      return
      end


