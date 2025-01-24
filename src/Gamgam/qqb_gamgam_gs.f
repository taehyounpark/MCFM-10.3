!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_gamgam_gs(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell                                            *
c     October, 2002.                                                   *
c     Modified by CW to Include photon fragmentation dipoles Feb 11    *
c    Matrix element SUBTRACTION squared averag'd over init'l colors    *
c    and spins                                                         *
c     f(-p1) + f(-p2) -->  gamma(p3) + gamma(p4) + parton(p5)          *
c***********************************************************************
c======= C. Williams edited March 2013 to be more efficient in dipoles
c======= for fragmentation pieces.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      integer:: j,k,nd
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub15_2v,sub25_1v,
     & msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & msq35_2(-nf:nf,-nf:nf),
     & msq45_2(-nf:nf,-nf:nf),
     & sub35_2,sub45_2
      external qqb_gamgam,donothing_gvec,qqb_dirgam

      if(frag) then
         ndmax=4
      else
         ndmax=2
      endif


c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_gamgam,donothing_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_gamgam,donothing_gvec)

c---- Extra photon dipoles
      if (frag) then
         call dipsfrag(3,p,3,5,2,sub35_2,msq35_2,qqb_dirgam)
         call dipsfrag(4,p,4,5,2,sub45_2,msq45_2,qqb_dirgam)
         do j=3,4
           phot_dip(j)=.true.
         enddo
      endif


      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if(frag) then
         if((j /= 0) .and. (k==0)) then
            msq(3,j,k)=Q(j)**2*msq35_2(j,k)*sub35_2*half
            msq(4,j,k)=Q(j)**2*msq45_2(j,k)*sub45_2*half
         elseif((j==0).and.(k /= 0)) then
            msq(3,j,k)=Q(k)**2*msq35_2(j,k)*sub35_2*half
            msq(4,j,k)=Q(k)**2*msq45_2(j,k)*sub45_2*half
          endif
       endif

       if ((j  /=  0) .and. (k  /=  0) .and. (j /= -k)) goto 19

c--- do only q-qb and qb-q cases
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
         msq(1,j,k)=2._dp*cf*sub15_2(qq)*msq15_2(j,k)
         msq(2,j,k)=2._dp*cf*sub25_1(qq)*msq25_1(j,k)
      elseif ((j  /=  0) .and. (k == 0)) then
         msq(2,j,k)=2._dp*tr*sub25_1(qg)*msq25_1(j,-j)
      elseif ((j == 0) .and. (k  /=  0)) then
         msq(1,j,k)=2._dp*tr*sub15_2(qg)*msq15_2(-k,k)
      elseif ((j == 0) .and. (k == 0)) then
c JC removing gg contribution
c--- Comment out the following four lines to remove gg contribution
c         msq(1,j,k)=2._dp*xn*(sub15_2(gg)*msq15_2(j,k)
c     &                     +sub15_2v*msq15_2v(j,k))
c         msq(2,j,k)=2._dp*xn*(sub25_1(gg)*msq25_1(j,k)
c     &                     +sub25_1v*msq25_1v(j,k))
      endif

 19   continue
      enddo
      enddo

      return
      end

