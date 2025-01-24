!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine WW_Hgaga(p,msq)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion by W-W exchange only
c---Matrix element squared averaged over initial colors and spins

c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
c---- Extension to photon decay contributed by Fabian Stoeckli
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: p(mxpart,4),fac
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_du,uub_ddb
      real(dp):: msqgamgam

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      msq(:,:)=0._dp

      call dotem(6,p,s)


      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

      fac=0.25_dp*gwsq**3*hdecay
c Color cancels, 0.25_dp is spin average

c q-q and qbar-qbar
c--- u(1)+d(2) -> d(5)+u(6)
c--- ub(1)+db(2) -> db(5)+ub(6)
      call msqpieces_gaga_ww(1,2,6,5,ud_du)

c q-qbar and qbar-q
c--- u(1)+ub(2) -> d(5)+db(6)
c--- ub(1)+d(2) -> db(5)+u(6)
      call msqpieces_gaga_ww(1,6,2,5,uub_ddb)

c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1
      msq(j,k)=0._dp
        if     ((j > 0) .and. (k < 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j < 0) .and. (k > 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j > 0) .and. (k > 0)) then
          if (pn(j)+pn(k) == +3) msq(j,k)=fac*ud_du
        elseif ((j < 0) .and. (k < 0)) then
          if (pn(j)+pn(k) == -3) msq(j,k)=fac*ud_du
        endif
      enddo
      enddo

      return
      end


      subroutine msqpieces_gaga_ww(i1,i2,i5,i6,wll)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'sprods_com.f'
      real(dp):: wll,htheta
      real(dp):: propw,x
      integer:: i1,i2,i5,i6
c--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propw(i1,i2)=sign(one,(s(i1,i2)-wmass**2))*sqrt(
     &((s(i1,i2)-wmass**2)**2+htheta(s(i1,i2))*(wmass*wwidth)**2)/wmass)
      wll=s(i1,i2)*s(i5,i6)/(propw(i1,i6)*propw(i2,i5))**2
      return
      end
