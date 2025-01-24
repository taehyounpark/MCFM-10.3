!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wbfromc(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
c----averaged over initial colours and spins
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + bbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ b(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qgWq,qbgWqb,gqbWqb,gqWq,w1cjet

      msq(:,:)=0._dp
c--setup s products through common block
      call dotem(5,p,s)

      fac=gwsq**2*gsq*V

      qgWq=0._dp
      gqWq=0._dp
      qbgWqb=0._dp
      gqbWqb=0._dp
      if (nwz == 1) then
      qgWq=  -aveqg*fac*w1cjet(1,5,3,4,2)
      gqWq=  -aveqg*fac*w1cjet(2,5,3,4,1)
      elseif (nwz == -1) then
      gqbWqb=-aveqg*fac*w1cjet(2,5,4,3,1)
      qbgWqb=-aveqg*fac*w1cjet(1,5,4,3,2)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if (((j == 2) .or. (j == 4)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,-5)*qgWq
      elseif ((j == 0) .and. ((k == +2).or.(k == +4))) then
          msq(j,k)=Vsq(-5,k)*gqWq
      elseif (((j == -2).or.(j == -4)) .and. (k == 0))then
          msq(j,k)=Vsq(j,+5)*qbgWqb
      elseif ((j == 0) .and. ((k == -2).or.(k == -4))) then
          msq(j,k)=Vsq(+5,k)*gqbWqb
      endif
      enddo
      enddo
      return
      end



