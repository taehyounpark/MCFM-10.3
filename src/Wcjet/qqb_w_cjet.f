!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_w_cjet(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
c----averaged over initial colours and spins
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + cbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ckm.f'
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

      qgWq=  -aveqg*fac*w1cjet(1,5,3,4,2)
      gqWq=  -aveqg*fac*w1cjet(2,5,3,4,1)
      gqbWqb=-aveqg*fac*w1cjet(2,5,4,3,1)
      qbgWqb=-aveqg*fac*w1cjet(1,5,4,3,2)

      do j=-nflav,nflav
      do k=-nflav,nflav
      if (((j == 1) .or. (j == 3)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,-4)*qgWq
      elseif ((j == 0) .and. ((k == +1).or.(k == +3))) then
          msq(j,k)=Vsq(-4,k)*gqWq
      elseif (((j == -1).or.(j == -3)) .and. (k == 0))then
          msq(j,k)=Vsq(j,+4)*qbgWqb
      elseif ((j == 0) .and. ((k == -1).or.(k == -3))) then
          msq(j,k)=Vsq(+4,k)*gqbWqb
      endif

      enddo
      enddo
      return
      end

      function w1cjet(j1,j2,j3,j4,j5)
      implicit none
      include 'types.f'
      real(dp):: w1cjet
c     Matrix element squared for s(1) cbar(2) -> e-(3) nu(4) g(5)
      include 'mxpart.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: j1,j2,j5,j3,j4
      real(dp):: prop,masssq

      prop=((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      masssq=s(j1,j3)+s(j1,j4)+s(j1,j5)+s(j3,j4)+s(j3,j5)+s(j4,j5)
      w1cjet=-2._dp*s(j4,j1)*masssq*(s(j3,j2)+s(j3,j5))/s(j2,j5)/s(j2,j5)
     & -s(j3,j2)/s(j2,j5)
     & *(s(j4,j2)-(s(j4,j1)+s(j4,j5))*(s(j2,j5)+s(j2,j1))/s(j1,j5))
     & -s(j4,j1)/s(j1,j5)
     & *(s(j3,j1)-(s(j3,j2)+s(j3,j5))*(s(j1,j5)+s(j2,j1))/s(j2,j5))
      w1cjet=w1cjet/prop
      return
      end

