!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8qid(j1,j2,j3,j4,j5,j6,j7,j8,C,jud1,jud2,ff,ampsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'srdiags.f'
      include 'nf.f'
      include 'zcouple_cms.f'
c     Performs the square and the sum over polarization of the diagrams
c     which have the W+W- radiated off the same external line.
c     Because C, which is t3456 dependent, is passed as an argument,
c     3,4,5,6 must remain the leptons although they can permuted
      integer j1,j2,j3,j4,j5,j6,j7,j8,jud1,jud2,pol17,pol28
c     jud1,(jud2) specify whether the upper(lower) lines
c     are down-type=1, or up-type=2
c     quarks in the final state are identical so there is interference.
      real(dp):: ampsq,ff
c      real(dp):: bit1(2,2),bit2(2,2)
      complex(dp):: C(2,2),uppera(2,2),upperb(2,2),amp(2),
     & lowera(2,2),lowerb(2,2)
      complex(dp):: uppesa(2,2),uppesb(2,2),
     & lowesa(2,2),lowesb(2,2),
     & uppersr1(2,2),uppersr2(2,2),lowersr1(2,2),lowersr2(2,2),
     & uppessr1(2,2),uppessr2(2,2),lowessr1(2,2),lowessr2(2,2)

      ampsq=0d0
      call a8new(j1,j2,j3,j4,j5,j6,j7,j8,uppera,upperb)
      call a8new(j2,j1,j3,j4,j5,j6,j8,j7,lowera,lowerb)
      call a8new(j1,j2,j3,j4,j5,j6,j8,j7,uppesa,uppesb)
      call a8new(j2,j1,j3,j4,j5,j6,j7,j8,lowesa,lowesb)

      if (srdiags) then
c Note: just use 3,4,5,6 here instead of j3,j4,j5,j6 because we do not
c need additional interchange that is used in call to facilitate u<->d
        call sr4q(j1,j2,5,6,4,3,j7,j8,jud1,.false.,uppersr1)
        call sr4q(j1,j2,3,4,6,5,j7,j8,jud1,.true.,uppersr2)
        call sr4q(j2,j1,5,6,4,3,j8,j7,jud2,.false.,lowersr1)
        call sr4q(j2,j1,3,4,6,5,j8,j7,jud2,.true.,lowersr2)
        call sr4q(j1,j2,5,6,4,3,j8,j7,jud1,.false.,uppessr1)
        call sr4q(j1,j2,3,4,6,5,j8,j7,jud1,.true.,uppessr2)
        call sr4q(j2,j1,5,6,4,3,j7,j8,jud2,.false.,lowessr1)
        call sr4q(j2,j1,3,4,6,5,j7,j8,jud2,.true.,lowessr2)
      else
        uppersr1(:,:)=czip
        uppersr2(:,:)=czip
        lowersr1(:,:)=czip
        lowersr2(:,:)=czip
        uppessr1(:,:)=czip
        uppessr2(:,:)=czip
        lowessr1(:,:)=czip
        lowessr2(:,:)=czip
      endif

      do pol17=1,2
      do pol28=1,2

      amp(1)=uppera(pol17,pol28)+C(pol17,jud1)*upperb(pol17,pol28)
     &      +lowera(pol28,pol17)+C(pol28,jud2)*lowerb(pol28,pol17)
     &      +(uppersr1(pol17,pol28)+uppersr2(pol17,pol28)
     &       +lowersr1(pol28,pol17)+lowersr2(pol28,pol17))*zxw
      amp(2)=uppesa(pol17,pol28)+C(pol17,jud1)*uppesb(pol17,pol28)
     &      +lowesa(pol28,pol17)+C(pol28,jud2)*lowesb(pol28,pol17)
     &      +(uppessr1(pol17,pol28)+uppessr2(pol17,pol28)
     &       +lowessr1(pol28,pol17)+lowessr2(pol28,pol17))*zxw

      if (pol17  ==  pol28) then
      ampsq=ampsq+ff*(cdabs(amp(1))**2+cdabs(amp(2))**2
     & +2d0/xn*dble(amp(1)*dconjg(amp(2))))
c      bit1(pol17,pol28)=ff*(cdabs(amp(1))**2+cdabs(amp(2))**2
c     & +2d0/xn*dble(amp(1)*dconjg(amp(2))))

c      write(6,*) pol17,pol28,bit1(pol17,pol28)
      else
      ampsq=ampsq+ff*(cdabs(amp(1))**2+cdabs(amp(2))**2)
c      bit1(pol17,pol28)=ff*cdabs(amp(1))**2
c      bit2(pol17,pol28)=ff*cdabs(amp(2))**2
c      write(6,*) pol17,pol28,bit1(pol17,pol28)
c      write(6,*) pol17,pol28,bit2(pol17,pol28)
      endif

c      write(6,*) 'a8qid:',pol17,pol28,bit(pol17,pol28)
      enddo
      enddo
      return
      end
