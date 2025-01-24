!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8fourqsq(j1,j2,j3,j4,j5,j6,j7,j8,C,ff,
     & ampsq_dd,ampsq_du,ampsq_ud,ampsq_uu)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'srdiags.f'
      include 'nf.f'
      include 'zcouple_cms.f'
c     Performs the square and the sum over polarization of the diagrams
c     which have the W+W- radiated off the external line.
c     The routine returns four squared matrix elements, corresponding
c     to all combinations of down and up-type quarks on each line.
c     Because C which is t3456 dependent, is passed as an argument,
c     3,4,5,6 must remain the leptons although they can permuted
      integer j1,j2,j3,j4,j5,j6,j7,j8,jud1,jud2,pol17,pol28
c     jud1,(jud2) specify whether the upper(lower) lines
c     are down-type=1, or up-type=2
      real(dp):: ampsq(2,2),ff,ampsq_dd,ampsq_du,ampsq_ud,ampsq_uu
      complex(dp):: C(2,2),
     & uppera_u(2,2),upperb_u(2,2),lowera_u(2,2),lowerb_u(2,2),
     & uppera_d(2,2),upperb_d(2,2),lowera_d(2,2),lowerb_d(2,2),
     & uppersr1_u(2,2),uppersr2_u(2,2),lowersr1_u(2,2),lowersr2_u(2,2),
     & uppersr1_d(2,2),uppersr2_d(2,2),lowersr1_d(2,2),lowersr2_d(2,2)

c--- call all basic amplitudes
      call a8new(j1,j2,j3,j4,j5,j6,j7,j8,uppera_u,upperb_u)
      call a8new(j1,j2,j5,j6,j3,j4,j7,j8,uppera_d,upperb_d)
      if (srdiags) then
        call sr4q(j1,j2,j5,j6,j4,j3,j7,j8,2,.false.,uppersr1_u)
        call sr4q(j1,j2,j3,j4,j6,j5,j7,j8,2,.true.,uppersr2_u)
        call sr4q(j1,j2,j5,j6,j4,j3,j7,j8,1,.false.,uppersr1_d)
        call sr4q(j1,j2,j3,j4,j6,j5,j7,j8,1,.true.,uppersr2_d)
      else
        uppersr1_u(:,:)=czip
        uppersr2_u(:,:)=czip
        uppersr1_d(:,:)=czip
        uppersr2_d(:,:)=czip
      endif

      call a8new(j2,j1,j3,j4,j5,j6,j8,j7,lowera_u,lowerb_u)
      call a8new(j2,j1,j5,j6,j3,j4,j8,j7,lowera_d,lowerb_d)
      if (srdiags) then
        call sr4q(j2,j1,j5,j6,j4,j3,j8,j7,2,.false.,lowersr1_u)
        call sr4q(j2,j1,j3,j4,j6,j5,j8,j7,2,.true.,lowersr2_u)
        call sr4q(j2,j1,j5,j6,j4,j3,j8,j7,1,.false.,lowersr1_d)
        call sr4q(j2,j1,j3,j4,j6,j5,j8,j7,1,.true.,lowersr2_d)
      else
        lowersr1_u(:,:)=czip
        lowersr2_u(:,:)=czip
        lowersr1_d(:,:)=czip
        lowersr2_d(:,:)=czip
      endif

c--- loop over all combinations and apply appropriate couplings
      do jud1=1,2
      do jud2=1,2
      ampsq(jud1,jud2)=0d0

      do pol17=1,2
      do pol28=1,2

      if     ((jud1  ==  1) .and. (jud2  ==  1)) then
        ampsq(jud1,jud2)=ampsq(jud1,jud2)+ff*cdabs(
     &   +uppera_d(pol17,pol28)+C(pol17,jud1)*upperb_d(pol17,pol28)
     &   +lowera_d(pol28,pol17)+C(pol28,jud2)*lowerb_d(pol28,pol17)
     &   +(uppersr1_d(pol17,pol28)+uppersr2_d(pol17,pol28)
     &    +lowersr1_d(pol28,pol17)+lowersr2_d(pol28,pol17))*zxw)**2
      elseif ((jud1  ==  1) .and. (jud2  ==  2)) then
        ampsq(jud1,jud2)=ampsq(jud1,jud2)+ff*cdabs(
     &   +uppera_d(pol17,pol28)+C(pol17,jud1)*upperb_d(pol17,pol28)
     &   +lowera_u(pol28,pol17)+C(pol28,jud2)*lowerb_u(pol28,pol17)
     &   +(uppersr1_d(pol17,pol28)+uppersr2_d(pol17,pol28)
     &    +lowersr1_u(pol28,pol17)+lowersr2_u(pol28,pol17))*zxw)**2
      elseif ((jud1  ==  2) .and. (jud2  ==  1)) then
        ampsq(jud1,jud2)=ampsq(jud1,jud2)+ff*cdabs(
     &   +uppera_u(pol17,pol28)+C(pol17,jud1)*upperb_u(pol17,pol28)
     &   +lowera_d(pol28,pol17)+C(pol28,jud2)*lowerb_d(pol28,pol17)
     &   +(uppersr1_u(pol17,pol28)+uppersr2_u(pol17,pol28)
     &    +lowersr1_d(pol28,pol17)+lowersr2_d(pol28,pol17))*zxw)**2
      elseif ((jud1  ==  2) .and. (jud2  ==  2)) then
        ampsq(jud1,jud2)=ampsq(jud1,jud2)+ff*cdabs(
     &   +uppera_u(pol17,pol28)+C(pol17,jud1)*upperb_u(pol17,pol28)
     &   +lowera_u(pol28,pol17)+C(pol28,jud2)*lowerb_u(pol28,pol17)
     &   +(uppersr1_u(pol17,pol28)+uppersr2_u(pol17,pol28)
     &    +lowersr1_u(pol28,pol17)+lowersr2_u(pol28,pol17))*zxw)**2
      endif
      enddo
      enddo

      enddo
      enddo

      ampsq_dd=ampsq(1,1)
      ampsq_du=ampsq(1,2)
      ampsq_ud=ampsq(2,1)
      ampsq_uu=ampsq(2,2)

      return
      end
