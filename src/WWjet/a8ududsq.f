!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8ududsq(j1,j2,j3,j4,j5,j6,j7,j8,C,jud1,jud2,ff,
     & ampsqid,ampsqnotid)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'srdiags.f'
      include 'nf.f'
      include 'zcouple_cms.f'
c     Performs the square and the sum over polarization of the diagrams
c     which have the W+W- radiated off each line, plus the diagrams
c     where the W+ is radiated from the upper and W- from the lower
c     (returned in ampsqid)
c     Diagrams with a W from each line only are also returned
c     separately in ampsqnotid
c     Because C which is t3456 dependent, is passed as an argument,
c     3,4,5,6 must remain the leptons although they can permuted
      integer j1,j2,j3,j4,j5,j6,j7,j8,jud1,jud2,pol17,pol28
c     jud1,(jud2) specify whether the upper(lower) lines
c     are down-type=1, or up-type=2
      real(dp):: ampsqid,ampsqnotid,ff
      complex(dp):: C(2,2),uppera(2,2),upperb(2,2),
     & lowera(2,2),lowerb(2,2),wpwmamp,extra,
     & uppersr1(2,2),uppersr2(2,2),lowersr1(2,2),lowersr2(2,2)

      if (jud1  ==  2) then
      call a8new(j1,j2,j3,j4,j5,j6,j7,j8,uppera,upperb)
      elseif (jud1  == 1) then
      call a8new(j1,j2,j5,j6,j3,j4,j7,j8,uppera,upperb)
      endif
      if (srdiags) then
        call sr4q(j1,j2,j5,j6,j4,j3,j7,j8,jud1,.false.,uppersr1)
        call sr4q(j1,j2,j3,j4,j6,j5,j7,j8,jud1,.true.,uppersr2)
      else
        uppersr1(:,:)=czip
        uppersr2(:,:)=czip
      endif

      if (jud2  ==  2) then
      call a8new(j2,j1,j3,j4,j5,j6,j8,j7,lowera,lowerb)
      elseif (jud2 == 1) then
      call a8new(j2,j1,j5,j6,j3,j4,j8,j7,lowera,lowerb)
      endif
      if (srdiags) then
        call sr4q(j2,j1,j5,j6,j4,j3,j8,j7,jud2,.false.,lowersr1)
        call sr4q(j2,j1,j3,j4,j6,j5,j8,j7,jud2,.true.,lowersr2)
      else
        lowersr1(:,:)=czip
        lowersr2(:,:)=czip
      endif

c--- the additional diagrams involving a W radiated off each line
c--- only contribute to the left-left polarization. The interference
c--- contribution is suppressed by an additional factor due to colour
c--- (down by -1/N). The sign difference should probably be in wpwmamp.f
      ampsqid=0d0
      do pol17=1,2
      do pol28=1,2
      if ((pol17  ==  1) .and. (pol28  ==  1)) then
        extra=wpwmamp(j1,j2,j3,j4,j5,j6,j8,j7)
        ampsqnotid=ff*cdabs(extra)**2
      else
        extra=czip
      endif
      ampsqid=ampsqid+ff*(
     & cdabs(+uppera(pol17,pol28)+C(pol17,jud1)*upperb(pol17,pol28)
     &       +lowera(pol28,pol17)+C(pol28,jud2)*lowerb(pol28,pol17)
     &       +(uppersr1(pol17,pol28)+uppersr2(pol17,pol28)
     &        +lowersr1(pol28,pol17)+lowersr2(pol28,pol17))*zxw
     &       )**2
     &+cdabs(extra)**2
     &+2d0/xn*dble(dconjg(extra)*
     &      (+uppera(pol17,pol28)+C(pol17,jud1)*upperb(pol17,pol28)
     &       +lowera(pol28,pol17)+C(pol28,jud2)*lowerb(pol28,pol17)
     &       +(uppersr1(pol17,pol28)+uppersr2(pol17,pol28)
     &        +lowersr1(pol28,pol17)+lowersr2(pol28,pol17))*zxw)))
      enddo
      enddo

      return
      end
