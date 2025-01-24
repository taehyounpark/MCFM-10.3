!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine wwggamps(j1,j2,j3,j4,j5,j6,j7,j8,fAB,fBA,Qid)
c     Amplitudes for the process
c     0 --> qbar(j1) + q(j2)+ve(j3)+e^+(p4)+mu-(p5)+vm~(p6)+g(p7)+g(p8)
c     only s34,s56,propagators included.
c     f(jtype,2,2,2) \equiv f(jtype,polg1,polg2,polq)
c     jtype=1 \equiv t-channel
c     jtype=2 \equiv s-channel, z,gamma
c     jtype=3 singly-resonant diagrams

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'srdiags.f'

      integer j1,j2,j3,j4,j5,j6,j7,j8,mplus,minus,polg1,polg2,polq,Qid
      complex(dp):: A8a(2,2),A8b(2,2),A8sr3465(2,2,2),A8sr5643(2,2,2)
      complex(dp):: fAB(3,2,2,2),fBA(3,2,2,2)
      parameter(minus=1,mplus=2)

      fAB(:,:,:,:)=czip
      fBA(:,:,:,:)=czip

c---Handle all cases with negative helicity quark line
c---Remember that the a8tree's are in DKS notation,
c-  so we swap (1<->2) and (5<-->6) to get to our standard notation
      call a8treea(j2,j1,j3,j4,j6,j5,j7,j8,za,zb,a8a)
      call a8treeb(j2,j1,j3,j4,j6,j5,j7,j8,za,zb,a8b)
      do polg1=1,2
      do polg2=1,2
      fAB(1,polg1,polg2,minus)=+a8a(polg1,polg2)
      fAB(2,polg1,polg2,minus)=-a8b(polg1,polg2)
c ---sign?
      enddo
      enddo

c---Other color order swap (7<-->8)
      call a8treea(j2,j1,j3,j4,j6,j5,j8,j7,za,zb,a8a)
      call a8treeb(j2,j1,j3,j4,j6,j5,j8,j7,za,zb,a8b)
      do polg1=1,2
      do polg2=1,2
      fBA(1,polg1,polg2,minus)=+a8a(polg2,polg1)
      fBA(2,polg1,polg2,minus)=-a8b(polg2,polg1)
c ---sign?
      enddo
      enddo

c-- opposite quark helicity by taking complex conjugation (za<->zb)
c-- and swapping all helicities
      call a8treeb(j2,j1,j4,j3,j5,j6,j7,j8,zb,za,a8b)
      do polg1=1,2
      do polg2=1,2
      fAB(1,polg1,polg2,mplus)=czip
      fAB(2,polg1,polg2,mplus)=+a8b(3-polg1,3-polg2)
c      write(6,*) +a8b(3-polg1,3-polg2)
c ---sign?
      enddo
      enddo

c      write(6,*) 'a8b org',a8b
c--- just switch j1 and j2 to change quark line polarization
c      call a8treeb(j1,j2,j3,j4,j6,j5,j8,j7,za,zb,a8b)
c      write(6,*) 'a8b alt',a8b
c      pause
c      do polg1=1,2
c      do polg2=1,2
c      fAB(1,polg1,polg2,mplus)=czip
c      fAB(2,polg1,polg2,mplus)=+a8b(polg2,polg1)
c      write(6,*) +a8b(polg2,polg1)
c ---sign?
c      enddo
c      enddo
c      pause
c---Other color order, swap (7<-->8)
      call a8treeb(j2,j1,j4,j3,j5,j6,j8,j7,zb,za,a8b)
      do polg1=1,2
      do polg2=1,2
      fBA(1,polg1,polg2,mplus)=czip
      fBA(2,polg1,polg2,mplus)=+a8b(3-polg2,3-polg1)
c      write(6,*) +a8b(3-polg2,3-polg1)
c ---sign?
      enddo
      enddo

c      call a8treeb(j1,j2,j3,j4,j6,j5,j7,j8,za,zb,a8b)
c      do polg1=1,2
c      do polg2=1,2
c      fBA(1,polg1,polg2,mplus)=czip
c      fBA(2,polg1,polg2,mplus)=+a8b(polg1,polg2)
c      write(6,*) +a8b(polg1,polg2)
c ---sign?
c      enddo
c      enddo
c      pause

      if (srdiags .eqv. .false.) return

c For down-type amplitudes, no need to permute 3456 -> 5634, so do some
c work to permute them back again;
c Also, add amplitude with 3<->5 and 4<->6, also swapping e and nu couplings
      if (Qid == 2) then
        call srgg(j2,j1,j5,j6,j4,j3,j7,j8,Qid,.false.,A8sr5643)
        call srgg(j2,j1,j3,j4,j6,j5,j7,j8,Qid,.true.,A8sr3465)
      else
        call srgg(j2,j1,j3,j4,j6,j5,j7,j8,Qid,.false.,A8sr5643)
        call srgg(j2,j1,j5,j6,j4,j3,j7,j8,Qid,.true.,A8sr3465)
      endif
      do polg1=1,2
      do polg2=1,2
      do polq=1,2
      fAB(3,polg1,polg2,polq)=A8sr3465(polq,polg1,polg2)
     &                       +A8sr5643(polq,polg1,polg2)
      enddo
      enddo
      enddo

c---Other color order swap (8<-->7)
      if (Qid == 2) then
        call srgg(j2,j1,j5,j6,j4,j3,j8,j7,Qid,.false.,A8sr5643)
        call srgg(j2,j1,j3,j4,j6,j5,j8,j7,Qid,.true.,A8sr3465)
      else
        call srgg(j2,j1,j3,j4,j6,j5,j8,j7,Qid,.false.,A8sr5643)
        call srgg(j2,j1,j5,j6,j4,j3,j8,j7,Qid,.true.,A8sr3465)
      endif
      do polg1=1,2
      do polg2=1,2
      do polq=1,2
      fBA(3,polg1,polg2,polq)=A8sr3465(polq,polg2,polg1)
     &                       +A8sr5643(polq,polg2,polg1)
      enddo
      enddo
      enddo

c Phase is different for jtype=1,2 amplitudes, where we have used za<->zb
c interchange to obtain 'mplus' from 'minus' casess
      fAB(3,:,:,mplus)=-fAB(3,:,:,mplus)
      fBA(3,:,:,mplus)=-fBA(3,:,:,mplus)

      return
      end

